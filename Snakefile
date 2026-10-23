import csv
import os


# ============================================================
# Configuration
# ============================================================

configfile: "config/config.yaml"

SAMPLES_FILE = config["samples"]

samples = {}

with open(SAMPLES_FILE) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        samples[row["sample"]] = row["pod5"]

SAMPLE_NAMES = list(samples.keys())


# ============================================================
# Default automated workflow endpoint
#
# The automated workflow produces:
#   1. Autocycler consensus assembly
#   2. Independent Plassembler reconstruction
#
# These outputs must be manually reviewed before final genome
# configuration. Chromosome rotation with dnaapler is therefore
# intentionally excluded from the default target.
# ============================================================

rule all:
    input:
        expand(
            "results/05_assembly/{sample}/autocycler/consensus_assembly.fasta",
            sample=SAMPLE_NAMES
        ),
        expand(
            "results/05_assembly/{sample}/plassembler",
            sample=SAMPLE_NAMES
        )


# ============================================================
# 1. Dorado SUP basecalling
# ============================================================

rule dorado_basecall:
    input:
        pod5=lambda wildcards: samples[wildcards.sample]

    output:
        directory("results/01_basecalling/{sample}")

    params:
        dorado=config["dorado"]["executable"],
        model=config["dorado"]["model"],
        models_dir=config["dorado"]["models_directory"],
        device=config["dorado"]["device"]

    threads: 8

    shell:
        r"""
        mkdir -p {output}

        {params.dorado} basecaller \
            {params.model} \
            {input.pod5} \
            --device {params.device} \
            --models-directory {params.models_dir} \
            --output-dir {output} \
            --emit-summary
        """


# ============================================================
# 2. Convert Dorado BAM output to FASTQ
# ============================================================

rule dorado_to_fastq:
    conda: "envs/read_qc.yaml"
    input:
        basecalling="results/01_basecalling/{sample}"

    output:
        "results/02_read_qc/fastq/{sample}.fastq.gz"

    threads: 4

    shell:
        r"""
        mkdir -p $(dirname {output})

        samtools fastq \
            -@ {threads} \
            results/01_basecalling/{wildcards.sample}/*.bam \
        | gzip -c > {output}
        """


# ============================================================
# 3. Filter reads
#
# Keep reads:
#   mean Q >= 15
#   length >= 300 bp
# ============================================================

rule filtlong:
    conda: "envs/read_qc.yaml"
    input:
        "results/02_read_qc/fastq/{sample}.fastq.gz"

    output:
        "results/02_read_qc/filtered/{sample}.Q15.L300.fastq.gz"

    params:
        min_length=config["filtlong"]["min_length"],
        min_q=config["filtlong"]["min_mean_q"]

    shell:
        r"""
        mkdir -p $(dirname {output})

        filtlong \
            --min_length {params.min_length} \
            --min_mean_q {params.min_q} \
            {input} \
        | gzip -c > {output}
        """


# ============================================================
# 4. Estimate genome size with Lrge
# ============================================================

rule lrge:
    conda: "envs/lrge.yaml"
    input:
        "results/02_read_qc/filtered/{sample}.Q15.L300.fastq.gz"

    output:
        "results/03_genome_size/{sample}/lrge.txt"

    params:
        platform=config["lrge"]["platform"]

    threads:
        config["lrge"]["threads"]

    shell:
        r"""
        mkdir -p $(dirname {output})

        lrge \
            --platform {params.platform} \
            --threads {threads} \
            --output {output} \
            {input}
        """


# ============================================================
# 5. Autocycler read subsampling
# ============================================================

rule autocycler_subsample:
    conda: "envs/autocycler.yaml"
    input:
        reads="results/02_read_qc/filtered/{sample}.Q15.L300.fastq.gz",
        genome_size="results/03_genome_size/{sample}/lrge.txt"

    output:
        yaml="results/04_subsampling/{sample}/subsample.yaml",
        reads=expand(
            "results/04_subsampling/{{sample}}/sample_{n:02d}.fastq",
            n=range(1, config["autocycler"]["subsamples"] + 1)
        )

    params:
        outdir="results/04_subsampling/{sample}",
        n_subsamples=config["autocycler"]["subsamples"],
        seed=config["autocycler"]["seed"],
        min_depth=config["autocycler"]["min_read_depth"]

    shell:
        r"""
        GENOME_SIZE=$(cat {input.genome_size})

        autocycler subsample \
            --reads {input.reads} \
            --out_dir {params.outdir} \
            --genome_size "$GENOME_SIZE" \
            --count {params.n_subsamples} \
            --min_read_depth {params.min_depth} \
            --seed {params.seed}
        """


# ============================================================
# 5. Candidate assemblies for Autocycler
# ============================================================

ASSEMBLERS = ["flye", "hifiasm", "raven", "miniasm"]
SUBSAMPLES = [
    f"sample_{n:02d}"
    for n in range(1, config["autocycler"]["subsamples"] + 1)
]



rule autocycler_helper:
    conda: "envs/autocycler.yaml"
    input:
        reads="results/04_subsampling/{sample}/{subsample}.fastq",
        genome_size="results/03_genome_size/{sample}/lrge.txt"

    output:
        status="results/05_assembly/{sample}/status/{subsample}_{assembler}.status.tsv"

    params:
        out_prefix="results/05_assembly/{sample}/assemblies/{sample}_{subsample}_{assembler}",
        work_dir="results/05_assembly/{sample}/work/{subsample}_{assembler}"

    threads: 8

    wildcard_constraints:
        assembler="flye|hifiasm|raven|miniasm",
        subsample="sample_[0-9]+"

    shell:
        r"""
        GENOME_SIZE=$(cat {input.genome_size})

        mkdir -p "$(dirname {params.out_prefix})"
        mkdir -p "{params.work_dir}"
        mkdir -p "$(dirname {output.status})"

        set +e

        "autocycler" helper \
            --reads "{input.reads}" \
            --out_prefix "{params.out_prefix}" \
            --genome_size "$GENOME_SIZE" \
            --threads {threads} \
            --read_type ont_r10 \
            --dir "{params.work_dir}" \
            {wildcards.assembler}

        EXIT_CODE=$?

        set -e

        FASTA="{params.out_prefix}.fasta"
        GFA="{params.out_prefix}.gfa"

        if [[ "$EXIT_CODE" -eq 0 && -s "$FASTA" ]]; then
            STATUS="success"
        else
            STATUS="failed"
            rm -f "$FASTA" "$GFA"
        fi

        printf "sample\tsubsample\tassembler\tstatus\texit_code\tfasta\n" > "{output.status}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "{wildcards.sample}" \
            "{wildcards.subsample}" \
            "{wildcards.assembler}" \
            "$STATUS" \
            "$EXIT_CODE" \
            "$FASTA" >> "{output.status}"
        """


def assembly_status_files(wildcards):
    return [
        f"results/05_assembly/{wildcards.sample}/status/{subsample}_{assembler}.status.tsv"
        for subsample in SUBSAMPLES
        for assembler in ASSEMBLERS
    ]


# ============================================================
# 6. Autocycler compression
# ============================================================

rule autocycler_compress:
    conda: "envs/autocycler.yaml"
    input:
        statuses=assembly_status_files

    output:
        gfa="results/05_assembly/{sample}/autocycler/input_assemblies.gfa",
        yaml="results/05_assembly/{sample}/autocycler/input_assemblies.yaml"

    params:
        assemblies_dir="results/05_assembly/{sample}/assemblies",
        autocycler_dir="results/05_assembly/{sample}/autocycler",
        kmer=config["autocycler_consensus"]["kmer"],
        max_contigs=config["autocycler_consensus"]["max_contigs"]

    threads:
        config["autocycler_consensus"]["compress_threads"]

    shell:
        r"""
        "autocycler" compress \
            --assemblies_dir "{params.assemblies_dir}" \
            --autocycler_dir "{params.autocycler_dir}" \
            --kmer {params.kmer} \
            --max_contigs {params.max_contigs} \
            --threads {threads}
        """


# ============================================================
# 7. Autocycler clustering
# ============================================================

rule autocycler_cluster:
    conda: "envs/autocycler.yaml"
    input:
        gfa="results/05_assembly/{sample}/autocycler/input_assemblies.gfa",
        yaml="results/05_assembly/{sample}/autocycler/input_assemblies.yaml"

    output:
        tsv="results/05_assembly/{sample}/autocycler/clustering/clustering.tsv",
        yaml="results/05_assembly/{sample}/autocycler/clustering/clustering.yaml"

    params:
        autocycler_dir="results/05_assembly/{sample}/autocycler"

    threads:
        config["autocycler_consensus"]["cluster_threads"]

    shell:
        r"""
        "autocycler" cluster \
            --autocycler_dir "{params.autocycler_dir}"
        """


# ============================================================
# 8. Trim selected Autocycler clusters
# ============================================================

def selected_clusters(wildcards):
    return config["autocycler_clusters"][wildcards.sample]


rule autocycler_trim:
    conda: "envs/autocycler.yaml"
    input:
        clustering="results/05_assembly/{sample}/autocycler/clustering/clustering.yaml",
        untrimmed="results/05_assembly/{sample}/autocycler/clustering/qc_pass/{cluster}/1_untrimmed.gfa"

    output:
        trimmed="results/05_assembly/{sample}/autocycler/clustering/qc_pass/{cluster}/2_trimmed.gfa"

    params:
        cluster_dir="results/05_assembly/{sample}/autocycler/clustering/qc_pass/{cluster}"

    threads: 8

    shell:
        r"""
        "autocycler" trim \
            --cluster_dir "{params.cluster_dir}" \
            --threads {threads}
        """


# ============================================================
# 9. Resolve selected Autocycler clusters
# ============================================================

rule autocycler_resolve:
    conda: "envs/autocycler.yaml"
    input:
        trimmed="results/05_assembly/{sample}/autocycler/clustering/qc_pass/{cluster}/2_trimmed.gfa"

    output:
        final="results/05_assembly/{sample}/autocycler/clustering/qc_pass/{cluster}/5_final.gfa"

    params:
        cluster_dir="results/05_assembly/{sample}/autocycler/clustering/qc_pass/{cluster}"

    shell:
        r"""
        "autocycler" resolve \
            --cluster_dir "{params.cluster_dir}"
        """


# ============================================================
# 10. Combine resolved clusters
# ============================================================

def resolved_cluster_gfas(wildcards):
    return [
        f"results/05_assembly/{wildcards.sample}/autocycler/clustering/qc_pass/{cluster}/5_final.gfa"
        for cluster in selected_clusters(wildcards)
    ]


rule autocycler_combine:
    conda: "envs/autocycler.yaml"
    input:
        gfas=resolved_cluster_gfas,
        reads="results/02_read_qc/filtered/{sample}.Q15.L300.fastq.gz"

    output:
        fasta="results/05_assembly/{sample}/autocycler/consensus_assembly.fasta",
        gfa="results/05_assembly/{sample}/autocycler/consensus_assembly.gfa",
        yaml="results/05_assembly/{sample}/autocycler/consensus_assembly.yaml"

    params:
        autocycler_dir="results/05_assembly/{sample}/autocycler"

    threads: 8

    shell:
        r"""
        "autocycler" combine \
            --autocycler_dir "{params.autocycler_dir}" \
            --in_gfas {input.gfas} \
            --reads "{input.reads}" \
            --threads {threads}
        """



# ============================================================
# 11. Independent plasmid validation with Plassembler
# ============================================================



rule plassembler:
    conda: "envs/plassembler.yaml"
    input:
        reads="results/02_read_qc/filtered/{sample}.Q15.L300.fastq.gz"

    output:
        directory("results/05_assembly/{sample}/plassembler")

    params:
        database=config["plassembler"]["database"],
        chromosome_length=config["plassembler"]["chromosome_length"]

    threads:
        config["plassembler"]["threads"]

    shell:
        r"""
        "plassembler" long \
            --database "{params.database}" \
            --longreads "{input.reads}" \
            --chromosome {params.chromosome_length} \
            --outdir "{output}" \
            --threads {threads} \
            --prefix "{wildcards.sample}" \
            --skip_qc \
            --keep_fastqs \
            --keep_chromosome
        """


# ============================================================
# 12. Chromosome reorientation
#
# IMPORTANT:
# chromosome.pre_rotation.fasta is a manually curated input.
# Review the Autocycler assembly graph together with the
# independent Plassembler results before creating this file.
# ============================================================

rule dnaapler_chromosome:
    conda: "envs/plassembler.yaml"
    input:
        chromosome="results/05_assembly/{sample}/final_preparation/chromosome.pre_rotation.fasta"

    output:
        fasta="results/05_assembly/{sample}/final_preparation/dnaapler_chromosome/{sample}_chromosome_reoriented.fasta"

    params:
        outdir="results/05_assembly/{sample}/final_preparation/dnaapler_chromosome",
        evalue=config["dnaapler"]["evalue"],
        autocomplete=config["dnaapler"]["autocomplete"],
        seed=config["dnaapler"]["seed"]

    threads:
        config["dnaapler"]["threads"]

    shell:
        r"""
        "dnaapler" chromosome \
            --input "{input.chromosome}" \
            --output "{params.outdir}" \
            --threads {threads} \
            --force \
            --prefix "{wildcards.sample}_chromosome" \
            --evalue {params.evalue} \
            --autocomplete {params.autocomplete} \
            --seed_value {params.seed}
        """
