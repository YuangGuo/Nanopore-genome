#!/usr/bin/env python3

import sys
from pathlib import Path

sample = sys.argv[1]
assemblies_dir = Path(sys.argv[2])
output = Path(sys.argv[3])

assemblers = ["flye", "hifiasm", "raven", "miniasm"]
subsamples = [f"sample_{i:02d}" for i in range(1, 5)]

rows = []

for subsample in subsamples:
    for assembler in assemblers:
        fasta = assemblies_dir / f"{sample}_{subsample}_{assembler}.fasta"

        if fasta.is_file() and fasta.stat().st_size > 0:
            status = "success"
        else:
            status = "failed_or_missing"

        rows.append((sample, subsample, assembler, status, str(fasta)))

output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w") as f:
    f.write("sample\tsubsample\tassembler\tstatus\tfasta\n")
    for row in rows:
        f.write("\t".join(row) + "\n")

successful = sum(row[3] == "success" for row in rows)

print(f"{sample}: {successful}/{len(rows)} successful candidate assemblies")

if successful == 0:
    raise SystemExit(f"ERROR: no successful candidate assemblies found for {sample}")
