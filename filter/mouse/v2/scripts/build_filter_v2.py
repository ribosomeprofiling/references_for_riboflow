#!/usr/bin/env python3

import gzip
import re
import subprocess
from pathlib import Path

RAW = Path("../raw")
OUT_BASE = Path("../mouse_rtRNA")


def resolve(name):
    for cand in (RAW / name, RAW / (name + ".gz")):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"{name}[.gz] not found in {RAW}")


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    header, seq = None, []

    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())

    if header is not None:
        yield header, "".join(seq)


records = []

for h, s in read_fasta(resolve("ncbi_rRNA.fa")):
    acc = h.split()[0]
    m = re.search(r"\(([^)]*)\)", h)            # gene symbol in parens, e.g. (Rn18s)
    sym = m.group(1) if m else acc              # fall back to the accession if unnamed
    records.append((f"ncbi_rRNA|{acc}|{sym}", s))

for h, s in read_fasta(resolve("ncbi_mt.fa")):
    loc, sym, _ = h.split()
    records.append((f"ncbi_mt|{loc}|{sym}", s))

for fname, prefix in (("mm39-tRNAs.fa", "gtrnadb_tRNA"),):
    for h, s in read_fasta(resolve(fname)):
        tid = re.search(r"ID: ([^)]*)\)", h).group(1)
        records.append((f"{prefix}|{h.split()[0].removeprefix('Mus_musculus_')}|{tid}", s))

names = [n for n, _ in records]
assert len(names) == len(set(names)), "duplicate sequence names"

with gzip.open(OUT_BASE.with_suffix(".fa.gz"), "wt") as out:
    for n, s in records:
        out.write(f">{n}\n")
        for i in range(0, len(s), 60):
            out.write(s[i:i + 60] + "\n")

counts = {}
for n in names:
    counts[n.split("|")[0]] = counts.get(n.split("|")[0], 0) + 1

for k, v in sorted(counts.items()):
    print(f"{k:16} {v}")
print(f"total: {len(records)}")

subprocess.run(["bowtie2-build", "-q", str(OUT_BASE.with_suffix(".fa.gz")), str(OUT_BASE)], check=True)
print(sorted(p.name for p in OUT_BASE.parent.glob("mouse_rtRNA*")))
