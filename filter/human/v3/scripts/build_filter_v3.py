#!/usr/bin/env python3
import gzip
import re
import subprocess
from pathlib import Path

RAW = Path("../raw")
OUT_BASE = Path("../human_rtRNA")
V2_FASTA = Path("../../v2/human_rtRNA.fa.gz")
V2_KEEP = re.compile(r"^hg38_(RNU4-46Px|RNU6-1330P|RN7SKP176|rmsk_AluY|rmsk_L1MB5)")


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

for h, s in read_fasta(RAW / "ncbi_rRNA.fa.gz"):
    records.append((f"ncbi_rRNA|{h.split()[0]}|{re.search(r'\(([^)]*)\)', h).group(1)}", s))

for h, s in read_fasta(RAW / "ncbi_mt.fa.gz"):
    loc, sym, _ = h.split()
    records.append((f"ncbi_mt|{loc}|{sym}", s))

for fname, prefix in (("gtrnadb_hg38-tRNAs.fa.gz", "gtrnadb_tRNA"),):
    for h, s in read_fasta(RAW / fname):
        tid = re.search(r"ID: ([^)]*)\)", h).group(1)
        records.append((f"{prefix}|{h.split()[0].removeprefix('Homo_sapiens_')}|{tid}", s))

v2 = [(h, s) for h, s in read_fasta(V2_FASTA) if V2_KEEP.match(h)]
assert len(v2) == 5, f"expected 5 v2 additions, found {len(v2)}"

for h, s in v2:
    locus = re.search(r"range=(\S+)", h).group(1)
    strand = re.search(r"strand=(\S)", h).group(1)
    records.append((f"v2|{h.split()[0]}|{locus}({strand})", s))

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
print(sorted(p.name for p in OUT_BASE.parent.glob("human_rtRNA*")))
