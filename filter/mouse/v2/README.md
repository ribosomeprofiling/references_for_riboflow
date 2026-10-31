# README

Filter sequences for the mouse v2 RiboFlow filter. RiboFlow aligns reads to these sequences before transcriptome alignment, and discards aligned reads.

`scripts/download_sources.R` fetches the NCBI RefSeq `.fa` files into `raw/`, and `scripts/build_filter_v2.py` builds `mouse_rtRNA.fa.gz` and its Bowtie2 index `mouse_rtRNA.*.bt2` from the `raw/` files. The GtRNAdb mm39 tRNA file (`raw/mm39-tRNAs.fa.gz`) is downloaded from GtRNAdb separately and is already present in `raw/` (the build accepts plain or gzipped).

 v2 uses NCBI RefSeq sequences for rRNAs, mitochondrial rRNAs and mitochondrial tRNAs. Nuclear rRNA uses the curated RefSeq **transcripts** (the transcribed rRNA)

Only GtRNAdb **high-confidence** nuclear tRNAs are included.

v1 used genomic rDNA repeat (`BK000964.3`) and the GtRNAdb high-confidence nuclear tRNAs. v2 replaces the rDNA repeat with the RefSeq rRNA transcripts and adds the mitochondrial rRNA/tRNA genes.

## Sources

| File                       | Source                                                                                                                                                                                        |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ncbi_rRNA.fa`           | NCBI RefSeq rRNA transcripts: `NR_046233.2` (45S pre-rRNA, Rn45s), `NR_003278.3` (18S, Rn18s), `NR_003279.1` (28S, Rn28s1), `NR_003280.2` (5.8S, Rs5-8s1), `NR_030686.1` (5S, Rn5s) |
| `ncbi_mt.fa`             | NCBI RefSeq mouse mitochondrion `NC_005089.1`, rRNA + tRNA features extracted by coordinate                                                                                                 |
| `mm39-tRNAs.fa.gz`       | GtRNAdb mm39 (GRCm39) high-confidence tRNA set, genomic sequence (introns included)                                                                                                           |

## Contents (436 sequences)

| Header           | N Sequences | Type                                                                       |
| ---------------- | ----------- | -------------------------------------------------------------------------- |
| `ncbi_rRNA`    | 5           | 45S pre-rRNA (Rn45s), 18S (Rn18s), 28S (Rn28s1), 5.8S (Rs5-8s1), 5S (Rn5s) |
| `ncbi_mt`      | 24          | mt-Rnr1 (12S), mt-Rnr2 (16S), and the 22 mitochondrial tRNAs               |
| `gtrnadb_tRNA` | 407         | nuclear tRNA genes, GtRNAdb mm39 high-confidence set                       |

## Differences from v1

v1 has 409 sequences (225 unique) and v2 has 436 (253 unique). They share 224 unique sequences (the nuclear tRNAs, identical between v1's mm10 build and v2's mm39 build). v1's genomic rDNA repeat is replaced by RefSeq rRNA transcripts.

| Set        | N (unique) | Content                                                                                |
| ---------- | ---------- | -------------------------------------------------------------------------------------- |
| in both    | 224        | GtRNAdb high-confidence nuclear tRNAs                                             |
| only in v1 | 1          | the genomic rDNA repeat `BK000964.3` (replaced in v2 by the RefSeq rRNA transcripts) |
| only in v2 | 29         | 5 RefSeq nuclear rRNA transcripts (45S/18S/28S/5.8S/5S) and 24 mitochondrial rRNA/tRNA |
