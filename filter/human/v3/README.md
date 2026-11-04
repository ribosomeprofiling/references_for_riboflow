# README

Filter sequences for the human v3 RiboFlow Filter. RiboFlow aligns reads to these sequences before transcriptome alignment, and discards aligned reads.

`scripts/build_filter_v3.py` builds `human_rtRNA.fa.gz` and its Bowtie2 index `human_rtRNA.*.bt2` from the files in `raw/`:

* `ncbi_rRNA.fa.gz`: curated human RefSeq rRNAs, from the NCBI nuccore search `txid9606[orgn:noexp] AND biomol_rrna[prop] AND srcdb_refseq_known[prop]`
* `ncbi_mt.fa.gz`: the 2 rRNA and 22 tRNA genes of the mitochondrial genome `NC_012920.1`, cut from the genome at their NCBI Gene coordinates
* `gtrnadb_hg38-tRNAs.fa.gz` and `gtrnadb_hg38-filtered-tRNAs.fa.gz`: from https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz

v2 did not contain mitochondrial rRNA/tRNA genes, whereas v3 does. We update the tRNA sequences from hg19 UCSC, to the GtRNAdb hg38 reference. We retain the three small-RNA pseudogenes and repeat elements from v2.

v3 uses NCBI RefSeq sequences for rRNAs, mitochondrial rRNAs and mitochondrial tRNAs.

GtRNAdb **high-confidence** nuclear tRNAs are included.

## Contents (498 sequences)

| Header           | N Sequences | Type                                                                                                           |
| ---------------- | ----------- | -------------------------------------------------------------------------------------------------------------- |
| `ncbi_rRNA`    | 37          | 45S, 18S, 28S, 5.8S for rDNA copies N1–N5 (20), 5S genes RNA5S1–17 (17)                                      |
| `ncbi_mt`      | 24          | MT-RNR1 (12S), MT-RNR2 (16S), the 22 mitochondrial tRNAs                                                       |
| `gtrnadb_tRNA` | 432         | nuclear tRNA genes, GtRNAdb high-confidence set, genomic sequence (introns included)                           |
| `v2`           | 5           | RNU4-46Px, RNU6-1330P, RN7SKP176 (U4/U6/7SK pseudogenes), RepeatMasker AluY and L1MB5 (see v2 README for loci) |

## Differences from v2

v2 has 635 sequences (469 unique) and v3 has 498 (313 unique). v3 is a rebuild, not a superset: the nuclear tRNAs moved from the hg19 UCSC set to GtRNAdb hg38, and the rRNAs come from NCBI RefSeq rather than v2's five hand-picked sequences. See `scripts/compare_v2_and_v3.ipynb` for the full comparison.

By unique sequence:

| Set        | N (unique) | Content                                                                                                |
| ---------- | ---------- | ------------------------------------------------------------------------------------------------------ |
| in both    | 261        | nuclear tRNAs identical between the hg19 and hg38 builds, plus the 5 snRNA-pseudogene/repeat additions |
| only in v2 | 208        | 204 hg19 tRNAs and 4 of v2's hand-picked rRNAs                                                         |
| only in v3 | 52         | 24 mitochondrial rRNA/tRNA (new), 15 RefSeq nuclear rRNAs, and 13 hg38 tRNAs new to the GtRNAdb set    |
