
# README

## Annotation Files
### Appris
- appris_data.principal.txt

Appris Version 2025_08.v50
Appris list downloaded from : https://appris.bioinfo.cnio.es/#/downloads
on April 1 2026

### GENCODE:
- gencode.v48.annotation.gtf.gz
- gencode.v48.pc_transcripts.fa.gz

All other files have been downloaded from gencode website: https://www.gencodegenes.org/human/
Version: Release 48 (GRCh38.p14)

## Filtering Transcripts

**Input filtering — PRINCIPAL transcripts only:**
- Reads the APPRIS annotation file and retains only transcripts with a `PRINCIPAL` category
- Parses the numeric level from `PRINCIPAL:N` (e.g., `PRINCIPAL:1` → level 1); non-numeric levels are assigned level 99
- Also tracks which transcripts are labeled `MANE_Select`

**Per-gene best-transcript selection (in priority order):**
1. `PRINCIPAL:1` + `MANE_Select` — ties broken by longest transcript
2. `PRINCIPAL:2` + `MANE_Select` — ties broken by longest transcript
3. `PRINCIPAL:1` only — ties broken by longest transcript
4. `PRINCIPAL:2` only — ties broken by longest transcript
5. `MANE_Select` at any PRINCIPAL level — ties broken by longest transcript
6. Longest PRINCIPAL transcript at any level (last resort)

**Post-selection filters:**
- Retain only transcripts annotated as `protein_coding` gene type (from GTF)
- Retain only transcripts with a valid CDS and both a proper start codon (ATG or single-nucleotide Hamming-1 variants: CTG, GTG, TTG, AAG, ACG, AGG, ATA, ATC, ATT) and a proper stop codon (TAG, TAA, TGA)

## Bowtie2 Indexes

`bowtie2-build ${SELECTED_FASTA} ${OUTPUT_BASE}_selected`

## Differences with transcriptome/human/v2

v3 is built with the APPRIS selection logic described above. To compare it
against the older v2 reference we processed
[GSM1606107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1606107) and
[GSM1606108](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1606108)
through RiboFlow against both references (read-length window 28–32 nt). The
resulting ribo files (`v2.ribo`, `v3.ribo` at the repo root) drive the Jupyter
notebook
[`v3/scripts/mapping_comparison_of_v2_and_v3.ipynb`](v3/scripts/mapping_comparison_of_v2_and_v3.ipynb).

**Headline findings from that notebook:**
- QC metrics (length distribution, start/stop metagenes, region splits) are visually indistinguishable between v2 and v3.
- Across the 14,973 genes present in both references, log2 CDS occupancy agrees with Pearson r ≈ 0.994 for both samples.
- v2 has 19,734 aliased CDS rows vs. 19,871 for v3 — 4,761 v2-only and 4,898 v3-only entries, reflecting APPRIS isoform-selection differences.
- Only ~55–70 genes per sample show |Δ| > 100 reads. The top divergent genes (`SYNCRIP-202`, `RPSA-201`, `PPP2CA-203`, `IMPDH2-201`, `SKP1-202`) are paralog/homolog mappability redistributions, not biological signal.
- Total CDS reads drop ~0.5% from v2 to v3.

In short: v3 is a safe drop-in replacement for v2; gene-level results that depend on one of the divergent genes above should be cross-checked against v2.