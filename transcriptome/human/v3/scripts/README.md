# HUMAN TRANSCRIPTOME REFERENCE AND ANNOTATION

## APPRIS Isoform Selection Logic (`v3/scripts/generate_appris_reference.py`)

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

**Outputs:**
- Selected and discarded transcript FASTA files (gzipped)
- BED file annotating UTR5, CDS, and UTR3 regions for each selected transcript
- TSV of transcript lengths for selected transcripts

## Driver script (`generate_riboflow_ref_and_annot.sh`)

Wrapper that runs the full v3 build end-to-end from this directory:

1. Calls `generate_appris_reference.py` with the v48 GENCODE FASTA + GTF and the APPRIS principal list from `../raw/sequence/` and `../raw/annotation/`.
2. Writes the selected/discarded FASTAs, BED, and length TSV into `../riboflow_annot_and_ref/` under the `appris_human_v3` basename.
3. Runs `bowtie2-build` on the selected FASTA to produce the `.bt2` index files RiboFlow consumes.

## Notebooks

- [`mapping_comparison_of_v2_and_v3.ipynb`](mapping_comparison_of_v2_and_v3.ipynb) — v2-vs-v3 ribo-seq comparison on GSM1606107 / GSM1606108. Consumes `v2.ribo` and `v3.ribo` at the repo root; full conclusions are rendered in the notebook.
- [`v3_appris_explore.ipynb`](v3_appris_explore.ipynb) — exploratory sanity checks on the APPRIS isoform-selection output (transcript counts, MANE coverage, divergence vs. v2).

**Environment note.** `ribopy` depends on the removed `pipes` stdlib module and will not import on Python ≥ 3.13. For these notebooks we use a dedicated Python 3.9 venv at the repo root (`.venv-ribo/`); launch with `.venv-ribo/bin/jupyter lab` from the repo root.
