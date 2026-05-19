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
