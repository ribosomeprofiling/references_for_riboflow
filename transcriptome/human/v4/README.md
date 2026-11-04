# README

Version 4 of the RiboFlow human transcriptome: one transcript per protein-coding gene, 19,336 genes.

## Sources

- `raw/annotation/appris_data.principal.txt`: APPRIS principal isoforms for GENCODE Release 50, from https://appris.bioinfo.cnio.es/#/downloads
- `raw/annotation/gencode.v50.annotation.gtf.gz` and `raw/sequence/gencode.v50.pc_transcripts.fa.gz`: GENCODE Release 50 (GRCh38.p14), from https://www.gencodegenes.org/human/

The two GENCODE files are over GitHub's 100 MB limit, so they are not in the repository. Download them into `raw/` before building.

## Building

Run `scripts/generate_riboflow_ref_and_annot.sh` from `scripts/`. It runs `generate_appris_reference.py` and builds the Bowtie2 index in `riboflow_annot_and_ref/`. It also writes `appris_human_v4_discarded.fa.gz`, the transcripts that were not selected. That file is not in the repository either.

## Changes from v3

Transcripts are picked the same way as in v3. Three new steps remove duplicated sequence, so fewer reads are discarded as multimappers:

1. Read-through genes are removed.
2. Transcripts that lie inside a longer transcript are removed.
3. Bicistronic transcripts (one RNA annotated as two genes) keep only the gene with the longer CDS.

See `scripts/v4_appris_explore.ipynb` for the selection step by step, and `scripts/mapping_comparison_of_v3_and_v4.ipynb` for the mapping comparison with v3.
