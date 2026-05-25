# Transcriptome reference, annotation and filter sequences for RiboFlow

In this repository, you can find RiboFlow reference and annotation files. These files can be used directly in RiboFlow without any additional preparation or processing.

## Transcriptome Reference and Annotation

Transcriptome reference and annotation files are in the subfolder 

**transcriptome/organism/v[0-9]/riboflow_annot_and_ref**

In particular, the contents are as follows.

   * Bowtie2 reference files (.bt2)
   * Transcriptome annotation (.bed)
   * Transcript Lengths (.tsv)


## Filter Reference

Before aligning the reads to the transcriptome, RiboFlow filters out sequences that map to the filter sequences. 
Typically, these sequences consist of ribosomal RNAs and tRNAs.
Bowtie2 references and sequence files can be found in the subfolder 

**filter/organism/v[0-9]/**

## Organisms

We provide reference and annotation for the following organisms. Depending on the demand, we may expand this list in the future.

  * Human — versions `v1`, `v2`, `v3` (current default: `v3`)
  * Mouse

## Comparison notebooks

The v2-vs-v3 human-reference comparison (ribo-seq QC and CDS-occupancy concordance on GSM1606107 / GSM1606108) lives at
[`transcriptome/human/v3/scripts/mapping_comparison_of_v2_and_v3.ipynb`](transcriptome/human/v3/scripts/mapping_comparison_of_v2_and_v3.ipynb).

