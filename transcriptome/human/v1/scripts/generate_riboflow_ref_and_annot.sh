#!/bin/#!/usr/bin/env bash

INPUT_FASTA=../raw/sequence/gencode.v29.transcripts.fa.gz
INPUT_GTF=../raw/annotation/gencode.v29.annotation.gtf.gz
INPUT_APPRIS_LIST=../raw/annotation/appris_data.principal.txt.gz

mkdir -p ../riboflow_annot_and_ref
OUTPUT_BASE=../riboflow_annot_and_ref/appris_human_v1
SELECTED_FASTA="${OUTPUT_BASE}_selected.fa.gz"


python generate_appris_reference.py -f ${INPUT_FASTA} \
                                    -g ${INPUT_GTF} \
                                    -a ${INPUT_APPRIS_LIST} \
                                    -o ${OUTPUT_BASE}

echo "#################################"
echo "Generating Bowtie2 Reference"
echo "#################################"

bowtie2-build ${SELECTED_FASTA} ${OUTPUT_BASE}_selected
