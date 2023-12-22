#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=../nextflow_test.config
#Options: 'PBS_apptainer','local_apptainer','local_singularity', 'PBS_singularity'
NFX_PROFILE='PBS_apptainer'
#Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
NFX_ENTRY='align_call_peaks'
#The output prefix on filenames for reports/logs
REPORT=${1:-"pipeline_report"}

# Nextflow run to execute the workflow 
PREFIX=${REPORT}_${DATE}
nextflow -c ${NFX_CONFIG}\
    -log artifacts/reports/${PREFIX}_nextflow.log \
    run main.nf \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report artifacts/reports/${PREFIX}.html \
    -with-dag artifacts/dag/${PREFIX}_dag.pdf \
    -cache TRUE \
    -resume \
    "$@"
