#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow_test.config
#Options: 'PBS_apptainer','local_apptainer','local_singularity', 'PBS_singularity'
NFX_PROFILE='local_apptainer'
#Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
NFX_ENTRY='align_call_peaks'

# REPORT is the output prefix on filenames for reports/logs
# if 1st positional argument is a command line flag for nextflow run, then use the default report name
ALL_ARGS="$@"
LEN="$#"
if [ $(echo "$ALL_ARGS" | grep -Ec "^-{1,2}[a-zA-Z]") == 0 ] && [ $LEN -le 1 ]
then
    REPORT=${1:-"pipeline_report"}
    ARGS=''
elif [ $(echo "$ALL_ARGS" | grep -Ec "^-{1,2}[a-zA-Z]") != 0 ]
then
    REPORT="pipeline_report"
    ARGS=$ALL_ARGS
else
    REPORT=${@:1:1}
    ARGS=${@:2:$LEN}
fi

# Nextflow run to execute the workflow 
PREFIX=${REPORT}_${DATE}
nextflow -c ${NFX_CONFIG}\
    -log artifacts/reports/${PREFIX}_nextflow.log \
    run main.nf \
    $ARGS \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report artifacts/reports/${PREFIX}.html \
    -with-dag artifacts/dag/${PREFIX}_dag.pdf \
    -cache FALSE
