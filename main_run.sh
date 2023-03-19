#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow.config
#Options: 'local_singularity', 'PBS_singularity'
NFX_PROFILE='PBS_singularity'
#Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
NFX_ENTRY='align_call_peaks'
#The output prefix on filenames for reports/logs
REPORT=${1:-"pipeline_report"}

# Load the modules 
if [[ $NFX_PROFILE =~ "singularity" ]]
then
    module load singularity
fi

# Nextflow run to execute the workflow 
# See https://github.com/nextflow-io/nextflow/commit/b3a4bf857cef0a8b16814baf9e13aa9296ca208a
# export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true
PREFIX=${REPORT}_${DATE}
nextflow -c ${NFX_CONFIG} -q \
    -log reports/${PREFIX}_nextflow.log \
    run main.nf \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report reports/${PREFIX}.html \
    -with-dag dag/${PREFIX}_dag.pdf \
    -cache TRUE \
    -resume
