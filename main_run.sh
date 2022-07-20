#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow.config
#Options: 'local_singularity', 'PBS_singularity', and 'local_docker'
NFX_PROFILE='PBS_singularity'
#Options: 'bowtie2_index', 'call_peaks'
NFX_ENTRY='call_peaks'
#The output prefix on filenames for reports/logs
REPORT=${1:-"pipeline_report"}

# Load the modules 
module load singularity/3.9.9

# Nextflow run to execute the workflow 
PREFIX=${REPORT}_${DATE}
nextflow -c ${NFX_CONFIG} -log reports/${PREFIX}_nextflow.log run main.nf \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report reports/${PREFIX}.html \
    -with-dag dag/${PREFIX}_dag.pdf \
    -cache TRUE \
    -resume
