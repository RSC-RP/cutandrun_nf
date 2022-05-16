#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow.config
#Options: local_singularity, PBS_singularity, and local_docker
NFX_PROFILE='PBS_singularity'
#Options: star_index or rnaseq_count
NFX_ENTRY='rnaseq_count'
#The output prefix on filenames for reports/logs
REPORT="kappe_s_multispecies_index_quant2"

# Load the modules 
module load singularity/3.9.9

# Nextflow run to execute the workflow 
nextflow -c ${NFX_CONFIG} -log reports/${REPORT}_nextflow.log run main.nf \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report reports/${REPORT}_${DATE}.html \
    -with-dag dag/${REPORT}_${DATE}_dag.pdf \
    -cache TRUE \
    -resume
