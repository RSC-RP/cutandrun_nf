#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow.config
#Options: 'PBS_apptainer','local_apptainer','local_singularity', 'PBS_singularity'
NFX_PROFILE='local_apptainer'
#Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
NFX_ENTRY='align_call_peaks'
#The output prefix on filenames for reports/logs
REPORT=${1:-"pipeline_report"}

# Load the modules 
if [[ $NFX_PROFILE =~ "singularity" ]]
then
    module load singularity
fi

if [[ $NFX_PROFILE =~ "apptainer" ]]
then
    module load apptainer
    apptainer --version
    apptainer cache clean --force
fi

# openssl s_client -brief -connect repo.anaconda.com:443 -servername cybertron
# openssl s_client -brief -connect quay.io:443 -servername cybertron > /dev/null
# URL="https://github.com/nextflow-io/nextflow/releases/download/v23.12.0-edge/nextflow-23.12.0-edge-all"
# NXF_EXEC=$(basename $URL)
# wget -O $NXF_EXEC $URL
# chmod +x ./$NXF_EXEC
# ./$NXF_EXEC -version
export SINGULARITY_MKSQUASHFS_PROCS=2
echo "set SINGULARITY_MKSQUASHFS_PROCS = $SINGULARITY_MKSQUASHFS_PROCS"
# Nextflow run to execute the workflow 
PREFIX=${REPORT}_${DATE}
nextflow -C ${NFX_CONFIG}\
    -log reports/${PREFIX}_nextflow.log \
    run main.nf \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report reports/${PREFIX}.html \
    -with-dag dag/${PREFIX}_dag.pdf \
    -cache TRUE \
    -resume \
    "$@"
