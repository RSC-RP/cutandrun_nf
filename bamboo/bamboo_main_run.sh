#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow.config
#Options: 'PBS_singularity', 'local_singularity'
NFX_PROFILE='PBS_singularity'
#Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
NFX_ENTRY='align_call_peaks'

# Load the modules. Set and clean cacheDirs for images.
echo "loading HPC modules for $NFX_PROFILE profile"
if [[ $NFX_PROFILE =~ "singularity" ]]
then
    module load singularity
    # print out the software version information 
    singularity --version
    echo "Setting NXF_SINGULARITY_CACHEDIR = $IMAGE_CACHE"
    export NXF_SINGULARITY_CACHEDIR=$IMAGE_CACHE
    export SINGULARITY_CACHEDIR=$IMAGE_CACHE
    echo "Clean the container image cache"
    singularity cache clean --force
elif [[ $NFX_PROFILE =~ "apptainer" ]]
then
    module load apptainer
    # print out the software version information 
    apptainer --version
    echo "Setting NXF_APPTAINER_CACHEDIR = $IMAGE_CACHE"
    export NXF_APPTAINER_CACHEDIR=$IMAGE_CACHE
    export APPTAINER_CACHEDIR=$IMAGE_CACHE
    echo "Clean the container image cache"
    apptainer cache clean --force
fi

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
nextflow -C ${NFX_CONFIG}\
    -log artifacts/reports/${PREFIX}_nextflow.log \
    run main.nf \
    $ARGS \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report artifacts/reports/${PREFIX}.html \
    -with-dag artifacts/dag/${PREFIX}_dag.pdf \
    -cache FALSE

# add permissions to reports and dags
chmod -R 644 artifacts/reports
chmod -R 644 artifacts/dag
