#!/bin/bash -l

set -eou pipefail

# bamboo automation assumes these environment variables are passed in from build_pipeline.sh called by build.sh
export WORK_DIR=$WORK_DIR
export IMAGE_CACHE=$IMAGE_CACHE
cd $TEMP_DIR

# Set-up nexflow conda image
DATE=$(date +%F)
IMAGE=nxf_${DATE}.sif
CONDA_ENV_NAME="nxf_temp"

# nextflow is installed in a non-standard directory in the container
# see https://apptainer.org/docs/user/main/environment_and_metadata.html#manipulating-path
# export APPTAINERENV_APPEND_PATH="/opt/conda/bin:/depot/apps/apptainer/1.1.9/bin:/opt/pbspro/2020.1/bin"
# export APPTAINERENV_APPEND_PATH="/opt/conda/bin"
# appt_version=$(ml apptainer; apptainer exec $IMAGE /bin/bash -c "nextflow -version")
# echo $appt_version
export SINGULARITYENV_APPEND_PATH="/opt/conda/bin"
sing_version=$(ml singularity; singularity exec $IMAGE /bin/bash -c "nextflow -version")
echo $sing_version

echo "create new mamba environment"
# print out the software version information. Assumes svc account has mamba installed in default location.
VER=$(mamba --version)
echo "mamba version:" $VER
mamba env create --quiet --force -f env/nextflow.yaml --name "$CONDA_ENV_NAME"

echo "activate nextflow mamba environment"
conda activate $CONDA_ENV_NAME
mamba env list

echo "create artifact dir"
mkdir -p artifacts
OUTDIR='./artifacts/results/mouse'

echo "create nextflow work directory"
WORK_DIR=$WORK_DIR/$(basename $TEMP_DIR)
mkdir -p $WORK_DIR
echo $WORK_DIR

echo "create artifacts"
# run the pipeline using the default parameters
REPORT1="default"
bash ./bamboo/bamboo_main_run.sh $REPORT1 --outdir "$OUTDIR/$REPORT1" -work-dir "$WORK_DIR"

# run the pipeline using the default parameters and building the bowtie indexes
REPORT2="build_index"
# bash ./bamboo/bamboo_main_run.sh $REPORT2 --outdir "$OUTDIR/$REPORT2" --build_index 'true' --build_spike_index 'true'

# run with different executors 

# run different entry points

# Deactivate and delete the environment 
echo "Deactivate and delete the environment"
conda deactivate
mamba env remove --name $CONDA_ENV_NAME --yes