#!/bin/bash -l

set -eou pipefail

# bamboo automation assumes workingdir is environment variable TEMP_DIR passed in from qsub -v
cd $TEMP_DIR

# Set-up nexflow conda image
DATE=$(date +%F)
IMAGE=nxf_${DATE}.sif
CONDA_ENV_NAME="nxf_temp"

# nextflow is installed in a non-standard directory in the container
# see https://apptainer.org/docs/user/main/environment_and_metadata.html#manipulating-path
# export APPTAINERENV_APPEND_PATH="/opt/conda/bin:/depot/apps/apptainer/1.1.9/bin:/opt/pbspro/2020.1/bin"
export APPTAINERENV_APPEND_PATH="/opt/conda/bin"
# apptainer exec $IMAGE /bin/bash -c "nextflow -version"
# apptainer exec $IMAGE -B /opt/pbspro/2020.1 -B /depot/apps/apptainer/1.1.9/ /bin/bash -c "nextflow -version"
# apptainer shell -B $PWD -B /opt/pbspro/2020.1 -B /depot/apps/apptainer/1.1.9/ $IMAGE

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
# export WORK_DIR=$WORK_DIR/$(basename $TEMP_DIR)

echo $WORK_DIR

echo "create artifacts"
cp bamboo/bamboo_main_run.sh ./artifacts
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
mamba env remove --name $CONDA_ENV_NAME