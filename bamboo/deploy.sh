#!/bin/bash

set -eou pipefail 

echo "run deploy script"
echo "the current working directory is:" $PWD
echo "the bamboo agent is on" $HOSTNAME

echo "bamboo project variables"
echo $bamboo_test_server
echo $bamboo_web_server
echo $bamboo_log_root
echo $bamboo_build_server
echo $bamboo_svc_user
echo $bamboo_svc_pass

echo "local variable assignment"
TEST_SERVER=$bamboo_test_server
WEB_SERVER=$bamboo_web_server
LOG_ROOT=$bamboo_log_root
BUILD_SERVER=$bamboo_build_server
SVC_USER=$bamboo_svc_user
SVC_PASS=$bamboo_svc_pass
# positional arguments defined in the bamboo deploy stage 
DEPLOY=$1
OUTNAME=$2

echo "test_server = $TEST_SERVER"
echo "web_server = $WEB_SERVER"
echo "log_root = $LOG_ROOT"
echo "build_server = $BUILD_SERVER"
echo "svc_user = $SVC_USER"
echo "svc_pass = $SVC_PASS"
echo "results directory outname = $OUTNAME"

echo "define artifact dirs"
ART_DIR=./artifacts
echo $ART_DIR
DEPLOY_DIR_LOG=$LOG_ROOT/RPDEV/cutandrun_nf/deploy/$OUTNAME
DEPLOY_DIR=/gpfs/assoc/rsc/nextflow_outs/$OUTNAME
echo $DEPLOY_DIR_LOG

echo "clean out deploy_dir_log"
if [[ -d $DEPLOY_DIR_LOG ]]
then
    rm -rf $DEPLOY_DIR_LOG/*
fi

echo "deploy artifacts logs"
mkdir -p $DEPLOY_DIR_LOG
cp -R $ART_DIR/* $DEPLOY_DIR_LOG/ || { echo "artifacts not found"; exit 1; }
chmod -R 775 $DEPLOY_DIR_LOG/*

echo "Passed argument: $DEPLOY"
echo "define the deployment directory"
if [ $DEPLOY == 'main' ] 
then
    OUTDIR=$DEPLOY_DIR
elif [ $DEPLOY == 'dev' ]
then
    OUTDIR=$DEPLOY_DIR_LOG
fi

echo "clean out deployment directory"
if [[ -d $OUTDIR ]]
then
    rm -rf $OUTDIR/*
fi

echo "copy artifacts to deploy machine"
# echo "if something needs to fail, have it exit with non-zero error code and Bamboo will fail the task"
mkdir -p $OUTDIR
cp -R $ART_DIR/* $OUTDIR/ || { echo "artifacts not found"; exit 1; }
chmod -R 775 $OUTDIR/*

echo "FINISHED DEPLOYMENT"
exit
