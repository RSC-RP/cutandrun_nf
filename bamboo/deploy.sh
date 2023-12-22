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
# positional arguments defined on the bamboo deploy stage 
DEPLOY=$1
OUTNAME=$2

echo "test_server = $TEST_SERVER"
echo "web_server = $WEB_SERVER"
echo "log_root = $LOG_ROOT"
echo "build_server = $BUILD_SERVER"
echo "svc_user = $SVC_USER"
echo "svc_pass = $SVC_PASS"
echo "results directory outname = $OUTNAME"

echo "define artifact dir"
ART_DIR=./artifacts
echo $ART_DIR
DEPLOY_DIR_LOG=$LOG_ROOT/RPDEV/cutandrun_nf/deploy
echo $DEPLOY_DIR_LOG

echo "clean out deploy_dir_log"
if [[ -d $DEPLOY_DIR_LOG ]]
then
    rm -rf $DEPLOY_DIR_LOG/*
fi

echo "define the deployment directory"
#/gpfs/assoc/rsc/nextflow_outs/
OUTDIR=$DEPLOY_DIR_LOG/$OUTNAME
echo "clean out deployment directory"
if [[ -d $OUTDIR ]]
then
    rm -rf $OUTDIR/*
fi

echo "deploy artifacts logs"
cp -R $ART_DIR/* $DEPLOY_DIR_LOG/ || { echo "artifacts not found"; exit 1; }
chmod -R 775 $DEPLOY_DIR_LOG/*

echo "Passed argument: $DEPLOY"
# echo "if something needs to fail, have it exit with non-zero error code and Bamboo will fail the task"
if [ $DEPLOY == 'main' ] || [ $DEPLOY == 'dev' ]; then
  echo "create output dir on deploy machine"
  mkdir -p $OUTDIR
#   TEMP_DIR=$(sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "mktemp -d -p ~/bamboo_tmp")
  echo "created $OUTDIR"

  echo "copy artifacts to deploy machine"
  cp -R $ART_DIR/* $OUTDIR/ || { echo "artifacts not found"; exit 1; }
  chmod -R 775 $OUTDIR/*
#   sshpass -f $SVC_PASS scp -r $ART_DIR/* $SVC_USER@$BUILD_SERVER:$OUTDIR
#   sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "ls -Ra $OUTDIR"

  if [ $DEPLOY == 'main' ]; then
    echo ""
  fi
  
#   echo "clean up deploy machine"
#   sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "rm -rf $TEMP_DIR"
fi

echo "FINISHED DEPLOYMENT"
exit
