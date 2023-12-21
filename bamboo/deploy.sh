#!/bin/bash

echo "run deploy script"
ls -R

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

echo "test_server = $TEST_SERVER"
echo "web_server = $WEB_SERVER"
echo "log_root = $LOG_ROOT"
echo "build_server = $BUILD_SERVER"
echo "svc_user = $SVC_USER"
echo "svc_pass = $SVC_PASS"

echo "define artifact dir"
ART_DIR=./artifacts
echo $ART_DIR
DEPLOY_DIR_LOG=$LOG_ROOT/RPDEV/cutandrun_nf/deploy
echo $DEPLOY_DIR_LOG

echo "clean out deploy_dir_log"
rm -rf $DEPLOY_DIR_LOG/*

echo "deploy artifacts"
cp -R $ART_DIR/* $DEPLOY_DIR_LOG/ || { echo "artifacts not found"; exit 1; }

echo "Passed argument: $1"
if [ $1 == 'main' ] || [ $1 == 'dev' ]; then
  echo "create working dir on deploy machine"
  TEMP_DIR=$(sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "mktemp -d -p ~/bamboo_tmp")
  echo "created $TEMP_DIR"

  echo "copy artifacts to deploy machine"
  sshpass -f $SVC_PASS scp -r $ART_DIR/* $SVC_USER@$BUILD_SERVER:$TEMP_DIR
  sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "ls -Ra $TEMP_DIR"

  echo "if something needs to fail, have it exit with non-zero error code and Bamboo will fail the task"

  if [ $1 == 'main' ]; then
    echo "do something specific to main branch"
  fi
  
  echo "clean up deploy machine"
  sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "rm -rf $TEMP_DIR"
fi

echo "FINISHED DEPLOYMENT"
exit
