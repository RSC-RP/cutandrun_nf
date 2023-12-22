#!/bin/bash

set -eou pipefail

# start Bamboo agent on rplbam01
# check on agent status
#   sudo systemctl status bamboo-agent
# manually start the service
#   sudo systemctl start bamboo-agent.service

# make this file executable in the repo: https://medium.com/@akash1233/change-file-permissions-when-working-with-git-repos-on-windows-ea22e34d5cee
# see current file permissions (last 3 digits of the first 6 digits) from the following command
#   git ls-files --stage
# set executable with the following command
#   git update-index --chmod=+x 'name-of-shell-script'

echo "build on HPC with user:"
echo $USER

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
DOCKER_SERVER=rdldoc01
SVC_USER=$bamboo_svc_user
SVC_PASS=$bamboo_svc_pass

echo "test_server = $TEST_SERVER"
echo "web_server = $WEB_SERVER"
echo "log_root = $LOG_ROOT"
echo "build_server = $BUILD_SERVER"
echo "docker_server = $DOCKER_SERVER"
echo "svc_user = $SVC_USER"
echo "svc_pass = $SVC_PASS"

echo "create working dir on build machine"
TEMP_DIR=$(sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "mktemp -d -p /home/$SVC_USER/bamboo_tmp")
echo "created $TEMP_DIR on $BUILD_SERVER"

# echo "copy repo to docker build machine"
# sshpass -f $SVC_PASS ssh $SVC_USER@$DOCKER_SERVER "mkdir -p $TEMP_DIR"
# sshpass -f $SVC_PASS scp -r * $SVC_USER@$DOCKER_SERVER:$TEMP_DIR
# echo "created $TEMP_DIR on $DOCKER_SERVER"

# echo "build nextflow docker image on docker build machine" 
# sshpass -f $SVC_PASS ssh $SVC_USER@$DOCKER_SERVER "cd $TEMP_DIR; ./bamboo/build_image.sh"

echo "copy repo with nextflow image to build machine tmp"
# sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "cd $(dirname $TEMP_DIR); scp -r $SVC_USER@$DOCKER_SERVER:$TEMP_DIR ."
sshpass -f $SVC_PASS scp -r * $SVC_USER@$BUILD_SERVER:$TEMP_DIR

echo "schedule the build remotely"
sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "$TEMP_DIR/bamboo/pbs_remote.sh $TEMP_DIR/bamboo/build.pbs $TEMP_DIR"
echo "remote job scheduled"
wait # wait for pbs jobs to finish running

echo "copy build output to bamboo machine"
if [[ -d "artifacts" ]]
then
    rm -rf artifacts
fi
sshpass -f $SVC_PASS scp -r $SVC_USER@$BUILD_SERVER:$TEMP_DIR/artifacts .

echo "clean up build machine"
sshpass -f $SVC_PASS ssh $SVC_USER@$BUILD_SERVER "rm -rf $TEMP_DIR"

echo "FINISHED BUILD"
exit