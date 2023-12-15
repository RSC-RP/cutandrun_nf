#!/bin/bash

# echo "env"
# env

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

echo "create working dir on build machine"
TEMP_DIR=$(sshpass -f $SVC_PASS ssh $SVC_USER@$TEST_SERVER "mktemp -d -p /home/$SVC_USER/bamboo_tmp")
echo "created $TEMP_DIR"

echo "copy repo to test machine tmp"
sshpass -f $SVC_PASS scp -r * $SVC_USER@$TEST_SERVER:$TEMP_DIR

echo "schedule the test remotely"
sshpass -f $SVC_PASS ssh $SVC_USER@$TEST_SERVER "$TEMP_DIR/bamboo/pbs_remote.sh $TEMP_DIR/bamboo/test.pbs $TEMP_DIR" &
echo "remote job(s) scheduled"
wait # wait for pbs job(s) to finish running

echo "copy test output to bamboo machine"
sshpass -f $SVC_PASS scp -r $SVC_USER@$TEST_SERVER:$TEMP_DIR/test-reports .

echo "clean up test machine"
sshpass -f $SVC_PASS ssh $SVC_USER@$TEST_SERVER "rm -rf $TEMP_DIR"

echo "FINISHED TEST"
exit