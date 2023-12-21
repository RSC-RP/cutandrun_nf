#!/bin/bash

set -eou pipefail 

DATE=$(date +%F)

echo "docker build with tag: nxf:$DATE"
# change directories up 1 level, since docker build must be executed in parent directory, not bamboo/ dir
echo "$PWD"
docker build --tag nxf:$DATE -f Docker/Dockerfile .

echo "apptainer build from local cache"
apptainer build nxf_${DATE}.sif docker-daemon:nxf:$DATE

# echo "remove docker and apptainer images"
# docker image rm nxf:$DATE
# apptainer delete nxf:$DATE
# rm nxf_${DATE}.sif