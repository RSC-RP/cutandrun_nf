#!/bin/bash

set -eou pipefail 

DATE=$(date +%F)

echo "docker build with tag: nxf:$DATE"
docker build --tag nxf:$DATE -f Docker/Dockerfile .

echo "apptainer build from local cache"
apptainer build nxf_${DATE}.sif docker-daemon:nxf:$DATE

# echo "remove docker and apptainer images"
# docker image rm nxf:$DATE
# apptainer delete nxf:$DATE
# rm nxf_${DATE}.sif