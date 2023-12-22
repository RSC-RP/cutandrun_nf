#!/bin/bash

set -eou pipefail 

DATE=$(date +%F)

echo "docker build with tag: nxf:$DATE"
docker build --tag nxf:$DATE -f Docker/Dockerfile .

echo "apptainer build from local cache"
# Getting warnings after the conversion - warn rootless{opt/conda/share/terminfo/a/altos-2} ignoring (usually) harmless EPERM on setxattr "user.rootlesscontainers"
# see https://github.com/EESSI/software-layer/issues/12
apptainer --silent build --force nxf_${DATE}.sif docker-daemon:nxf:$DATE

# echo "remove docker and apptainer images"
# docker image rm nxf:$DATE
# apptainer delete nxf:$DATE
# rm nxf_${DATE}.sif