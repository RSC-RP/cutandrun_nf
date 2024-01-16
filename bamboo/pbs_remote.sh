#!/bin/bash
module load pbspro/2020.1

# query loop logic from: https://stackoverflow.com/questions/25674240/how-to-know-when-pbs-batch-jobs-are-complete
numJobs=0
# jobid=$(~/tmp/pbs_wrapper.sh ~/tmp/knitRForBiols.pbs)

if [ ! -f "$1" ]; then
    echo "PBS script file does not exist: $1"
    exit 1
fi

if [ ! -d "$2" ]; then
    echo "TEMP_DIR does not exist: $2"
    exit 1
fi

jobid=$(qsub -v TEMP_DIR="$2" $1)
echo "job: $jobid"
echo "PBS Script: $1"
echo "TEMP_DIR: $2"
(( numJobs++ ))
numDone=0
while [ $numDone -lt $numJobs ]; do
  if qstat -xf $jobid | grep 'Exit_status = '
  then
    (( numDone++ ))
    echo "pbs job returned: $jobid"
    echo "PBS Script: $1"
    echo "TEMP_DIR: $2"
  else
    echo "sleeping while waiting on job to finish: $jobid"
    sleep 10
  fi
done
