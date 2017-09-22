#!/bin/bash
#echo -e "------------------- START --------------------"

cmsenv_dir=$1
cwd=$2
to_run=`echo ${@:3}`

echo -e "---------------- Environments ----------------"

echo -e "\n[0] cd $cmsenv_dir/src"
cd $cmsenv_dir/src

echo -e "\n[1] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo -e "\n[2] cmsenv"
eval `scramv1 runtime -sh`

echo -e "\n[3] cd $cwd"
cd $cwd

echo -e "\n------------------ Combine ------------------"

i=4
# Multiple commands can be run, arguments should be separated by commas (",")
# and white spaces replaced by 3 underscores ("___")
if (( `echo $to_run | grep "," | wc -l`)); then
    for job in `echo $to_run | sed "s;,; ;g;"`; do
	job=`echo $job | sed "s;___; ;g"`
	echo -e "\n[$i] time $job"
	time $job
	i=$((i+1))
    done
else
    echo -e "\n[$i] time $to_run"
    time $to_run
fi

echo -e "\n"
echo -e "-------------------- END ---------------------\n"

