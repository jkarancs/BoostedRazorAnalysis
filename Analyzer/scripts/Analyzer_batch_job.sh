#!/bin/bash
cwd=$1
log=$2
out=$4
log_tmp=$PWD/`basename $log`
out_tmp=$PWD/`basename $out`
to_run=`echo ${@:3} | sed "s;$out;$out_tmp;1"`

echo -e "------------------- START --------------------"
echo -e "---------------- Environments ----------------"

echo -e "\n[0] cd $cwd"
cd $cwd

echo -e "\n[1] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo -e "\n[2] cmsenv"
eval `scramv1 runtime -sh`

echo -e "\n------------------ Analyzer ------------------"

if (( `echo $cwd | grep '^/eos' | wc -l` )); then
    echo -e "\n[3] time $to_run > $log 2>&1"
    time $to_run > $log 2>&1
    
    echo -e "\n[4] mv $log_tmp $cwd/$log"
    mv $log_tmp $cwd/$log
    
    echo -e "\n[5] mv $out_tmp $cwd/$out"
    mv $out_tmp $cwd/$out
else
    echo -e "\n[3] time $to_run "
    time $to_run > $log 2>&1
    
    echo -e "\n[4] mv $out_tmp $cwd/$out"
    mv $out_tmp $cwd/$out
fi

echo -e "\n"
echo -e "-------------------- END ---------------------\n"

