#!/bin/bash
cwd=$1
to_run=${@:2}



echo -e "------------------- START --------------------"
echo -e "---------------- Environments ----------------"

echo -e "\n[0] cd $cwd"
cd $cwd

echo -e "\n[1] source /afs/cern.ch/cms/cmsset_default.sh"
source /afs/cern.ch/cms/cmsset_default.sh

echo -e "\n[2] cmsenv"
eval `scramv1 runtime -sh`

echo -e "\n------------------ Analyzer ------------------"

echo -e "\n[3] $to_run"
$to_run

echo -e "\n"
echo -e "-------------------- END ---------------------\n"


