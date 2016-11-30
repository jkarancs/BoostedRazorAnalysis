#!/bin/bash
cwd=$1
out=$3
out_tmp=$PWD/`basename $3`
to_run=`echo ${@:2} | sed "s;$out;$out_tmp;1"`

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

echo -e "\n[4] cp -p $out_tmp $out"
mv $out_tmp $cwd/$out

echo -e "\n"
echo -e "-------------------- END ---------------------\n"


