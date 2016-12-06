#!/bin/bash
cwd=$1
out=$3
out_tmp=$PWD/`basename $3`
to_run=`echo ${@:2} | sed "s;$out;$out_tmp;1"`
USE_EOS_MOUNT=0

echo -e "------------------- START --------------------"
echo -e "---------------- Environments ----------------"

echo -e "\n[0] cd $cwd"
cd $cwd

echo -e "\n[1] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo -e "\n[2] cmsenv"
eval `scramv1 runtime -sh`

if (( $USE_EOS_MOUNT )); then
    echo -e "\n[3] eosmount eos_mount_dir"
    /afs/cern.ch/project/eos/installation/cms/bin/eos.select -b fuse mount eos_mount_dir
fi

echo -e "\n------------------ Analyzer ------------------"

if (( $USE_EOS_MOUNT )); then
    echo -e "\n[4] time $to_run && mv $out_tmp $cwd/$out"
else
    echo -e "\n[3] time $to_run && mv $out_tmp $cwd/$out"
fi
time $to_run && mv $out_tmp $cwd/$out

echo -e "\n"
echo -e "-------------------- END ---------------------\n"


