#!/bin/bash

#datasets="/data/jkarancs/CMSSW/ntuple/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_Oct24_input.txt"
datasets="/data/jkarancs/CMSSW/ntuple/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_Jan12_input.txt"
xsecfile="/data/jkarancs/CMSSW/ntuple/Analysis/B2GTTrees/test/crab3/cross_sections.txt"
ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Jan12"

echo -n "" > BackGroundXSec.txt
while read short dataset; do
    if (( `echo $dataset | grep MINIAODSIM | grep -v "SMS" | wc -l` )); then
	primary_dataset=`echo $dataset | sed "s;/; ;g" | awk '{ print $1 }'`
	xsec=`grep '^'$primary_dataset $xsecfile | tail -1 | awk '{ printf "%f", $2*$3 }'`
	totweight=`root -l 'scripts/get_totweight.C("'$ntuple/$short/*.root'")' | grep "totweight:" | awk '{ print $2 }'`
	echo $short $xsec $totweight
	echo $short $primary_dataset $xsec $totweight | awk '{ printf "%-50s %-85s %f %f\n", $1, $2, $3, $4  }' >> BackGroundXSec.txt
    fi
done < $datasets



ls -ltr BackGroundXSec.txt
