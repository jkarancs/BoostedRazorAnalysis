#!/bin/bash

xsecfile="/data/jkarancs/CMSSW/ntuple/Analysis/B2GTTrees/test/crab3/cross_sections.txt"

#datasets="/data/jkarancs/CMSSW/ntuple/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_Oct24_input.txt"
#datasets="/data/jkarancs/CMSSW/ntuple/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_Jan12_input.txt"
#ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Jan12"

#datasets="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_8_0_26_patch2/src/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_May10_part4_input.txt"
#ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/May10"
#datasets="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_8_0_26_patch2/src/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_May10_part5_input.txt"
#ntuple="/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/May10_part5"
#datasets="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_8_0_26_patch2/src/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_May10_part6_input.txt"
#ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/May10_part6"
#datasets="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_8_0_26_patch2/src/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_May10_part7_input.txt"
#ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/May10_part7"
#datasets="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_8_0_26_patch2/src/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_Nov30_part*_input.txt"
#ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Nov30_part1"
datasets="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_8_0_26_patch2/src/Analysis/B2GTTrees/test/crab3/MINIAODv2_80X_Nov30_part3_input.txt"
ntuple="/data_6tb/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Nov30_part3"

cat $datasets | grep -v '^\#' | grep MINIAODSIM | grep -v "SMS" > tmp_datasets.txt
echo -n "" > BackGroundXSec.txt
while read short dataset; do
    primary_dataset=`echo $dataset | sed "s;/; ;g" | awk '{ print $1 }'`
    xsec=`grep '^'$primary_dataset $xsecfile | tail -1 | awk '{ printf "%f", $2*$3 }'`
    totweight=`root -l 'scripts/get_totweight.C("'$ntuple/$short/*.root'")' | grep "totweight:" | awk '{ print $2 }'`
    echo $totweight
    # look for recovery tasks and add them if they exist
    if [ -d $ntuple/$short'_recovery' ]; then
	totweight=`echo $totweight | sed -e 's/[eE]+*/\\*10\\^/'`
	totweight2=`root -l 'scripts/get_totweight.C("'$ntuple/$short'_recovery'/*.root'")' | grep "totweight:" | awk '{ print $2 }' | sed -e 's/[eE]+*/\\*10\\^/'`
	echo $totweight2
	totweight=`echo $totweight + $totweight2 | bc`
	echo $totweight
    fi
    #echo $short $xsec $totweight
    echo $short $primary_dataset $xsec $totweight | awk '{ printf "%-50s %-85s %f %f\n", $1, $2, $3, $4  }'
    echo $short $primary_dataset $xsec $totweight | awk '{ printf "%-50s %-85s %f %f\n", $1, $2, $3, $4  }' >> BackGroundXSec.txt
done < tmp_datasets.txt
rm tmp_datasets.txt

# In a second step merge weights for extension/backup etc datasets
echo -n "" > BackGroundXSec_temp.txt
while read short primary_dataset xsec totweight; do
    if (( `echo $short | grep -E "_backup|_ext1|_ext2|_ext3|_ext4|_2" | wc -l` == 0 )); then
	totweight=`grep -E '^'$short' '\|'^'$short'_backup '\|'^'$short'_ext1 '\|'^'$short'_ext2 '\|'^'$short'_ext3 '\|'^'$short'_ext4 '\|'^'$short'_2 ' BackGroundXSec.txt | awk '{ sum+=$4 }END{ printf "%f", sum }'`
	echo $short $primary_dataset $xsec $totweight | awk '{ printf "%-50s %-85s %f %f\n", $1, $2, $3, $4  }' >> BackGroundXSec_temp.txt
    fi
done < BackGroundXSec.txt

mv BackGroundXSec_temp.txt BackGroundXSec.txt

ls -ltr BackGroundXSec.txt
