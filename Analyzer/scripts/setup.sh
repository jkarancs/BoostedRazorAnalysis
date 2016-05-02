#!/bin/bash

# 74X
#LATEST_NTUPLE_EOS="Skim_Feb22_1AK8JetPt300"
#LATEST_NTUPLE_GRID18="Skim_Feb22_1AK8JetPt300"

# 76X
LATEST_NTUPLE_EOS="Skim_Apr28_1AK8JetPt300"
#LATEST_NTUPLE_GRID18="Apr13_edm_Apr01"
LATEST_NTUPLE_GRID18="Skim_Apr28_1AK8JetPt300"

cmsenv

ANA_BASE=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer
DIR="$ANA_BASE/ntuple/Latest"

if (( `echo $HOSTNAME | grep lxplus | wc -l` )); then
    # Mount EOS if not mounted already
    if [ ! -d $ANA_BASE/eos_mount_dir ]; then mkdir -p $ANA_BASE/eos_mount_dir; chmod 444 $ANA_BASE/eos_mount_dir; fi
    if (( ! `ls $ANA_BASE/eos_mount_dir | wc -l` )); then eosmount $ANA_BASE/eos_mount_dir; fi
    # Remake softlink to point to latest ntuple location - on eos
    if [ -L $ANA_BASE/ntuple/Latest ]; then rm $ANA_BASE/ntuple/Latest; fi
    ln -s eos/$LATEST_NTUPLE_EOS $ANA_BASE/ntuple/Latest
elif (( `echo $HOSTNAME | grep "grid18\.kfki\.hu" | wc -l` )); then
    # Remake softlink to point to latest ntuple location - on grid18
    if [ -L $ANA_BASE/ntuple/Latest ]; then rm $ANA_BASE/ntuple/Latest; fi
    ln -s grid18/$LATEST_NTUPLE_GRID18 $ANA_BASE/ntuple/Latest
else
    echo "Error, not on lxplus or grid18 (Budapest)"
    exit
fi

# (Re)make all file lists based on the latest available files
mkdir -p $ANA_BASE/filelists/backgrounds $ANA_BASE/filelists/signals $ANA_BASE/filelists/data
if (( `ls $ANA_BASE/filelists/* | grep "\.txt" | wc -l` )); then rm $ANA_BASE/filelists/*/*.txt; fi

# Data: Merge anything that has eras in their name
# Signals: Simplified models, like T1tttt, T2tt, T5ttcc etc.
# Background: All the rest
ls -l $DIR/ | grep "^d" | awk '{ print $9 }' > dirs_temp.txt
for data_era in `grep "20[1-2][0-9][A-G]" dirs_temp.txt`; do ls -tr $DIR/$data_era/*.root > $ANA_BASE/filelists/data/$data_era.txt; done
for sig in `grep "T[1-9][t,b,c,q][t,b,c,q]" dirs_temp.txt`; do ls -tr $DIR/$sig/*.root > $ANA_BASE/filelists/signals/$sig.txt; done
for bkg in `grep -Ev '20[1-2][0-9][A-G]|T[1-9][t,b,c,q][t,b,c,q]|_ext$' dirs_temp.txt`; do
    ls -tr $DIR/$bkg/*.root > $ANA_BASE/filelists/backgrounds/$bkg.txt
    # Add also extension datasets if available ( "_ext" at the end of same directory )
    if [ -d $DIR/"$bkg"_ext ]; then ls -tr $DIR/"$bkg"_ext/*.root >> $ANA_BASE/filelists/backgrounds/$bkg.txt; fi
done
rm dirs_temp.txt
