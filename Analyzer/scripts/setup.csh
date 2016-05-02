# 74X
set LATEST_NTUPLE_EOS="Skim_Feb22_1AK8JetPt300"
#set LATEST_NTUPLE_GRID18="Skim_Feb22_1AK8JetPt300"

# 76X
set LATEST_NTUPLE_GRID18="Apr13_edm_Apr01"

cmsenv

set ANA_BASE=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer
set DIR="$ANA_BASE/ntuple/Latest"

if ( `echo $HOSTNAME | grep lxplus | wc -l` ) then
    # Mount EOS if not mounted already
    if ( ! -d $ANA_BASE/eos_mount_dir ) then
	mkdir -p $ANA_BASE/eos_mount_dir
	chmod 444 $ANA_BASE/eos_mount_dir
    endif
    if ( ! `ls $ANA_BASE/eos_mount_dir | wc -l` ) then
	eosmount $ANA_BASE/eos_mount_dir
    endif
    # Remake softlink to point to latest ntuple location - on eos
    if ( -l $ANA_BASE/ntuple/Latest ) rm $ANA_BASE/ntuple/Latest
    ln -s eos/$LATEST_NTUPLE_EOS $ANA_BASE/ntuple/Latest
else if ( `echo $HOSTNAME | grep "grid18\.kfki\.hu" | wc -l` ) then
    # Remake softlink to point to latest ntuple location - on grid18
    if ( -l $ANA_BASE/ntuple/Latest ) rm $ANA_BASE/ntuple/Latest
    ln -s grid18/$LATEST_NTUPLE_GRID18 $ANA_BASE/ntuple/Latest
else
    echo "Error, not on lxplus or grid18 (Budapest)"
    exit
endif

# Remake all file lists based on the latest available files
mkdir -p $ANA_BASE/filelists/backgrounds $ANA_BASE/filelists/signals $ANA_BASE/filelists/data
if ( `ls $ANA_BASE/filelists/* | grep "\.txt" | wc -l` ) rm $ANA_BASE/filelists/*/*.txt

# Data: Merge anything that has eras in their name
# Signals: Simplified models, like T1tttt, T2tt, T5ttcc etc.
# Background: All the rest
ls -l $DIR/ | grep "^d" | awk '{ print $9 }' > dirs_temp.txt
ls $DIR | grep -Ev "\.root|log" >! dirs_temp.txt
foreach data_era ( `grep "20[1-2][0-9][A-G]" dirs_temp.txt` )
    ls -tr $DIR/$data_era/*.root >! filelists/data/$data_era.txt
end
foreach sig ( `grep "T[1-9][t,b,c,q][t,b,c,q]" dirs_temp.txt` )
    ls -tr $DIR/$sig/*.root >! filelists/signals/$sig.txt
end
foreach bkg ( `grep -Ev '20[1-2][0-9][A-G]|T[1-9][t,b,c,q][t,b,c,q]|_ext$' dirs_temp.txt` )
    ls -tr $DIR/$bkg/*.root >! filelists/backgrounds/$bkg.txt
    # Add also extension datasets if available ( "_ext" at the end of same directory )
    if ( -d $DIR/"$bkg"_ext ) ls -tr $DIR/"$bkg"_ext/*.root >> filelists/backgrounds/$bkg.txt
end
rm dirs_temp.txt
