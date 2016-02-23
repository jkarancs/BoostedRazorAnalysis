set LATEST_NTUPLE_EOS="Dec06"
set LATEST_NTUPLE_GRID18="Jan12_edm_Jan06"
#set LATEST_NTUPLE_GRID18="Skim_Jan26_1AK8JetPt350"

cmsenv
cd $CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer

if ( `echo $HOSTNAME | grep lxplus | wc -l` ) then
    # Mount EOS if not mounted already
    if ( ! -d eos_mount_dir ) then
	mkdir -p eos_mount_dir
	chmod 444 eos_mount_dir
    endif
    if ( ! `ls eos_mount_dir | wc -l` ) then
	eosmount eos_mount_dir
    endif
    # Remake softlink to point to latest ntuple location - on eos
    if ( -l ntuple/Latest ) rm ntuple/Latest
    ln -s eos/$LATEST_NTUPLE_EOS ntuple/Latest
else if ( `echo $HOSTNAME | grep "grid18\.kfki\.hu" | wc -l` ) then
    # Remake softlink to point to latest ntuple location - on grid18
    if ( -l ntuple/Latest ) rm ntuple/Latest
    ln -s grid18/$LATEST_NTUPLE_GRID18 ntuple/Latest
else
    echo "Error, not on lxplus or grid18 (Budapest)"
    exit
endif

# Remake all file lists based on the latest available files
if ( `ls filelists/* | grep "\.txt" | wc -l` ) rm filelists/*/*.txt
source scripts/make_file_lists.csh

cd -
