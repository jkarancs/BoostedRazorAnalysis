set dir="$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/ntuple/Latest"

mkdir -p filelists/backgrounds
mkdir -p filelists/signals
mkdir -p filelists/data

# Make file lists for Data, Backgrounds and Signals
# Since the files are already merged, these are mostly single files
# (except for data)
foreach data_era ( `ls $dir | grep -Ev "\.root|log" | grep "_25ns" | sed "s;_2015.*;;" | sort -u` )
    ls -tr $dir/"$data_era"*/*.root >! filelists/data/$data_era.txt
end


foreach bkg ( `ls $dir | grep -Ev '\.root|log|_ext$|_25ns|T1tttt' ` )
    ls -tr $dir/$bkg/*.root >! filelists/backgrounds/$bkg.txt
    # Add also extension datasets if available ( "_ext" at the end of same directory )
    if ( -d $dir/"$bkg"_ext ) ls -tr $dir/"$bkg"_ext/*.root >> filelists/backgrounds/$bkg.txt
end

ls -tr $dir/T1tttt_*/*.root | grep FullSim >! filelists/signals/T1tttt_FullSim.txt
ls -tr $dir/T1tttt_*/*.root | grep -v FullSim >! filelists/signals/T1tttt_FastSim.txt
