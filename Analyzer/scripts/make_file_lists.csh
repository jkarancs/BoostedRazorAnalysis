set dir="$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/ntuple/Latest"

mkdir -p filelists/backgrounds
mkdir -p filelists/signals
mkdir -p filelists/data

# Make file lists for Data, Backgrounds and Signals
ls $dir | grep -Ev "\.root|log" >! dirs_temp.txt

# Data: Merge anything that has eras in their name
foreach data_era ( `grep "_25ns" dirs_temp.txt` )
#foreach data_era ( `grep "20[1-2][0-9][A-G]" dirs_temp.txt` )
    ls -tr $dir/$data_era/*.root >! filelists/data/$data_era.txt
end

# Signals: Simplified models, like T1tttt, T2tt, T5ttcc etc.
foreach sig ( `grep "T[1-9][t,b,c,q][t,b,c,q]" dirs_temp.txt` )
    ls -tr $dir/$sig/*.root >! filelists/signals/$sig.txt
end

# Background: All the rest
foreach bkg ( `grep -Ev '_25ns|T[1-9][t,b,c,q][t,b,c,q]|_ext$' dirs_temp.txt` )
#foreach bkg ( `grep -Ev '20[1-2][0-9][A-G]|T[1-9][t,b,c,q][t,b,c,q]|_ext$' dirs_temp.txt` )
    ls -tr $dir/$bkg/*.root >! filelists/backgrounds/$bkg.txt
    # Add also extension datasets if available ( "_ext" at the end of same directory )
    if ( -d $dir/"$bkg"_ext ) ls -tr $dir/"$bkg"_ext/*.root >> filelists/backgrounds/$bkg.txt
end

rm dirs_temp.txt
