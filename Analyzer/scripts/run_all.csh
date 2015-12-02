set IntLumi=2500 # pb^-1

set DATE=`date | sed "s; ;_;g;s;:;h;1;s;:;m;1"`

mkdir -p results/log

# Print commands for backgrounds
foreach bkg ( `ls ../data/filelists/backgrounds | grep "\.txt" | sed "s;\.txt;;"` )
    echo "./Analyzer ../data/filelists/backgrounds/$bkg.txt 0 0 $IntLumi results/Bkg_$bkg.root Bkg_$bkg Pileup_True >! results/log/Bkg_$bkg"_"$DATE.log"
end

# Print commands for data
foreach data ( `ls ../data/filelists/data | grep "\.txt" | sed "s;\.txt;;"` )
    echo "./Analyzer ../data/filelists/data/$data.txt -1 -1 -1 results/Data_$data.root Data_$data >! results/log/Data_$data"_"$DATE.log"
end
