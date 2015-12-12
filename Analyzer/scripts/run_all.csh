set DATE=`date | sed "s; ;_;g;s;:;h;1;s;:;m;1"`

mkdir -p results/log

# Print commands for backgrounds
foreach bkg ( `ls filelists/backgrounds | grep "\.txt" | sed "s;\.txt;;"` )
    echo "./Analyzer results/Bkg_$bkg.root filelists/backgrounds/$bkg.txt >! results/log/Bkg_$bkg"_"$DATE.log"
end

# Print commands for data
foreach data ( `ls filelists/data | grep "\.txt" | sed "s;\.txt;;"` )
    echo "./Analyzer results/Data_$data.root filelists/data/$data.txt >! results/log/Data_$data"_"$DATE.log"
end
