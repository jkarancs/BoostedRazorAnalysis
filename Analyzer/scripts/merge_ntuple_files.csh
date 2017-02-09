set dir=$1

# Merge files - one file for each dataset
# check if directory is not currently written to (min 5 minutes since modified)
if ( -e files.txt ) rm files.txt
foreach sample ( `find $dir/* -maxdepth 0 -type d -mmin +5` )
    # and it has files in it (eg. did not already merge files)
    if ( `ls $sample | grep "\.root" | wc -l` > 0 ) then
        echo "hadd $sample.root $sample/*.root"
        echo "rm $sample/*.root"
	echo $sample.root >>! files.txt
    endif
end

# Merge extension datasets to a single file
foreach ext ( `cat files.txt | grep "_ext\.root"` )
    set base=`echo $ext | sed "s;_ext;;"`
    set merged=`echo $ext | sed "s;_ext;_merged;"`
    echo "hadd $merged $base $ext"
    echo "rm $base $ext"
end

# Make file lists for data, backgrounds and signals
# Since the files are already merged, these are mostly single files
# (except for data)
foreach data ( `ls $dir | grep "\.root" | grep ns_2015 | sed "s;_2015.*;;" | sort -u` )
    ls -tr $dir/"$data"_*.root >! ../data/filelists/data/$data.txt
end

foreach bkg ( `ls $dir | grep "\.root" | grep -v ns_2015 | grep -v T1tttt | sed "s;\.root;;"` )
    echo $dir/$bkg.root >! ../data/filelists/backgrounds/$bkg.txt
end

ls -tr $dir/T1tttt_mGluino-*.root >! ../data/filelists/signals/T1tttt_FastSim.txt
echo $dir/T1tttt_FullSim_mGluino-1500_mLSP-100.root >! ../data/filelists/signals/T1tttt_FullSim.txt
