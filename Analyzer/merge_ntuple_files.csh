set dir="Nov18_ntuple"
# Merge files - one file for each dataset
# check if directory is not currently written to (min 5 minutes since modified)
foreach sample ( `find $dir/* -maxdepth 0 -type d -mmin +5` )
    # and it has files in it (eg. did not already merge files)
    if ( `ls $sample | grep "\.root" | wc -l` > 0 ) then
        echo "hadd $sample.root $sample/*.root"
        echo "rm $sample/*.root"
    endif
end
