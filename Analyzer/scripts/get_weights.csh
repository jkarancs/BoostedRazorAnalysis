if ( -e weights.txt ) rm weights.txt
foreach list ( `ls ../data/filelists/backgrounds/*.txt` )
    foreach file ( `cat $list` )
	root 'get_weight.C("'$file'")' | grep "SumWeights" | tee -a weights.txt
    end
end

echo "\n\nFiles with wrong weights:\n\n"

set N=`wc -l weights.txt | awk '{ print $1 }'`
foreach i ( `seq 1 $N` )
    set line=`sed -e $i'p;d' weights.txt`
    set weight=`echo $line | awk '{ print $NF }'`
    set weight_file=`echo $line | awk '{ print $(NF-2) }'`
    if ( `echo $weight $weight_file | awk '{ printf "%f\n", $1-$2  }'` != "0.000000" ) echo $line
end
