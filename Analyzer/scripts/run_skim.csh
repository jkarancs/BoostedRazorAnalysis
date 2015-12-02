set DATE=`date | sed "s; ;_;g;s;:;h;1;s;:;m;1"`

mkdir -p results/skim/log

# Print commands for backgrounds
foreach skim ( `ls results | grep "\.root" | sed "s;\.root;;" ` )
    set logfile=`ls -tr results/log/$skim*.log | tail -1`
    set xsect=`grep '^xsect' $logfile | awk '{ print $NF }'`
    set totweight=`grep '^totweight' $logfile | awk '{ print $NF }'`
    set lumi=`grep '^lumi' $logfile | awk '{ print $NF }'`
    echo "./Analyzer results/$skim.root $xsect $totweight $lumi results/skim/$skim.root $skim Pileup_True >! results/skim/log/$skim"_"$DATE.log"
end
