set DATE=`date | cut -f2- -d" " | sed "s; ;_;g;s;:;h;1;s;:;m;1"`
set OUTDIR="results/run_$DATE"
set SKIMOUTDIR="ntuple/grid18/Skim_Feb22_1AK8JetPt300"
set SCRIPT="run_$DATE.csh"

set args=""
foreach n ( `seq 1 $#argv` )
    set args="$args "$argv[$n]
end


# options
set opt_run=0
set opt_skim=0
set opt_test=0
set opt_replot=0
set no_opt=0
set N=0
if ( `echo "$args" | grep "\-\-run" | wc -l` ) then
    set opt_run=1
    alias eval_or_echo   'echo \!* >! temp_$DATE.csh; source temp_$DATE.csh; rm temp_$DATE.csh'

else if ( `echo "$args" | grep "\-\-skim" | wc -l` ) then
    set opt_skim=1
    set SCRIPT="skim_$DATE.csh"
    set OUTDIR="results/skim_$DATE"
    alias eval_or_echo   'echo \!* >! temp_$DATE.csh; source temp_$DATE.csh; rm temp_$DATE.csh'

else if ( `echo "$args" | grep "\-\-test" | wc -l` ) then
    set opt_test=1
    alias eval_or_echo 'echo "echo \\[$N\\] "\!*"\necho" >>! temp_$DATE.csh; set N=`expr $N + 1`; echo \!* >>! temp_$DATE.csh'

else if ( `echo "$args" | grep "\-\-replot" | wc -l` ) then
    set opt_replot=1
    alias eval_or_echo 'echo "echo \\[$N\\] "\!*"\necho" >>! temp_$DATE.csh; set N=`expr $N + 1`; echo \!* >>! temp_$DATE.csh'

else
    set no_opt=1
    alias eval_or_echo   'echo "[$N] "\!*; echo; set N=`expr $N + 1`; echo \!* >>! temp_$DATE.csh'
    echo "The script will run the following if --run added:\n"

endif

if ( $opt_run || $opt_skim || $no_opt ) then
    eval_or_echo "make clean"
    eval_or_echo "make Analyzer ||  echo '\\nCompilation failed. Exiting...' && setenv makefailed && exit"
    eval_or_echo "make Plotter ||  echo '\\nCompilation failed. Exiting...' && setenv makefailed && exit"
    eval_or_echo "mkdir -p $OUTDIR/log"
    eval_or_echo "mkdir -p $OUTDIR/backup/common"
    eval_or_echo "cp -rp *.h *.cc Makefile* $OUTDIR/backup/"
    eval_or_echo "cp -rp common/*.h common/*.cc $OUTDIR/backup/common/"
    
    set all_output=""
    # Print commands for backgrounds
    echo -n "" >! $SCRIPT
    foreach bkg ( `ls filelists/backgrounds | grep "\.txt" | sed "s;\.txt;;"` )
	if ( $opt_skim) then
	    eval_or_echo "mkdir -p $SKIMOUTDIR/$bkg"
	    echo "./Analyzer $SKIMOUTDIR/$bkg/Skim.root filelists/backgrounds/$bkg.txt >! $OUTDIR/log/Bkg_$bkg"_"$DATE.log" >> $SCRIPT
	    set all_output="$all_output $SKIMOUTDIR/$bkg/Skim.root"
	else
	    echo "./Analyzer $OUTDIR/Bkg_$bkg.root filelists/backgrounds/$bkg.txt >! $OUTDIR/log/Bkg_$bkg"_"$DATE.log" >> $SCRIPT
	    set all_output="$all_output $OUTDIR/Bkg_$bkg.root"
	endif
    end
    
    # Print commands for data
    foreach data ( `ls filelists/data | grep "\.txt" | sed "s;\.txt;;"` )
	if ( $opt_skim) then
            eval_or_echo "mkdir -p $SKIMOUTDIR/$data"
            echo "./Analyzer $SKIMOUTDIR/$data/Skim.root filelists/data/$data.txt >! $OUTDIR/log/Data_$data"_"$DATE.log" >> $SCRIPT
            set all_output="$all_output $SKIMOUTDIR/$data/Skim.root"
	else
	    echo "./Analyzer $OUTDIR/Data_$data.root filelists/data/$data.txt >! $OUTDIR/log/Data_$data"_"$DATE.log" >> $SCRIPT
	    set all_output="$all_output $OUTDIR/Data_$data.root"
	endif
    end
    
    # Print commands for signal
    foreach signal ( `ls filelists/signals | grep "\.txt" | sed "s;\.txt;;"` )
	if ( $opt_skim ) then
            eval_or_echo "mkdir -p $SKIMOUTDIR/$signal"
            echo "./Analyzer $SKIMOUTDIR/$signal/Skim.root filelists/signals/$signal.txt >! $OUTDIR/log/Signal_$signal"_"$DATE.log" >> $SCRIPT
            set all_output="$all_output $SKIMOUTDIR/$signal/Skim.root"
	else
	    #echo "./Analyzer $OUTDIR/Signal_$signal.root filelists/signals/$signal.txt >! $OUTDIR/log/Signal_$signal"_"$DATE.log" >> $SCRIPT
	    #set all_output="$all_output $OUTDIR/Signal_$signal.root"
	endif
    end
    
    # Run scripts (Analyzers in parallel, then Plotter)
    if ( $opt_skim ) then
	eval_or_echo "source $SCRIPT"
    else
	eval_or_echo "par_source $SCRIPT 4"
    endif
    eval_or_echo "sleep 1s"
    eval_or_echo "mv $SCRIPT $OUTDIR/backup/"
    eval_or_echo "./Plotter results/Plotter_out_$DATE.root $all_output"
    eval_or_echo "if ( -e Latest_plots.root ) rm Latest_plots.root"
    eval_or_echo "ln -s results/Plotter_out_$DATE.root Latest_plots.root"

else if ( $opt_replot ) then
    set infile=`ls -tr results/Plotter_out_* | grep -v replot | tail -1`
    set outfile=`echo $infile | sed "s;\.root;_replot\.root;"`
    set orig_outdir=`echo $infile | sed "s;Plotter_out_;;;s;\.root;;"`
    eval_or_echo "make clean"
    eval_or_echo "make Plotter ||  echo '\\nCompilation failed. Exiting...' && setenv makefailed && exit"
    eval_or_echo "mkdir -p $orig_outdir/backup_replot/common"
    eval_or_echo "cp -rp *.h *.cc Makefile* $orig_outdir/backup_replot/"
    eval_or_echo "cp -rp common/*.h common/*.cc $orig_outdir/backup_replot/common/"
    eval_or_echo "./Plotter $outfile $infile"
    eval_or_echo "if ( -e Latest_plots.root ) rm Latest_plots.root"
    eval_or_echo "ln -s $outfile Latest_plots.root"
    echo "\nRemaking the latest plots...\n"
    source temp_$DATE.csh
    if ( ! $?makefailed ) then
        echo "\nDone.\nShowing the result..."
	root 'scripts/show_result.C("'Latest_plots.root'")'
    else
	unsetenv makefailed
    endif
    rm temp_$DATE.csh

else if ( $opt_test ) then
    eval_or_echo "make clean"
    eval_or_echo "make Analyzer ||  echo '\\nCompilation failed. Exiting...' && setenv makefailed && exit"
    eval_or_echo "make Plotter ||  echo '\\nCompilation failed. Exiting...' && setenv makefailed && exit"
    eval_or_echo "./Analyzer quickTest=1 results/Analyzer_test_Data.root filelists/data/JetHT_25ns.txt"
    eval_or_echo "./Analyzer quickTest=1 results/Analyzer_test_QCD.root filelists/backgrounds/QCD_HT500to700.txt"
    eval_or_echo "./Analyzer quickTest=1 results/Analyzer_test_TTbar.root filelists/backgrounds/TTJets_HT-600to800.txt"
    eval_or_echo "./Plotter results/Plotter_test.root results/Analyzer_test_Data.root results/Analyzer_test_QCD.root results/Analyzer_test_TTbar.root"
    echo "\nRunning the following test...\n"
    source temp_$DATE.csh
    if ( ! $?makefailed ) then
        echo "\nDone.\nShowing the result..."
        root 'scripts/show_result.C("'results/Plotter_test.root'")'
    else
	unsetenv makefailed
    endif
    rm temp_$DATE.csh

endif

if ( $no_opt ) then
    echo "------------------------------------\n"
    echo "Possible options:"
    echo "  --run       runs above lines"
    echo "  --skim      creates skimmed ntuples in $SKIMOUTDIR"
    echo "  --replot    remake plots from the output of last run"
    echo "  --test      runs a quick test instead\n"
    echo "None of the above specified, do you want to run?"
    echo -n "Please type your answer (yes/no): "
    set answer=$<
    echo
    if ( $answer == "yes" ) then
        echo "\nOk. Running...\n"
        source temp_$DATE.csh
        if ( ! $?makefailed ) then
            echo "\nDone.\nShowing the result..."
            root 'scripts/show_result.C("'Latest_plots.root'")'
        else
            unsetenv makefailed
        endif
    else
	echo "OK. Exiting...\n"
	rm $SCRIPT
    endif
    rm temp_$DATE.csh
endif
