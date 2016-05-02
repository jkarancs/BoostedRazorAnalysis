import os, sys, glob, time, multiprocessing, subprocess

# Parse command line arguments
cmdargs = sys.argv
options = cmdargs[1:]

# ----------------------  Settings -----------------------
DATE = time.strftime("%Y_%m_%d_%Hh%Mm%S", time.localtime())
OUTDIR = "results/run_"+DATE # log files, backup files, output files for non-skims
SKIMNAME="Skim_Apr28_1AK8JetPt300"
EOSDIR = "srm://srm-eoscms.cern.ch/eos/cms/store/caf/user/jkarancs/B2GTTreeNtuple/"
NPROC = 4 # Number of processors to use for Analyzer jobs

if "--skim" in options:
    OUTDIR = OUTDIR.replace("run_", "skim_")

# Analyzer:
# Each element supplies 3 arguments for each Analyzer job:
# [output filename, input file list, output log]
# For skimming/full running, all datasets are used
# for testing a selected few

ana_arguments_test = [
    ["results/Analyzer_test_data.root",  "filelists/data/JetHT_25ns_2015C.txt",           "results/Analyzer_test_data.log"],
    ["results/Analyzer_test_ttbar.root", "filelists/backgrounds/TTJets_amcatnloFXFX.txt", "results/Analyzer_test_ttbar.log"],
    ["results/Analyzer_test_qcd.root",   "filelists/backgrounds/QCD_HT700to1000.txt",     "results/Analyzer_test_qcd.log"],
]

ana_arguments_full = []
ana_arguments_skim = []
for input_file in glob.glob("filelists/data/*.txt") + glob.glob("filelists/backgrounds/*.txt") + glob.glob("filelists/signals/*.txt"):
    output_file = OUTDIR+"/"+input_file.split("/")[-1].replace("txt", "root")
    skim_output_file = "ntuple/grid18/"+SKIMNAME+"/"+input_file.split("/")[-1].replace(".txt","/Skim.root")
    log_file = OUTDIR+"/log/"+input_file.split("/")[-1].replace("txt", "log")
    ana_arguments_full.append([output_file,      input_file, log_file])
    ana_arguments_skim.append([skim_output_file, input_file, log_file])

#ana_arguments_skim = [ana_arguments_skim[0], ana_arguments_skim[22], ana_arguments_skim[36]]

# Plotter:
# It gets input automatically from the output of Analyzer, output will be on screen
plotter_output_file = "results/Plotter_out_"+DATE+".root"

# ---------------------- Cmd Line  -----------------------

opt_dry = int(not "--run" in options)
opt_plot = int("--plot" in options)
opt_quick = 0
opt_test = 0
opt_transfer = 0
if (opt_dry): print "--run option not specified, doing a dry run (only printing out commands)"

if "--full" in options:
    print "Running with option: --full"
    ana_arguments = ana_arguments_full
elif "--skim" in options:
    print "Running with option: --skim"
    ana_arguments = ana_arguments_skim
    opt_transfer = 1
elif "--replot" in options:
    print "Running with option: --replot"
    ORIGOUTDIR = max(glob.glob(OUTDIR.replace(DATE,"*")), key=os.path.getmtime)
elif "--test" in options:
    opt_test = 1
else:
    opt_test = 1
    opt_quick = 1
if "--quick" in options:
    opt_quick = 1

if opt_test:
    print "Running with option: --test (few files)"
    ana_arguments = ana_arguments_test
    plotter_output_file = "results/Plotter_test.root"
if opt_quick:
    print "Running with option: --quick (1/100 statistics)"
if opt_plot:
    print "Running with option: --plot (will produce plots)"

# --------------------- Functions ------------------------
# Show and run command with stdout on screen
icommand=0
def special_call(cmd, verbose=1):
    global icommand, opt_dry
    if verbose:
        if opt_dry:
            print("[dry] "),
        else:
            print("[%d] " % icommand),
        for i in xrange(len(cmd)): print cmd[i]+" ",
        print ""
    if not opt_dry:
        if subprocess.call(cmd):
            print "Process failed, exiting ..."
            sys.exit()
        if verbose: print ""
    sys.stdout.flush()
    icommand+=1

# Run command with stdout/stderr saved to logfile
def logged_call(cmd, logfile):
    global opt_dry
    if not os.path.exists(os.path.dirname(logfile)):
        special_call(["mkdir", "-p", os.path.dirname(logfile)], 0)
    if opt_dry:
        subprocess.call(["echo", "[dry]"]+cmd+[">", logfile])
    else:
        log = open(logfile, 'a')
        proc = subprocess.Popen(cmd, stdout=log, stderr=log)
        proc.wait()
        log.close()

# Compile programs
def compile(Ana = 1, Plotter = 1):
    print "Compiling ..."
    print
    special_call(["make", "clean"])
    if Ana: special_call(["make", "Analyzer"])
    if Plotter: special_call(["make", "Plotter"])
    print "Compilation successful."
    print

# backup files for bookkeeping
def backup_files(backup_dir):
    print "Backing up files in: "+backup_dir
    print
    special_call(["mkdir", "-p", backup_dir])
    special_call(["cp", "-rp", "pileup", "systematics", "filelists", "common", "scripts"] + glob.glob("*.h") + glob.glob("*.cc") + glob.glob("Makefile*") + [backup_dir+"/"])
    print

# Run a single Analyzer instance (on a single input list, i.e. one dataset)
def analyzer_job((output_file, input_list, output_log)):
    global opt_dry, opt_quick, opt_transfer, opt_plot
    if not opt_dry:
        print "Start Analyzing: "+input_list
    if not os.path.exists(os.path.dirname(output_file)):
        special_call(["mkdir", "-p", os.path.dirname(output_file)], 0)
    cmd = ["./Analyzer", output_file, input_list]
    if not opt_plot: cmd.append("noPlots=1")
    if opt_quick: cmd.append("quickTest=1")
    logged_call(cmd, output_log)
    if opt_transfer:
        outpath = output_file.split("/")[-3]+"/"+output_file.split("/")[-2]+"/"+output_file.split("/")[-1]
        logged_call(["lcg-cp", "-v", output_file, EOSDIR+outpath], output_log)
    return output_file

# Run all Analyzer jobs in parallel
def analysis(ana_arguments, nproc):
    njob = len(ana_arguments)
    if njob<nproc: nproc = njob
    print "Running "+str(njob)+" instances of Analyzer jobs:"
    print
    workers = multiprocessing.Pool(processes=nproc)
    output_files = workers.map(analyzer_job, ana_arguments, chunksize=1)
    print "All Analyzer jobs finished."
    print
    return output_files

# Run Plotter, output of Analyzer is input for this code
def plotter(input_files, output_file):
    print "Start plotting from output files"
    print
    special_call(["./Plotter", output_file] + input_files)
    print "Plotting finished."
    print

def show_result(plotter_out):
    print "Showing the result in root: "
    print
    special_call(["root", "-l", 'scripts/show_result.C("'+plotter_out+'")'])

# ---------------------- Running -------------------------

if "--replot" in options:
    compile(0)
    backup_files(ORIGOUTDIR+"/backup_replot")
    #plotter_input_files = glob.glob(ORIGOUTDIR+"/*.root")
    plotter_input_file = max(glob.glob("results/Plotter_out_*.root"), key=os.path.getmtime).replace("_replot","")
    plotter_output_file = plotter_input_file.replace(".root","_replot.root")
    plotter([plotter_input_file], plotter_output_file)
    show_result(plotter_output_file)
else:
    compile(1, "--plot" in options)
    if not opt_test: backup_files(OUTDIR+"/backup")
    plotter_input_files = analysis(ana_arguments, NPROC)
    if "--plot" in options:
        plotter(plotter_input_files, plotter_output_file)
        show_result(plotter_output_file)

print "Done."


"""

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
    eval_or_echo "mkdir -p $OUTDIR/log $OUTDIR/backup/pileup $OUTDIR/backup/systematics"
    eval_or_echo "mkdir -p $OUTDIR/backup/filelists $OUTDIR/backup/common $OUTDIR/backup/scripts"
    eval_or_echo "cp -rp *.h *.cc Makefile* $OUTDIR/backup/"
    eval_or_echo "cp -rp pileup/* $OUTDIR/backup/pileup/"
    eval_or_echo "cp -rp systematics/* $OUTDIR/backup/systematics/"
    eval_or_echo "cp -rp filelists/* $OUTDIR/backup/filelists/"
    eval_or_echo "cp -rp common/*.h common/*.cc $OUTDIR/backup/common/"
    eval_or_echo "cp -rp scripts/*.C scripts/*.csh scripts/*.py $OUTDIR/backup/scripts/"
    
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
	eval_or_echo "par_source $SCRIPT 2"
    endif
    eval_or_echo "sleep 1s"
    eval_or_echo "mv $SCRIPT $OUTDIR/backup/"
    eval_or_echo "./Plotter results/Plotter_out_$DATE.root $all_output"
    eval_or_echo "if ( -e Latest_plots.root ) rm Latest_plots.root"
    eval_or_echo "ln -sf results/Plotter_out_$DATE.root Latest_plots.root"

else if ( $opt_replot ) then
    set infile=`ls -tr results/Plotter_out_* | grep -v replot | tail -1`
    set outfile=`echo $infile | sed "s;\.root;_replot\.root;"`
    set orig_outdir=`echo $infile | sed "s;Plotter_out_;;;s;\.root;;"`
    eval_or_echo "make clean"
    eval_or_echo "make Plotter ||  echo '\\nCompilation failed. Exiting...' && setenv makefailed && exit"
    eval_or_echo "mkdir -p $OUTDIR/backup_replot/pileup $OUTDIR/backup_replot/systematics"
    eval_or_echo "mkdir -p $OUTDIR/backup_replot/filelists $OUTDIR/backup_replot/common $OUTDIR/backup_replot/scripts"
    eval_or_echo "cp -rp *.h *.cc Makefile* $OUTDIR/backup_replot/"
    eval_or_echo "cp -rp pileup/* $OUTDIR/backup_replot/pileup/"
    eval_or_echo "cp -rp systematics/* $OUTDIR/backup_replot/systematics/"
    eval_or_echo "cp -rp filelists/* $OUTDIR/backup_replot/filelists/"
    eval_or_echo "cp -rp common/*.h common/*.cc $OUTDIR/backup_replot/common/"
    eval_or_echo "cp -rp scripts/*.C scripts/*.csh scripts/*.py $OUTDIR/backup_replot/scripts/"
    eval_or_echo "./Plotter $outfile $infile"
    eval_or_echo "if ( -e Latest_plots.root ) rm Latest_plots.root"
    eval_or_echo "ln -sf $outfile Latest_plots.root"
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
"""
