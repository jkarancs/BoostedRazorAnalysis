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

ana_arguments_full = []
ana_arguments_skim = []
for input_file in glob.glob("filelists/data/*.txt") + glob.glob("filelists/backgrounds/*.txt") + glob.glob("filelists/signals/*.txt"):
    output_file = OUTDIR+"/"+input_file.split("/")[-1].replace("txt", "root")
    skim_output_file = "ntuple/grid18/"+SKIMNAME+"/"+input_file.split("/")[-1].replace(".txt","/Skim.root")
    log_file = OUTDIR+"/log/"+input_file.split("/")[-1].replace("txt", "log")
    ana_arguments_full.append([output_file,      input_file, log_file])
    ana_arguments_skim.append([skim_output_file, input_file, log_file])

#ana_arguments_skim = [ana_arguments_skim[0], ana_arguments_skim[22], ana_arguments_skim[36]]

def subset(list, match):
    result = []
    for i in range(0, len(list)):
        if match in list[i][0]: result.append(list[i])
    return result

jetht                = subset(ana_arguments_full, "JetHT")
met                  = subset(ana_arguments_full, "MET")
single_ele           = subset(ana_arguments_full, "SingleElectron")
single_mu            = subset(ana_arguments_full, "SingleMuon")
dy_z                 = subset(ana_arguments_full, "DYJets")
dy_z                += subset(ana_arguments_full, "ZJets")
qcd                  = subset(ana_arguments_full, "QCD_HT")
#qcd                += subset(ana_arguments_full, "QCD_GenJets5_HT")
singletop            = subset(ana_arguments_full, "ST")
ttx                  = subset(ana_arguments_full, "TTG")
ttx                 += subset(ana_arguments_full, "TTW")
ttx                 += subset(ana_arguments_full, "TTZ")
ttx                 += subset(ana_arguments_full, "TTTT")
tt_madgraph_ht       = subset(ana_arguments_full, "TTJets_HT")
tt_mcatnlo           = subset(ana_arguments_full, "TTJets_amcatnloFXFX")
tt_madgraph          = subset(ana_arguments_full, "TTJets_madgraph")
tt_mcatnlo_her       = subset(ana_arguments_full, "TT_amcatnlo_pythia8")
tt_mcatnlo_py8       = subset(ana_arguments_full, "TT_amcatnlo_pythia8")
tt_powheg_her        = subset(ana_arguments_full, "TT_powheg_herwig")
tt_powheg_py8        = subset(ana_arguments_full, "TT_powheg_pythia8.txt")
tt_powheg_py8_mpiOff = subset(ana_arguments_full, "TT_powheg_pythia8_mpiOFF")
tt_powheg_py8_noCR   = subset(ana_arguments_full, "TT_powheg_pythia8_noCR")
wjets                = subset(ana_arguments_full, "WJets")
diboson              = subset(ana_arguments_full, "WW")
diboson             += subset(ana_arguments_full, "WZ")
diboson             += subset(ana_arguments_full, "ZZ")
diboson             += subset(ana_arguments_full, "ZH")
T1tttt               = subset(ana_arguments_full, "SMS-T1tttt_mGluino-1500_mLSP-100_FullSim")
#T1tttt              += subset(ana_arguments_full, "SMS-T1tttt_mGluino-1200_mLSP-800_FullSim")

#ana_arguments_test = [
#    ["results/Analyzer_test_data.root",  "filelists/data/JetHT_25ns_2015C.txt",           "results/Analyzer_test_data.log"],
#    ["results/Analyzer_test_ttbar.root", "filelists/backgrounds/TTJets_amcatnloFXFX.txt", "results/Analyzer_test_ttbar.log"],
#    ["results/Analyzer_test_qcd.root",   "filelists/backgrounds/QCD_HT700to1000.txt",     "results/Analyzer_test_qcd.log"],
#]
ana_arguments_test = jetht + tt_powheg_py8 + qcd


# Plotter:
# It gets input automatically from the output of Analyzer, output will be on screen
plotter_output_file = "results/Plotter_out_"+DATE+".root"

# ---------------------- Cmd Line  -----------------------

opt_dry = int(not "--run" in options)
opt_plot = int("--plot" in options)
opt_quick = 0
opt_test = 0
opt_skim = 0
if (opt_dry): print "--run option not specified, doing a dry run (only printing out commands)"

if "--full" in options:
    print "Running with option: --full"
    ana_arguments = ana_arguments_full
elif "--skim" in options:
    print "Running with option: --skim"
    ana_arguments = ana_arguments_skim
    opt_skim = 1
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
    global opt_dry, opt_quick, opt_plot, opt_skim
    if not opt_dry:
        print "Start Analyzing: "+input_list
    if not os.path.exists(os.path.dirname(output_file)):
        special_call(["mkdir", "-p", os.path.dirname(output_file)], 0)
    cmd = ["./Analyzer", output_file, input_list]
    if opt_skim and not opt_plot: cmd.append("noPlots=1")
    if opt_quick: cmd.append("quickTest=1")
    logged_call(cmd, output_log)
    if opt_skim:
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
