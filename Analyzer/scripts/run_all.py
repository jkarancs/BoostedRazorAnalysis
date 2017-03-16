import os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, ROOT
from optparse import OptionParser

# ---------------------- Cmd Line  -----------------------

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option("--run",     dest="run",     action="store_true", default=False,   help="Without this option, script only prints cmds it would otherwise excecute")
parser.add_option("--full",    dest="full",    action="store_true", default=False,   help="Run on all datasets found in filelists directory")
parser.add_option("--test",    dest="test",    action="store_true", default=False,   help="Run only on some test files (jetht, ttbar, qcd, T5ttcc)")
parser.add_option("--batch",   dest="batch",   action="store_true", default=False,   help="Send the jobs to batch")
parser.add_option("--queue",   dest="QUEUE",   type="string",       default="1nh",   help="Specify which batch queue to use (Default=1nh)")
parser.add_option("--quick",   dest="NQUICK",  type="int",          default=0,       help="Run only on a subset of events (1/NQUICK)")
parser.add_option("--nevt",    dest="NEVT",    type="int",          default=-1,      help="Tells how many event to run as a maximum in a single job (Default=-1 all)")
parser.add_option("--nfile",   dest="NFILE",   type="int",          default=-1,      help="Tells how many input files to run in a single job (Default=-1 all)")
parser.add_option("--sleep",   dest="SLEEP",   type="int",          default=3,       help="Wait for this number of seconds between submitting each batch job (Default 3s)")
parser.add_option("--useprev", dest="useprev", action="store_true", default=False,   help="Use previously created temporary filelists")
parser.add_option("--nproc",   dest="NPROC",   type="int",          default=1,       help="Tells how many parallel interactive jobs to start (Default=3)")
parser.add_option("--outdir",  dest="OUTDIR",  type="string",       default="",      help="Output directory (Default: results/run_[DATE])")
parser.add_option("--skimout", dest="SKIMOUT", type="string",       default="",      help="Output directory for skimming")
parser.add_option("--skim",    dest="skim",    action="store_true", default=False,   help="Skim output to --skimout directory (change in script)")
parser.add_option("--skimopt", dest="skimopt", action="store_true", default=False,   help="Optimize skim output size based on measured event ratios")
parser.add_option("--mirror",  dest="mirror",  action="store_true", default=False,   help="Also copy skim output to EOS")
parser.add_option("--plot",    dest="plot",    action="store_true", default=False,   help="Make plots after running using Plotter (Janos)")
parser.add_option("--replot",  dest="replot",  action="store_true", default=False,   help="Remake latest set of plots using Plotter (Janos)")
parser.add_option("--recover", dest="recover", action="store_true", default=False,   help="Recover stopped task (eg. due to some error)")
(opt,args) = parser.parse_args()

# ----------------------  Settings -----------------------
# Some further (usually) fixed settings, should edit them in this file

# Output directories/files
DATE = time.strftime("%Y_%m_%d_%Hh%Mm%S", time.localtime())
if opt.OUTDIR == "" and not opt.skim and not opt.replot:
    opt.OUTDIR = "results/run_"+DATE # log files, backup files, output files for non-skims

if opt.skim:
    COPYSCRIPT = ""
    if opt.OUTDIR == "":
        opt.OUTDIR = "results/skim_"+DATE # log files, backup files, output files for non-skims
    # Mirror also here
    if opt.SKIMOUT == "":
        print "ERROR: Give a suitable --skimout argument, eg. --skimout ntuple/grid18/Skim_Oct31_2Jet_1JetAK8"
        sys.exit()
    if opt.NFILE == -1 and opt.NEVT == -1 and not opt.useprev:
        print "ERROR: Give a suitable --nfile or --nevt argument, otherwise output might become too large!"
        sys.exit()
    if opt.NQUICK>1:
        if opt.mirror:
            print "ERROR: Please, don't mirror stuff to EOS, when testing!"
            sys.exit()
    else:
        if opt.mirror:
            # --mirror copies here
            EOS_JANOS  = "gsiftp://eoscmsftp.cern.ch//eos/cms/store/caf/user/jkarancs/B2GTTreeNtuple/"
        else:
            # If not, then makes a script for Viktor
            EOS_VIKTOR = "gsiftp://eoscmsftp.cern.ch//eos/cms/store/caf/user/veszpv/B2GTTreeNtuple/"
            COPYSCRIPT = opt.SKIMOUT.replace(opt.SKIMOUT.split("/")[-1],"")+"mirror_to_Viktors_EOS_"+opt.SKIMOUT.split("/")[-1]+".sh"
            print "Warning: Don't you want to mirror to EOS? Add: --mirror option!"
            print "         If not, ignore this message!"
            print "         Creating a copy script for Viktor: "+COPYSCRIPT
if opt.batch and opt.NEVT == -1 and not opt.useprev:
    print "ERROR: Give a suitable --nevt argument, otherwise some jobs will run too long on batch!"
    print "       Recommended option for 1nh queue: --nevt=2000000 (~30min/job, ~1000 evt/s is usual)"
    print "       Or use --useprev option to run on previously created temporary filelists"
    sys.exit()
if opt.plot:
    #PLOTTER_IN automatic
    if opt.test:
        PLOTTER_OUT = "results/Plotter_test.root"
    else:
        PLOTTER_OUT = opt.OUTDIR.replace("run_","Plotter_out_")+".root"
if opt.replot:
    if opt.OUTDIR == "":
        # Find last working directory automatically and find output files there
        opt.OUTDIR  = max(glob.glob("results/run_*"), key=os.path.getmtime)
        PLOTTER_IN  = max(glob.glob("results/Plotter_out_*.root"), key=os.path.getmtime).replace("_replot","")
        PLOTTER_OUT = PLOTTER_IN.replace(".root","_replot.root")
        PLOTTER_IN  = [PLOTTER_IN]
    else:
        PLOTTER_IN  = glob.glob(opt.OUTDIR+"/*.root")
        PLOTTER_OUT = opt.OUTDIR.replace("run_","Plotter_out_")+"_replot.root"

# Working directory, during running we cd here and back
if opt.test:
    EXEC_PATH = os.getcwd()
elif opt.replot:
    EXEC_PATH = opt.OUTDIR+"/backup_replot"
elif not opt.test:
    EXEC_PATH = opt.OUTDIR+"/backup"

# Print some options for logging
if not opt.run:
    print "--run option not specified, doing a dry run (only printing out commands)"

if opt.full:
    print "Running with option: --full"
if opt.skim:
    print "Running with option: --skim"
elif opt.replot:
    print "Running with option: --replot"
    opt.plot = 0 # for safety
elif opt.test:
    print "Running with option: --test (few files)"

if opt.plot:
    print "Running with option: --plot (will produce plots with Plotter)"

if opt.NQUICK>1:
    print "Running with option: --quick "+str(opt.NQUICK)+" (1/"+str(opt.NQUICK)+" statistics)"

# Some automatic filelists
if (opt.full):
    input_filelists  = glob.glob("filelists/data/*.txt")
    input_filelists += glob.glob("filelists/backgrounds/*.txt")
    input_filelists += glob.glob("filelists/signals/*.txt")
elif opt.test:
    input_filelists  = glob.glob("filelists/data/JetHT*.txt")
    input_filelists += glob.glob("filelists/backgrounds/QCD_HT*.txt")
    input_filelists += glob.glob("filelists/backgrounds/TT_powheg-pythia8_ext4*.txt")
    input_filelists += glob.glob("filelists/signals/FastSim_SMS-T5ttcc*.txt")
elif not opt.replot and len(args) < 1:
    print "Always tell me what filelists to run over (except when using --full or --test options)!"
    print "For more help, run as python %s -h" % (sys.argv[0])
    sys.exit()
else:
    input_filelists = args

# ----------------- Analyzer Arguments -------------------
# Analyzer (see below in functions):
# Each element supplies 3 arguments for each Analyzer job:
# [output filename, input file list, output log]
# For skimming/full running, all datasets are used
# for testing a selected few

if opt.useprev:
    print "Reusing previously created temporary filelists for split jobs (eg. --batch) in filelists_tmp/:"
elif (opt.NFILE != -1 or opt.NEVT != -1):
    print "Start creating new temporary filelists for split jobs (eg. batch) in filelists_tmp/:"
    for tmp_txtfile in glob.glob('filelists_tmp/*/*.txt'): os.remove(tmp_txtfile)

ana_arguments = []
# Loop over all filelists
if opt.NEVT != -1: bad_files = open("bad_files_found.txt", "w")
for filelist in input_filelists:
    # Will put all files into the OUTDIR and its subdirectories
    log_file = opt.OUTDIR+"/log/"+filelist.split("/")[-1].replace("txt", "log")
    if opt.skim:
        # Except for skim, where we send the large output files to a different directory
        # keeping subdirectory structure (suitable for a future input )
        output_file = opt.SKIMOUT+"/"+filelist.split("/")[-1].replace(".txt","/Skim.root")
    else:
        output_file = opt.OUTDIR +"/"+filelist.split("/")[-1].replace("txt", "root")
    # Now let's make the argument list for the Analyzer jobs
    options = []
    if opt.NQUICK>1: options.append("quickTest="+str(opt.NQUICK))
    if opt.skim and not opt.plot: options.append("noPlots=1")
    # Temporary filelists
    if opt.useprev:
        # Use previously created lists
        if not opt.skim: options.append("fullFileList="+EXEC_PATH+"/"+filelist) # Need full ntuple to correctly normalize weights
        prev_lists = glob.glob(filelist.replace("filelists","filelists_tmp").replace(".txt","_[0-9]*.txt"))
        for jobnum in range(1, len(prev_lists)+1):
            tmp_filelist = prev_lists[jobnum-1]
            args = [output_file.replace(".root","_"+str(jobnum)+".root"), [EXEC_PATH+"/"+tmp_filelist], options, log_file.replace(".log","_"+str(jobnum)+".log")]
            ana_arguments.append(args)
    elif opt.NEVT != -1:
        # SPLIT MODE (recommended for batch): Each jobs runs on max opt.NEVT
        JOB_NEVT = opt.NEVT
        # Further optimize this number for skimming
        # based on measured unskimmed to skimmed ratios (found in skim_ratios.txt)
        if opt.skimopt:
            samplename = filelist.split("/")[-1][:-4]
            with open("skim_ratios.txt") as ratios:
                for line in ratios:
                    column = line.split()
                    if samplename == column[2]:
                        JOB_NEVT *= int(column[0])
            #print str(JOB_NEVT)+" "+samplename
        # Need full ntuple to correctly normalize weights
        if not opt.skim: options.append("fullFileList="+EXEC_PATH+"/"+filelist) # Need full ntuple to correctly normalize weights
        # loop on file lists and split to tmp_filists for nevt < JOB_NEVT
        with open(filelist) as f:
            files = f.read().splitlines()
            jobnum = 0
            totalevt = 0
            for i in range(0, len(files)):
                # First get the number of events in the file
                f = ROOT.TFile.Open(files[i])
                if not f:
                    print files[i]+" is not a root file"
                    print>>bad_files, files[i]+" bad file"
                    continue
                tree = f.Get("B2GTree")
                if not tree: tree = f.Get("B2GTTreeMaker/B2GTree")
                nevt = tree.GetEntries()
                #print str(nevt)+" "+files[i] # comment in if need event number per files (for skim_ratios.txt)
                if nevt == 0:
                    print files[i]+" has 0 entries"
                    print>>bad_files, files[i]+" 0 entry"
                # Create a new list after every JOB_NEVT
                if i==0 or (totalevt + nevt > JOB_NEVT):
                    jobnum += 1
                    tmp_filelist = filelist.replace("filelists","filelists_tmp").replace(".txt","_"+str(jobnum)+".txt")
                    args = [output_file.replace(".root","_"+str(jobnum)+".root"), [EXEC_PATH+"/"+tmp_filelist], options, log_file.replace(".log","_"+str(jobnum)+".log")]
                    ana_arguments.append(args)
                    totalevt = 0
                totalevt += nevt
                with open(tmp_filelist, "a") as job_filelist:
                    print>>job_filelist, files[i]
        print "  "+filelist.replace("filelists","filelists_tmp").replace(".txt","_*.txt")+" created"
    elif opt.NFILE != -1:
        # SPLIT MODE: Each jobs runs on max opt.NFILE
        if not opt.skim: options.append("fullFileList="+EXEC_PATH+"/"+filelist) # Need full ntuple to correctly normalize weights
        with open(filelist) as f:
            files = f.read().splitlines()
            for n in range(1, (len(files)-1)/opt.NFILE+2):
                tmp_filelist = filelist.replace("filelists","filelists_tmp").replace(".txt","_"+str(n)+".txt")
                with open(tmp_filelist, "w") as job_filelist:
                    for i in range((n-1)*opt.NFILE, min(n*opt.NFILE,len(files))):
                        print>>job_filelist, files[i]
                args = [output_file.replace(".root","_"+str(n)+".root"), [EXEC_PATH+"/"+tmp_filelist], options, log_file.replace(".log","_"+str(n)+".log")]
                ana_arguments.append(args)
    else:
        # In case of a single job/dataset
        ana_arguments.append([output_file, [EXEC_PATH+"/"+filelist], options, log_file])


if opt.NEVT != -1: bad_files.close()

if opt.NFILE != -1 or opt.NEVT != -1 and not opt.useprev:
    print "All temporary filelist ready."

# --------------------- Functions ------------------------
# Show and run command with stdout on screen
icommand=0
def special_call(cmd, verbose=1):
    global icommand, opt
    if verbose:
        if opt.run:
            print("[%d]" % icommand),
        else:
            print("[dry]"),
        for i in xrange(len(cmd)): print cmd[i],
        print ""
    if opt.run:
        if subprocess.call(cmd):
            print "ERROR: Problem executing command:"
            print("[%d]" % icommand)
            for i in xrange(len(cmd)): print cmd[i],
            print ""
            print "exiting."
            sys.exit()
        if verbose: print ""
    sys.stdout.flush()
    icommand+=1

# Run command with stdout/stderr saved to logfile
def logged_call(cmd, logfile):
    global opt
    dirname = os.path.dirname(logfile)
    if dirname != "" and not os.path.exists(dirname):
        special_call(["mkdir", "-p", os.path.dirname(logfile)], 0)
    if opt.run:
        with open(logfile, "w") as log:
            proc = subprocess.Popen(cmd, stdout=log, stderr=log, close_fds=True)
            proc.wait()
    else:
        proc = subprocess.call(["echo", "[dry]"]+cmd+[">", logfile])

# Compile programs
def compile(Ana = 1, Plotter = 1):
    global opt, EXEC_PATH
    print "Compiling ..."
    print
    saved_path = os.getcwd()
    if opt.run: os.chdir(EXEC_PATH)
    special_call(["make", "clean"])
    if Ana: special_call(["make", "Analyzer"])
    if Plotter: special_call(["make", "Plotter"])
    if opt.run: os.chdir(saved_path)
    print "Compilation successful."
    print

# backup files for bookkeeping
def backup_files(backup_dir):
    print "Backing up files in: "+backup_dir
    print
    special_call(["mkdir", "-p", backup_dir])
    #special_call(["cp", "-rp", "btag_eff", "pileup", "scale_factors", "trigger_eff", "common", "filelists", "filelists_tmp", "scripts", "systematics"] + glob.glob("*.h") + glob.glob("*.cc") + glob.glob("Makefile*") + [backup_dir+"/"])
    special_call(["cp", "-rp", "pileup", "scale_factors", "common", "filelists", "filelists_tmp", "scripts", "systematics"] + glob.glob("*.h") + glob.glob("*.cc") + glob.glob("Makefile*") + [backup_dir+"/"])
    print



# Run a single Analyzer instance (on a single input list, i.e. one dataset)
def analyzer_job((jobindex)):
    global ana_arguments, opt, EXEC_PATH, COPYSCRIPT
    output_file = ana_arguments[jobindex][0]
    input_list  = ana_arguments[jobindex][1]
    options     = ana_arguments[jobindex][2]
    output_log  = ana_arguments[jobindex][3]
    if opt.run:
        if opt.batch:
            print "Sending job to batch (queue: "+opt.QUEUE+"), expected output: "+output_file
        else:
            print "Start Analyzing, expected output: "+output_file
    if not os.path.exists(os.path.dirname(output_file)):
        special_call(["mkdir", "-p", os.path.dirname(output_file)], 0)
    cmd = [EXEC_PATH+"/Analyzer", output_file] + options + input_list
    if opt.batch:
        logdirname = os.path.dirname(output_log)
        if logdirname != "" and not os.path.exists(logdirname): special_call(["mkdir", "-p", logdirname], 0)
        cmd = shlex.split('bsub -q '+opt.QUEUE+' -J '+DATE+'_'+str(jobindex)+' -oo '+output_log+' -L /bin/bash '+os.getcwd()+'/scripts/Analyzer_batch_job.sh '+os.getcwd())+cmd
        special_call(cmd, not opt.run)
    else:
        if opt.NPROC>3: cmd = ['nice']+cmd
        logged_call(cmd, output_log)
    # Mirror output (copy to EOS)
    if opt.batch: time.sleep(opt.SLEEP)
    elif opt.skim:
        outpath = output_file.split("/")[-3]+"/"+output_file.split("/")[-2]+"/"+output_file.split("/")[-1]
        if opt.mirror: logged_call(["env", "--unset=LD_LIBRARY_PATH", "gfal-copy", "-r", output_file, EOS_JANOS+outpath], output_log)
        elif opt.run and opt.NQUICK==0:
            with open(COPYSCRIPT, "a") as script:
                print>>script, 'env --unset=LD_LIBRARY_PATH gfal-copy -r '+output_file+' '+EOS_VIKTOR+outpath
                #print>>script, 'srm-set-permissions -type=CHANGE -group=RW '+EOS_VIKTOR+outpath
    return output_file


# Run all Analyzer jobs in parallel
def analysis(ana_arguments, nproc):
    global opt
    njob = len(ana_arguments)
    if njob<nproc: nproc = njob
    print "Running "+str(njob)+" instances of Analyzer jobs:"
    print
    saved_path = os.getcwd()
    if opt.batch:
        # Running on the batch and babysitting all jobs until completion
        output_files = []
        njob = len(ana_arguments)
        if opt.run:
            last_known_status = [-1] * njob # -1: start, 0: finished, time(): last submit/status check
            finished = 0
            # Loop until all jobs are finished
            while finished != njob:
                if finished != 0: time.sleep(30)
                finished = 0
                for jobindex in range(0, njob):
                    output_file = ana_arguments[jobindex][0]
                    #output_log  = ana_arguments[jobindex][3]
                    if last_known_status[jobindex] == -1:
                        if opt.recover and os.path.isfile(output_file):
                            finished += 1
                            output_files.append(output_file)
                            last_known_status[jobindex] = 0
                        else:
                            # Submit all jobs
                            analyzer_job(jobindex)
                            last_known_status[jobindex] = time.time()
                    elif last_known_status[jobindex] == 0:
                        finished += 1
                    else:
                        # First check if output file exists
                        if os.path.isfile(output_file):
                            finished += 1
                            output_files.append(output_file)
                            last_known_status[jobindex] = 0
                        # If the last submission/check is older than 10 minutes check job status with bjobs
                        elif time.time() - last_known_status[jobindex] > (600 if opt.NQUICK<2 else 600/opt.NQUICK):
                            jobname = DATE+'_'+str(jobindex)
                            logged_call(shlex.split('bjobs -J '+jobname), 'jobstatus_'+jobname+'.txt')
                            with open('jobstatus_'+jobname+'.txt') as jobstatus:
                                lines = jobstatus.readlines()
                                if 'Job <'+jobname+'> is not found' in lines[0]:
                                    # Job is missing, so resubmit
                                    print "Job "+jobname+" is missing, resubmitting ...         "
                                    analyzer_job(jobindex)
                                #else:
                                #    status = lines[1].split()[2]
                                #    if status == 'PEND' or status == 'RUN':
                                #        print "Job "+jobname+" is pending/running"
                            os.remove('jobstatus_'+jobname+'.txt')
                            last_known_status[jobindex] = time.time()
                # Finally print status
                print "Analyzer jobs on batch (Done/All): "+str(finished)+"/"+str(njob)+"   \r",
                sys.stdout.flush()
            print "\nAll batch jobs finished."
        else:
            for jobindex in range(0, njob):
                output_files.append(analyzer_job(jobindex))
    else:
        # Use the N CPUs in parallel on the current computer to analyze all jobs
        workers = multiprocessing.Pool(processes=nproc)
        njob = len(ana_arguments)
        output_files = workers.map(analyzer_job, range(njob), chunksize=1)
        workers.close()
        workers.join()
        print "All Analyzer jobs finished."
    print
    return output_files

# Run Plotter, output of Analyzer is input for this code
def plotter(input_files, output_file):
    global opt, EXEC_PATH
    print "Start plotting from output files"
    print
    special_call([EXEC_PATH+"/Plotter", output_file] + input_files)
    print "Plotting finished."
    print

def show_result(plotter_out):
    print "Showing the result in root: "
    print
    special_call(["root", "-l", 'scripts/show_result.C("'+plotter_out+'")'])

# ---------------------- Running -------------------------

# renew token before running
if opt.replot:
    if not opt.recover:
        backup_files(EXEC_PATH)
        compile(0)
    plotter(PLOTTER_IN, PLOTTER_OUT)
    #if not 'lxplus' in socket.gethostname():
    #    show_result(PLOTTER_OUT)
else:
    if not opt.recover:
        if not opt.test:
            backup_files(EXEC_PATH)
        compile(1, opt.plot)
    plotter_input_files = analysis(ana_arguments, opt.NPROC)
    if opt.plot:
        plotter(plotter_input_files, PLOTTER_OUT)
        #if not 'lxplus' in socket.gethostname():
        #    show_result(PLOTTER_OUT)

print "Done."
