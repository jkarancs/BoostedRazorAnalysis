import re, os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, getpass, ROOT
from optparse import OptionParser
from common_functions import *

# ---------------------- Cmd Line  -----------------------

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option("--run",         dest="run",         action="store_true", default=False,   help="Without this option, script only prints cmds it would otherwise excecute")
parser.add_option("--full",        dest="full",        action="store_true", default=False,   help="Run on all datasets found in filelists directory")
parser.add_option("--test",        dest="test",        action="store_true", default=False,   help="Run only on some test files (jetht, ttbar, qcd, T5ttcc)")
parser.add_option("--batch",       dest="batch",       action="store_true", default=False,   help="Send the jobs to batch")
parser.add_option("--queue",       dest="QUEUE",       type="string",       default="1nh",   help="Specify which batch queue to use (Default=1nh)")
parser.add_option("--optim",       dest="optim",       action="store_true", default=False,   help="Optimize job event number based log files in --prevdir, or measured skim ratios")
parser.add_option("--prevdir",     dest="PREVDIR",     type="string",       default="",      help="Previous running directory used to optimize jobs (default=last dir in results/)")
parser.add_option("--jobtime",     dest="JOBTIME",     type="int",          default=1500,    help="Desired job running time in s (default=1500)")
parser.add_option("--quick",       dest="NQUICK",      type="int",          default=0,       help="Run only on a subset of events (1/NQUICK)")
parser.add_option("--nevt",        dest="NEVT",        type="int",          default=-1,      help="Tells how many event to run as a maximum in a single job (Default=-1 all)")
parser.add_option("--nfile",       dest="NFILE",       type="int",          default=-1,      help="Tells how many input files to run in a single job (Default=-1 all)")
parser.add_option("--sleep",       dest="SLEEP",       type="int",          default=3,       help="Wait for this number of seconds between submitting each batch job (Default 3s)")
parser.add_option("--useprev",     dest="useprev",     action="store_true", default=False,   help="Use previously created temporary filelists")
parser.add_option("--nproc",       dest="NPROC",       type="int",          default=1,       help="Tells how many parallel interactive jobs to start (Default=3)")
parser.add_option("--outdir",      dest="OUTDIR",      type="string",       default="",      help="Output directory (Default: results/run_[SUBTIME])")
parser.add_option("--skimout",     dest="SKIMOUT",     type="string",       default="",      help="Output directory for skimming")
parser.add_option("--skim",        dest="skim",        action="store_true", default=False,   help="Skim output to --skimout directory (change in script)")
parser.add_option("--mirror",      dest="mirror",      action="store_true", default=False,   help="Also copy skim output to EOS")
parser.add_option("--mirror_user", dest="mirror_user", action="store_true", default=False,   help="Also copy skim output to Janos' EOS")
parser.add_option("--plot",        dest="plot",        action="store_true", default=False,   help="Make plots after running using Plotter (Janos)")
parser.add_option("--replot",      dest="replot",      action="store_true", default=False,   help="Remake latest set of plots using Plotter (Janos)")
parser.add_option("--recover",     dest="recover",     action="store_true", default=False,   help="Recover stopped task (eg. due to some error)")
parser.add_option("--nohadd",      dest="nohadd",      action="store_true", default=False,   help="Disable hadding output files")
parser.add_option("--nocheck",     dest="nocheck",     action="store_true", default=False,   help="Disable nevent check for output files")
parser.add_option("--haddonly",    dest="haddonly",    action="store_true", default=False,   help="Do not submit any jobs, only merge output")
(opt,args) = parser.parse_args()

# ----------------------  Settings -----------------------
# Some further (usually) fixed settings, should edit them in this file

# Make sure we start with a fresh token (lasts 24h)
currtime = time.time() # used for tokens/recovery jobs
special_call(["kinit", "-R"], 0)

# Output directories/files
SUBTIME = time.strftime("%Y_%m_%d_%Hh%Mm%S", time.localtime())
TMPDIR = "/tmp/"+getpass.getuser()+"/"
if opt.OUTDIR == "" and not opt.skim and not opt.replot:
    opt.OUTDIR = "results/run_"+SUBTIME # log files, backup files, output files for non-skims

if opt.skim:
    COPYSCRIPT = ""
    if opt.OUTDIR == "":
        opt.OUTDIR = "results/skim_"+SUBTIME # log files, backup files, output files for non-skims
    # Mirror also here
    if opt.SKIMOUT == "":
        print "ERROR: Give a suitable --skimout argument, eg. --skimout ntuple/grid18/Skim_Oct31_2Jet_1JetAK8"
        sys.exit()
    if opt.NFILE == -1 and opt.NEVT == -1 and not opt.optim and not opt.useprev:
        print "ERROR: Give a suitable --nfile, --nevt or --optim argument, otherwise output might become too large!"
        sys.exit()
    if opt.NQUICK>1:
        if opt.mirror or opt.mirror_user:
            print "ERROR: Please, don't mirror stuff to EOS, when testing!"
            sys.exit()
    else:
        if opt.mirror:
            # --mirror copies here
            EOS_JANOS  = "gsiftp://eoscmsftp.cern.ch//eos/cms/store/caf/user/jkarancs/B2GTTreeNtuple/"
        elif opt.mirror_user:
            # --mirror copies here
            EOS_JANOS  = "root://eosuser.cern.ch//eos/user/j/jkarancs/B2GTTreeNtuple/"
        else:
            # If not, then makes a script for Viktor
            EOS_VIKTOR = "gsiftp://eoscmsftp.cern.ch//eos/cms/store/caf/user/veszpv/B2GTTreeNtuple/"
            COPYSCRIPT = opt.SKIMOUT.replace(opt.SKIMOUT.split("/")[-1],"")+"mirror_to_Viktors_EOS_"+opt.SKIMOUT.split("/")[-1]+".sh"
            print "Warning: Don't you want to mirror to EOS? Add: --mirror option!"
            print "         If not, ignore this message!"
            print "         Creating a copy script for Viktor: "+COPYSCRIPT
if opt.batch and opt.NEVT == -1 and not opt.optim and not opt.useprev:
    print "ERROR: Give a suitable --nevt or --optim argument, otherwise some jobs will run too long on batch!"
    print "       Recommended to start with: --queue=8nh --nevt=1000000"
    print "       Then next time, you can use default 1nh queue with ~1500s jobs using: --optim"
    print "       If too many jobs fail, try lowering job runtime eg: --optim --jobtime=1200"
    print "       Or use --useprev option to run on previously created temporary filelists"
    sys.exit()
#if opt.optim and not opt.skim and opt.PREVDIR == "":
#    sorted_dirs = sorted(glob.glob("results/*/"), key=os.path.getmtime)
#    if len(sorted_dirs):
#        opt.PREVDIR = sorted_dirs[-1][:-1]
#        print "Using --optim option, but no --prevdir argument was given, using previous working directory: "+opt.PREVDIR
#    else:
#        print "ERROR: Using --optim option, but no --prevdir argument was given, and cannot find previous working directory in results/*"
#        sys.exit()
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
else:
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
if opt.recover:
    saved_path = os.getcwd()
    os.chdir(EXEC_PATH)
    input_filelists  = glob.glob("filelists/data/*.txt")
    input_filelists += glob.glob("filelists/backgrounds/*.txt")
    input_filelists += glob.glob("filelists/signals/*.txt")    
    ####    input_filelists = glob.glob("filelists/signals/FastSim_SMS-T5ttcc_*.txt")    
    os.chdir(saved_path)
elif opt.full:
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

# Read some options from included settings_*.h file
with open("Analyzer.cc") as ana:
    for line in ana:
        if '#include "settings_' in line:
            settings_file = line.split()[1].replace('"','')
vary_syst = False
with open(settings_file) as settings:
    for line in settings:
        if "define SYST" in line and not "0" in line:
            vary_syst = True


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
if opt.NEVT != -1 or opt.optim: bad_files = open("bad_files_found.txt", "w")
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
        #print filelist
        #print len(ana_arguments)
        # Use previously created lists
        if not opt.skim: options.append("fullFileList="+filelist) # Need full ntuple to correctly normalize weights
        if opt.recover:
            saved_path = os.getcwd()
            os.chdir(EXEC_PATH)
        prev_lists = glob.glob(filelist.replace("filelists","filelists_tmp").replace(".txt","_[0-9]*.txt"))
        new_lists = []
        for jobnum in range(1,1+len(prev_lists)):
            new_lists.append(filelist.replace("filelists","filelists_tmp").replace(".txt","_"+str(jobnum)+".txt"))
        #prev_lists = new_lists
        if opt.recover:
            os.chdir(saved_path)
        for i in range(0, len(new_lists)):
            tmp_filelist = new_lists[i]
            #jobnum = tmp_filelist.replace(".txt","").split("_")[-1]
            jobnum = i+1
            job_args = [output_file.replace(".root","_"+str(jobnum)+".root"), [EXEC_PATH+"/"+tmp_filelist], options, log_file.replace(".log","_"+str(jobnum)+".log")]
            ana_arguments.append(job_args)        
    elif opt.NEVT != -1 or opt.optim:
        # SPLIT MODE (recommended for batch): Each jobs runs on max opt.NEVT
        JOB_NEVT = opt.NEVT
        # Further optimize this number
        # based on measured unskimmed to skimmed ratios (found in skim_ratios.txt)
        optim_found = False
        if opt.optim:
            samplename = filelist.split("/")[-1][:-4]
            # Skimming uses a text file (skim_ratios.txt)
            if opt.skim:
                with open("skim_ratios.txt") as ratios:
                    for line in ratios:
                        column = line.split()
                        #if samplename == column[2]:
                        if column[1] in samplename and not optim_found:
                            JOB_NEVT = int(float(column[0])*JOB_NEVT)
                            optim_found = True
            # Other jobs use the results from a previous run
            else:
                # Previous txt file method (job_ratios.txt)
                with open("systjob_ratios.txt" if vary_syst else "job_ratios.txt") as ratios:
                    for line in ratios:
                        column = line.split()
                        if samplename == column[1]:
                            JOB_NEVT *= float(column[0]) * 0.8
                            optim_found = True
                # if no optimization found, use 0.2
                if not optim_found:
                    JOB_NEVT *= 0.2
                ##  # Get information from log files in opt.PREVDIR
                ##  min_nps  = float(1e6)
                ##  tot_nevt = float(0)
                ##  tot_time = float(0)
                ##  for log_file in sorted(glob.glob(opt.PREVDIR+"/log/"+samplename+"*.log")):
                ##      with open (log_file) as log:
                ##          for line in log:
                ##              if re.search("JobMonitoringReport", line):
                ##                  column  = line.split()
                ##                  runtime = column[2]
                ##                  nevt    = column[4]
                ##                  nps     = column[6]
                ##                  #print samplename+" "+nevt+"/"+runtime+" = "+nps
                ##                  if (float(nps)<min_nps): min_nps = float(nps)
                ##                  tot_nevt += float(nevt)
                ##                  tot_time += float(runtime)
                ##  max_nevt = int(opt.JOBTIME * min_nps)
                ##  #print ("%-50s   MIN=%7.1f  AVG=%7.1f   MAX_NEVT=%10d" % (samplename, min_nps, (tot_nevt/tot_time), max_nevt))
                ##  JOB_NEVT = max_nevt
        # Need full ntuple to correctly normalize weights
        if not opt.skim: options.append("fullFileList="+EXEC_PATH+"/"+filelist) # Need full ntuple to correctly normalize weights
        # loop on file lists and split to tmp_filists for nevt < JOB_NEVT
        with open(filelist) as f:
            files = f.read().splitlines()
            jobnum = 0
            totalevt = 0
            for i in range(0, len(files)):
                # First get the number of events in the file
                #print files[i]
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
                    job_args = [output_file.replace(".root","_"+str(jobnum)+".root"), [EXEC_PATH+"/"+tmp_filelist], options, log_file.replace(".log","_"+str(jobnum)+".log")]
                    ana_arguments.append(job_args)
                    totalevt = 0
                totalevt += nevt
                ntry = 0
                while True:
                    try:
                        with open(tmp_filelist, "a") as job_filelist:
                            print>>job_filelist, files[i]
                    except:
                        print "Warning: Could not write to disk (IOError), wait 10s and continue"
                        time.sleep(10)
                        ntry += 1
                        if ntry == 20: sys.exit()
                        continue
                    break
        if opt.optim:
            if optim_found:
                print "  "+filelist.replace("filelists","filelists_tmp").replace(".txt","_*.txt")+" created (MAX NEVT (optim) = "+str(JOB_NEVT)+")"
            else:
                print "  "+filelist.replace("filelists","filelists_tmp").replace(".txt","_*.txt")+" created (MAX NEVT (optim) = "+str(JOB_NEVT)+") - sample not in job_ratios.txt"
        else:
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
                job_args = [output_file.replace(".root","_"+str(n)+".root"), [EXEC_PATH+"/"+tmp_filelist], options, log_file.replace(".log","_"+str(n)+".log")]
                ana_arguments.append(job_args)
    else:
        # In case of a single job/dataset
        ana_arguments.append([output_file, [EXEC_PATH+"/"+filelist], options, log_file])

#print "Number of jobs: "+str(len(ana_arguments))
#sys.exit()

# Recover info about previously submitted jobs
last_known_status = [-1] * len(ana_arguments) # -1: start, 0: finished, 1: recovered, time(): last submit/status check
if opt.recover and opt.batch:
    # Check if the submission time is available
    if os.path.exists(EXEC_PATH+"/submission_time.txt"):
        with open(EXEC_PATH+"/submission_time.txt") as sub_time:
            SUBTIME=sub_time.readline().replace('\n','')
    print "Recovering previous jobs with submission time: "+SUBTIME
    # Check status of all running/pending jobs on batch
    nrecov = 0
    logged_call(shlex.split('bjobs -W -noheader'), TMPDIR+'batchstatus_'+SUBTIME+'.txt', opt.run)
    with open(TMPDIR+'batchstatus_'+SUBTIME+'.txt') as batchstatus:
        lines = batchstatus.readlines()
        for line in lines:
            if not "job found" in line: # output can be: "No unfinished job found"
                jobname = line.split()[6]
                if jobname.startswith(SUBTIME):
                    jobindex = int(jobname.split("_")[-1])
                    last_known_status[jobindex] = currtime
                    nrecov += 1
    os.remove(TMPDIR+'batchstatus_'+SUBTIME+'.txt')
    print "Successfully recovered "+str(nrecov)+" jobs from previous submission"

if opt.NEVT != -1:
    ntry = 0
    while True:
        try:
            bad_files.close()
        except:
            print "Warning: Could not close file (IOError), wait 10s and continue"
            time.sleep(10)
            ntry += 1
            if ntry == 20: sys.exit()
            continue
        break

if opt.NFILE != -1 or opt.NEVT != -1 and not opt.useprev:
    print "All temporary filelist ready."

# --------------------- Functions ------------------------

# Compile programs
def compile(Ana = 1, Plotter = 1):
    global opt, EXEC_PATH
    print "Compiling ..."
    print
    saved_path = os.getcwd()
    if opt.run: os.chdir(EXEC_PATH)
    special_call(["make", "clean"], opt.run)
    if Ana:
        special_call(["make", "-j8", "Analyzer"], opt.run)
        special_call(["chmod", "777", "Analyzer"], opt.run)
    if Plotter:
        special_call(["make", "-j8", "Plotter"], opt.run)
        special_call(["chmod", "777", "Plotter"], opt.run)
    if opt.run: os.chdir(saved_path)
    print "Compilation successful."
    print

# backup files for bookkeeping
def backup_files(backup_dir, submission_time):
    print "Backing up files in: "+backup_dir
    print
    special_call(["mkdir", "-p", backup_dir], opt.run)
    special_call(["cp", "-rp", "btag_eff", "pileup", "scale_factors", "trigger_eff", "common", "filelists", "filelists_tmp", "scripts", "systematics"] + glob.glob("*.h") + glob.glob("*.cc") + glob.glob("Makefile*") + [backup_dir+"/"], opt.run)
    with open(backup_dir+"/submission_time.txt","w") as sub_time:
        print>>sub_time, submission_time
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
        special_call(["mkdir", "-p", os.path.dirname(output_file)], opt.run, 0)
    cmd = [EXEC_PATH+"/Analyzer", output_file] + options + input_list
    if opt.batch:
        logdirname = os.path.dirname(output_log)
        if logdirname != "" and not os.path.exists(logdirname): special_call(["mkdir", "-p", logdirname], opt.run, 0)
        if os.getcwd().startswith("/afs"):
            cmd = shlex.split('bsub -q '+opt.QUEUE+' -J '+SUBTIME+'_'+str(jobindex)+' -oo '+output_log+' -L /bin/bash '+os.getcwd()+'/scripts/Analyzer_batch_job.sh '+os.getcwd()+' '+output_log)+cmd
        else:
            # Currently bsub cannot send the log file to EOS, so in order to avoid annoying e-mails and LSFJOB directories,
            # we send the output to a dummy afs file. The log output will be copied inside the script instead
            #job_log = "/tmp/"+getpass.getuser()+"/"+os.path.dirname(opt.OUTDIR+"/").split("/")[-1]+"/"+os.path.basename(output_log)
            job_log = "/tmp/"+os.path.basename(output_log)
            #if not os.path.exists(os.path.dirname(job_log)):
            #    special_call(["mkdir", "-p", os.path.dirname(job_log)], opt.run, 0)
            #    special_call(['chmod', '-R', '777', "/tmp/"+getpass.getuser()], opt.run, 0)
            #cmd = shlex.split('bsub -q '+opt.QUEUE+' -J '+SUBTIME+'_'+str(jobindex)+' -oo '+job_log+' -L /bin/bash '+os.getcwd()+'/scripts/Analyzer_batch_job.sh '+os.getcwd()+' '+output_log)+cmd
            if vary_syst:
                cmd = shlex.split('bsub -M 3500000 -v 5000000 -q '+opt.QUEUE+' -J '+SUBTIME+'_'+str(jobindex)+' -oo '+job_log+' -L /bin/bash '+os.getcwd()+'/scripts/Analyzer_batch_job.sh '+os.getcwd()+' '+output_log)+cmd
            else:
                cmd = shlex.split('bsub -M 1500000 -v 3000000 -q '+opt.QUEUE+' -J '+SUBTIME+'_'+str(jobindex)+' -oo '+job_log+' -L /bin/bash '+os.getcwd()+'/scripts/Analyzer_batch_job.sh '+os.getcwd()+' '+output_log)+cmd                
        special_call(cmd, opt.run, not opt.run)
    else:
        if opt.NPROC>3: cmd = ['nice']+cmd
        logged_call(cmd, output_log, opt.run)
    # Mirror output (copy to EOS)
    if opt.batch: time.sleep(opt.SLEEP)
    elif opt.skim:
        outpath = output_file.split("/")[-3]+"/"+output_file.split("/")[-2]+"/"+output_file.split("/")[-1]
        if opt.mirror or opt.mirror_user: logged_call(["env", "--unset=LD_LIBRARY_PATH", "gfal-copy", "-f", "-r", output_file, EOS_JANOS+outpath], output_log, opt.run)
        elif opt.run and opt.NQUICK==0:
            with open(COPYSCRIPT, "a") as script:
                print>>script, 'env --unset=LD_LIBRARY_PATH gfal-copy -r '+output_file+' '+EOS_VIKTOR+outpath
                #print>>script, 'srm-set-permissions -type=CHANGE -group=RW '+EOS_VIKTOR+outpath
    return output_file

def merge_output(ana_arguments, last_known_status):
    # Check list of files ready to be merged (with hadd)
    if not os.path.exists(opt.OUTDIR+"/hadd"): special_call(["mkdir", "-p", opt.OUTDIR+"/hadd"], opt.run, 0)
    prev_sample = ""
    mergeables = []
    all_mergeables = []
    njob = len(ana_arguments)
    for i in range(0, njob):
        job_index = ana_arguments[i][0][:-5].split("_")[-1]
        sample = ana_arguments[i][0][:-(6+len(job_index))]
        if sample != prev_sample:
            prev_sample = sample
            ready_to_merge = True
            if len(mergeables)>1: all_mergeables.append(mergeables)
            mergeables = [sample.rsplit("/",1)[0]+"/hadd/"+sample.rsplit("/",1)[1]+".root"]
        if ready_to_merge:
            if last_known_status[i]==0:
                mergeables.append(ana_arguments[i][0])
            else:
                ready_to_merge = False
                mergeables = []
        #print ("%4d - %d - missing=%d - %s" % (i, ready_to_merge, last_known_status[i]==-1, ana_arguments[i][0]))
    if ready_to_merge: all_mergeables.append(mergeables)
    # Merge them if they are ready
    for i in range(0, len(all_mergeables)):
        output   = all_mergeables[i][0]
        log      = output.rsplit("/",1)[0]+"/log/"+output.rsplit("/",1)[1].replace(".root",".log")
        allinput = all_mergeables[i][1:]
        if not os.path.exists(output):
            if len(allinput)==1:
                # Single files can simply be copied
                print "File for "+all_mergeables[i][0]+" is ready"
                special_call(["cp","-p", allinput[0], output], opt.run, 0)
            else:
                # Multiple files will be hadded
                print str(len(allinput))+" files for "+output+" are ready to be merged"
                while not os.path.exists(output):
                    #logged_call(["hadd", "-f", "-v", "-n", "200"]+all_mergeables[i], all_mergeables[i][0].rsplit("/",1)[0]+"/log/"+all_mergeables[i][0].rsplit("/",1)[1].replace(".root",".log"), opt.run)
                    # hadd produces problems when merging too many files
                    # so we merge files in chunks of 100 files each
                    # problems happen typically over a few hundred input files
                    Nmerge = 100
                    alltmp = []
                    if len(allinput)<Nmerge:
                        # Simple hadding all output files
                        logged_call(["hadd", "-f", "-v", output]+allinput, log, opt.run)
                    else:
                        # Two staged hadding:
                        # - First merge every Nmerge files to temp files
                        for n in range(1, (len(allinput)-1)/Nmerge+2):
                            tmplist = []
                            for i in range((n-1)*Nmerge, min(n*Nmerge,len(allinput))): tmplist.append(allinput[i])
                            tmpoutput = output.replace(".root","_"+str(n)+".root")
                            tmplog    = tmpoutput.rsplit("/",1)[0]+"/log/"+tmpoutput.rsplit("/",1)[1].replace(".root",".log")
                            alltmp.append(tmpoutput)
                            print "- Merging into temp file: "+tmpoutput
                            logged_call(["hadd", "-f", "-v", tmpoutput]+tmplist, tmplog, opt.run)
                        # - then merge the resulting temp files into a single file
                        #   and remove the temporary files
                        print "- Merging temp files into: "+output
                        logged_call(["hadd", "-f", "-v", output]+alltmp, log, opt.run)
                        for tmpfile in alltmp:
                            if os.path.isfile(tmpfile):
                                os.remove(tmpfile)
                            else:
                                print "Something went wrong with the hadding of tmp file: "+tmpfile
                                sys.exit()
                    #  Check that the result has the right size (if not, delete)
                    if os.path.isfile(output):
                        if os.path.getsize(output) < 1000: os.remove(output)

# Run all Analyzer jobs in parallel
def analysis(ana_arguments, last_known_status, nproc):
    global opt, currtime
    njob = len(ana_arguments)
    if not opt.batch and njob<nproc: nproc = njob
    output_files = []
    saved_path = os.getcwd()
    if opt.haddonly:
        for jobindex in range(0, njob):
            # Check if output file exists and size is larger than 1000 bytes
            output_file = ana_arguments[jobindex][0]
            file_size = 0
            if os.path.isfile(output_file):
                file_size = os.path.getsize(output_file)
                if file_size > 1000:
                    output_files.append(output_file)
                    last_known_status[jobindex] = 0
        merge_output(ana_arguments, last_known_status)
    elif not opt.batch:
        print "Running "+str(njob)+" instances of Analyzer jobs:"
        print
        # Use the N CPUs in parallel on the current computer to analyze all jobs
        workers = multiprocessing.Pool(processes=nproc)
        output_files = workers.map(analyzer_job, range(njob), chunksize=1)
        workers.close()
        workers.join()
        last_known_status = [0] * njob
        print "All Analyzer jobs finished."
        if not opt.nohadd:
            merge_output(ana_arguments, last_known_status)
    else:
        print "Running "+str(njob)+" instances of Batch jobs:"
        print
        # Running on the batch and babysitting all jobs until completion
        output_files = []
        njob = len(ana_arguments)
        if opt.run:
            finished = 0
            # Loop until all jobs are finished
            while finished != njob:
                # Renew token every ~10 hours
                if (time.time()-currtime)>36000:
                    special_call(["kinit", "-R"], 0)
                    currtime = time.time()
                if finished != 0: time.sleep(30)
                finished = 0
                for jobindex in range(0, njob):
                    output_file = ana_arguments[jobindex][0]
                    input_txtfile = ana_arguments[jobindex][1][0]
                    #output_log  = ana_arguments[jobindex][3]
                    #if jobindex==3: sys.exit()
                    if last_known_status[jobindex] == -1:
                        # Initial step (can be also a recovery task)
                        # First check if output file exists, has a size larger than 1000 bytes
                        # it can be opened, not corrupt and the counts histo
                        # returns the same number of events as the total of all input files
                        file_size = 0
                        if os.path.isfile(output_file):
                            file_size = os.path.getsize(output_file)
                            #print output_file+"="+str(file_size)
                        if file_size <= 1000 or not opt.recover:
                            # Submit all jobs
                            analyzer_job(jobindex)
                            last_known_status[jobindex] = time.time()
                        else:
                            # 
                            # Check if processed events exactly match input event counts
                            pass_nevent_check = True
                            if not opt.nocheck:
                                input_count = 0
                                output_count = 0
                                with open(input_txtfile) as txt:
                                    for infile in txt:
                                        infile = infile.replace('\n','')
                                        fin = ROOT.TFile.Open(infile)
                                        tree = fin.Get("B2GTree")
                                        input_count += tree.GetEntries()
                                        fin.Close()
                                if os.path.isfile(output_file):
                                    with suppress_stdout_stderr():
                                        fout = ROOT.TFile.Open(output_file)
                                        if fout and not fout.TestBit(ROOT.TFile.kRecovered):
                                            h_counts = fout.Get("counts")
                                            if h_counts: output_count += h_counts.GetBinContent(1)
                                        if fout: fout.Close()
                                pass_nevent_check = (input_count == output_count)
                                #print input_txtfile+"="+str(input_count)
                                #print output_file+"="+str(output_count)
                            if pass_nevent_check:
                                #print output_file+" - OK!"
                                finished += 1
                                output_files.append(output_file)
                                last_known_status[jobindex] = 0
                            else:
                                if os.path.isfile(output_file): os.remove(output_file)
                                # Resubmit job
                                print "Job["+str(jobindex)+"] "+output_file+" failed event check - input_count="+str(input_count)+" output_count="+str(output_count)
                                analyzer_job(jobindex)
                                last_known_status[jobindex] = time.time()
                    elif last_known_status[jobindex] == 0:
                        # Already completed
                        finished += 1
                    else:
                        # Repeated step
                        # First check if output file exists, has a size larger than 1000 bytes
                        # it can be opened, not corrupt and the counts histo
                        # returns the same number of events as the total of all input files
                        file_size = 0
                        if os.path.isfile(output_file): file_size = os.path.getsize(output_file)
                        if file_size > 1000:
                            # Check if processed events exactly match input event counts
                            pass_nevent_check = True
                            if not opt.nocheck:
                                input_count = 0
                                output_count = 0
                                with open(input_txtfile) as txt:
                                    for infile in txt:
                                        infile = infile.replace('\n','')
                                        fin = ROOT.TFile.Open(infile)
                                        tree = fin.Get("B2GTree")
                                        input_count += tree.GetEntries()
                                        fin.Close()
                                with suppress_stdout_stderr():
                                    fout = ROOT.TFile.Open(output_file)
                                    if fout and not fout.TestBit(ROOT.TFile.kRecovered):
                                        h_counts = fout.Get("counts")
                                        if h_counts: output_count += h_counts.GetBinContent(1)
                                    if fout: fout.Close()
                                pass_nevent_check = (input_count == output_count)
                            if pass_nevent_check:
                                finished += 1
                                output_files.append(output_file)
                                last_known_status[jobindex] = 0
                            else:
                                os.remove(output_file)
                                # Resubmit job
                                print "Job["+str(jobindex)+"] "+output_file+" failed event check - input_count="+str(input_count)+" output_count="+str(output_count)
                                analyzer_job(jobindex)
                                last_known_status[jobindex] = time.time()                           
                        # If the last submission/check is older than 10 minutes check job status with bjobs
                        elif time.time() - last_known_status[jobindex] > (600 if opt.NQUICK<2 else 600/opt.NQUICK):
                            jobname = SUBTIME+'_'+str(jobindex)
                            logged_call(shlex.split('bjobs -J '+jobname), TMPDIR+'jobstatus_'+jobname+'.txt', opt.run)
                            with open(TMPDIR+'jobstatus_'+jobname+'.txt') as jobstatus:
                                lines = jobstatus.readlines()
                                if 'Job <'+jobname+'> is not found' in lines[0]:
                                    # Job is missing, so resubmit
                                    print "Job["+str(jobindex)+"] "+jobname+" is missing, resubmitting ...         "
                                    analyzer_job(jobindex)
                                #else:
                                #    status = lines[1].split()[2]
                                #    if status == 'PEND' or status == 'RUN':
                                #        print "Job "+jobname+" is pending/running"
                            os.remove(TMPDIR+'jobstatus_'+jobname+'.txt')
                            last_known_status[jobindex] = time.time()
                # Merge output files
                if not opt.nohadd:
                    merge_output(ana_arguments, last_known_status)
                # Print status
                print "Analyzer jobs on batch (Done/All): "+str(finished)+"/"+str(njob)+"   \r",
                sys.stdout.flush()
            print "\nAll batch jobs finished."
        else:
            for jobindex in range(0, njob):
                output_files.append(analyzer_job(jobindex))
    print
    return output_files

# Run Plotter, output of Analyzer is input for this code
def plotter(input_files, output_file):
    global opt, EXEC_PATH
    print "Start plotting from output files"
    print
    special_call([EXEC_PATH+"/Plotter", output_file] + input_files, opt.run)
    print "Plotting finished."
    print

def show_result(plotter_out):
    print "Showing the result in root: "
    print
    special_call(["root", "-l", 'scripts/show_result.C("'+plotter_out+'")'], opt.run)

# ---------------------- Running -------------------------

# renew token before running
if opt.replot:
    if not opt.recover:
        backup_files(EXEC_PATH, SUBTIME)
        compile(0)
    plotter(PLOTTER_IN, PLOTTER_OUT)
    #if not 'lxplus' in socket.gethostname():
    #    show_result(PLOTTER_OUT)
else:
    if not opt.recover:
        if not opt.test:
            backup_files(EXEC_PATH, SUBTIME)
        compile(1, opt.plot)
    plotter_input_files = analysis(ana_arguments, last_known_status, opt.NPROC)
    if opt.plot:
        if not opt.nohadd: plotter_input_files = glob.glob(opt.OUTDIR+"/hadd/*.root")
        plotter(plotter_input_files, PLOTTER_OUT)
        #if not 'lxplus' in socket.gethostname():
        #    show_result(PLOTTER_OUT)

print "Done."
