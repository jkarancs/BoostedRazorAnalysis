import re, os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, getpass, ROOT
from optparse import OptionParser


# ---------------------- Cmd Line  -----------------------

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option('-d','--dir',    dest="dir",         type="string",       default="",         help="Input/output directory (use output of Analyzer)")
parser.add_option('-m','--model',  dest="model",       type="string",       default="T5ttcc",   help='Signal model (default="T5ttcc")')
parser.add_option('--nproc',       dest="NPROC",       type="int",          default=6,          help="Tells how many parallel combine to start (Default=6)")
(opt,args) = parser.parse_args()

boxes = args

if len(boxes)<2:
    print "Error: give to additional arguments for the boxes to merge (eg. python scripts/combine_cards.py -d <inputdir> -m <model> WAna_nj35 WAna_nj6"
    sys.exit()

# --------------------- Functions ------------------------
# Show and run command with stdout on screen
icommand=0
def special_call(cmd, verbose=1):
    global icommand, opt
    if verbose:
        print("[%d]" % icommand),
        for i in xrange(len(cmd)): print cmd[i],
        print ""
    ntry = 0
    while True:
        try:
            if subprocess.call(cmd):
                print "ERROR: Problem executing command:"
                print("[%d]" % icommand)
                for i in xrange(len(cmd)): print cmd[i],
                print ""
                print "exiting."
                sys.exit()
        except:
            print "Could not excecute command: "
            print("[%d]" % icommand)
            for i in xrange(len(cmd)): print cmd[i],
            print ""
            print "Wait 10s and continue"
            time.sleep(10)
            ntry += 1
            if ntry == 20: sys.exit()
            continue
        break
    if verbose: print ""
    sys.stdout.flush()
    icommand+=1

# Run command with stdout/stderr saved to logfile
def logged_call(cmd, logfile, append=False):
    global opt
    dirname = os.path.dirname(logfile)
    if dirname != "" and not os.path.exists(dirname):
        special_call(["mkdir", "-p", os.path.dirname(logfile)], 0)
    ntry = 0
    while True:
        try:
            with open(logfile, "a" if append else "w") as log:
                proc = subprocess.Popen(cmd, stdout=log, stderr=log, close_fds=True)
                proc.wait()
        except:
            print "Could not write to disk (IOError), wait 10s and continue"
            time.sleep(10)
            ntry += 1
            if ntry == 20: sys.exit()
            continue
        break

# Run a single Combine instance (on a single input list, i.e. one dataset)
combine_cmds = []
def combine_job((jobindex)):
    global combine_cmds, opt
    cmd = combine_cmds[jobindex][0]
    log = combine_cmds[jobindex][1]
    print "Starting Combine, expected output: "+log
    if not os.path.exists(os.path.dirname(log)):
        special_call(["mkdir", "-p", os.path.dirname(log)], 0)
    #if opt.NPROC>3: cmd = ['nice']+cmd
    logged_call(cmd, log)
    return log

# Run all Combine jobs in parallel
def run_combine(combine_cmds, nproc):
    global opt
    njob = len(combine_cmds)
    if njob<nproc: nproc = njob
    print "Running "+str(njob)+" instances of Combine jobs:"
    print
    # Use the N CPUs in parallel on the current computer to run combine for all jobs
    workers = multiprocessing.Pool(processes=nproc)
    output_files = workers.map(combine_job, range(njob), chunksize=1)
    workers.close()
    workers.join()
    print "All Combine jobs finished."
    print
    return output_files

# --------------- Run the Combination --------------------

print "Looping on Signal points and merging data cards"
cards = []
for card in glob.glob(opt.dir+"/cards/RazorBoost_"+boxes[0]+"_"+opt.model+"_*.txt"):
    mg   = card.replace(".txt","").split("_")[-2:][0]
    mchi = card.replace(".txt","").split("_")[-2:][1]
    cards_to_merge = []
    for box in boxes:
        cards_to_merge.append(opt.dir+"/cards/RazorBoost_"+box+"_"+opt.model+"_"+mg+"_"+mchi+".txt")
    cmd = ["combineCards.py"]
    for i in range(len(boxes)):
        cmd.append(boxes[i]+"="+opt.dir+"/cards/RazorBoost_"+boxes[i]+"_"+opt.model+"_"+mg+"_"+mchi+".txt")
    combined_card = opt.dir+"/cards/RazorBoost_combined_"+opt.model+"_"+mg+"_"+mchi+".txt"
    logged_call(cmd, combined_card)
    cards.append(combined_card)


time.sleep(10)

print "All combined data cards ("+str(len(cards))+") ready, running combine"

results = []
combine_cmds = []
for card in cards:
    combine_out_filename = card.replace("cards/RazorBoost","combine/RazorBoost").replace(".txt",".log")
    results.append(combine_out_filename)
    combine_cmds.append((["combine", "-M", "AsymptoticLimits", "-d", card], combine_out_filename))

run_combine(combine_cmds, opt.NPROC)

print "Done."

