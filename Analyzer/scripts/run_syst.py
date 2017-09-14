import re, os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, getpass, ROOT
from optparse import OptionParser

# ---------------------- Cmd Line  -----------------------

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option('-d','--dir',    dest="dir",         type="string",       default="",         help="Input/output directory (use output of Analyzer)")
parser.add_option('-b','--box',    dest="box",         type="string",       default="WAna_nj6", help='Analysis box, eg. TopAna (default="WAna_nj6")')
parser.add_option('-m','--model',  dest="model",       type="string",       default="T5ttcc",   help='Signal model (default="T5ttcc")')
parser.add_option('--nohadd',      dest="nohadd",      action="store_true", default=False,      help='Do not merge input files (default=merge them)')
parser.add_option('--nocards',     dest="nocards",     action="store_true", default=False,      help='Do not create data cards, i.e. run on existing ones (default=create them)')
parser.add_option('--nocombine',   dest="nocombine",   action="store_true", default=False,      help='Do not rerun combine, i.e. run on existing results (default=run combine)')
parser.add_option('--test',        dest="TEST",        type="int",          default=0,          help="Run only on a N signal points (default=0 - all)")
parser.add_option('--nproc',       dest="NPROC",       type="int",          default=6,          help="Tells how many parallel combine to start (Default=6)")
(opt,args) = parser.parse_args()

BIN = ""
if "_nj35" in opt.box: BIN = "_nj35"
if "_nj6"  in opt.box: BIN = "_nj6"

# ---------------------- Settings ------------------------

# List of input files
data = [
    "JetHT_Run2016B_03Feb2017_v2",
    "JetHT_Run2016C_03Feb2017",
    "JetHT_Run2016D_03Feb2017",
    "JetHT_Run2016E_03Feb2017",
    "JetHT_Run2016F_03Feb2017",
    "JetHT_Run2016G_03Feb2017",
    "JetHT_Run2016H_03Feb2017_v2",
    "JetHT_Run2016H_03Feb2017_v3"
]

signal = []
if opt.model == "T1tttt":
    signal.append("FastSim_SMS-T1tttt")
elif opt.model == "T1ttbb":
    signal.append("FastSim_SMS-T1ttbb")
elif opt.model == "T2tt":
    signal.append("FastSim_SMS-T2tt_mStop-150to250")
    signal.append("FastSim_SMS-T2tt_mStop-250to350")
    signal.append("FastSim_SMS-T2tt_mStop-350to400")
    signal.append("FastSim_SMS-T2tt_mStop-400to1200")
elif opt.model == "T5tttt":
    signal.append("FastSim_SMS-T5tttt")
elif opt.model == "T5ttcc":
    signal.append("FastSim_SMS-T5ttcc")
    signal.append("FastSim_SMS-T5ttcc_mGluino1750to2300")
else:
    "Error: unknown signal model: "+opt.model
    sys.exit()

top = [
    # single top
    "ST_s-channel_4f_InclusiveDecays",
    "ST_t-channel_antitop_4f_inclusiveDecays",
    "ST_t-channel_top_4f_inclusiveDecays",
    "ST_tW_antitop_5f_inclusiveDecays",
    "ST_tW_top_5f_inclusiveDecays",
    # ttbar
    "TT_powheg-pythia8"
]
multijet = [
    "QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
    # V --> QQ
    "DYJetsToQQ_HT180",
    "WJetsToQQ_HT180",
    "ZJetsToQQ_HT600toInf",
    # VV --> QQQQ
    "WWTo4Q",
    "ZZTo4Q"
]
wjets = [
    "WJetsToLNu_Wpt-0To50",
    "WJetsToLNu_Wpt-100to200",
    "WJetsToLNu_Wpt-200toInf",
    "WJetsToLNu_Wpt-50To100"
]
ztoinv = [
    "ZJetsToNuNu_HT-1200To2500",
    "ZJetsToNuNu_HT-200To400",
    "ZJetsToNuNu_HT-2500ToInf",
    "ZJetsToNuNu_HT-400To600",
    "ZJetsToNuNu_HT-600To800",
    "ZJetsToNuNu_HT-800To1200"
]
other = [
    # g*/Z --> ll
    "DYJetsToLL_M-5to50_HT-100to200",
    "DYJetsToLL_M-5to50_HT-200to400",
    "DYJetsToLL_M-5to50_HT-400to600",
    "DYJetsToLL_M-5to50_HT-600toInf",
    "DYJetsToLL_M-50_HT-200to400",
    "DYJetsToLL_M-50_HT-2500toInf",
    "DYJetsToLL_M-50_HT-400to600",
    "DYJetsToLL_M-50_HT-600to800",
    "DYJetsToLL_M-50_HT-800to1200",
    "DYJetsToLL_M-50_HT-1200to2500",
    # gamma
    "GJets_HT-40To100",
    "GJets_HT-100To200",
    "GJets_HT-200To400",
    "GJets_HT-400To600",
    "GJets_HT-600ToInf",
    # ttbar + X
    "TTGJets",
    "TTWJetsToLNu",
    "TTWJetsToQQ",
    "TTZToLLNuNu",
    "TTZToQQ",
    "TTTT",
    # VV
    "WWTo2L2Nu",
    "WWToLNuQQ",
    "WZTo1L1Nu2Q",
    "WZTo1L3Nu",
    "WZTo2L2Q",
    "WZTo2Q2Nu",
    "WZTo3LNu",
    "ZZTo2L2Nu",
    "ZZTo2L2Q",
    "ZZTo2Q2Nu",
    "ZZTo4L",
    # VVV
    "WWW",
    "WWZ",
    "WZZ",
    "ZZZ"
]

# List of systematics to consider
systematics = [
    "",
#    "_lumiUp",
#    "_lumiDown",
    "_pileupUp",
    "_pileupDown",
    "_alphasUp",
    "_alphasDown",
    "_facscaleUp",
    "_facscaleDown",
    "_renscaleUp",
    "_renscaleDown",
    "_facrenscaleUp", 
    "_facrenscaleDown", 
    "_triggerUp",
    "_triggerDown",
    "_jesUp",
    "_jesDown",
    "_jerUp",
    "_jerDown",
    "_metUp", 
    "_metDown", 
    "_elerecoUp",
    "_elerecoDown",
    "_eleidUp",
    "_eleidDown",
    "_eleisoUp",
    "_eleisoDown",
    "_elefastsimUp",
    "_elefastsimDown",
    "_muontrkUp",
    "_muontrkDown",
    "_muonidisoUp",
    "_muonidisoDown",
    "_muonfastsimUp",
    "_muonfastsimDown",
    "_btagUp",
    "_btagDown",
    "_btagfastsimUp",
    "_btagfastsimDown",
    "_wtagUp",
    "_wtagDown",
    "_wtagfastsimUp",
    "_wtagfastsimDown",
    "_toptagUp",
    "_toptagDown",
    "_toptagfastsimUp",
    "_toptagfastsimDown",
]

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

def load(f, name, pf=""):
    h = f.Get(name)
    h_new = ROOT.TH1D(h.GetName()+pf,h.GetTitle(),h.GetNbinsX(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
    for i in range (0, h.GetNbinsX()+2):
        h_new.SetBinContent(i,h.GetBinContent(i));
        h_new.SetBinError(i,h.GetBinError(i));
    h_new.SetEntries(h.GetEntries())
    h_new.SetDirectory(0)
    return h_new

def bg_est(name, data, sub, mult, div):
    est   = data.Clone(name)
    mult2 = mult.Clone(name+"_mult")
    div2  = div .Clone(name+"_div")
    for hist in sub:
        est.Add(hist, -1)
    # Zero bins with negative counts
    # Coming from NLO generator scale variations (flipped sign)
    # See thread: https://hypernews.cern.ch/HyperNews/CMS/get/generators/3510.html
    # TODO: Should we zero the event weight during filling of a histogram if original event weight was positive?
    for binx in range(1, est.GetNbinsX()+1):
        if est.GetBinContent(binx)<0:
            est.SetBinContent(binx,0)
            est.SetBinError(binx,0)
        if mult2.GetBinContent(binx)<0:
            mult2.SetBinContent(binx,0)
            mult2.SetBinError(binx,0)
        if div2.GetBinContent(binx)<0:
            div2.SetBinContent(binx,0)
            div2.SetBinError(binx,0)
    est.Multiply(mult2)
    est.Divide(div2)
    #if "Top_MJ_est" in name:
    for binx in range(1, est.GetNbinsX()+1):
        if est.GetBinContent(binx)<0:
            print "2nd pass, negative bin content: "+name+" binx="+str(binx)+" cont="+str(est.GetBinContent(binx))
            print " - mult: "+str(mult.GetBinContent(binx))
            print " - div : "+str(div .GetBinContent(binx))
    return est

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

# ----------------- Merge histograms --------------------

if not opt.nohadd:
    print "Merging histograms from directory: "+opt.dir
    
    data_files = []
    for name in data: data_files.append(opt.dir+"/hadd/"+name+".root")
    
    signal_files = []
    for name in signal: signal_files.append(opt.dir+"/hadd/"+name+".root")
    
    data_files = []
    for name in data: data_files.append(opt.dir+"/hadd/"+name+".root")
    signal_files = []
    for name in signal: signal_files.append(opt.dir+"/hadd/"+name+".root")
    multijet_files = []
    for name in multijet: multijet_files.append(opt.dir+"/hadd/"+name+".root")
    top_files = []
    for name in top: top_files.append(opt.dir+"/hadd/"+name+".root")
    wjets_files = []
    for name in wjets: wjets_files.append(opt.dir+"/hadd/"+name+".root")
    ztoinv_files = []
    for name in ztoinv: ztoinv_files.append(opt.dir+"/hadd/"+name+".root")
    other_files = []
    for name in other: other_files.append(opt.dir+"/hadd/"+name+".root")
    
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/data.root"]+data_files,         "syst_"+opt.dir+"/hadd/log/data.log")
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/signal.root"]+signal_files,     "syst_"+opt.dir+"/hadd/log/signal.log")
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/multijet.root"]+multijet_files, "syst_"+opt.dir+"/hadd/log/multijet.log")
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/top.root"]+top_files,           "syst_"+opt.dir+"/hadd/log/top.log")
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/wjets.root"]+wjets_files,       "syst_"+opt.dir+"/hadd/log/wjets.log")
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/ztoinv.root"]+ztoinv_files,     "syst_"+opt.dir+"/hadd/log/ztoinv.log")
    logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/other.root"]+other_files,       "syst_"+opt.dir+"/hadd/log/other.log")

# ----------------- Harvest histograms -------------------

# Load:
# Q_data, Q_TT, Q_MJ, Q_WJ, Q_ZI, Q_OT
# W_data, W_TT, W_MJ, W_WJ, W_ZI, W_OT
# T_data, T_TT, T_MJ, T_WJ, T_ZI, T_OT
# S_data, S_TT, S_MJ, S_WJ, S_ZI, S_OT

print "Loading histograms"

# Data
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/data.root")
Q_data = load(f,"MRR2_Q_data"+BIN,"_data")
W_data = load(f,"MRR2_W_data"+BIN,"_data")
T_data = load(f,"MRR2_T_data"+BIN,"_data")
S_data = load(f,"MRR2_S_data"+BIN,"_data")

# Signal
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal.root")
S_signal = []
counter = 0
for ikey in range(0, f.GetListOfKeys().GetEntries()):
    name = f.GetListOfKeys().At(ikey).GetName()
    if name.startswith("MRR2_S_signal") and not "Up" in name and not "Down" in name:
        if not "_nj35" in name and not "_nj6" in name:
            counter+=1
            S_syst = []
            for syst in systematics:
                S_syst.append(load(f, name+BIN+syst, "_sig"))
            S_signal.append(S_syst)
    if opt.TEST>0:
        if counter==opt.TEST:
            break

# Background
# top
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/top.root")
Q_TT = []
W_TT = []
T_TT = []
S_TT = []
for syst in systematics:
    Q_TT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_TT"))
    W_TT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_TT"))
    T_TT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_TT"))
    S_TT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_TT"))
# multijet
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/multijet.root")
Q_MJ = []
W_MJ = []
T_MJ = []
S_MJ = []
for syst in systematics:
    Q_MJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_MJ"))
    W_MJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_MJ"))
    T_MJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_MJ"))
    S_MJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_MJ"))
# wjets
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/wjets.root")
Q_WJ = []
W_WJ = []
T_WJ = []
S_WJ = []
for syst in systematics:
    Q_WJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_WJ"))
    W_WJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_WJ"))
    T_WJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_WJ"))
    S_WJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_WJ"))
# ztoinv
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ztoinv.root")
Q_ZI = []
W_ZI = []
T_ZI = []
S_ZI = []
for syst in systematics:
    Q_ZI.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_ZI"))
    W_ZI.append(load(f,"MRR2_W_bkg"+BIN+syst,"_ZI"))
    T_ZI.append(load(f,"MRR2_T_bkg"+BIN+syst,"_ZI"))
    S_ZI.append(load(f,"MRR2_S_bkg"+BIN+syst,"_ZI"))
# other
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/other.root")
Q_OT = []
W_OT = []
T_OT = []
S_OT = []
for syst in systematics:
    Q_OT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_OT"))
    W_OT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_OT"))
    T_OT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_OT"))
    S_OT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_OT"))


# --------------- Background Estimation ------------------

# Formulas for bkg estimate:
# S_MJ_est = (Q_data - Q_notMJ) * S_MJ / Q_MJ
# T_MJ_est = (Q_data - Q_notMJ) * T_MJ / Q_MJ
# S_WJ_est = (W_data - W_notWJ) * S_WJ / W_WJ
# S_TT_est = (T_data - T_MJ_est - T_notTTorMJ) * S_TT / T_TT

# Estimate:
#T_MJ_est = (Q_data - Q_TT     - Q_WJ - Q_ZI - Q_OT) * T_MJ / Q_MJ
#S_MJ_est = (Q_data - Q_TT     - Q_WJ - Q_ZI - Q_OT) * S_MJ / Q_MJ
#S_WJ_est = (W_data - W_TT     - W_MJ - W_ZI - W_OT) * S_WJ / W_WJ
#S_TT_est = (T_data - T_MJ_est - T_WJ - T_ZI - T_OT) * S_TT / T_TT
#S_ZI_est = S_ZI
#S_OT_est = S_OT
#S_est    = S_MJ_est + S_WJ_est + S_TT_est + S_ZI_est + S_OT_est

# Caluclating background estimates
print "Calculating background estimates"

Top_est      = []
MultiJet_est = []
WJets_est    = []
ZInv_est     = []
Other_est    = []
for i in range(0, len(systematics)):
    T_MJ_est = bg_est("Top_MJ_est",                  Q_data, [Q_TT[i],            Q_WJ[i], Q_ZI[i], Q_OT[i]], T_MJ[i], Q_MJ[i])
    S_TT_est = bg_est("Top"         +systematics[i], T_data, [          T_MJ_est, T_WJ[i], T_ZI[i], T_OT[i]], S_TT[i], T_TT[i])
    S_MJ_est = bg_est("MultiJet"    +systematics[i], Q_data, [Q_TT[i],            Q_WJ[i], Q_ZI[i], Q_OT[i]], S_MJ[i], Q_MJ[i])
    S_WJ_est = bg_est("WJets"       +systematics[i], W_data, [W_TT[i],  W_MJ[i],           W_ZI[i], W_OT[i]], S_WJ[i], W_WJ[i])
    S_ZI_est = S_ZI[i].Clone("ZInv" +systematics[i])
    S_OT_est = S_OT[i].Clone("Other"+systematics[i])
    Top_est     .append(S_TT_est)
    MultiJet_est.append(S_MJ_est)
    WJets_est   .append(S_WJ_est)
    ZInv_est    .append(S_ZI_est)
    Other_est   .append(S_OT_est)

# Now save a different root file for each signal point
if opt.nocards:
    print "Reusing previous data cards"
else:
    print "Looping on Signal points and creating data cards"
    if not os.path.exists("syst_"+opt.dir+"/cards"):
        special_call(["mkdir", "-p", "syst_"+opt.dir+"/cards"], 0)

cards = []
for signal_syst in S_signal:
    scan_point = signal_syst[0].GetName()[:-4].replace("MRR2_S_signal_","").replace(BIN,"")
    root_filename = "syst_"+opt.dir+"/cards/RazorBoost_"+opt.box+"_"+opt.model+"_"+scan_point+".root"
    if not opt.nocards:
        fout = ROOT.TFile.Open(root_filename,"recreate")
        print "  Creating root file: "+root_filename
        # Add signal systematics
        for i in range(0, len(systematics)):
            signal_syst[i].Write("Signal"+systematics[i])
        # Add background estimates
        for bkg in [Top_est, MultiJet_est, WJets_est, ZInv_est, Other_est]:
            for syst_var in bkg:
                syst_var.Write()
        # Add data counts
        S_data.Write("data_obs")
    card_filename = root_filename.replace(".root",".txt")
    cards.append(card_filename)
    if not opt.nocards:
        print "  Creating data card: "+card_filename
        card=open(card_filename, 'w+')
        card.write(
    '''    imax 1 number of channels
    jmax 5 number of backgrounds
    kmax * number of nuisance parameters
    ------------------------------------------------------------
    observation	''')
        card.write(str(S_data.Integral()))
        card.write(
    '''
    ------------------------------------------------------------
    shapes * * ''')
        card.write(root_filename)
        card.write(
    ''' $PROCESS $PROCESS_$SYSTEMATIC
    ------------------------------------------------------------
    bin		''')
        card.write("%s\t%s\t%s\t%s\t%s\t%s" % (opt.box, opt.box, opt.box, opt.box, opt.box, opt.box))
        card.write(
    '''
    process		Signal	Top	MultiJet	WJets	ZInv	Other
    process		0	1	2	3	4	5
    rate		''')
        card.write("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (signal_syst[0].Integral(), Top_est[0].Integral(), MultiJet_est[0].Integral(), WJets_est[0].Integral(), ZInv_est[0].Integral(), Other_est[0].Integral()) )
        card.write(
    '''
    ------------------------------------------------------------
    lumi		lnN	1.026	1.026	1.026	1.026	1.026	1.026
    pileup		shape	-	1.0	1.0	1.0	1.0	1.0
    jes		shape	1.0	1.0	1.0	1.0	1.0	1.0
    jer 		shape	1.0	1.0	1.0	1.0	1.0	1.0
    met		shape	1.0	1.0	1.0	1.0	1.0	1.0
    trigger		shape	1.0	1.0	1.0	1.0	1.0	1.0
    facscale	shape	1.0	1.0	1.0	1.0	1.0	1.0
    renscale	shape	1.0	1.0	1.0	1.0	1.0	1.0
    facrenscale	shape	1.0	1.0	1.0	1.0	1.0	1.0
    alphas		shape	1.0	1.0	1.0	1.0	1.0	1.0
    elereco		shape	1.0	1.0	1.0	1.0	1.0	1.0
    eleid		shape	1.0	1.0	1.0	1.0	1.0	1.0
    eleiso		shape	1.0	1.0	1.0	1.0	1.0	1.0
    elefastsim	shape	1.0	-	-	-	-	-
    muontrk		shape	1.0	1.0	1.0	1.0	1.0	1.0
    muonidiso	shape	1.0	1.0	1.0	1.0	1.0	1.0
    muonfastsim	shape	1.0	-	-	-	-	-
    btag		shape	1.0	1.0	1.0	1.0	1.0	1.0
    btagfastsim	shape	1.0	-	-	-	-	-
    wtag		shape	1.0	1.0	1.0	1.0	1.0	1.0
    wtagfastsim	shape	1.0	-	-	-	-	-
    toptag		shape	1.0	1.0	1.0	1.0	1.0	1.0
    toptagfastsim	shape	1.0	-	-	-	-	-
    '''
        )
        card.close()

time.sleep(10)

if opt.nocombine:
    print "Reusing previous combine results"
else:
    print "All data cards ready, running combine"
    
results = []
combine_cmds = []
for card in cards:
    combine_out_filename = card.replace("cards/RazorBoost","combine/RazorBoost").replace(".txt",".log")
    results.append(combine_out_filename)
    if  not opt.nocombine:
        combine_cmds.append((["combine", "-M", "AsymptoticLimits", "-d", card], combine_out_filename))

run_combine(combine_cmds, opt.NPROC)

##    print "Creating summary plots"
##    
##    
##    if not os.path.exists("syst_"+opt.dir+"/results"):
##        special_call(["mkdir", "-p", "syst_"+opt.dir+"/results"], 0)
##    fout = ROOT.TFile.Open("syst_"+opt.dir+"/results/limits.root","recreate")
##    # limit histo(s)
##    obs_limit_T5ttcc       = ROOT.TH2D("obs_limit_T5ttcc",      "T5ttcc;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Observed Limit",  201,-12.5,5012.5, 201,-12.5,5012.5)
##    exp_limit_2Down_T5ttcc = ROOT.TH2D("exp_limit_2Down_T5ttcc","T5ttcc;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Observed Limit",  201,-12.5,5012.5, 201,-12.5,5012.5)
##    exp_limit_1Down_T5ttcc = ROOT.TH2D("exp_limit_1Down_T5ttcc","T5ttcc;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Observed Limit",  201,-12.5,5012.5, 201,-12.5,5012.5)
##    exp_limit_T5ttcc       = ROOT.TH2D("exp_limit_T5ttcc",      "T5ttcc;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Observed Limit",  201,-12.5,5012.5, 201,-12.5,5012.5)
##    exp_limit_1Up_T5ttcc   = ROOT.TH2D("exp_limit_1Up_T5ttcc",  "T5ttcc;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Observed Limit",  201,-12.5,5012.5, 201,-12.5,5012.5)
##    exp_limit_2Up_T5ttcc   = ROOT.TH2D("exp_limit_2Up_T5ttcc",  "T5ttcc;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Observed Limit",  201,-12.5,5012.5, 201,-12.5,5012.5)
##    
##    for result in results:
##        with open(result) as log:
##            mMother = result[:-4].split("_")[-2:][0]
##            mLSP    = result[:-4].split("_")[-2:][1]
##            for line in log:
##                if "Observed Limit:" in line:
##                    obs_limit = line.split()[4]
##                    print "Limit for "+mMother+", "+mLSP+": "+obs_limit
##                    obs_limit_T5ttcc.Fill(float(mMother),float(mLSP),float(obs_limit))
##                elif "Expected  2.5%:" in line:
##                    exp_limit_2Down = line.split()[4]
##                    exp_limit_2Down_T5ttcc.Fill(float(mMother),float(mLSP),float(exp_limit_2Down))
##                elif "Expected 16.0%:" in line:
##                    exp_limit_1Down = line.split()[4]
##                    exp_limit_1Down_T5ttcc.Fill(float(mMother),float(mLSP),float(exp_limit_1Down))
##                elif "Expected 50.0%:" in line:
##                    exp_limit = line.split()[4]
##                    exp_limit_T5ttcc.Fill(float(mMother),float(mLSP),float(exp_limit))
##                elif "Expected 84.0%:" in line:
##                    exp_limit_1Up = line.split()[4]
##                    exp_limit_1Up_T5ttcc.Fill(float(mMother),float(mLSP),float(exp_limit_1Up))
##                elif "Expected 97.5%:" in line:
##                    exp_limit_2Up = line.split()[4]
##                    exp_limit_2Up_T5ttcc.Fill(float(mMother),float(mLSP),float(exp_limit_2Up))
##    
##    obs_limit_T5ttcc.Write()
##    exp_limit_2Down_T5ttcc.Write()
##    exp_limit_1Down_T5ttcc.Write()
##    exp_limit_T5ttcc.Write()
##    exp_limit_1Up_T5ttcc.Write()
##    exp_limit_2Up_T5ttcc.Write()
##    
##    time.sleep(120)
print "Done."
