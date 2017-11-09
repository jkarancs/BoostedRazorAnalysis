import re, os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, getpass, ROOT
from array import array
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
#parser.add_option('--nocombine',   dest="nocombine",   action="store_true", default=False,      help='Do not rerun combine, i.e. run on existing results (default=run combine)')
parser.add_option('--test',        dest="TEST",        type="int",          default=0,          help="Run only on a N signal points (default=0 - all)")
parser.add_option('--nproc',       dest="NPROC",       type="int",          default=6,          help="Tells how many parallel combine to start (Default=6)")
(opt,args) = parser.parse_args()

BIN = ""
if "_nj35" in opt.box: BIN = "_nj35"
if "_nj6"  in opt.box: BIN = "_nj6"

# ---------------------- Settings ------------------------

lumi = 35867 # /pb
ntuple = "ntuple/Latest"
combine_bins = False

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
elif opt.model == "T5qqqqVV":
    signal.append("FastSim_SMS-T5qqqqVV")
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
    "_topptUp",
    "_topptDown",
    "_isrUp", 
    "_isrDown", 
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

def load(f, name, pf="", combine = False):
    h = f.Get(name)
    if combine:
        h_new = combinebins(h, pf)
    else:
        h_new = ROOT.TH1D(h.GetName()+pf,h.GetTitle(),h.GetNbinsX(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
        for i in range (0, h.GetNbinsX()+2):
            h_new.SetBinContent(i,h.GetBinContent(i));
            h_new.SetBinError(i,h.GetBinError(i));
        h_new.SetEntries(h.GetEntries())
        h_new.SetDirectory(0)
    return h_new

def combinebins(h, pf=""):
    binmap = {
        1:   1,
        2:   2,
        3:   3,
        4:   4,
        5:   5,
        6:   6,
        7:   7,
        8:   8,
        9:   9,
        10: 10,
        11: 11,
        12: 12,
        13: 13,
        14: 14, # Merge 2 MR [1200, 1600], R2 [0.24, 0.5, 1.0]
        15: 14,
        16: 15,
        17: 16,
        18: 17, # Merge 3 MR [1600, 2000], R2 [0.16, 0.24, 0.5, 1.0]
        19: 17,
        20: 17,
        21: 18,
        22: 19, # Merge 4 MR [2000, 4000], R2 [0.12, 0.16, 0.24, 0.5, 1.0]
        23: 19,
        24: 19,
        25: 19 }
    h_new = ROOT.TH1D(h.GetName()+pf,h.GetTitle(), 19,0,19)
    for i in range (1, h.GetNbinsX()+1):
        h_new.SetBinContent(binmap[i], h_new.GetBinContent(binmap[i]) + h.GetBinContent(i))
        h_new.SetBinError  (binmap[i], (h_new.GetBinError(binmap[i])**2 + h.GetBinError(i)**2)**0.5)
    h_new.SetEntries(h.GetEntries())
    h_new.SetDirectory(0)
    return h_new

def loadclone(f, name, pf=""):
    h = f.Get(name)
    h_new = h.Clone(name+pf)
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
        #if mult2.GetBinContent(binx)<0:
        #    mult2.SetBinContent(binx,0)
        #    mult2.SetBinError(binx,0)
        #if div2.GetBinContent(binx)<0:
        #    div2.SetBinContent(binx,0)
        #    div2.SetBinError(binx,0)
    # bin-by-bin k factor
    #est.Multiply(mult2)
    #est.Divide(div2)
    # common k factor
    est.Scale(mult.Integral()/div.Integral())
    #if "Top_MJ_est" in name:
    #for binx in range(1, est.GetNbinsX()+1):
    #    if est.GetBinContent(binx)<0:
    #        print "2nd pass, negative bin content: "+name+" binx="+str(binx)+" cont="+str(est.GetBinContent(binx))
    #        print " - mult: "+str(mult.GetBinContent(binx))
    #        print " - div : "+str(div .GetBinContent(binx))
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

    # data
    if not os.path.exists("syst_"+opt.dir+"/hadd/data.root"):
        data_files = []
        for name in data: data_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/data.root"]+data_files,               "syst_"+opt.dir+"/hadd/log/data.log")

    # signal
    if not os.path.exists("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root"):
        signal_files = []
        for name in signal: signal_files.append(opt.dir+"/hadd/"+name+".root")
        if len(signal_files)>1:
            logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/signal_"+opt.model+".root"]+signal_files, "syst_"+opt.dir+"/hadd/log/signal_"+opt.model+".log")
        else:
            #print (["cp", "-p"]+signal_files+["syst_"+opt.dir+"/hadd/signal_"+opt.model+".root"])
            logged_call(["cp", "-p"]+signal_files+["syst_"+opt.dir+"/hadd/signal_"+opt.model+".root"], "syst_"+opt.dir+"/hadd/log/signal_"+opt.model+".log")
        # Get unskimmed number of events in low/high npv region (for signal acceptance)
        input_files = []
        for sig in signal: input_files += glob.glob(ntuple+"/"+sig+"/*.root")
        for i in range(len(input_files)):
            f = ROOT.TFile.Open(input_files[i])
            if "T2tt" in opt.model:
                h = f.Get("npvLowHigh_T2tt")
            else:
                h = f.Get("npvLowHigh_T1tttt")
            if i==0:
                npvLowHighHist_allevt = h.Clone(h.GetName()+"_allevt")
                npvLowHighHist_allevt.SetDirectory(0)
            else:
                npvLowHighHist_allevt.Add(h)
            f.Close()
        f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root","UPDATE")
        npvLowHighHist_allevt.Write()
        f.Close()

    # background
    if not os.path.exists("syst_"+opt.dir+"/hadd/multijet.root"):
        multijet_files = []
        for name in multijet: multijet_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/multijet.root"]+multijet_files,       "syst_"+opt.dir+"/hadd/log/multijet.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/top.root"):
        top_files = []
        for name in top: top_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/top.root"]+top_files,                 "syst_"+opt.dir+"/hadd/log/top.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/wjets.root"):
        wjets_files = []
        for name in wjets: wjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/wjets.root"]+wjets_files,             "syst_"+opt.dir+"/hadd/log/wjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/ztoinv.root"):
        ztoinv_files = []
        for name in ztoinv: ztoinv_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/ztoinv.root"]+ztoinv_files,           "syst_"+opt.dir+"/hadd/log/ztoinv.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/other.root"):
        other_files = []
        for name in other: other_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/other.root"]+other_files,             "syst_"+opt.dir+"/hadd/log/other.log")

# ----------------- Harvest histograms -------------------

# Load:
# BG estimate
# Q_data, Q_TT, Q_MJ, Q_WJ, Q_ZI, Q_OT
# W_data, W_TT, W_MJ, W_WJ, W_ZI, W_OT
# T_data, T_TT, T_MJ, T_WJ, T_ZI, T_OT
# S_data, S_TT, S_MJ, S_WJ, S_ZI, S_OT
# For kappa uncertainty:
# q_data, q_TT, q_MJ, q_WJ, q_ZI, q_OT

print "Loading histograms"

# Data
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/data.root")
Q_data = load(f,"MRR2_Q_data"+BIN,"_data", combine_bins)
W_data = load(f,"MRR2_W_data"+BIN,"_data", combine_bins)
T_data = load(f,"MRR2_T_data"+BIN,"_data", combine_bins)
S_data = load(f,"MRR2_S_data"+BIN,"_data", combine_bins)
q_data = load(f,"MRR2_q_data"+BIN,"_data", combine_bins)
npvHist = load(f,"nvtx","_data")
npvHist.Scale(1/npvHist.Integral())

# Signal
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root")
S_signal = []
counter = 0
for ikey in range(0, f.GetListOfKeys().GetEntries()):
    name = f.GetListOfKeys().At(ikey).GetName()
    if name.startswith("MRR2_S_signal") and not "Up" in name and not "Down" in name:
        if not "_nj35" in name and not "_nj6" in name:
            counter+=1
            S_syst = []
            for syst in systematics:
                S_syst.append(load(f, name+BIN+syst, "_sig", combine_bins))
            S_signal.append(S_syst)
    if opt.TEST>0:
        if counter==opt.TEST:
            break
# Histos for pileup acceptance systematic
if "T2tt" in opt.model:
    npvLowHighHist        = loadclone(f,"npvLowHigh_T2tt","_sig")
    npvLowHighHist_allevt = loadclone(f,"npvLowHigh_T2tt_allevt","_sig")
else:
    npvLowHighHist        = loadclone(f,"npvLowHigh_T1tttt","_sig")
    npvLowHighHist_allevt = loadclone(f,"npvLowHigh_T1tttt_allevt","_sig")

# Background
# top
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/top.root")
Q_TT = []
W_TT = []
T_TT = []
S_TT = []
q_TT = []
for syst in systematics:
    Q_TT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_TT", combine_bins))
    W_TT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_TT", combine_bins))
    T_TT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_TT", combine_bins))
    S_TT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_TT", combine_bins))
    q_TT.append(load(f,"MRR2_q_bkg"+BIN+syst,"_TT", combine_bins))
# multijet
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/multijet.root")
Q_MJ = []
W_MJ = []
T_MJ = []
S_MJ = []
q_MJ = []
for syst in systematics:
    Q_MJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_MJ", combine_bins))
    W_MJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_MJ", combine_bins))
    T_MJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_MJ", combine_bins))
    S_MJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_MJ", combine_bins))
    q_MJ.append(load(f,"MRR2_q_bkg"+BIN+syst,"_MJ", combine_bins))
# wjets
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/wjets.root")
Q_WJ = []
W_WJ = []
T_WJ = []
S_WJ = []
q_WJ = []
for syst in systematics:
    Q_WJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_WJ", combine_bins))
    W_WJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_WJ", combine_bins))
    T_WJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_WJ", combine_bins))
    S_WJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_WJ", combine_bins))
    q_WJ.append(load(f,"MRR2_q_bkg"+BIN+syst,"_WJ", combine_bins))
# ztoinv
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ztoinv.root")
Q_ZI = []
W_ZI = []
T_ZI = []
S_ZI = []
q_ZI = []
for syst in systematics:
    Q_ZI.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_ZI", combine_bins))
    W_ZI.append(load(f,"MRR2_W_bkg"+BIN+syst,"_ZI", combine_bins))
    T_ZI.append(load(f,"MRR2_T_bkg"+BIN+syst,"_ZI", combine_bins))
    S_ZI.append(load(f,"MRR2_S_bkg"+BIN+syst,"_ZI", combine_bins))
    q_ZI.append(load(f,"MRR2_q_bkg"+BIN+syst,"_ZI", combine_bins))
# other
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/other.root")
Q_OT = []
W_OT = []
T_OT = []
S_OT = []
q_OT = []
for syst in systematics:
    Q_OT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_OT", combine_bins))
    W_OT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_OT", combine_bins))
    T_OT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_OT", combine_bins))
    S_OT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_OT", combine_bins))
    q_OT.append(load(f,"MRR2_q_bkg"+BIN+syst,"_OT", combine_bins))

# ------------- Signal MET systematics ------------------

# Implementing method:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17#Special_treatment_of_MET_uncerta
# --> https://github.com/RazorCMS/RazorAnalyzer/blob/master/python/SMSTemplates.py#L140-L170

print "Correcting yields using gen-MET vs PF MET comparison"
i_metUp = -1
for i in range(len(S_signal[0])):
    if "metUp" in S_signal[0][i].GetName(): i_metUp = i

for i in range(len(S_signal)):
    # pfmet:  S_signal[i][0]       nominal
    # genmet: S_signal[i][i_metUp] "metUp"
    for ibin in range(1, S_signal[i][0].GetNbinsX()+1):
        pfmetYield = S_signal[i][0].GetBinContent(ibin)
        genmetYield = S_signal[i][i_metUp].GetBinContent(ibin)        
        central = (pfmetYield+genmetYield)/2.0
        unc = (pfmetYield-genmetYield)/2.0
        S_signal[i][i_metUp]  .SetBinContent(ibin, central + unc)
        S_signal[i][i_metUp+1].SetBinContent(ibin, central - unc)        
        if pfmetYield <= 0:
            continue
        #print "Multiplying bin content by",central/pfmetYield
        for hist in S_signal[i]:
            if '_metUp' in hist.GetName() or '_metDown' in hist.GetName(): 
                continue
            hist.SetBinContent(ibin, hist.GetBinContent(ibin) * central/pfmetYield)

# ------------- Signal pile-up systematic ---------------

# Implementing method:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsMoriond17#Pileup_lumi
# --> https://twiki.cern.ch/twiki/pub/CMS/SUSRecommendationsMoriond17/pileup_acceptance_systematic.pdf
# Similar implementation to: https://github.com/RazorCMS/RazorAnalyzer/blob/master/python/SMSTemplates.py#L60-L138

print "Performing NPV extrapolation procedure for signal MC"
i_pileupUp = -1
for i in range(len(S_signal[0])):
    if "pileupUp" in S_signal[0][i].GetName(): i_pileupUp = i

#get theory cross sections
xsecs = {}
for line in open('./data/stop13TeV.txt' if ('T2' in opt.model) else 'data/gluino13TeV.txt','r'):
    line = line.replace('\n','')
    xsecs[float(line.split(',')[0])]=float(line.split(',')[1]) #pb

w_pileup_nom = ROOT.TH2D("w_pileup_nominal","",npvLowHighHist_allevt.GetNbinsX(),npvLowHighHist_allevt.GetXaxis().GetXmin(),npvLowHighHist_allevt.GetXaxis().GetXmax(),
                         npvLowHighHist_allevt.GetNbinsY(),npvLowHighHist_allevt.GetYaxis().GetXmin(),npvLowHighHist_allevt.GetYaxis().GetXmax())
w_pileup_err = ROOT.TH2D("w_pileup_error","",npvLowHighHist_allevt.GetNbinsX(),npvLowHighHist_allevt.GetXaxis().GetXmin(),npvLowHighHist_allevt.GetXaxis().GetXmax(),
                         npvLowHighHist_allevt.GetNbinsY(),npvLowHighHist_allevt.GetYaxis().GetXmin(),npvLowHighHist_allevt.GetYaxis().GetXmax())

for i in range(len(S_signal)):
    histLowNPV  = S_signal[i][i_pileupUp+1]
    histHighNPV = S_signal[i][i_pileupUp]
    
    # Get gluino/stop and lsp bins for npvHighLow (signal)
    scan_point = S_signal[i][0].GetName()[:-4].replace("MRR2_S_signal_","").replace(BIN,"")
    mg_bin  = npvLowHighHist.GetXaxis().FindBin(float(scan_point.split("_")[0]))
    mch_bin = npvLowHighHist.GetXaxis().FindBin(float(scan_point.split("_")[1]))
    
    # Bin centers of low/high vertex distributions in data - TODO: Ask Dustin
    npvHist.GetXaxis().SetRangeUser(npvLowHighHist.GetZaxis().GetBinLowEdge(1),
                                    npvLowHighHist.GetZaxis().GetBinLowEdge(1)+npvLowHighHist.GetZaxis().GetBinWidth(1)-1)
    xmean_low = npvHist.GetMean()
    npvHist.GetXaxis().SetRangeUser(npvLowHighHist.GetZaxis().GetBinLowEdge(2),
                                    npvLowHighHist.GetZaxis().GetBinLowEdge(2)+npvLowHighHist.GetZaxis().GetBinWidth(2))
    xmean_high = npvHist.GetMean()
    x = array('d', [xmean_low, xmean_high]) # bin centers of npvLowHighHist
    ex = array('d', [0,0])
    pL = 0
    pH = 0
    epL = 0
    epH = 0
    for ibin in range(1, histLowNPV.GetNbinsX()+1):
        # Calculate relative acceptances
        pL += histLowNPV.GetBinContent(ibin)
        pH += histHighNPV.GetBinContent(ibin)
        epL = (epL*epL + histLowNPV.GetBinError(ibin)*histLowNPV.GetBinError(ibin)) ** 0.5
        epH = (epH*epH + histHighNPV.GetBinError(ibin)*histHighNPV.GetBinError(ibin)) ** 0.5
    
    p = pL+pH
    NL = npvLowHighHist_allevt.GetBinContent(mg_bin,mch_bin,1)
    NH = npvLowHighHist_allevt.GetBinContent(mg_bin,mch_bin,2)
    N = NL+NH
    evt_weight = (lumi * xsecs[float(scan_point.split("_")[0])] / N)
    if not epL: epL = 1.83 * evt_weight
    if not epH: epH = 1.83 * evt_weight
    #print "p (L,H): "+str(pL)+", "+str(pH)+"   N (L,H): "+str(NL)+", "+str(NH)
    relAccLow  = pL/(p if p>0 else 1)  * (N/NL)
    relAccHigh = pH/(p if p>0 else 1)  * (N/NH)
    errRelAccLow  = epL/(p if p>0 else 1)  * (N/NL)
    errRelAccHigh = epH/(p if p>0 else 1)  * (N/NH)
    #print "relAcc (low,high): "+str(relAccLow)+", "+str(relAccHigh)
    
    # Perform linear fit to 2-bin NPV distribution
    y = array('d', [relAccLow, relAccHigh])
    ey = array('d', [errRelAccLow, errRelAccHigh])
    graph = ROOT.TGraphErrors(2, x, y, ex, ey)
    fitResult = graph.Fit("pol1", "SMFQ")
    p0 = fitResult.Parameter(0)
    p1 = fitResult.Parameter(1)
    fitCov = fitResult.GetCovarianceMatrix()
    
    # Convolve linear fit with NPV distribution in data
    averageAcceptance = 0
    averageAcceptanceErr = 0
    for npvBin in range(1, npvHist.GetNbinsX()+1):
        npv = npvHist.GetXaxis().GetBinCenter(npvBin)
        npvWeight = npvHist.GetBinContent(npvBin)
        fitPred = p0 + npv * p1
        fitError = ( npv*npv * fitCov(1,1) + 
                     fitCov(0,0) + 2*npv * fitCov(0,1) )**(0.5)
        averageAcceptance += npvWeight * fitPred
        averageAcceptanceErr += npvWeight * fitError
    averageAcceptance = max(0.01, averageAcceptance)
    w_pileup_nom.SetBinContent(mg_bin, mch_bin, averageAcceptance)
    w_pileup_err.SetBinContent(mg_bin, mch_bin, averageAcceptanceErr)
    
    # Put the up/down errors from this procedure 
    # into the npvextrap histograms
    for ibin in range(1, histLowNPV.GetNbinsX()+1):
        nominal = S_signal[i][0].GetBinContent(ibin)
        S_signal[i][i_pileupUp]  .SetBinContent(ibin, (averageAcceptance+averageAcceptanceErr)*nominal)
        S_signal[i][i_pileupUp+1].SetBinContent(ibin, max(0.01,averageAcceptance-averageAcceptanceErr)*nominal)
    
    # Scale all the rest of the signal histograms based on these weights
    #print "pL,pH: "+str(pL)+", "+str(pH)+" acceptance factors ",
    #print "(nominal, up, down): "+str(averageAcceptance)+",",
    #print " "+str(averageAcceptance+averageAcceptanceErr)+",",
    #print " "+str(max(0,averageAcceptance-averageAcceptanceErr))
    weightedOverNominal = averageAcceptance
    #print "In bin %d, weighted yield is %.3f of nominal"%(ibin,weightedOverNominal)
    for hist in S_signal[i]:
        if '_pileupUp' in hist.GetName() or '_pileupDown' in hist.GetName(): 
            continue
        for ibin in range(1, histLowNPV.GetNbinsX()+1):
            hist.SetBinContent(ibin, hist.GetBinContent(ibin) * weightedOverNominal)

f_pu = ROOT.TFile.Open("signal_pileup_syst_"+opt.model+"_"+opt.box+".root","RECREATE")
w_pileup_nom.Write()
w_pileup_err.Write()
f_pu.Close()

# --------------- kappa uncertainty ---------------------

##  # Determined from the Closure of Q --> Q'
##  
##  # Merging MC yields in Q and Q'
##  Q_MC = Q_TT[0].Clone("Q_MC")
##  Q_MC.Add(Q_MJ[0])
##  Q_MC.Add(Q_WJ[0])
##  Q_MC.Add(Q_ZI[0])
##  Q_MC.Add(Q_OT[0])
##  q_MC = q_TT[0].Clone("q_MC")
##  q_MC.Add(q_MJ[0])
##  q_MC.Add(q_WJ[0])
##  q_MC.Add(q_ZI[0])
##  q_MC.Add(q_OT[0])
##  
##  # Calculating bin-by-bin Q' estimate
##  q_est = Q_data.Clone("q_est")
##  q_est.Multiply(q_MC)
##  q_est.Divide(Q_MC)
##  
##  # Divide data counts by this estimate
##  kappa_unc = q_data.Clone("kappa_unc")
##  kappa_unc.Divide(q_est)
##  
##  # Convolute gaussians with area of event counts in data
##  functions = []
##  fmax = 0
##  all_func = ""
##  for xbin in range(1, kappa_unc.GetNbinsX()+1):
##      #print "Drawing gaus with parameters: norm="+str(q_data.GetBinContent(xbin))+", mean="+str(kappa_unc.GetBinContent(xbin))+", sigma="+str(kappa_unc.GetBinError(xbin))
##      gaus = ROOT.TF1("gaus"+str(xbin),"[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))", -10,10)
##      #gaus.SetParameters(q_data.GetBinContent(xbin)**0.5,kappa_unc.GetBinContent(xbin),kappa_unc.GetBinError(xbin))
##      gaus.SetParameters(q_data.GetBinContent(xbin),kappa_unc.GetBinContent(xbin),kappa_unc.GetBinError(xbin))
##      all_func += "gaus("+str((xbin-1)*3)+")+"
##      gaus.SetLineColor(605+xbin)
##      functions.append(gaus)
##      if gaus.GetMaximum() > fmax: fmax = gaus.GetMaximum()
##  
##  kappa_uncertainty = 0.6
##  plotmax = 500
##  plot = ROOT.TCanvas("kappa_unc_"+opt.box)
##  hist = ROOT.TH1D("h",";Data/Prediction;A.U.",20,-10,10)
##  hist.GetXaxis().SetRangeUser(0,2)
##  hist.GetYaxis().SetRangeUser(0,fmax*1.1)
##  hist.SetStats(0)
##  hist.Draw()
##  conv = ROOT.TF1("conv", all_func[:-1], -10,10)
##  conv.SetLineColor(1)
##  for i in range(0, len(functions)):
##      functions[i].Draw("SAME")
##      conv.SetParameter(i*3,  functions[i].GetParameter(0))
##      conv.SetParameter(i*3+1,functions[i].GetParameter(1))
##      conv.SetParameter(i*3+2,functions[i].GetParameter(2))
##  conv.Draw("SAME")
##  plotmax = 1.2 * conv.GetMaximum()
##  hist.GetYaxis().SetRangeUser(0,plotmax)
##  l1 = ROOT.TLine(1-kappa_uncertainty,0, 1-kappa_uncertainty,plotmax)
##  l1.Draw()
##  l2 = ROOT.TLine(1+kappa_uncertainty,0, 1+kappa_uncertainty,plotmax)
##  l2.Draw()
##  lat = ROOT.TLatex(0.6, 0.9*plotmax, "Integral in range = %.3f" % (conv.Integral(1-kappa_uncertainty, 1+kappa_uncertainty)/conv.Integral(-10, 10)))
##  lat.Draw()
##  
##  plot.SaveAs(plot.GetName()+".png")
##  
##  fout = ROOT.TFile("kappa_unc_"+opt.box+".root","RECREATE")
##  Q_MC.Write()
##  q_MC.Write()
##  q_est.Write()
##  Q_data.Write("Q_data")
##  q_data.Write("q_data")
##  kappa_unc.Write()
##  plot.Write()
##  fout.Close()

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
T_MJ_est = bg_est("Top_MJ_est",         Q_data, [Q_TT[0],            Q_WJ[0], Q_ZI[0], Q_OT[0]], T_MJ[0], Q_MJ[0])
for i in range(0, len(systematics)):
    S_TT_est = bg_est("Top"         +systematics[i], T_data, [          T_MJ_est, T_WJ[0], T_ZI[0], T_OT[0]], S_TT[i], T_TT[i])
    S_MJ_est = bg_est("MultiJet"    +systematics[i], Q_data, [Q_TT[0],            Q_WJ[0], Q_ZI[0], Q_OT[0]], S_MJ[i], Q_MJ[i])
    S_WJ_est = bg_est("WJets"       +systematics[i], W_data, [W_TT[0],  W_MJ[0],           W_ZI[0], W_OT[0]], S_WJ[i], W_WJ[i])
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
    #root_filename = "syst_"+opt.dir+"/cards/RazorBoost_"+opt.box+"_"+opt.model+"_"+scan_point+".root"
    root_filename = "syst_"+opt.dir+"/cards/RazorBoost_SMS-"+opt.model+"_"+scan_point+"_"+opt.box+".root"
    if not opt.nocards:
        fout = ROOT.TFile.Open(root_filename,"recreate")
        #print "  Creating root file: "+root_filename
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
        #print "  Creating data card: "+card_filename
        card=open(card_filename, 'w+')
        card.write(
'''imax 1 number of channels
jmax 5 number of backgrounds
kmax * number of nuisance parameters
------------------------------------------------------------
observation	'''
            )
        card.write(str(S_data.Integral()))
        card.write(
'''
------------------------------------------------------------
shapes * * '''
            )
        card.write(root_filename)
        card.write(
''' $PROCESS $PROCESS_$SYSTEMATIC
------------------------------------------------------------
bin		'''
            )
        card.write("%s\t%s\t%s\t%s\t%s\t%s" % (opt.box, opt.box, opt.box, opt.box, opt.box, opt.box))
        card.write(
'''
process		Signal	Top	MultiJet	WJets	ZInv	Other
process		0	1	2	3	4	5
rate		'''
            )
        card.write("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (signal_syst[0].Integral(), Top_est[0].Integral(), MultiJet_est[0].Integral(), WJets_est[0].Integral(), ZInv_est[0].Integral(), Other_est[0].Integral()) )
        card.write(
'''
------------------------------------------------------------
lumi		lnN	1.025	1.025	1.025	1.025	1.025	1.025
pileup		shape	1.0	1.0	1.0	1.0	1.0	1.0
toppt		shape	-	1.0	1.0	1.0	1.0	1.0
isr		shape	1.0	-	-	-	-	-
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

print "All data cards ready"

#results = []
#combine_cmds = []
#for card in cards:
#    combine_out_filename = card.replace("cards/RazorBoost","combine/RazorBoost").replace(".txt",".log")
#    results.append(combine_out_filename)
#    if  not opt.nocombine:
#        combine_cmds.append((["combine", "-M", "AsymptoticLimits", "-d", card], combine_out_filename))
#
#run_combine(combine_cmds, opt.NPROC)

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
#print "Done."
