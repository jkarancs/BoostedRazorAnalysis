import re, os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, getpass, ROOT, io
from array import array
from optparse import OptionParser
from common_functions import *
import tdrstyle

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
if   "_nj35" in opt.box: BIN = "_nj35"
elif "_nj45" in opt.box: BIN = "_nj45"
elif "_nj46" in opt.box: BIN = "_nj46"
elif "_nj6"  in opt.box: BIN = "_nj6"
elif "_nj7"  in opt.box: BIN = "_nj7"
DATE = "_".join(opt.dir.split("/")[-1].split("_")[1:4])

# ---------------------- Settings ------------------------


lumi = 35867 # /pb
ntuple = "ntuple/Latest"
combine_bins = True
binned_k = False

mrbins = [ 800, 1000, 1200, 1600, 2000, 4000 ]
r2bins = [ 0.08, 0.12, 0.16, 0.24, 0.4, 2.0 ]
nrazorbin = (len(mrbins)-1)*(len(r2bins)-1)

tdrstyle.setTDRStyle()

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
elif opt.model == "T1ttbb_dM5to25":
    signal.append("FastSim_SMS-T1ttbb_deltaM5to25")
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
    #"ST_tW_antitop_5f_inclusiveDecays",
    #"ST_tW_top_5f_inclusiveDecays",
    "ST_tW_antitop_5f_NoFullyHadronicDecays",
    "ST_tW_top_5f_NoFullyHadronicDecays"
]
ttbar = [
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

gjets = [
    # gamma
    "GJets_HT-40To100",
    "GJets_HT-100To200",
    "GJets_HT-200To400",
    "GJets_HT-400To600",
    "GJets_HT-600ToInf",
]

dyjets = [
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
    "DYJetsToLL_M-50_HT-1200to2500"
]

rest = [
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

nonwjets  = top + ttbar + multijet         + ztoinv + gjets + dyjets + rest
nongjets  = top + ttbar + multijet + wjets + ztoinv         + dyjets + rest
nondyjets = top + ttbar + multijet + wjets + ztoinv + gjets          + rest
other = rest + gjets + dyjets

bkg = ttbar + top + wjets + multijet + ztoinv + other

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
    "_ak8scaleUp",
    "_ak8scaleDown",
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
    "_wmistagUp",
    "_wmistagDown",
    "_wmistagfastsimUp",
    "_wmistagfastsimDown",
    "_wmasstagUp",
    "_wmasstagDown",
    "_wantitagUp",
    "_wantitagDown",
    "_toptagUp",
    "_toptagDown",
    "_toptagfastsimUp",
    "_toptagfastsimDown",
    "_topmistagUp",
    "_topmistagDown",
    "_topmistagfastsimUp",
    "_topmistagfastsimDown",
    "_top0bmasstagUp",
    "_top0bmasstagDown",
    "_topmasstagUp",
    "_topmasstagDown",
    "_topantitagUp",
    "_topantitagDown",
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
    if dirname != "" and not os.path.exists(dirname): os.makedirs(dirname)
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
        h_new = combinebins(h, h.GetName()+pf)
    else:
        h_new = ROOT.TH1D(h.GetName()+pf,h.GetTitle(),h.GetNbinsX(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
        for i in range (0, h.GetNbinsX()+2):
            h_new.SetBinContent(i,h.GetBinContent(i));
            h_new.SetBinError(i,h.GetBinError(i));
        h_new.SetEntries(h.GetEntries())
        h_new.SetDirectory(0)
    return h_new

def combinebins(h, name=""):
    h_new = ROOT.TH1D(name,h.GetTitle(), 22,0,22)
    h_new.SetDirectory(0)
    h_new.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
    h_new.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
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
        14: 14,
        15: 15,
        16: 16,
        17: 17,
        18: 18,
        19: 19, # Merge 2 MR [1600, 2000], R2 [0.24, 0.5, 1.0]
        20: 19,
        21: 20,
        22: 21,
        23: 22, # Merge 3 MR [2000, 4000], R2 [0.16, 0.24, 0.5, 1.0]
        24: 22,
        25: 22 }
    for i in range (1, h.GetNbinsX()+1):
        h_new.SetBinContent(binmap[i], h_new.GetBinContent(binmap[i]) + h.GetBinContent(i))
        h_new.SetBinError  (binmap[i], (h_new.GetBinError(binmap[i])**2 + h.GetBinError(i)**2)**0.5)
    h_new.SetEntries(h.GetEntries())
    add_bin_labels(h_new, 1, mrbins, r2bins)
    return h_new

def loadclone(f, name, pf=""):
    h = f.Get(name)
    h_new = h.Clone(name+pf)
    h_new.SetDirectory(0)
    return h_new

def calc_factorized_kappa(sr_2d, cr_2d, name, title):
    nbin = sr_2d.GetNbinsX()*sr_2d.GetNbinsY()
    h_new = ROOT.TH1D(name, ";;"+title+" transfer factor", nbin,0,nbin)
    # Calculate the transfer factor in bins of MR
    # Simply ratio of projections along x (MR)
    h_kappa_MR = sr_2d.ProjectionX()
    h_kappa_MR.Divide(cr_2d.ProjectionX())
    # And also in bins of R2 (and divide by average, in order to factorize)
    # projection along y (R2)
    # divided by the average factor
    h_kappa_R2 = sr_2d.ProjectionY()
    h_kappa_R2.Divide(cr_2d.ProjectionY())
    h_kappa_R2.Scale(cr_2d.Integral()/sr_2d.Integral())
    # Multiply factors and their errors in unrolled format
    unrolled_bin = 0
    for binx in range(1, sr_2d.GetNbinsX()+1):
        for biny in range(1, sr_2d.GetNbinsY()+1):
            unrolled_bin += 1
            tf1     = h_kappa_MR.GetBinContent(binx)
            tf1_err = h_kappa_MR.GetBinError  (binx)
            tf2     = h_kappa_R2.GetBinContent(biny)
            tf2_err = h_kappa_R2.GetBinError  (biny)
            tf     = tf1*tf2
            tf_err = (tf1_err*tf1_err*tf2*tf2 + tf2_err*tf2_err*tf1*tf1) ** 0.5
            h_new.SetBinContent(unrolled_bin, tf)
            h_new.SetBinError  (unrolled_bin, tf_err)
    return h_new

def bg_est(name, data, sub, sr, cr, combine_bins, write=False, binned_k = True):
    if not ("Up" in name or "Down" in name) and write:
        sr.Write()
        cr.Write()
    # Subtract MC from data counts
    if combine_bins and binned_k:
        est = combinebins(data, name)
    else:
        if combine_bins:
            est   = data.Clone(name+"_nocomb")
        else:
            est   = data.Clone(name)            
    for hist in sub:
        if combine_bins and binned_k:
            est.Add(combinebins(hist, hist.GetName()+"_combined"), -1)
        else:
            est.Add(hist, -1)
    if not ("Up" in name or "Down" in name) and write: est.Write(name+"_nonzeroed")
    # Zero bins with negative counts
    # Coming from NLO generator scale variations (flipped sign)
    # See thread: https://hypernews.cern.ch/HyperNews/CMS/get/generators/3510.html
    # TODO: Should we zero the event weight during filling of a histogram if original event weight was positive?
    for binx in range(1, est.GetNbinsX()+1):
        if est.GetBinContent(binx)<0:
            est.SetBinContent(binx,0)
            est.SetBinError(binx,0)
    if not ("Up" in name or "Down" in name) and write: est.Write(name+"_zeroed")
    # Calculate the transfer factors
    # Binned version
    processname = name.split("_")[0]
    srname = sr.GetName().split("_")[1].replace("s","S'").replace("q","Q'")
    crname = cr.GetName().split("_")[1]
    h_sr = sr.Clone()
    h_cr = cr.Clone()
    if combine_bins:
        h_sr = combinebins(h_sr, sr.GetName()+"_combined")
        h_cr = combinebins(h_cr, cr.GetName()+"_combined")
    nbin = h_sr.GetNbinsX()
    h_kappa_binned        = ROOT.TH1D("kappa_binned_"+srname+"_"+crname, ";;"+srname+"/"+crname+" transfer factor", nbin,0,nbin)
    add_bin_labels(h_kappa_binned, combine_bins)
    h_kappa_binned.Divide(h_sr, h_cr)
    nbin_nocomb = sr.GetNbinsX()
    h_kappa_binned_nocomb = ROOT.TH1D("kappa_binned_nocomb_"+srname+"_"+crname, ";;"+srname+"/"+crname+" transfer factor", nbin_nocomb,0,nbin_nocomb)
    add_bin_labels(h_kappa_binned_nocomb, 0)
    h_kappa_binned_nocomb.Divide(sr, cr)
    # Zero also negative kappas (for binned version only)
    for binx in range(1, h_sr.GetNbinsX()+1):
        if h_kappa_binned.GetBinContent(binx)<0:
            h_kappa_binned.SetBinContent(binx,0)
            h_kappa_binned.SetBinError(binx,0)
    for binx in range(1, sr.GetNbinsX()+1):
        if h_kappa_binned_nocomb.GetBinContent(binx)<0:
            h_kappa_binned_nocomb.SetBinContent(binx,0)
            h_kappa_binned_nocomb.SetBinError(binx,0)
    # Factorized form
    h_kappa_fact = calc_factorized_kappa(rollback(sr), rollback(cr), "kappa_factorized_"+srname+"_"+crname, srname+"/"+crname)
    if not ("Up" in name or "Down" in name) and write:
        h_kappa_binned.Write()
        h_kappa_binned_nocomb.Write()
        h_kappa_fact.Write()
    # Final estimate
    # with the kappa of choice
    if binned_k:
        est.Multiply(h_kappa_binned)
    else:
        est.Multiply(h_kappa_fact)
    # For the factorized kappa version
    # Merge only the bins after the estimate
    if (not binned_k) and combine_bins:
        est = combinebins(est, name)
    # Draw kappa comparison plot
    if not ("Up" in name or "Down" in name):
        if write: est.Write(name+"_est")
        can_tf = custom_can(h_kappa_binned, "kappa_"+processname+"_"+srname+"_"+crname+"_"+opt.box, 0,0)
        #can_tf.SetLogy(1)
        h_kappa_binned_nocomb.SetMarkerColor(633)
        h_kappa_binned_nocomb.SetLineColor  (633)
        #h_kappa_binned_nocomb.GetYaxis().SetRangeUser(1e-2,1e3)
        h_kappa_binned_nocomb.Draw("PE1")
        h_kappa_fact.Draw("SAME PE1")
        leg = ROOT.TLegend(0.3,0.7,0.6,0.85,"")
        leg.AddEntry(h_kappa_fact,   "#color[1]{Factorized form}", "pl")
        leg.AddEntry(h_kappa_binned, "#color[633]{Binned}",        "pl")
        leg.SetTextSize(0.04)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw("SAME")
        save_plot(can_tf, "", "Plots/kappa/"+DATE+"/"+can_tf.GetName())
    ##print  if name == "Top":
    ##print      for binx in range(1, est.GetNbinsX()+1):
    ##print          print ("Top: bin=%-2d, est=%5.2f" % (binx, est.GetBinContent(binx)))
    ##print  if name == "Top":
    ##print      for binx in range(1, est.GetNbinsX()+1):
    ##print          print ("Top: bin=%-2d, sr=%5.2f +- %5.2f" % (binx, sr2.GetBinContent(binx), sr2.GetBinError(binx)))
    ##print  if name == "Top":
    ##print      for binx in range(1, est.GetNbinsX()+1):
    ##print          print ("Top: bin=%-2d, cr=%5.2f" % (binx, cr2.GetBinContent(binx)))
    ##print  if name == "Top":
    ##print      for binx in range(1, est.GetNbinsX()+1):
    ##print          print ("Top: bin=%-2d, trans=%5.2f" % (binx, sr2.GetBinContent(binx)/max(0.001,cr2.GetBinContent(binx))))
    ##print      print ("Top: avg trans=%5.2f" % (sr2.Integral()/max(0.001,cr2.Integral())))
    # common k factor
    ##common  est.Scale(sr.Integral()/cr.Integral())
    #if "Top_MJ_est" in name:
    #for binx in range(1, est.GetNbinsX()+1):
    #    if est.GetBinContent(binx)<0:
    #        print "2nd pass, negative bin content: "+name+" binx="+str(binx)+" cont="+str(est.GetBinContent(binx))
    #        print " - sr: "+str(sr.GetBinContent(binx))
    #        print " - cr : "+str(cr .GetBinContent(binx))
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

# Silence stdout/stderr
class suppress_stdout_stderr(object):
    def __init__(self):
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        self.save_fds = [os.dup(1), os.dup(2)]
    def __enter__(self):
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)
    def __exit__(self, *_):
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

def save_plot(can, name, plotname, write=True):
    # Check if the directory exists (if not create it first)
    dirname = os.path.dirname(plotname)
    if dirname != "" and not os.path.exists(dirname):
        os.makedirs(dirname)
    with suppress_stdout_stderr():
        can.SaveAs(plotname+".png")
    if write:
        if name != "":
            can.Write(name)
        else:
            can.Write()

def fit_fraction(chiso, temp1, temp2, plotdir, bin):
    temp1.Write()
    temp2.Write()
    EBEE = ("EB" if "EB" in chiso.GetName() else "EE")
    data    = chiso.Clone("data_" +EBEE+"_"+bin)
    can = ROOT.TCanvas("CHiso_"+EBEE+"_fit_"+bin)
    data.Draw("hist")
    prompt_val = ROOT.Double(0)
    prompt_err = ROOT.Double(0)
    fake_val = ROOT.Double(0)
    fake_err = ROOT.Double(0)
    if data.Integral()>0.:
        mctemp1 = temp1.Clone("temp1_"+EBEE+"_"+bin)
        mctemp2 = temp2.Clone("temp2_"+EBEE+"_"+bin)
        templates = ROOT.TObjArray(2)
        templates.Add(mctemp1)
        templates.Add(mctemp2)
        fitter = ROOT.TFractionFitter(data,templates)
        fitter.Constrain(0, 0,1)
        fitter.Constrain(1, 0,1)
        with suppress_stdout_stderr(): status = fitter.Fit()
        result = fitter.GetPlot();
        fitter.GetResult(0, prompt_val, prompt_err)
        fitter.GetResult(1, fake_val, fake_err)
        #print "prompt: "+str(prompt_val)+" fake_value: "+str(fake_val)
        result.SetLineColor(2)
        result.Draw("same")
        leg = ROOT.TLegend(0.28,0.70,0.9,0.9)
        leg.SetTextSize(0.03)
        leg.SetHeader("Purity : "+str("%4.3f" % prompt_val)+" in "+bin+", "+EBEE)
        leg.AddEntry(data, "Data", "LPE")
        leg.AddEntry(result, "Total Fit", "L")
        leg.Draw("SAME")
        save_plot(can, "", plotdir+"/Fit_"+bin+"_"+EBEE+"_"+opt.box)
    return prompt_val, prompt_err

def get_zslice(input, name, binx1, binx2, biny1, biny2):
    out = ROOT.TH1D(name,";"+input.GetXaxis().GetTitle(), input.GetNbinsZ(),
                    input.GetZaxis().GetXmin(),input.GetZaxis().GetXmax())
    for binz in range(1, input.GetNbinsZ()+1):
        for binx in range(binx1, binx2+1):
            for biny in range(biny1, biny2+1):
                out.SetBinContent(binz, out.GetBinContent(binz)+input.GetBinContent(binx,biny,binz))
                out.SetBinError  (binz, (out.GetBinError(binz)**2+input.GetBinError(binx,biny,binz)**2)**0.5)
    return out

def calc_syst_err(vh, name):
    nominal = vh[0].Clone(name)
    syst_up = vh[0].Clone(name+"_up")
    syst_dn = vh[0].Clone(name+"_dn")
    for binx in range(1, nominal.GetNbinsX()+1):
        err_up = nominal.GetBinError(binx) ** 2
        err_dn = nominal.GetBinError(binx) ** 2
        for i in range(1, len(vh)):
            syst_err = vh[i].GetBinContent(binx)-nominal.GetBinContent(binx)
            if syst_err>0:
                err_up += syst_err ** 2
            else:
                err_dn += syst_err ** 2
        err_up = err_up ** 0.5
        err_dn = err_dn ** 0.5
        syst_up.SetBinContent(binx, err_up)
        syst_dn.SetBinContent(binx, err_dn)
    return [nominal, syst_up, syst_dn]

def unroll(h):
    nbin = h.GetNbinsX()*h.GetNbinsY()
    h_new = ROOT.TH1D(h.GetName()+"_unrolled", "", nbin,0,nbin)
    h_new.SetDirectory(0)
    for binx in range(1, h.GetNbinsX()+1):
        for biny in range(1, h.GetNbinsY()+1):
            unrolled_bin = (binx-1)*h.GetNbinsY()+biny
            h_new.SetBinContent(unrolled_bin, h.GetBinContent(binx, biny))
            h_new.SetBinError  (unrolled_bin, h.GetBinError  (binx, biny))
    return h_new

def rollback(h):
    h_new = ROOT.TH2D(h.GetName()+"_2D", ";M_{R} (GeV);R^{2}", len(mrbins)-1,0,len(mrbins)-1, len(r2bins)-1,0,len(r2bins)-1)
    h_new.SetDirectory(0)
    for i in range(h.GetNbinsX()):
        bin_mr = (i/(len(r2bins)-1))+1
        bin_r2 = (i%(len(mrbins)-1))+1
        h_new.SetBinContent(bin_mr, bin_r2, h.GetBinContent(i+1))
        h_new.SetBinError  (bin_mr, bin_r2, h.GetBinError  (i+1))
    return h_new


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
        # TODO: Remove this when new ntuple is ready
        ##  f = ROOT.TFile.Open("syst_results/run_2017_11_13_syst/hadd/signal_"+opt.model+".root")
        ##  if "T2tt" in opt.model: npvLowHighHist_allevt = f.Get("npvLowHigh_T2tt_allevt")
        ##  else:                   npvLowHighHist_allevt = f.Get("npvLowHigh_T1tttt_allevt")
        ##  npvLowHighHist_allevt = npvLowHighHist_allevt.Clone()
        ##  npvLowHighHist_allevt.SetDirectory(0)
        ##  f.Close()
        f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root","UPDATE")
        npvLowHighHist_allevt.Write()
        f.Close()

    # background
    if not os.path.exists("syst_"+opt.dir+"/hadd/multijet.root"):
        multijet_files = []
        for name in multijet: multijet_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/multijet.root"]+multijet_files,       "syst_"+opt.dir+"/hadd/log/multijet.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/ttbar.root"):
        ttbar_files = []
        for name in ttbar: ttbar_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/ttbar.root"]+ttbar_files,             "syst_"+opt.dir+"/hadd/log/ttbar.log")
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
    
    # all bkg summed
    if not os.path.exists("syst_"+opt.dir+"/hadd/bkg.root"):
        bkg_files = []
        for name in bkg: bkg_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/bkg.root"]+bkg_files,                 "syst_"+opt.dir+"/hadd/log/bkg.log")

    # other (needed for Z(nunu)
    if not os.path.exists("syst_"+opt.dir+"/hadd/nonwjets.root"):
        nonwjets_files = []
        for name in nonwjets: nonwjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/nonwjets.root"]+nonwjets_files,       "syst_"+opt.dir+"/hadd/log/nonwjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/gjets.root"):
        gjets_files = []
        for name in gjets: gjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/gjets.root"]+gjets_files,             "syst_"+opt.dir+"/hadd/log/gjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/nongjets.root"):
        nongjets_files = []
        for name in nongjets: nongjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/nongjets.root"]+nongjets_files,       "syst_"+opt.dir+"/hadd/log/nongjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/dyjets.root"):
        dyjets_files = []
        for name in dyjets: dyjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/dyjets.root"]+dyjets_files,           "syst_"+opt.dir+"/hadd/log/dyjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/nondyjets.root"):
        nondyjets_files = []
        for name in nondyjets: nondyjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/nondyjets.root"]+nondyjets_files,     "syst_"+opt.dir+"/hadd/log/nondyjets.log")

#sys.exit()

# ----------------- Harvest histograms -------------------

# Load:
# BG estimate
# Q_data, Q_TT, Q_MJ, Q_WJ, Q_ZI, Q_OT
# W_data, W_TT, W_MJ, W_WJ, W_ZI, W_OT
# T_data, T_TT, T_MJ, T_WJ, T_ZI, T_OT
# S_data, S_TT, S_MJ, S_WJ, S_ZI, S_OT
# For kappa uncertainty:
# q_data, q_TT, q_MJ, q_WJ, q_ZI, q_OT
# For closure tests
# L_data, L_TT, L_MJ, L_WJ, L_ZI, L_OT

print "Loading histograms"

# Data
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/data.root")
Q_data = load(f,"MRR2_Q_data"+BIN,"_data", 0)
W_data = load(f,"MRR2_W_data"+BIN,"_data", 0)
L_data = load(f,"MRR2_L_data"+BIN,"_data", 0)
T_data = load(f,"MRR2_T_data"+BIN,"_data", 0)
S_data = load(f,"MRR2_S_data"+BIN,"_data", combine_bins)
s_data = load(f,"MRR2_s_data"+BIN,"_data", combine_bins)
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
        if not "_nj" in name:
            counter+=1
            S_syst = []
            for syst in systematics:
                S_syst.append(load(f, name+BIN+syst, "_sig", 0))
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
# Merge statistics in Mglu/Mstop vs Mlsp
npvLowHighHist       .Rebin3D(4,4,1)
npvLowHighHist_allevt.Rebin3D(4,4,1)

# Background
# top + ttbar
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ttbar.root")
Q_TT = []
W_TT = []
L_TT = []
T_TT = []
S_TT = []
s_TT = []
q_TT = []
for syst in systematics:
    Q_TT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_TT", 0))
    W_TT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_TT", 0))
    L_TT.append(load(f,"MRR2_L_bkg"+BIN+syst,"_TT", 0))
    T_TT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_TT", 0))
    S_TT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_TT", 0))
    s_TT.append(load(f,"MRR2_s_bkg"+BIN+syst,"_TT", 0))
    q_TT.append(load(f,"MRR2_q_bkg"+BIN+syst,"_TT", 0))
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/top.root")
for i in range(len(systematics)):
    # Fix problem with nonexistent scale weights for single top
    syst = systematics[i]
    if "scale" in syst: syst = ""
    Q_TT[i].Add(load(f,"MRR2_Q_bkg"+BIN+syst,"_T", 0))
    W_TT[i].Add(load(f,"MRR2_W_bkg"+BIN+syst,"_T", 0))
    L_TT[i].Add(load(f,"MRR2_L_bkg"+BIN+syst,"_T", 0))
    T_TT[i].Add(load(f,"MRR2_T_bkg"+BIN+syst,"_T", 0))
    S_TT[i].Add(load(f,"MRR2_S_bkg"+BIN+syst,"_T", 0))
    s_TT[i].Add(load(f,"MRR2_s_bkg"+BIN+syst,"_T", 0))
    q_TT[i].Add(load(f,"MRR2_q_bkg"+BIN+syst,"_T", 0))
# multijet
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/multijet.root")
Q_MJ = []
W_MJ = []
L_MJ = []
T_MJ = []
S_MJ = []
s_MJ = []
q_MJ = []
for syst in systematics:
    Q_MJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_MJ", 0))
    W_MJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_MJ", 0))
    L_MJ.append(load(f,"MRR2_L_bkg"+BIN+syst,"_MJ", 0))
    T_MJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_MJ", 0))
    S_MJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_MJ", 0))
    s_MJ.append(load(f,"MRR2_s_bkg"+BIN+syst,"_MJ", 0))
    q_MJ.append(load(f,"MRR2_q_bkg"+BIN+syst,"_MJ", 0))
# wjets
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/wjets.root")
Q_WJ = []
W_WJ = []
L_WJ = []
T_WJ = []
S_WJ = []
L_WJ = []
s_WJ = []
q_WJ = []
for syst in systematics:
    Q_WJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_WJ", 0))
    W_WJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_WJ", 0))
    L_WJ.append(load(f,"MRR2_L_bkg"+BIN+syst,"_WJ", 0))
    T_WJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_WJ", 0))
    S_WJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_WJ", 0))
    s_WJ.append(load(f,"MRR2_s_bkg"+BIN+syst,"_WJ", 0))
    q_WJ.append(load(f,"MRR2_q_bkg"+BIN+syst,"_WJ", 0))
# ztoinv
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ztoinv.root")
Q_ZI = []
W_ZI = []
L_ZI = []
T_ZI = []
S_ZI = []
s_ZI = []
q_ZI = []
for syst in systematics:
    Q_ZI.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_ZI", 0))
    W_ZI.append(load(f,"MRR2_W_bkg"+BIN+syst,"_ZI", 0))
    L_ZI.append(load(f,"MRR2_L_bkg"+BIN+syst,"_ZI", 0))
    T_ZI.append(load(f,"MRR2_T_bkg"+BIN+syst,"_ZI", 0))
    S_ZI.append(load(f,"MRR2_S_bkg"+BIN+syst,"_ZI", 0))
    s_ZI.append(load(f,"MRR2_s_bkg"+BIN+syst,"_ZI", 0))
    q_ZI.append(load(f,"MRR2_q_bkg"+BIN+syst,"_ZI", 0))
# other
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/other.root")
Q_OT = []
W_OT = []
L_OT = []
T_OT = []
S_OT = []
s_OT = []
q_OT = []
for syst in systematics:
    Q_OT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_OT", 0))
    W_OT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_OT", 0))
    L_OT.append(load(f,"MRR2_L_bkg"+BIN+syst,"_OT", 0))
    T_OT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_OT", 0))
    S_OT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_OT", 0))
    s_OT.append(load(f,"MRR2_s_bkg"+BIN+syst,"_OT", 0))
    q_OT.append(load(f,"MRR2_q_bkg"+BIN+syst,"_OT", 0))

# ---------------- Z(nunu) estimate ---------------------

purity_in_Gm1 = False

# Loading plots
#ZInv_dir = "syst_"+opt.dir
ZInv_dir1 = "syst_"+(opt.dir.replace("2018_01_08_syst","2018_01_11_CHIsoNoEACorr"))
ZInv_dir2 = "syst_"+(opt.dir.replace("2018_01_08_syst","2018_01_03_syst"))


# Needed for 1-lepton estimate
##  f = ROOT.TFile.Open(ZInv_dir1+"/hadd/data.root")
##  L_data_2D   = loadclone(f,"MR_R2_L"+BIN,"_data")
##  f = ROOT.TFile.Open(ZInv_dir1+"/hadd/wjets.root")
##  L_WJ_2D     = loadclone(f,"MR_R2_L"+BIN,"_WJ")
##  f = ROOT.TFile.Open(ZInv_dir1+"/hadd/nonwjets.root")
##  L_NONWJ  = loadclone(f,"MR_R2_L"+BIN,"_NONWJ")

# Z(nunu) estimate - using SR/CR(G) transfer factors
#f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/bkg.root")
f = ROOT.TFile.Open(ZInv_dir2+"/hadd/bkg.root")
IsDirect_G_EB = loadclone(f, "MR_R2_IsDirect_G_EB", "_MC")
IsDirect_G_EE = loadclone(f, "MR_R2_IsDirect_G_EE", "_MC")
GDirectPrompt_2D = loadclone(f, "MR_R2_G_DirectPrompt", "_MC")
Z_MC      = loadclone(f, "MR_R2_Z",    "_MC")
G_MC      = loadclone(f, "MR_R2_G",    "_MC")
#G_MC = []
#for syst in systematics:
#    G_MC_SYST.append(rollback(load(f, "MRR2_G_bkg"+syst,"_TMP", 0), GDirectPrompt_2D, "MRR2_G"+syst))

#f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/gjets.root")
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/gjets.root")
CHIsoTemplate_Prompt_EB = loadclone(f, "CHIsoTemplate_Prompt_g_EB", "_MC")
CHIsoTemplate_Prompt_EE = loadclone(f, "CHIsoTemplate_Prompt_g_EE", "_MC")

#f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/data.root")
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/data.root")
G_data_EB = loadclone(f, "MR_R2_G_EB", "_data")
G_data_EE = loadclone(f, "MR_R2_G_EB", "_data")
CHIsoTemplate_Fake_EB = loadclone(f, "CHIsoTemplate_Fake_g_EB", "_data")
CHIsoTemplate_Fake_EE = loadclone(f, "CHIsoTemplate_Fake_g_EE", "_data")
if purity_in_Gm1:
    CHIso_GNoIso_EB = loadclone(f, "MR_R2_CHIso_gNoIso_EB", "_data")
    CHIso_GNoIso_EE = loadclone(f, "MR_R2_CHIso_gNoIso_EE", "_data")
else:
    CHIso_GNoIso_EB = loadclone(f, "MR_R2_CHIso_GNoIso_EB", "_data")
    CHIso_GNoIso_EE = loadclone(f, "MR_R2_CHIso_GNoIso_EE", "_data")    
Z_data = loadclone(f, "MR_R2_Z", "_data")
G_data = loadclone(f, "MR_R2_G", "_data")

# Subtracting backgrounds for double ratio
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/gjets.root")
G_GJ    = loadclone(f,"MR_R2_G", "_GJ")
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/nongjets.root")
G_NONGJ = loadclone(f,"MR_R2_G", "_NONGJ")
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/dyjets.root")
Z_DJ    = loadclone(f,"MR_R2_Z", "_DJ")
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/nondyjets.root")
Z_NONDJ = loadclone(f,"MR_R2_Z", "_NONDJ")

#f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ztoinv.root")
f = ROOT.TFile.Open(ZInv_dir1+"/hadd/ztoinv.root")
S_ZI_2D      = loadclone(f, "MR_R2_S"+BIN,  "_ZI")

# --------------------- Lepton Estimate ---------------------------

if purity_in_Gm1:
    f_zinv_est = ROOT.TFile.Open("zinv_est_"+opt.box+"_Gm1.root","recreate")
    plotdir = "Plots/z_inv_est/"+DATE+"/Gm1"
else:
    f_zinv_est = ROOT.TFile.Open("zinv_est_"+opt.box+".root","recreate")
    plotdir = "Plots/z_inv_est/"+DATE

print "Calculate 1-lepton based Z(nunu) estimate"
##  ZInv_lepton_est = L_data_2D.Clone("ZInv_lepton_est")
##  ZInv_lepton_est.Add(L_NONWJ, -1)
##  multiply = S_ZI_2D.Clone()
##  divide   = L_WJ_2D.Clone()
##  for binx in range(1, ZInv_lepton_est.GetNbinsX()+1):
##      for biny in range(1, ZInv_lepton_est.GetNbinsX()+1):
##          if ZInv_lepton_est.GetBinContent(binx, biny)<0:
##              ZInv_lepton_est.SetBinContent(binx,biny,0)
##              ZInv_lepton_est.SetBinError  (binx,biny,0)
##          if multiply.GetBinContent(binx,biny)<0:
##              multiply.SetBinContent(binx,biny,0)
##              multiply.SetBinError(binx,biny,0)
##          if divide.GetBinContent(binx,biny)<0:
##              divide.SetBinContent(binx,biny,0)
##              divide.SetBinError(binx,biny,0)
##  ZInv_lepton_est.Multiply(multiply)
##  ZInv_lepton_est.Divide(divide)
##  # unroll
##  for binx in range(1, ZInv_lepton_est.GetNbinsX()+1):
##      for biny in range(1, ZInv_lepton_est.GetNbinsX()+1):
##          unrolled_bin = (binx-1)*5+biny
##          ztonunu_lepton_est.SetBinContent(unrolled_bin, ZInv_lepton_est.GetBinContent(binx, biny))
##          ztonunu_lepton_est.SetBinError  (unrolled_bin, ZInv_lepton_est.GetBinError  (binx, biny))

ztonunu_lepton_est = bg_est("ZInv_S_L", L_data, [L_TT[0],  L_MJ[0],           L_ZI[0], L_OT[0]], S_WJ[i], L_WJ[i], combine_bins, 1, binned_k)
ztonunu_lepton_est.GetXaxis().SetTitle("Z(#nu#nu) 1-lepton estimate")

# --------------------- Photon Estimate ---------------------------

print "Calculate photon based Z(nunu) estimate"

# --------------------------------
#             Purity

print "- Fitting purity"

# In bins of MR
first = True
purity_MR_EB = ROOT.TH1D("purity_MR_EB",";MR bin;Photon purity",5,0.5,5.5)
purity_MR_EE = ROOT.TH1D("purity_MR_EE",";MR bin;Photon purity",5,0.5,5.5)
for binx in range(1, CHIso_GNoIso_EB.GetNbinsX()+1):
    binx1 = binx
    binx2 = binx
    biny1 = 1
    biny2 = CHIso_GNoIso_EB.GetNbinsY()
    binname= "".join(("MR_"+str(int(CHIso_GNoIso_EB.GetXaxis().GetBinLowEdge(binx1)))+"to",
                      str(int(CHIso_GNoIso_EB.GetXaxis().GetBinLowEdge(binx2)+CHIso_GNoIso_EB.GetXaxis().GetBinWidth(binx2)))+"_",
                      "R2_"+str(CHIso_GNoIso_EB.GetYaxis().GetBinLowEdge(biny1)).replace(".","p")+"to",
                      str(CHIso_GNoIso_EB.GetYaxis().GetBinLowEdge(biny2)+CHIso_GNoIso_EB.GetYaxis().GetBinWidth(biny2)).replace(".","p")))
    chiso_EB = get_zslice(CHIso_GNoIso_EB, "CHIso_EB_fit_"+binname, binx1,binx2,biny1,biny2)
    chiso_EE = get_zslice(CHIso_GNoIso_EE, "CHIso_EE_fit_"+binname, binx1,binx2,biny1,biny2)
    chiso_EB.GetXaxis().SetTitle("Photon Charged Isolation (GeV)")
    chiso_EB.GetYaxis().SetTitle("Photons (w/o ch. iso. cut) - Barrel")
    chiso_EE.GetXaxis().SetTitle("Photon Charged Isolation (GeV)")
    chiso_EE.GetYaxis().SetTitle("Photons (w/o ch. iso. cut) - Endcap")
    temp1_EB = get_zslice(CHIsoTemplate_Prompt_EB, "PrompTemplate_EB_"+binname, min(3,binx1),binx2,biny1,biny2)
    temp1_EE = get_zslice(CHIsoTemplate_Prompt_EE, "PrompTemplate_EE_"+binname, min(3,binx1),binx2,biny1,biny2)
    temp2_EB = get_zslice(CHIsoTemplate_Fake_EB,   "FakeTemplate_EB_"+binname,  min(3,binx1),binx2,biny1,biny2)
    temp2_EE = get_zslice(CHIsoTemplate_Fake_EE,   "FakeTemplate_EE_"+binname,  min(3,binx1),binx2,biny1,biny2)
    pur_EB, pur_EB_err = fit_fraction(chiso_EB, temp1_EB, temp2_EB, plotdir, binname)
    pur_EE, pur_EE_err = fit_fraction(chiso_EE, temp1_EE, temp2_EE, plotdir, binname)
    purity_MR_EB.SetBinContent(binx, pur_EB)
    purity_MR_EB.SetBinError  (binx, pur_EB_err)
    purity_MR_EE.SetBinContent(binx, pur_EE)
    purity_MR_EE.SetBinError  (binx, pur_EE_err)

pur_MR = ROOT.TCanvas("purity_MR")
purity_MR_EB.GetYaxis().SetRangeUser(0,2)
purity_MR_EB.Draw("PE")
purity_MR_EE.SetLineColor(2)
purity_MR_EE.Draw("SAME PE")
leg = ROOT.TLegend(0.28,0.75,0.9,0.9, "")
leg.SetTextSize(0.04)
leg.AddEntry(purity_MR_EB, "Barrel", "LPE")
leg.AddEntry(purity_MR_EE, "Endcap", "LPE")
leg.Draw("SAME")
save_plot(pur_MR, "", plotdir+"/Purity_vs_MR_"+opt.box)

# In bins of R2
first = True
purity_R2_EB = ROOT.TH1D("purity_R2_EB",";R2 bin;Photon purity",5,0.5,5.5)
purity_R2_EE = ROOT.TH1D("purity_R2_EE",";R2 bin;Photon purity",5,0.5,5.5)
for biny in range(1, CHIso_GNoIso_EB.GetNbinsY()+1):
    binx1 = 1
    binx2 = CHIso_GNoIso_EB.GetNbinsX()
    biny1 = biny
    biny2 = biny
    binname= "".join(("MR_"+str(int(CHIso_GNoIso_EB.GetXaxis().GetBinLowEdge(binx1)))+"to",
                      str(int(CHIso_GNoIso_EB.GetXaxis().GetBinLowEdge(binx2)+CHIso_GNoIso_EB.GetXaxis().GetBinWidth(binx2)))+"_",
                      "R2_"+str(CHIso_GNoIso_EB.GetYaxis().GetBinLowEdge(biny1)).replace(".","p")+"to",
                      str(CHIso_GNoIso_EB.GetYaxis().GetBinLowEdge(biny2)+CHIso_GNoIso_EB.GetYaxis().GetBinWidth(biny2)).replace(".","p")))
    chiso_EB = get_zslice(CHIso_GNoIso_EB, "CHIso_EB_fit_"+binname, binx1,binx2,biny1,biny2)
    chiso_EE = get_zslice(CHIso_GNoIso_EE, "CHIso_EE_fit_"+binname, binx1,binx2,biny1,biny2)
    chiso_EB.GetXaxis().SetTitle("Photon Charged Isolation (GeV)")
    chiso_EB.GetYaxis().SetTitle("Photons (w/o ch. iso. cut) - Barrel")
    chiso_EE.GetXaxis().SetTitle("Photon Charged Isolation (GeV)")
    chiso_EE.GetYaxis().SetTitle("Photons (w/o ch. iso. cut) - Endcap")
    temp1_EB = get_zslice(CHIsoTemplate_Prompt_EB, "PrompTemplate_EB_"+binname, binx1,binx2,min(3,biny1),biny2)
    temp1_EE = get_zslice(CHIsoTemplate_Prompt_EE, "PrompTemplate_EE_"+binname, binx1,binx2,min(3,biny1),biny2)
    temp2_EB = get_zslice(CHIsoTemplate_Fake_EB,   "FakeTemplate_EB_"+binname,  binx1,binx2,min(3,biny1),biny2)
    temp2_EE = get_zslice(CHIsoTemplate_Fake_EE,   "FakeTemplate_EE_"+binname,  binx1,binx2,min(3,biny1),biny2)
    pur_EB, pur_EB_err = fit_fraction(chiso_EB, temp1_EB, temp2_EB, plotdir, binname)
    pur_EE, pur_EE_err = fit_fraction(chiso_EE, temp1_EE, temp2_EE, plotdir, binname)
    purity_R2_EB.SetBinContent(biny, pur_EB)
    purity_R2_EB.SetBinError  (biny, pur_EB_err)
    purity_R2_EE.SetBinContent(biny, pur_EE)
    purity_R2_EE.SetBinError  (biny, pur_EE_err)

pur_R2 = ROOT.TCanvas("purity_R2")
purity_R2_EB.GetYaxis().SetRangeUser(0,2)
purity_R2_EB.Draw("PE")
purity_R2_EE.SetLineColor(2)
purity_R2_EE.Draw("SAME PE")
leg = ROOT.TLegend(0.28,0.75,0.9,0.9, "")
leg.SetTextSize(0.04)
leg.AddEntry(purity_R2_EB, "Barrel", "LPE")
leg.AddEntry(purity_R2_EE, "Endcap", "LPE")
leg.Draw("SAME")
save_plot(pur_R2, "", plotdir+"/Purity_vs_R2_"+opt.box)

# Measure average for EB/EE
first = True
purity_EB = ROOT.TH1D("purity_EB",";;Average photon purity",1,0,1)
purity_EE = ROOT.TH1D("purity_EE",";;Average photon purity",1,0,1)
binx1 = 1
binx2 = CHIso_GNoIso_EB.GetNbinsX()
biny1 = 1
biny2 = CHIso_GNoIso_EB.GetNbinsY()
binname= "avg"
chiso_EB = get_zslice(CHIso_GNoIso_EB,         "CHIso_EB_fit_avg",     binx1,binx2,biny1,biny2)
chiso_EE = get_zslice(CHIso_GNoIso_EE,         "CHIso_EE_fit_avg",     binx1,binx2,biny1,biny2)
chiso_EB.GetXaxis().SetTitle("Photon Charged Isolation (GeV)")
chiso_EB.GetYaxis().SetTitle("Photons (w/o ch. iso. cut) - Barrel")
chiso_EE.GetXaxis().SetTitle("Photon Charged Isolation (GeV)")
chiso_EE.GetYaxis().SetTitle("Photons (w/o ch. iso. cut) - Endcap")
temp1_EB = get_zslice(CHIsoTemplate_Prompt_EB, "PrompTemplate_EB_avg", binx1,binx2,biny1,biny2)
temp1_EE = get_zslice(CHIsoTemplate_Prompt_EE, "PrompTemplate_EE_avg", binx1,binx2,biny1,biny2)
temp2_EB = get_zslice(CHIsoTemplate_Fake_EB,   "FakeTemplate_EB_avg",  binx1,binx2,biny1,biny2)
temp2_EE = get_zslice(CHIsoTemplate_Fake_EE,   "FakeTemplate_EE_avg",  binx1,binx2,biny1,biny2)
pur_EB, pur_EB_err = fit_fraction(chiso_EB, temp1_EB, temp2_EB, plotdir, binname)
pur_EE, pur_EE_err = fit_fraction(chiso_EE, temp1_EE, temp2_EE, plotdir, binname)
purity_EB.SetBinContent(1, pur_EB)
purity_EB.SetBinError  (1, pur_EB_err)
purity_EE.SetBinContent(1, pur_EE)
purity_EE.SetBinError  (1, pur_EE_err)

pur = ROOT.TCanvas("purity_avg")
purity_EB.GetYaxis().SetRangeUser(0,2)
purity_EB.Draw("PE")
purity_EE.SetLineColor(2)
purity_EE.Draw("SAME PE")
leg = ROOT.TLegend(0.28,0.75,0.9,0.9, "")
leg.SetTextSize(0.04)
leg.AddEntry(purity_EB, "Barrel", "LPE")
leg.AddEntry(purity_EE, "Endcap", "LPE")
leg.Draw("SAME")
save_plot(pur, "", plotdir+"/Purity_Average_"+opt.box)

# Measure average purity (to subtract QCD etc. in data)
nevt_G_prompt = 0
nevt_G = 0
for binx in range(1, G_data_EB.GetNbinsX()+1):
    for biny in range(1, G_data_EB.GetNbinsY()+1):
        nevt_G        += G_data_EB.GetBinContent(binx,biny) + G_data_EE.GetBinContent(binx,biny)
        nevt_G_prompt += G_data_EB.GetBinContent(binx,biny) * pur_EB + G_data_EE.GetBinContent(binx,biny) * pur_EE
avg_purity_data = nevt_G_prompt / nevt_G
print "Average Purity = "+str(avg_purity_data)

# --------------------------------
#      Direct photon fraction
print "- Calculate direct photon fraction"

# Use an average direct photon fraction for EB/EE (It is ~0.9)
direct_frac_EB = IsDirect_G_EB.Project3D("z").GetMean()
direct_frac_EE = IsDirect_G_EE.Project3D("z").GetMean()

# --------------------------------
#         Transfer factors

print "- Calculate transfer factors"
# Plot the transfer factors
# Binned version
h_transfer_factor_binned = ROOT.TH1D("transfer_factor_binned", ";Bin;S_{Z(#nu#nu)}/G transfer factor", nrazorbin,0,nrazorbin)
can_tf = custom_can(h_transfer_factor_binned, "kappa_ZInv_S_G_"+opt.box)
can_tf.SetLogy(1)
S_ZI = unroll(S_ZI_2D)
G_DP = unroll(GDirectPrompt_2D)
S_ZI.Divide(G_DP)
for binx in range(1, S_ZI.GetNbinsX()+1):
    h_transfer_factor_binned.SetBinContent(binx, S_ZI.GetBinContent(binx))
    h_transfer_factor_binned.SetBinError  (binx, S_ZI.GetBinError  (binx))
h_transfer_factor_binned.SetMarkerColor(633)
h_transfer_factor_binned.SetLineColor  (633)
h_transfer_factor_binned.GetYaxis().SetRangeUser(1e-2,1e3)
h_transfer_factor_binned.Draw("PE1")
h_transfer_factor_binned.Write()
# Factorized form
h_transfer_factor = calc_factorized_kappa(S_ZI_2D, GDirectPrompt_2D, "transfer_factor_G_fact", "S_{Z(#nu#nu)}/G")
h_transfer_factor.Draw("SAME PE1")
h_transfer_factor.Write()

leg = ROOT.TLegend(0.3,0.7,0.6,0.85,"")
leg.AddEntry(h_transfer_factor,        "#color[1]{Factorized form}", "pl")
leg.AddEntry(h_transfer_factor_binned, "#color[633]{Binned}",   "pl")
leg.SetTextSize(0.04)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")

save_plot(can_tf, "", "Plots/kappa/"+DATE+"/"+can_tf.GetName())
#sys.exit()

# --------------------------------
#     Double ratio = k_Z / k_G

print "- Calculate double ratio"
##  k_Z = Z_data.Integral()/Z_MC.Integral()
##  k_G = G_data.Integral()/G_MC[0].Integral()
# For DY in Z region, use MC
Z_data_DY_est = Z_data.Clone("Z_data_DY_est")
Z_data_DY_est.Add(Z_NONDJ, -1)
k_Z = Z_data_DY_est.Integral()/Z_DJ.Integral()
# For GJets in G, instead use data measurement
G_data_G_est = G_data.Clone("G_data_G_est")
##  G_data_G_est.Add(G_NONGJ, -1)
G_data_G_est.Scale(avg_purity_data)
k_G = G_data_G_est.Integral()/G_GJ.Integral()
double_ratio = k_Z / k_G
print "Z_data:  "+str(Z_data.Integral())
print "Z_nonDY: "+str(Z_NONDJ.Integral())
print "Z_DY:    "+str(Z_DJ.Integral())
print "G_data:  "+str(G_data.Integral())
print "G_nonGJ: "+str(G_NONGJ.Integral())
print "G_pure:  "+str(G_data.Integral()*avg_purity_data)
print "G_GJ:    "+str(G_GJ.Integral())
print ("k_Z = %4.3f (%4.3f / %4.3f), k_G = %4.3f (%4.3f / %4.3f), double ratio = %4.3f" %
##     (k_Z, Z_data.Integral(), Z_MC.Integral(), k_G, G_data.Integral(), G_MC[0].Integral(), double_ratio))
       (k_Z, Z_data_DY_est.Integral(), Z_DJ.Integral(), k_G, G_data_G_est.Integral(), G_GJ.Integral(), double_ratio))

# Estimate background
ztonunu_photon_est               = ROOT.TH1D("ztonunu_photon_est",               ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
ztonunu_photon_est_purityUp      = ROOT.TH1D("ztonunu_photon_est_purityUp",      ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
ztonunu_photon_est_purityDown    = ROOT.TH1D("ztonunu_photon_est_purityDown",    ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
ztonunu_photon_est_dirfracUp     = ROOT.TH1D("ztonunu_photon_est_dirfracUp",     ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
ztonunu_photon_est_dirfracDown   = ROOT.TH1D("ztonunu_photon_est_dirfracDown",   ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
ztonunu_photon_est_leptonestUp   = ROOT.TH1D("ztonunu_photon_est_leptonestUp",   ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
ztonunu_photon_est_leptonestDown = ROOT.TH1D("ztonunu_photon_est_leptonestDown", ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
# Use average purities (pur_EB, pur_EE)
for binx in range(1, G_data_EB.GetNbinsX()+1):
    ##  # Use purity in bins of MR
    ##  pur_EB = purity_MR_EB.GetBinContent(binx)
    ##  pur_EE = purity_MR_EE.GetBinContent(binx)
    # And also the transfer factors
    for biny in range(1, G_data_EB.GetNbinsY()+1):
        # transfer factor is factorized in MR and R2
        unrolled_bin = (binx-1)*G_data_EB.GetNbinsY()+biny
        transfer_factor = h_transfer_factor.GetBinContent(unrolled_bin)
        npromptdirect  = G_data_EB.GetBinContent(binx,biny) * pur_EB * direct_frac_EB
        npromptdirect += G_data_EE.GetBinContent(binx,biny) * pur_EE * direct_frac_EE        
        npromptdirect_purityUp     = G_data_EB.GetBinContent(binx,biny) * (pur_EB + 0.1) * direct_frac_EB
        npromptdirect_purityUp    += G_data_EE.GetBinContent(binx,biny) * (pur_EE + 0.1) * direct_frac_EE
        npromptdirect_purityDown   = G_data_EB.GetBinContent(binx,biny) * (pur_EB - 0.1) * direct_frac_EB
        npromptdirect_purityDown  += G_data_EE.GetBinContent(binx,biny) * (pur_EE - 0.1) * direct_frac_EE
        npromptdirect_dirfracUp    = G_data_EB.GetBinContent(binx,biny) * pur_EB * (direct_frac_EB + 0.1)
        npromptdirect_dirfracUp   += G_data_EE.GetBinContent(binx,biny) * pur_EE * (direct_frac_EE + 0.1)
        npromptdirect_dirfracDown  = G_data_EB.GetBinContent(binx,biny) * pur_EB * (direct_frac_EB - 0.1)
        npromptdirect_dirfracDown += G_data_EE.GetBinContent(binx,biny) * pur_EE * (direct_frac_EE - 0.1)
        # Calculate errors
        # add statistical error
        npromptdirect_err              = (G_data_EB.GetBinError(binx,biny) * pur_EB * direct_frac_EB) ** 2
        npromptdirect_err             += (G_data_EE.GetBinError(binx,biny) * pur_EE * direct_frac_EE) ** 2
        npromptdirect_err_purityUp     = (G_data_EB.GetBinError(binx,biny) * (pur_EB + 0.1) * direct_frac_EB) ** 2
        npromptdirect_err_purityUp    += (G_data_EE.GetBinError(binx,biny) * (pur_EE + 0.1) * direct_frac_EE) ** 2
        npromptdirect_err_purityDown   = (G_data_EB.GetBinError(binx,biny) * (pur_EB - 0.1) * direct_frac_EB) ** 2
        npromptdirect_err_purityDown  += (G_data_EE.GetBinError(binx,biny) * (pur_EE - 0.1) * direct_frac_EE) ** 2
        npromptdirect_err_dirfracUp    = (G_data_EB.GetBinError(binx,biny) * pur_EB * (direct_frac_EB + 0.1)) ** 2
        npromptdirect_err_dirfracUp   += (G_data_EE.GetBinError(binx,biny) * pur_EE * (direct_frac_EE + 0.1)) ** 2
        npromptdirect_err_dirfracDown  = (G_data_EB.GetBinError(binx,biny) * pur_EB * (direct_frac_EB - 0.1)) ** 2
        npromptdirect_err_dirfracDown += (G_data_EE.GetBinError(binx,biny) * pur_EE * (direct_frac_EE - 0.1)) ** 2
        npromptdirect_err              = npromptdirect_err             ** 0.5
        npromptdirect_err_purityUp     = npromptdirect_err_purityUp    ** 0.5
        npromptdirect_err_purityDown   = npromptdirect_err_purityDown  ** 0.5
        npromptdirect_err_dirfracUp    = npromptdirect_err_dirfracUp   ** 0.5
        npromptdirect_err_dirfracDown  = npromptdirect_err_dirfracDown ** 0.5
        # fill measured direct prompt photon counts
        ztonunu_photon_est            .SetBinContent(unrolled_bin, npromptdirect)
        ztonunu_photon_est_purityUp   .SetBinContent(unrolled_bin, npromptdirect_purityUp)
        ztonunu_photon_est_purityDown .SetBinContent(unrolled_bin, npromptdirect_purityDown)
        ztonunu_photon_est_dirfracUp  .SetBinContent(unrolled_bin, npromptdirect_dirfracUp)
        ztonunu_photon_est_dirfracDown.SetBinContent(unrolled_bin, npromptdirect_dirfracDown)
        ztonunu_photon_est            .SetBinError  (unrolled_bin, npromptdirect_err)
        ztonunu_photon_est_purityUp   .SetBinError  (unrolled_bin, npromptdirect_err_purityUp)
        ztonunu_photon_est_purityDown .SetBinError  (unrolled_bin, npromptdirect_err_purityDown)
        ztonunu_photon_est_dirfracUp  .SetBinError  (unrolled_bin, npromptdirect_err_dirfracUp)
        ztonunu_photon_est_dirfracDown.SetBinError  (unrolled_bin, npromptdirect_err_dirfracDown)
        # Use Z(nunu) 1-lepton invisible estimate as a systematic
        zinv_est     = npromptdirect     * transfer_factor * double_ratio
        zinv_est_err = npromptdirect_err * transfer_factor * double_ratio
        lep_est      = ztonunu_lepton_est.GetBinContent(unrolled_bin)
        lep_est_syst = (max(zinv_est, lep_est) - min(zinv_est, lep_est)) / 2.0
        ztonunu_photon_est_leptonestUp   .SetBinContent(unrolled_bin, max(0, zinv_est + lep_est_syst))
        ztonunu_photon_est_leptonestUp   .SetBinError  (unrolled_bin, zinv_est_err)
        ztonunu_photon_est_leptonestDown .SetBinContent(unrolled_bin, max(0, zinv_est - lep_est_syst))
        ztonunu_photon_est_leptonestDown .SetBinError  (unrolled_bin, zinv_est_err)

# Finally, multiply with the transfer factors (factorized) and double ratio
ztonunu_photon_est            .Multiply(h_transfer_factor)
ztonunu_photon_est_purityUp   .Multiply(h_transfer_factor)
ztonunu_photon_est_purityDown .Multiply(h_transfer_factor)
ztonunu_photon_est_dirfracUp  .Multiply(h_transfer_factor)
ztonunu_photon_est_dirfracDown.Multiply(h_transfer_factor)
ztonunu_photon_est            .Scale(double_ratio)
ztonunu_photon_est_purityUp   .Scale(double_ratio)
ztonunu_photon_est_purityDown .Scale(double_ratio)
ztonunu_photon_est_dirfracUp  .Scale(double_ratio)
ztonunu_photon_est_dirfracDown.Scale(double_ratio)

# combine bins
if combine_bins:
    ztonunu_photon_est             = combinebins(ztonunu_photon_est,            "ZInv_combined")
    ztonunu_photon_est_purityUp    = combinebins(ztonunu_photon_est_purityUp,   "ZInv_purityUp_combined")
    ztonunu_photon_est_purityDown  = combinebins(ztonunu_photon_est_purityDown, "ZInv_purityDown_combined")
    ztonunu_photon_est_dirfracUp   = combinebins(ztonunu_photon_est_dirfracUp,  "ZInv_dirfracUp_combined")
    ztonunu_photon_est_dirfracDown = combinebins(ztonunu_photon_est_dirfracDown,"ZInv_dirfracDown_combined")
    ztonunu_photon_est_dirfracUp   = combinebins(ztonunu_photon_est_dirfracUp,  "ZInv_leptonestUp_combined")
    ztonunu_photon_est_dirfracDown = combinebins(ztonunu_photon_est_dirfracDown,"ZInv_leptonestDown_combined")

# Add them to vector (used in cards later)
ZInv_est = []
ZInv_est.append(ztonunu_photon_est              .Clone("ZInv"))
ZInv_est.append(ztonunu_photon_est_purityUp     .Clone("ZInv_purityUp"))
ZInv_est.append(ztonunu_photon_est_purityDown   .Clone("ZInv_purityDown"))
ZInv_est.append(ztonunu_photon_est_dirfracUp    .Clone("ZInv_dirfracUp"))
ZInv_est.append(ztonunu_photon_est_dirfracDown  .Clone("ZInv_dirfracDown"))
ZInv_est.append(ztonunu_photon_est_leptonestUp  .Clone("ZInv_leptonestUp"))
ZInv_est.append(ztonunu_photon_est_leptonestDown.Clone("ZInv_leptonestDown"))
for hist in ZInv_est: hist.SetDirectory(0)



# Finally scale by the relative number of events in the njet boxes (found in MC)
# Unroll also the ZToNuNu MC plot
##  ztonunu_mc    = ROOT.TH1D("ztonunu_mc",    ";Bin;Z(#nu#nu) estimate",nrazorbin,0,nrazorbin)
##  for binx in range(1, G_data_EB.GetNbinsX()+1):
##      for biny in range(1, G_data_EB.GetNbinsY()+1):
##          unrolled_bin = (binx-1)*5+biny
##          ztonunu_mc .SetBinContent(unrolled_bin, S_ZI_2D.GetBinContent(binx,biny))
##          ztonunu_mc .SetBinError  (unrolled_bin, S_ZI_2D.GetBinError  (binx,biny))
ztonunu_mc = unroll(S_ZI_2D)
if combine_bins:
    ztonunu_mc = combinebins(unroll(S_ZI_2D), "ztonunu_mc_combined")

# Make plot
can = ROOT.TCanvas("ztonunu_pred")
can.Divide(1,2)
pad = can.cd(1)
pad.SetPad(0,0.3,1,1)
pad.SetBottomMargin(0.02)
ztonunu_mc.GetXaxis().SetLabelSize(0)
ztonunu_mc.GetXaxis().SetLabelColor(0)
ymax = {"WAna_nj35": 20, "WAna_nj45": 20, "WAna_nj46": 20, "WAna_nj6": 8, "WAna_nj7": 8, "TopAna" : 4}
ztonunu_mc.GetYaxis().SetRangeUser(0, ymax[opt.box])
ztonunu_mc.GetYaxis().SetTitleSize(0.07)
ztonunu_mc.GetYaxis().SetTitleOffset(1.1)
ztonunu_mc.SetMarkerStyle(0)
ztonunu_mc.SetLineWidth(2)
ztonunu_mc.SetLineColor(633)
ztonunu_mc.Draw("HISTE")
ztonunu_photon_est.SetMarkerStyle(20)
ztonunu_photon_est.Draw("SAME PE1")
ztonunu_lepton_est.SetMarkerStyle(21)
ztonunu_lepton_est.SetMarkerColor(418)
ztonunu_lepton_est.SetLineColor(418)
ztonunu_lepton_est.Draw("SAME LE1")
title = {
    "WAna_nj35": "W analysis, 3#seqN_{jet}#seq5",
    "WAna_nj45": "W analysis, 4#seqN_{jet}#seq5",
    "WAna_nj46": "W analysis, 4#seqN_{jet}#seq6",
    "WAna_nj6": "W analysis, N_{jet}#seq6",
    "WAna_nj7": "W analysis, N_{jet}#seq7",
    "TopAna": "top analysis"
    }
leg = ROOT.TLegend(0.5,0.75,0.9,0.9, "")
leg.SetTextSize(0.04)
#leg.AddEntry(ztonunu_mc,         "ZToNuNu MC",       "LE")
#leg.AddEntry(ztonunu_photon_est, "#gamma estimate",  "PE")
#leg.AddEntry(ztonunu_lepton_est, "W(l#nu) estimate", "PE")
leg.SetNColumns(2)
leg.AddEntry(ztonunu_mc,         "Z(#nu#nu) MC",       "LE")
leg.AddEntry(0,                  ("#font[82]{%7.1f}" % ztonunu_mc.Integral()), "")
leg.AddEntry(ztonunu_photon_est, "#gamma estimate",  "PE")
leg.AddEntry(0,                  ("#font[82]{%7.1f}" % ztonunu_photon_est.Integral()), "")
leg.AddEntry(ztonunu_lepton_est, "W(l#nu) estimate", "PE")
leg.AddEntry(0,                  ("#font[82]{%7.1f}" % ztonunu_lepton_est.Integral()), "")
leg.Draw("SAME")
pad = can.cd(2)
pad.SetPad(0,0,1,0.3)
pad.SetTopMargin(0.04)
pad.SetBottomMargin(0.4)
pad.SetGridy(1)
# Add ratio plots
ztonunu_photon_ratio = ROOT.TH1D("ztonunu_photon_ratio", ";;Estmate/MC", ztonunu_mc.GetNbinsX(),0,ztonunu_mc.GetNbinsX())
ztonunu_photon_ratio.Divide(ztonunu_photon_est, ztonunu_mc)
ztonunu_photon_ratio.GetXaxis().SetLabelSize(0.125)
ztonunu_photon_ratio.GetYaxis().SetLabelSize(0.125)
ztonunu_photon_ratio.GetYaxis().SetNdivisions(502)
ztonunu_photon_ratio.GetXaxis().SetTitleSize(0.15)
ztonunu_photon_ratio.GetXaxis().SetTitleOffset(1.0)
ztonunu_photon_ratio.GetYaxis().SetTitleSize(0.15)
ztonunu_photon_ratio.GetYaxis().SetTitleOffset(0.5)
ztonunu_photon_ratio.GetYaxis().SetRangeUser(0,2)
ztonunu_photon_ratio.SetMarkerStyle(20)
add_bin_labels(ztonunu_photon_ratio, combine_bins, mrbins, r2bins)
ztonunu_photon_ratio.GetXaxis().LabelsOption("v")
ztonunu_photon_ratio.Draw("PE1")
ztonunu_lepton_ratio = ROOT.TH1D("ztonunu_lepton_ratio", ";;Estmate/MC", ztonunu_mc.GetNbinsX(),0,ztonunu_mc.GetNbinsX())
ztonunu_lepton_ratio.Divide(ztonunu_lepton_est, ztonunu_mc)
ztonunu_lepton_ratio.SetMarkerStyle(21)
ztonunu_lepton_ratio.SetMarkerColor(418)
ztonunu_lepton_ratio.SetLineColor(418)
ztonunu_lepton_ratio.Draw("SAME PE1")
save_plot(can, "", plotdir+"/ZInv_Estimate_"+opt.box)

ztonunu_photon_est.Write()
ztonunu_lepton_est.Write()
ztonunu_mc.Write()
f_zinv_est.Close()
#sys.exit()

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
    
    # Bin centers of low/high vertex distributions in data
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
#sys.exit()

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

# ------------------- Clousre Tests ----------------------

print "Doing closure tests"
if binned_k:
    fout = ROOT.TFile("closure_test_"+opt.box+"_binned_k.root","RECREATE")
else:
    fout = ROOT.TFile("closure_test_"+opt.box+"_factorized_k.root","RECREATE")    

# S'
# s_MJ_est = (Q_data - Q_notMJ) * s_MJ / Q_MJ
# s_WJ_est = (W_data - W_notWJ) * s_WJ / W_WJ
# s_TT_est = (T_data - T_notTT) * s_TT / T_TT
# s_ZI_est = (L_data - L_notWJ) * s_ZI / L_ZI
print "- in S'"

s_Top_est      = []
s_MultiJet_est = []
s_WJets_est    = []
s_ZInv_est     = []
s_Other_est    = []
#T_MJ_est = bg_est("MultiJet_T_Q_est",         Q_data, [Q_TT[0],            Q_WJ[0], Q_ZI[0], Q_OT[0]], T_MJ[0], Q_MJ[0], combine_bins, 1, binned_k)
for i in range(0, len(systematics)):
    #s_TT_est = bg_est("Top_s_T"       +systematics[i], T_data, [          T_MJ_est, T_WJ[0], T_ZI[0], T_OT[0]], s_TT[i], T_TT[i], combine_bins, 1, binned_k)
    s_TT_est = bg_est("Top_s_T"       +systematics[i], T_data, [          T_MJ[0], T_WJ[0], T_ZI[0], T_OT[0]], s_TT[i], T_TT[i], combine_bins, 1, binned_k)
    s_MJ_est = bg_est("MultiJet_s_Q"  +systematics[i], Q_data, [Q_TT[0],           Q_WJ[0], Q_ZI[0], Q_OT[0]], s_MJ[i], Q_MJ[i], combine_bins, 1, binned_k)
    s_WJ_est = bg_est("WJets_s_W"     +systematics[i], W_data, [W_TT[0],  W_MJ[0],          W_ZI[0], W_OT[0]], s_WJ[i], W_WJ[i], combine_bins, 1, binned_k)
    s_ZI_est = bg_est("ZInv_s_L"      +systematics[i], L_data, [L_TT[0],  L_MJ[0],          L_ZI[0], L_OT[0]], s_WJ[i], L_WJ[i], combine_bins, 1, binned_k)
    if combine_bins:
        s_OT_est = combinebins(s_OT[i], "Other_s"+systematics[i])        
    else:
        s_OT_est = s_OT[i].Clone("Other_s"+systematics[i])
    # Sometimes MC sum has negative counts (due to NLO MCs)
    for binx in range(1,s_OT_est.GetNbinsX()+1):
        if s_OT_est.GetBinContent(binx)<0:
            s_OT_est.SetBinContent(binx,0)
            s_OT_est.SetBinError(binx,0)
    s_Top_est     .append(s_TT_est)
    s_MultiJet_est.append(s_MJ_est)
    s_WJets_est   .append(s_WJ_est)
    s_ZInv_est    .append(s_ZI_est)
    s_Other_est   .append(s_OT_est)

s_TT_nom, s_TT_up, s_TT_dn = calc_syst_err(s_Top_est,     "s_Top_syst")
s_MJ_nom, s_MJ_up, s_MJ_dn = calc_syst_err(s_MultiJet_est,"s_MultiJet_syst")
s_WJ_nom, s_WJ_up, s_WJ_dn = calc_syst_err(s_WJets_est,   "s_WJets_syst")
s_ZI_nom, s_ZI_up, s_ZI_dn = calc_syst_err(s_ZInv_est,    "s_ZInv_syst")
s_OT_nom, s_OT_up, s_OT_dn = calc_syst_err(s_Other_est,   "s_Other_syst")
s_syst_err = s_TT_nom.Clone("s_Total")
s_syst_err.SetFillColor(1)
s_syst_err.SetFillStyle(3002)
s_syst_err.SetMarkerStyle(0)
for binx in range(1, s_syst_err.GetNbinsX()+1):
    err_up = (s_TT_up.GetBinContent(binx)**2 + s_MJ_up.GetBinContent(binx)**2 + s_WJ_up.GetBinContent(binx)**2 + s_ZI_up.GetBinContent(binx)**2 + s_OT_up.GetBinContent(binx)**2) ** 0.5
    err_dn = (s_TT_dn.GetBinContent(binx)**2 + s_MJ_dn.GetBinContent(binx)**2 + s_WJ_dn.GetBinContent(binx)**2 + s_ZI_dn.GetBinContent(binx)**2 + s_OT_dn.GetBinContent(binx)**2) ** 0.5
    nom = s_TT_nom.GetBinContent(binx) + s_MJ_nom.GetBinContent(binx) + s_WJ_nom.GetBinContent(binx) + s_ZI_nom.GetBinContent(binx) + s_OT_nom.GetBinContent(binx)
    s_syst_err.SetBinContent(binx, nom + (err_up - err_dn)/2)
    s_syst_err.SetBinError(binx, (err_up+err_dn)/2)

can = custom_can(s_data, "Closure_test_s_"+opt.box)
can.SetLogy(1)
s_data.GetYaxis().SetRangeUser(1.01e-1,1e4)
s_data.GetYaxis().SetTitle("Events/Bins")
s_data.Draw("PE1")
s_TT_nom.SetLineColor(633)
s_MJ_nom.SetLineColor(619)
s_WJ_nom.SetLineColor(418)
s_ZI_nom.SetLineColor(401)
s_OT_nom.SetLineColor(803)
s_TT_nom.SetLineWidth(2)
s_MJ_nom.SetLineWidth(2)
s_WJ_nom.SetLineWidth(2)
s_ZI_nom.SetLineWidth(2)
s_OT_nom.SetLineWidth(2)
s_total = ROOT.THStack("s_TotBkgEst","")
s_total.Add(s_OT_nom)
s_total.Add(s_ZI_nom)
s_total.Add(s_WJ_nom)
s_total.Add(s_MJ_nom)
s_total.Add(s_TT_nom)
s_total.Draw("SAME HIST")
s_syst_err.Draw("SAME E2")
s_data.Draw("SAMEPE1")
leg = ROOT.TLegend(0.65,0.55,0.95,0.85,"")
leg.AddEntry(s_data,   "#color[1]{Data}",                         "pl")
leg.AddEntry(s_TT_nom, "#color[633]{t#bar{t} + single top est.}", "l")
leg.AddEntry(s_MJ_nom, "#color[619]{Multijet est.}",              "l")
leg.AddEntry(s_WJ_nom, "#color[418]{W(#rightarrowl#nu) est.}",    "l")
leg.AddEntry(s_ZI_nom, "#color[401]{Z(#rightarrow#nu#nu) est.}",  "l")
leg.AddEntry(s_OT_nom, "#color[803]{Other (MC)}",                 "l")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
ROOT.gPad.Update()
can.Write()
ratio = add_stack_ratio_plot(can, 0,s_data.GetNbinsX(), 1, combine_bins, mrbins, r2bins)

# Q'
# q_MJ_est = (Q_data - Q_notMJ) * q_MJ / Q_MJ
# q_WJ_est = (W_data - W_notWJ) * q_WJ / W_WJ
# q_TT_est = (T_data - T_notTT) * q_TT / T_TT
# q_ZI_est = (L_data - L_notWJ) * q_ZI / L_ZI
print "- in Q'"

q_Top_est      = []
q_MultiJet_est = []
q_WJets_est    = []
q_ZInv_est     = []
q_Other_est    = []
for i in range(0, len(systematics)):
    #q_TT_est = bg_est("Top_q_T"         +systematics[i], T_data, [          T_MJ_est, T_WJ[0], T_ZI[0], T_OT[0]], q_TT[i], T_TT[i], combine_bins, 1, binned_k)
    q_TT_est = bg_est("Top_q_T"         +systematics[i], T_data, [          T_MJ[0], T_WJ[0], T_ZI[0], T_OT[0]], q_TT[i], T_TT[i], combine_bins, 1, binned_k)
    q_MJ_est = bg_est("MultiJet_q_Q"    +systematics[i], Q_data, [Q_TT[0],           Q_WJ[0], Q_ZI[0], Q_OT[0]], q_MJ[i], Q_MJ[i], combine_bins, 1, binned_k)
    q_WJ_est = bg_est("WJets_q_W"       +systematics[i], W_data, [W_TT[0],  W_MJ[0],          W_ZI[0], W_OT[0]], q_WJ[i], W_WJ[i], combine_bins, 1, binned_k)
    q_ZI_est = bg_est("ZInv_q_L"        +systematics[i], L_data, [L_TT[0],  L_MJ[0],          L_ZI[0], L_OT[0]], q_WJ[i], L_WJ[i], combine_bins, 1, binned_k)
    if combine_bins:
        q_OT_est = combinebins(q_OT[i], "Other_q"+systematics[i])        
    else:
        q_OT_est = q_OT[i].Clone("Other_q"+systematics[i])
    # Sometimes MC sum has negative counts (due to NLO MCs)
    for binx in range(1,q_OT_est.GetNbinsX()+1):
        if q_OT_est.GetBinContent(binx)<0:
            q_OT_est.SetBinContent(binx,0)
            q_OT_est.SetBinError(binx,0)
    q_Top_est     .append(q_TT_est)
    q_MultiJet_est.append(q_MJ_est)
    q_WJets_est   .append(q_WJ_est)
    q_ZInv_est    .append(q_ZI_est)
    q_Other_est   .append(q_OT_est)

q_TT_nom, q_TT_up, q_TT_dn = calc_syst_err(q_Top_est,     "q_Top_syst")
q_MJ_nom, q_MJ_up, q_MJ_dn = calc_syst_err(q_MultiJet_est,"q_MultiJet_syst")
q_WJ_nom, q_WJ_up, q_WJ_dn = calc_syst_err(q_WJets_est,   "q_WJets_syst")
q_ZI_nom, q_ZI_up, q_ZI_dn = calc_syst_err(q_ZInv_est,    "q_ZInv_syst")
q_OT_nom, q_OT_up, q_OT_dn = calc_syst_err(q_Other_est,   "q_Other_syst")
q_syst_err = q_TT_nom.Clone("q_Total")
q_syst_err.SetFillColor(1)
q_syst_err.SetFillStyle(3002)
q_syst_err.SetMarkerStyle(0)
for binx in range(1, q_syst_err.GetNbinsX()+1):
    err_up = (q_TT_up.GetBinContent(binx)**2 + q_MJ_up.GetBinContent(binx)**2 + q_WJ_up.GetBinContent(binx)**2 + q_ZI_up.GetBinContent(binx)**2 + q_OT_up.GetBinContent(binx)**2) ** 0.5
    err_dn = (q_TT_dn.GetBinContent(binx)**2 + q_MJ_dn.GetBinContent(binx)**2 + q_WJ_dn.GetBinContent(binx)**2 + q_ZI_dn.GetBinContent(binx)**2 + q_OT_dn.GetBinContent(binx)**2) ** 0.5
    nom = q_TT_nom.GetBinContent(binx) + q_MJ_nom.GetBinContent(binx) + q_WJ_nom.GetBinContent(binx) + q_ZI_nom.GetBinContent(binx) + q_OT_nom.GetBinContent(binx)
    q_syst_err.SetBinContent(binx, nom + (err_up - err_dn)/2)
    q_syst_err.SetBinError(binx, (err_up+err_dn)/2)

can = custom_can(q_data, "Closure_test_q_"+opt.box)
can.SetLogy(1)
q_data.GetYaxis().SetRangeUser(1.01e-1,1e4)
q_data.GetYaxis().SetTitle("Events/Bins")
q_data.Draw("PE1")
q_TT_nom.SetLineColor(633)
q_MJ_nom.SetLineColor(619)
q_WJ_nom.SetLineColor(418)
q_ZI_nom.SetLineColor(401)
q_OT_nom.SetLineColor(803)
q_TT_nom.SetLineWidth(2)
q_MJ_nom.SetLineWidth(2)
q_WJ_nom.SetLineWidth(2)
q_ZI_nom.SetLineWidth(2)
q_OT_nom.SetLineWidth(2)
q_total = ROOT.THStack("q_TotBkgEst","")
q_total.Add(q_OT_nom)
q_total.Add(q_ZI_nom)
q_total.Add(q_WJ_nom)
q_total.Add(q_MJ_nom)
q_total.Add(q_TT_nom)
q_total.Draw("SAME HIST")
q_syst_err.Draw("SAME E2")
q_data.Draw("SAMEPE1")
leg = ROOT.TLegend(0.65,0.55,0.95,0.85,"")
leg.AddEntry(q_data,   "#color[1]{Data}",                         "pl")
leg.AddEntry(q_TT_nom, "#color[633]{t#bar{t} + single top est.}", "l")
leg.AddEntry(q_MJ_nom, "#color[619]{Multijet est.}",              "l")
leg.AddEntry(q_WJ_nom, "#color[418]{W(#rightarrowl#nu) est.}",    "l")
leg.AddEntry(q_ZI_nom, "#color[401]{Z(#rightarrow#nu#nu) est.}",  "l")
leg.AddEntry(q_OT_nom, "#color[803]{Other (MC)}",                 "l")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
ROOT.gPad.Update()
can.Write()
ratio = add_stack_ratio_plot(can, 0,q_data.GetNbinsX(), 1, combine_bins, mrbins, r2bins)


fout.Close()

# Save the plots, but load them fresh from file
# because objects were lost in memory in the function add_stack_ratio_plot()
f = ROOT.TFile("closure_test_"+opt.box+".root")
can = f.Get("Closure_test_s_"+opt.box+"_Ratio")
if binned_k:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/binned_k/"+can.GetName(), 0)
else:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/factorized_k/"+can.GetName(), 0)
can = f.Get("Closure_test_q_"+opt.box+"_Ratio")
if binned_k:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/binned_k/"+can.GetName(), 0)
else:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/factorized_k/"+can.GetName(), 0)
f.Close()
sys.exit()

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

fout = ROOT.TFile.Open("bkg_estimate_"+opt.box+".root","recreate")
Top_est      = []
MultiJet_est = []
WJets_est    = []
####    ZInv_est     = []
Other_est    = []
#T_MJ_est = bg_est("MultiJet_T_Q_est",         Q_data, [Q_TT[0],            Q_WJ[0], Q_ZI[0], Q_OT[0]], T_MJ[0], Q_MJ[0], combine_bins, 1, binned_k)
for i in range(0, len(systematics)):
    #S_TT_est = bg_est("Top"         +systematics[i], T_data, [          T_MJ_est, T_WJ[0], T_ZI[0], T_OT[0]], S_TT[i], T_TT[i], combine_bins, 1, binned_k)
    S_TT_est = bg_est("Top"         +systematics[i], T_data, [          T_MJ[0], T_WJ[0], T_ZI[0], T_OT[0]], S_TT[i], T_TT[i], combine_bins, 1, binned_k)
    S_MJ_est = bg_est("MultiJet"    +systematics[i], Q_data, [Q_TT[0],           Q_WJ[0], Q_ZI[0], Q_OT[0]], S_MJ[i], Q_MJ[i], combine_bins, 1, binned_k)
    S_WJ_est = bg_est("WJets"       +systematics[i], W_data, [W_TT[0],  W_MJ[0],          W_ZI[0], W_OT[0]], S_WJ[i], W_WJ[i], combine_bins, 1, binned_k)
    ####    S_ZI_est = S_ZI[i].Clone("ZInv" +systematics[i])
    if combine_bins:
        S_OT_est = combinebins(S_OT[i], "Other"+systematics[i])
    else:
        S_OT_est = S_OT[i].Clone("Other"+systematics[i])
    # Sometimes MC sum has negative counts (due to NLO MCs)
    for binx in range(1,S_OT_est.GetNbinsX()+1):
        if S_OT_est.GetBinContent(binx)<0:
            S_OT_est.SetBinContent(binx,0)
            S_OT_est.SetBinError  (binx,0)
    Top_est     .append(S_TT_est)
    MultiJet_est.append(S_MJ_est)
    WJets_est   .append(S_WJ_est)
    ####    ZInv_est    .append(S_ZI_est)
    Other_est   .append(S_OT_est)
fout.Close()
sys.exit()

# Now save a different root file for each signal point
if not opt.nocards:
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
####    '''
####    ------------------------------------------------------------
####    lumi		lnN	1.025	1.025	1.025	1.025	1.025	1.025
####    pileup		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    toppt		shape	-	1.0	1.0	1.0	1.0	1.0
####    isr		shape	1.0	-	-	-	-	-
####    jes		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    jer 		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    met		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    trigger	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    facscale	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    renscale	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    facrenscale	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    alphas		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    elereco	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    eleid		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    eleiso		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    elefastsim	shape	1.0	-	-	-	-	-
####    muontrk	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    muonidiso	shape	1.0	1.0	1.0	1.0	1.0	1.0
####    muonfastsim	shape	1.0	-	-	-	-	-
####    btag		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    btagfastsim	shape	1.0	-	-	-	-	-
####    wtag		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    wtagfastsim	shape	1.0	-	-	-	-	-
####    toptag		shape	1.0	1.0	1.0	1.0	1.0	1.0
####    toptagfastsim	shape	1.0	-	-	-	-	-
####    '''
'''
------------------------------------------------------------
lumi			lnN	1.025	1.025	1.025	1.025	1.025	1.025
pileup			shape	1.0	1.0	1.0	1.0	-	1.0
toppt			shape	-	1.0	1.0	1.0	-	1.0
isr			shape	1.0	-	-	-	-	-
jes			shape	1.0	1.0	1.0	1.0	-	1.0
jer 			shape	1.0	1.0	1.0	1.0	-	1.0
met			shape	1.0	1.0	1.0	1.0	-	1.0
trigger			shape	1.0	1.0	1.0	1.0	-	1.0
facrenscale		shape	1.0	1.0	1.0	1.0	-	1.0
alphas			shape	1.0	1.0	1.0	1.0	-	1.0
ak8scale		shape	-	1.0	1.0	1.0	-	1.0
elereco			shape	1.0	1.0	1.0	1.0	-	1.0
eleid			shape	1.0	1.0	1.0	1.0	-	1.0
eleiso			shape	1.0	1.0	1.0	1.0	-	1.0
elefastsim		shape	1.0	-	-	-	-	-
muontrk			shape	1.0	1.0	1.0	1.0	-	1.0
muonidiso		shape	1.0	1.0	1.0	1.0	-	1.0
muonfastsim		shape	1.0	-	-	-	-	-
btag			shape	1.0	1.0	1.0	1.0	-	1.0
btagfastsim		shape	1.0	-	-	-	-	-
wtag			shape	1.0	1.0	1.0	1.0	-	1.0
wtagfastsim		shape	1.0	-	-	-	-	-
wmistag			shape	1.0	1.0	1.0	1.0	-	1.0
wmistagfastsim		shape	1.0	-	-	-	-	-
wmasstag		shape	-	-	-	1.0	-	-
wantitag		shape	-	-	1.0	-	-	-
toptag			shape	1.0	1.0	-	-	-	1.0
toptagfastsim		shape	1.0	-	-	-	-	-
topmistag		shape	1.0	1.0	1.0	1.0	-	1.0
topmistagfastsim	shape	1.0	-	-	-	-	-
top0bmasstag		shape	1.0	-	-	1.0	-	-
topmasstag		shape	1.0	1.0	1.0	1.0	-	-
topantitag		shape	1.0	-	1.0	-	-	-
purity			shape	-	-	-	-	1.0	-
dirfrac			shape	-	-	-	-	1.0	-
leptonest		shape	-	-	-	-	1.0	-
'''
            )
#facscale	shape	1.0	1.0	1.0	1.0	-	1.0
#renscale	shape	1.0	1.0	1.0	1.0	-	1.0
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
