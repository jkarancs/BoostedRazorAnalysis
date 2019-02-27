import re, os, sys, glob, time, logging, multiprocessing, socket, subprocess, shlex, getpass, ROOT, io
from array import array
from optparse import OptionParser
from HL_common_functions import *
import tdrstyle

keep = [] # Keep ROOT objects in memory

# ---------------------- Cmd Line  -----------------------

# Read options from command line
usage = "Usage: python %prog filelists [options]"
parser = OptionParser(usage=usage)
parser.add_option('-d','--dir',       dest="dir",         type="string",       default="",         help="Input/output directory (use output of Analyzer)")
parser.add_option('-b','--box',       dest="box",         type="string",       default="WAna_nj6", help='Analysis box, eg. TopAna (default="WAna_nj6")')
parser.add_option('-m','--model',     dest="model",       type="string",       default="T5ttcc",   help='Signal model (default="T5ttcc")')
parser.add_option('--nohadd',         dest="nohadd",      action="store_true", default=False,      help='Do not merge input files (default=merge them)')
parser.add_option('--nocards',        dest="nocards",     action="store_true", default=False,      help='Do not create data cards, i.e. run on existing ones (default=create them)')
#parser.add_option('--nocombine',      dest="nocombine",   action="store_true", default=False,      help='Do not rerun combine, i.e. run on existing results (default=run combine)')
parser.add_option('--test',           dest="TEST",        type="int",          default=0,          help="Run only on a N signal points (default=0 - all)")
parser.add_option('--nproc',          dest="NPROC",       type="int",          default=6,          help="Tells how many parallel combine to start (Default=6)")
parser.add_option('-s', '--scenario', dest="scenario",    type="int",          default=0,          help="0: Stat only, 1: YR2018, 2: Run2 systematics")
parser.add_option('-l', '--lumi',     dest="lumi",        type="float",        default=3000,       help="Total integrated luminosity (in fb^-1) to project")
parser.add_option('-e', '--energy',   dest="energy",      type="int",          default=14,         help="Energy (in TeV) to consider (HL: 14 or HE: 27 TeV)")
(opt,args) = parser.parse_args()

BIN = ""
if   "_nj45" in opt.box: BIN = "_nj45"
elif "_nj6"  in opt.box: BIN = "_nj6"
DATE = "_".join(opt.dir.split("/")[-1].split("_")[1:4])

# ---------------------- Settings ------------------------

lumi = opt.lumi # /fb
lumi_Run2 = 35.867
ntuple = "ntuple/Latest"
combine_bins = True
binned_k = True
#use_G = ("WAna" in opt.box)
use_G = True
# QCD: 24% W, 13% top
# DYToLL: 29% W, 19% top
# run script to get above numbers: python scripts/calc_systematics.py
qcd_syst = 0.24
dy_syst  = 0.29
purity_err  = 0.1
dirfrac_err = 0.1
if "TopAna" in opt.box:
    qcd_syst = 0.13
    dy_syst  = 0.19
prefix = ""
if opt.scenario == 0:
    prefix = "Stat. only, "
elif opt.scenario == 1:
    #prefix = "YR18 syst., "
    prefix = ""
    qcd_syst = qcd_syst / 2.0
    dy_syst  = dy_syst  / 2.0
    purity_err  = purity_err  / 10.0
    dirfrac_err = dirfrac_err / 10.0
elif opt.scenario == 2:
    prefix = "Run2 syst., "

mrbins = [ 800, 1000, 1200, 1600, 2000, 4000 ]
mrbins_TeV = [ 0.8, 1.0, 1.2, 1.6, 2.0, 4.0 ]
r2bins = [ 0.08, 0.12, 0.16, 0.24, 0.4, 1.5 ]
nrazorbin = (len(mrbins)-1)*(len(r2bins)-1)
ERA = 65
if opt.energy==27:
    ERA = 75

category = { "WAna_nj45": "#font[52]{W 4-5 jet category}", "WAna_nj6": "#font[52]{W 6 jet category}", "TopAna": "#font[52]{Top category}" }
BOX = category[opt.box]

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

example_signals = {
    "T1tttt" : ["FastSim_SMS-T1tttt"],
    "T2tt"   : ["FastSim_SMS-T2tt_mStop-150to250","FastSim_SMS-T2tt_mStop-250to350","FastSim_SMS-T2tt_mStop-350to400","FastSim_SMS-T2tt_mStop-400to1200"],
    "T5ttcc" : ["FastSim_SMS-T5ttcc", "FastSim_SMS-T5ttcc_mGluino1750to2300"]
}
    
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
    #"_facscaleUp",
    #"_facscaleDown",
    #"_renscaleUp",
    #"_renscaleDown",
    "_facrenscaleUp", 
    "_facrenscaleDown", 
    "_lostlepUp",
    "_lostlepDown",
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

def load_and_scale(vin, rel_scales, keep, name, pf="", combine = False):
    Run2_dir = "results/run_2018_06_08_syst_TopAna"
    if not "TopAna" in opt.box: Run2_dir = "results/run_2018_08_08_syst"
    h_new = 0
    for iFile in range(len(vin)):
        f = ROOT.TFile.Open(Run2_dir+"/hadd/"+vin[iFile]+".root")
        h = f.Get(name)
        if h:
            if h_new == 0:
                if combine:
                    h_new = h.Clone(name+pf+"_tmp")
                    h_new.Scale(rel_scales[vin[iFile]])
                else:
                    h_new = h.Clone(name+pf)                    
                    h_new.Scale(rel_scales[vin[iFile]])
                    h_new.SetDirectory(0)
                    keep.append(h_new)
            else:
                h_new.Add(h, rel_scales[vin[iFile]])
        f.Close()
    if combine:
        h_new = combinebins(h_new, name+pf)
        h_new.SetDirectory(0)
        keep.append(h_new)
    return h_new

def load_and_scale_signal_1d(f, scale, keep, name, pf="", combine = False):
    # signal contains one single xsec for all histos
    # so it is fine to open a single file to speed things up
    h_new = 0
    h = f.Get(name)
    if h:
        if combine:
            h_new = combinebins(h, name+pf)
        else:
            h_new = h.Clone(name+pf)                    
        h_new.Scale(scale)
        h_new.SetDirectory(0)
        keep.append(h_new)
    return h_new

def load_and_scale_signal_3d(vin, rel_scales_signal, keep, name, pf=""):
    Run2_dir = "results/run_2018_06_08_syst_TopAna"
    if not "TopAna" in opt.box: Run2_dir = "results/run_2018_08_08_syst"
    h_new = 0
    for iFile in range(len(vin)):
        fin = ROOT.TFile.Open(Run2_dir+"/hadd/"+vin[iFile]+".root")
        h = fin.Get(name)
        if h:
            h_tmp = h.Clone(name+pf)
            # scaling the counts by the relative xsec and luminosity
            for binx in range(1,h_tmp.GetNbinsX()+1):
                mass = int(h_tmp.GetXaxis().GetBinCenter(binx))
                if mass in rel_scales_signal.keys():
                    scale = rel_scales_signal[mass]
                    for biny in range(1,h_tmp.GetNbinsY()+1):
                        for binz in range(1,h_tmp.GetNbinsZ()+1):
                            c = h_tmp.GetBinContent(binx, biny, binz)
                            e = h_tmp.GetBinError  (binx, biny, binz)
                            h_tmp.SetBinContent(binx, biny, binz, c*scale)
                            h_tmp.SetBinError  (binx, biny, binz, e*scale)
            if h_new == 0:
                h_new = h_tmp
                h_new.SetDirectory(0)
                keep.append(h_tmp)
            else:
                h_new.Add(h_tmp)
        fin.Close()
    return h_new

def scale_down_syst(vh, scale_down_factors, scale_stat = False):
    for i in range(1, len(vh)):
        for syst in scale_down_factors.keys():
            if ("_"+syst+"Up") in vh[i].GetName() or ("_"+syst+"Down") in vh[i].GetName():
                for binx in range(1, vh[i].GetNbinsX()+1):
                    vh[i].SetBinContent(binx, vh[0].GetBinContent(binx)+(vh[i].GetBinContent(binx)-vh[0].GetBinContent(binx))*scale_down_factors[syst])
    if scale_stat:
        for i in range(len(vh)):
            for binx in range(1, vh[i].GetNbinsX()+1):
                vh[i].SetBinError(binx, vh[i].GetBinError(binx) * ((lumi_Run2/lumi) ** 0.5) )

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

def loadclone(f, name, pf="", add_postfix=True):
    h = f.Get(name)
    if add_postfix:
        h_new = h.Clone(name+pf)
    else:
        h_new = h.Clone(pf)        
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
            #if name == "transfer_factor_G_fact":
            #    print str(binx)+" "+str(biny)+" "+str(tf1)+" "+str(tf2)
            tf     = tf1*tf2
            tf_err = (tf1_err*tf1_err*tf2*tf2 + tf2_err*tf2_err*tf1*tf1) ** 0.5
            h_new.SetBinContent(unrolled_bin, tf)
            h_new.SetBinError  (unrolled_bin, tf_err)
    return h_new

##  def find_interval(func, symmetric = False, integ_range=64.0):
def find_interval(func, symmetric = False, integ_range=16.0):
    xmedian = 0
    xup    = 0
    xdown  = 0
    if symmetric:
        xmedian = 1.0
        # Symmetric interval centered on 1.0
        uncertainty = 0.0
        temp = integ_range
        for i in range(30):
            # Approximate both down/up intervals to cover at least 68.27/2
            integ_down = func.Integral(xmedian-uncertainty-temp, xmedian)/func.Integral(-integ_range, integ_range)
            integ_up   = func.Integral(xmedian, xmedian+uncertainty+temp)/func.Integral(-integ_range, integ_range)
            if integ_down < 0.6827/2 and integ_up < 0.6827/2:
                uncertainty += temp
            temp /= 2.0
        xdown = xmedian-uncertainty
        xup   = xmedian+uncertainty
    else:
        # Asymmetric interval
        # Find median (semi-integral) between -integ_range and integ_range
        temp = integ_range
        xmedian = -integ_range
        for i in range(30):
            integ = func.Integral(-integ_range, xmedian+temp)/func.Integral(-integ_range, integ_range)
            #if "MRR2_S_bkg_TT_combined_conv1" in func.GetName():
            #    print str(xmedian+temp)
            #    print func.Integral(-integ_range, xmedian+temp)
            #    print func.Integral(-integ_range, integ_range)
            #    print ("%d - %3.3f, %3.3f, %3.3f - %s" % (i, temp, integ, xmedian, func.GetName()))
            if integ<0.5: xmedian += temp
            temp /= 2.0
        #print "Median is: "+str(xmedian)+" MaxX is: "+str(func.GetMaximumX())
        # Find up/down intervals corresponding to 68.27% / 2 area each
        temp = integ_range/2.0
        up_interval   = 0.0
        down_interval = 0.0
        for i in range(25):
            integ_down = func.Integral(xmedian-down_interval-temp, xmedian)/func.Integral(-integ_range, integ_range)
            integ_up   = func.Integral(xmedian,   xmedian+up_interval+temp)/func.Integral(-integ_range, integ_range)
            if integ_down < 0.6827/2.0:
                down_interval += temp
            if integ_up < 0.6827/2.0:
                up_interval += temp
            temp /= 2.0
            #if "MRR2_S_bkg_TT_combined_conv1" in func.GetName():
            #    print ("%3.3f, %3.3f - %3.3f, %3.3f - %s" % (integ_down, integ_up, down_interval, up_interval, func.GetName()))
        down_interval = min(xmedian, down_interval)
        xdown = xmedian-down_interval
        xup   = xmedian+up_interval
    return [xmedian, xup, xdown]

saved_intervals = {}
saved_nominals  = {}
def fix_low_stat_bins(fout, h_S, h_S_LWP, extrap=0):
    legnames = { "TT":"t#bar{t} or single t", "MJ": "Multijet", "WJ": "W(l#nu)+jets", "ZI":"Z(#nu#nu)+jets", "OT":"Other" }
    names = { "TT":"Top", "MJ": "MultiJet", "WJ": "WJets", "ZI":"ZToNunu", "OT":"Other" }
    sample = names[h_S_LWP.GetName().replace("_combined","").split("_")[-1]]
    legname = legnames[h_S_LWP.GetName().replace("_combined","").split("_")[-1]]
    nosyst = not ("Up" in h_S.GetName() or "Down" in h_S.GetName())
    syst = ""
    if "Up" in h_S.GetName() or "Down" in h_S.GetName():
        if "Other" in h_S.GetName():
            syst = "_"+h_S.GetName().split("_")[-1]
        else:
            syst = "_"+h_S.GetName().split("_")[-3]
    # Count the number of missing bins
    nmissing1 = 0
    nmissing2 = 0
    if nosyst:
        h_added_LWP = h_S_LWP.Clone(h_S_LWP.GetName()+"_added")
        h_S_nominal = h_S.Clone(h_S.GetName()+"_nofix")
        saved_nominals[sample] = h_S_nominal
    else:
        h_S_nominal = saved_nominals[sample]
    for binx in range(1, h_S.GetNbinsX()+1):
        #c1 = h_S.GetBinContent(binx)
        c2 = h_S_LWP.GetBinContent(binx)
        # Check if we need extrapolation from the nominal histogram
        c1_nom = h_S_nominal.GetBinContent(binx)
        if c2>0: avgw2 = h_S_LWP.GetBinError(binx) ** 2 / c2
        if c1_nom <= 0:
            nmissing1 += 1
        if c2 <= 0:
            nmissing2 += 1
            w = 1.83 * prev_w
            #if i == 0: print "Counts for bin "+str(binx)+" is missing in loose region for "+sample+" setting content and error: "+str(w)
            h_S_LWP.SetBinContent(binx, 0)
            h_S_LWP.SetBinError  (binx, w)
            if nosyst:
                #h_added_LWP.SetBinContent(binx, 0)
                h_added_LWP.SetBinContent(binx, 1.01e-3) # To allow showing it on the plot
                h_added_LWP.SetBinError  (binx, w) 
        elif nosyst:
            h_added_LWP.SetBinContent(binx, 0)
            h_added_LWP.SetBinError  (binx, 0)
        prev_w = avgw2
    if nosyst and nmissing1>0:
        print ("- Fixed %2d missing bins for %s" % (nmissing1, sample) )
    # Calculate event ratio between S and LooseWP region
    sum1 = 0
    sum2 = 0
    for binx in range(1, h_S.GetNbinsX()+1):
        if h_S.GetBinContent(binx)>0 and h_S_LWP.GetBinContent(binx)>0:
            sum1 += h_S.GetBinContent(binx)
            sum2 += h_S_LWP.GetBinContent(binx)
    ratio = sum1/sum2
    # Plot the event ratios
    can = fout.Get("Loose_S_Region_"+sample+"_"+opt.box+"_Ratio")
    if not can and nosyst:
        can = custom_can(h_S, "Loose_S_Region_"+sample+"_"+opt.box, 0, 0)
        can.SetLogy(1)
        h_S.GetYaxis().SetTitle("Events / bin")
        h_S.GetYaxis().SetRangeUser(1.01e-3,1e4)
        h_S.SetMarkerStyle(20)
        h_S.SetMarkerColor(1)
        h_S.Draw("PE0")
        h_S_LWP.SetLineColor(1)
        h_S_LWP.SetLineWidth(2)
        h_S_LWP.SetFillColor(0)
        h_S_LWP.Draw("SAME PE0")
        if nmissing2>0:
            h_added_LWP.SetLineColor(2)
            h_added_LWP.SetLineWidth(2)
            h_added_LWP.SetFillColor(0)
            h_added_LWP.SetMarkerStyle(0)
            h_added_LWP.Draw("SAME PE0")
        leg = ROOT.TLegend(0.35,0.65,0.75,0.87, legname+", "+BOX)
        if "WAna" in opt.box:
            leg.AddEntry(h_S_LWP, "#color[1]{Signal region w/o tau_{21} req.}"+(" #color[2]{+Added}" if nmissing2>0 else ""), "l")
        else:
            leg.AddEntry(h_S_LWP, "#color[1]{Signal region w/ loosest top WP}"+(" #color[2]{+Added}" if nmissing2>0 else ""), "l")
        leg.AddEntry(h_S, "#color[1]{Signal region} #color[3]{+Extrapolated}", "pe")
        leg.SetTextSize(0.04)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw("SAME")
        keep.append(leg)
        draw_mr_bins([h_S], 1.01e-3,1e4, combine_bins, keep, mrbins_TeV, r2bins)
        can = add_ratio_plot(can, 0,h_S.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins, ratio)
        add_cms_era(can, ERA+1 if "WAna" in opt.box else ERA+2, keep, opt.energy, lumi, prefix)
        ROOT.gPad.Update()
        can.Write()
    pad2 = can.GetListOfPrimitives().At(1)
    h_ratio = pad2.GetListOfPrimitives().At(1)
    
    # Calculate the confidence interval corresponding to at least 68.3%
    # Convolute gaussians
    # the means are the ratio between Loose and S region
    # the errors are the sigmas
    # Normalize by the average ratio
    # Find the symmetric interval around 1.0 which covers 68.3%/2 on both sides
    # Calculate the interval for the nominal value and use for all other systematics
    if sample in saved_intervals:
        median = saved_intervals[sample][0]
        up     = saved_intervals[sample][1]
        down   = saved_intervals[sample][2]
    else:
        ##  integ_range = 64
        integ_range = 16
        functions = []
        fmax = 0
        all_func = ""
        i = 0
        for xbin in range(1, h_ratio.GetNbinsX()+1):
            if h_ratio.GetBinContent(xbin)>0:
                i = i+1
                #print "Drawing gaus with parameters: norm="+str(1)+", mean="+str(h_ratio.GetBinContent(xbin))+", sigma="+str(h_ratio.GetBinError(xbin))
                gaus = ROOT.TF1("gaus"+str(xbin),"[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))", -integ_range,integ_range)
                ##  gaus.SetParameters(1, h_ratio.GetBinContent(xbin)/ratio, h_ratio.GetBinError(xbin)/ratio)
                gaus.SetParameters(1, h_ratio.GetBinContent(xbin), h_ratio.GetBinError(xbin))
                all_func += "gaus("+str((i-1)*3)+")+"
                gaus.SetLineColor(605+xbin)
                functions.append(gaus)
                if gaus.GetMaximum() > fmax: fmax = gaus.GetMaximum()
        h_ratio_unc = ROOT.TH1D("ratio_unc_"+sample+"_"+opt.box, ";Signal/loose ratio;A. U.", 4*integ_range,-integ_range,integ_range)
        h_ratio_unc.GetYaxis().SetRangeUser(0,fmax*1.3)
        h_ratio_unc.SetStats(0)
        if nosyst:
            can_unc = custom_can(h_ratio_unc, "Ratio_uncertainty_"+sample+"_"+opt.box, 0,0, 500,500, 90,20,50)
            h_ratio_unc.Draw()
        # Total convolution
        conv = ROOT.TF1(h_S.GetName()+"_conv", all_func[:-1], -integ_range,integ_range)
        conv.SetLineColor(1)
        conv.SetLineWidth(2)
        for i in range(0, len(functions)):
            conv.SetParameter(i*3,  functions[i].GetParameter(0))
            conv.SetParameter(i*3+1,functions[i].GetParameter(1))
            conv.SetParameter(i*3+2,functions[i].GetParameter(2))
            if nosyst:
                functions[i].Draw("SAME")
        if nosyst:
            conv.Draw("SAME")
        # Plot also 5 separate convolutions for the 5 MR bins
        all_conv = []
        for mrbin in range(1,6):
            func = ""
            nfilled = 0
            params = []
            for r2bin in range(1, 6-(mrbin>3)-(mrbin>4)):
                xbin = (mrbin-1)*5-(mrbin>4)+r2bin
                if h_ratio.GetBinContent(xbin)>0:
                    func += "gaus("+str(nfilled*3)+")+"
                    params.append(1)
                    ##  params.append(h_ratio.GetBinContent(xbin)/ratio)
                    ##  params.append(h_ratio.GetBinError(xbin)/ratio)
                    params.append(h_ratio.GetBinContent(xbin))
                    params.append(h_ratio.GetBinError(xbin))
                    nfilled += 1
            binned_conv = ROOT.TF1(h_S.GetName()+"_conv"+str(mrbin), func[:-1], -integ_range,integ_range)
            binned_conv.SetLineColor(ROOT.kGreen-3+mrbin)
            for i in range(len(params)):
                binned_conv.SetParameter(i, params[i])
            all_conv.append(binned_conv)
            if nosyst:
                binned_conv.Draw("SAME")
        plotmax = 1.4 * conv.GetMaximum()
        h_ratio_unc.GetYaxis().SetRangeUser(0,plotmax)
        # Find assymmetric interval corresponding to 68.3% area
        # about the median of the convoluted distribution
        median, up, down = find_interval(conv)
        if nosyst:
            integral = conv.Integral(down, up)/conv.Integral(-integ_range, integ_range)
            half_range= 4*max(up-median,median-down)
            h_ratio_unc.GetXaxis().SetRangeUser(median-half_range, median+half_range)
            #print ("%3.3f +%3.3f -%3.3f, integral=%3.3f, %s" % (median, up-median, median-down, integral, sample))
            l0 = ROOT.TLine(median, 0, median, conv.GetMaximum())
            l1 = ROOT.TLine(down,   0, down,   conv.GetMaximum())
            l2 = ROOT.TLine(up,     0, up,     conv.GetMaximum())
            l1.SetLineColor(4)
            l2.SetLineColor(2)
            l0.SetLineStyle(7)
            l1.SetLineStyle(7)
            l2.SetLineStyle(7)
            l0.Draw("SAME")
            l1.Draw("SAME")
            l2.Draw("SAME")
            # Draw also binned limits
            # Calculate average intervals
            avg_median = 0
            avg_up     = 0
            avg_down   = 0
            for i in range(5):
                bin_median, bin_up, bin_down = find_interval(all_conv[i])
                #print str(i)+" "+str(bin_median)+" "+str(bin_up)+" "+str(bin_down)
                avg_median += bin_median/5
                avg_up     += bin_up/5
                avg_down   += bin_down/5
            #print ("%3.3f, %3.3f - %s" % (median, avg_median, h_S.GetName()))
            # Use these for the ucertainty
            median = avg_median
            up     = avg_up
            down   = avg_down
            saved_intervals[sample] = [median, up, down] # Save the interval for systematics variations
            avg_l0 = ROOT.TLine(avg_median, 0, avg_median, conv.GetMaximum())
            avg_l1 = ROOT.TLine(avg_down,   0, avg_down,   conv.GetMaximum())
            avg_l2 = ROOT.TLine(avg_up,     0, avg_up,     conv.GetMaximum())
            avg_l1.SetLineColor(ROOT.kBlue)
            avg_l2.SetLineColor(ROOT.kRed)
            avg_l0.SetLineWidth(2)
            avg_l1.SetLineWidth(2)
            avg_l2.SetLineWidth(2)
            avg_l0.Draw("SAME")
            avg_l1.Draw("SAME")
            avg_l2.Draw("SAME")
            #lat = ROOT.TLatex(median, 0.95*plotmax, "Integral = %.3f" % (integral))
            #lat.SetTextAlign(22)
            #lat.Draw()
            #lat2 = ROOT.TLatex(median, 0.87*plotmax, "Median/Down/Up =")
            #lat2.SetTextAlign(22)
            #lat2.Draw()
            #lat3 = ROOT.TLatex(median, 0.80*plotmax, "%.2f -%.2f +%.2f" % (median, median-down, up-median))
            #lat3.SetTextAlign(22)
            #lat3.Draw()
            leg = ROOT.TLegend(0.15,0.4,0.45,0.55)
            leg.AddEntry(conv,          "Conv - Total",     "l")
            leg.AddEntry(all_conv[2],   "Conv - M_{R} bins","l")
            leg.AddEntry(functions[10], "Gaus - All bins", "l")
            leg.SetTextSize(0.035)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw("SAME")
            leg2 = ROOT.TLegend(0.175,0.7,0.875,0.9, legname+", "+BOX)
            leg2.SetNColumns(2)
            leg2.AddEntry(l1,     "Total - lower limit",  "l")
            leg2.AddEntry(avg_l1, "Avg. M_{R} bins",      "l")
            leg2.AddEntry(l0,     "Total - median",       "l")
            leg2.AddEntry(avg_l0, "Avg. M_{R} bins",      "l")
            leg2.AddEntry(l2,     "Total - upper limit",  "l")
            leg2.AddEntry(avg_l2, "Avg. M_{R} bins",      "l")
            leg2.SetTextSize(0.035)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.Draw("SAME")
            h_ratio_unc.Write()
            add_cms_era(can_unc, ERA, keep, opt.energy, lumi, prefix)
            if nmissing1>0:
                save_plot(can_unc, "", "Plots/lowstat_bins/Extrap_unc_"+sample+"_"+opt.box, 1)
    
    # Set counts for bins where there are counts in the loose region
    h_extrap = h_S.Clone(h_S.GetName()+"_extrap")
    for binx in range(1, h_S.GetNbinsX()+1):
        c1 = h_S.GetBinContent(binx)
        e1 = h_S.GetBinError(binx)
        c2 = h_S_LWP.GetBinContent(binx)
        e2 = h_S_LWP.GetBinError(binx)
        # Check if we need extrapolation from the nominal histogram
        # and use it consitently for all systematics
        c1_nom = h_S_nominal.GetBinContent(binx)
        #if h_S.GetName() == "MRR2_S_bkg_MJ_combined" or h_S.GetName() == "MRR2_S_bkg_triggerUp_MJ_combined" or h_S.GetName() == "MRR2_S_bkg_triggerDown_MJ_combined":
        #    print ("bin=%2d  SR(syst)=%4.2f SR(nominal)=%4.2f - %s" % (binx, h_S.GetBinContent(binx), h_S_nominal.GetBinContent(binx), h_S.GetName()) )
        if c1_nom<=0 and (c2>0 or e2>0):
            # Extrapolate points from the looser region
            #h_S.SetBinContent           (binx, c2*ratio)
            #h_S.SetBinError             (binx, c2*ratio*uncertainty)
            if extrap==1:
                ##  h_S.SetBinContent(binx, c2*ratio*up)
                ##  h_S.SetBinError  (binx, c2*ratio*(up-median))
                h_S.SetBinContent(binx, c2*up)
                h_S.SetBinError  (binx, c2*(up-median))
            elif extrap==-1:
                ##  h_S.SetBinContent(binx, c2*ratio*down)
                ##  h_S.SetBinError  (binx, c2*ratio*(median-down))
                h_S.SetBinContent(binx, c2*down)
                h_S.SetBinError  (binx, c2*(median-down))
            else:
                ##  h_S.SetBinContent(binx, c2*ratio*median)
                ##  comb_err = ( (e2*ratio*median) ** 2 + (c2*ratio*(up-down)/2.0) ** 2 ) ** 0.5
                h_S.SetBinContent(binx, c2*median)
                comb_err = ( (e2*median) ** 2 + (c2*(up-down)/2.0) ** 2 ) ** 0.5
                h_S.SetBinError  (binx, comb_err)
                #if h_S.GetName() == "MRR2_S_bkg_MJ_combined" or h_S.GetName() == "MRR2_S_bkg_triggerUp_MJ_combined" or h_S.GetName() == "MRR2_S_bkg_triggerDown_MJ_combined":
                #    print ("Missing bin=%2d  SR(MC)=%4.2f+-%4.2f - c2=%4.2f, ratio=%4.2f, median=%4.2f - %s" % (binx, h_S.GetBinContent(binx), h_S.GetBinError(binx), c2, ratio, median, h_S.GetName()) )
            # For plotting
            if nosyst:
                if c2>0:
                    # Extrapolated counts
                    #h_extrap.SetBinContent(binx, c2*ratio)
                    #h_extrap.SetBinError  (binx, e2*ratio)             # Old version (without the extrapolation uncertainty)
                    #h_extrap.SetBinError  (binx, c2*ratio*uncertainty) # Newer version, using symmetric extrapolation uncertainty
                    # Latest version, using average asymmetric error (only for plotting) + added stat error from Loose region
                    ##  h_extrap.SetBinContent(binx, c2*ratio*median)
                    ##  comb_err = ( (e2*ratio*median) ** 2 + (c2*ratio*(up-down)/2.0) ** 2 ) ** 0.5
                    h_extrap.SetBinContent(binx, c2*median)
                    comb_err = ( (e2*median) ** 2 + (c2*(up-down)/2.0) ** 2 ) ** 0.5
                    h_extrap.SetBinError  (binx, comb_err)
                else:
                    # Zero counts with 1.83 error
                    h_extrap.SetBinContent(binx, 1.01e-3)         # To allow showing it on the plot
                    #h_extrap.SetBinError  (binx, e2*ratio)       # Old version (without the extrapolation uncertainty)
                    ##  h_extrap.SetBinError  (binx, e2*ratio*median) # Latest version, using average asymmetric error
                    h_extrap.SetBinError  (binx, e2*median) # Latest version, using average asymmetric error
        elif (nmissing1 and nosyst):
            h_extrap.SetBinContent(binx, 0)
            h_extrap.SetBinError  (binx, 0)
            #h_extrap.SetBinContent(binx, 1.01e-3) # To allow showing it on the plot
    if nosyst and nmissing1>0:
        print can.GetName()
        pad1 = can.GetListOfPrimitives().At(0)
        pad1.cd()
        h_extrap.SetMarkerColor(3)
        h_extrap.SetLineColor(3)
        h_extrap.SetLineWidth(2)
        h_extrap.Draw("SAME PE0")
        pad2 = can.GetListOfPrimitives().At(1)
        pad2.cd()
        h_ratio = pad2.GetListOfPrimitives().At(1)
        extrap_ratio = pad2.GetListOfPrimitives().At(1).Clone(h_S.GetName()+"_extrap_ratio")
        # Green extrapolated points in the ratio
        for binx in range(1, h_ratio.GetNbinsX()+1):
            if h_ratio.GetBinContent(binx)>0 or h_S_LWP.GetBinContent(binx)==0:
                extrap_ratio.SetBinContent(binx, 0)
                extrap_ratio.SetBinError  (binx, 0)
            else:
                # New version
                ##  extrap_ratio.SetBinContent(binx, ratio*median)
                extrap_ratio.SetBinContent(binx, median)
                # add also the extrapolation error to the statistical (which is given automatically by the ratio function)
                ##  comb_err = (extrap_ratio.GetBinError(binx) ** 2 + (ratio*(up-down)/2.0) ** 2) ** 0.5
                comb_err = (extrap_ratio.GetBinError(binx) ** 2 + ((up-down)/2.0) ** 2) ** 0.5
                #extrap_ratio.SetBinError  (binx, ratio*(up-down)/2.0)
                extrap_ratio.SetBinError  (binx, comb_err)
        extrap_ratio.SetMarkerStyle(20)
        extrap_ratio.SetMarkerColor(3)
        extrap_ratio.SetLineColor(3)
        extrap_ratio.SetLineWidth(1)
        extrap_ratio.Draw("SAME P")
        can.SetName(can.GetName()+"_PointsAdded")
        can.Write()
        can = fout.Get(can.GetName())
        save_plot(can, "", "Plots/lowstat_bins/Loose_S_Region_"+sample+"_"+opt.box, 1)


def div_err(b1, e1, b2, e2):
    # Taken from ROOT
    # https://root.cern.ch/root/html606/TH1_8cxx_source.html
    return ( (e1*e1*b2*b2 + e2*e2*b1*b1) / (b2*b2*b2*b2) ) ** 0.5

def mult_err(b1, e1, b2, e2):
    # Taken from ROOT
    # https://root.cern.ch/root/html606/TH1_8cxx_source.html
    return ( e1*e1*b2*b2 + e2*e2*b1*b1 ) ** 0.5

def bg_est(name, data, sub, sr, fout, cr, combine_bins, binned_k = True, scale=1.0, merge_CR_bins=True):
    # Create 1 and 2D versions of the MR-R2 plots
    if sr[0].InheritsFrom("TH2"):
        sr_1d = unroll(sr[0])
        sr_2d = sr[0]
    else:
        sr_1d = sr[0]
        sr_2d = rollback(sr[0])
    if len(sr)>1:
        if sr[1].InheritsFrom("TH2"):
            sr_lwp_1d = unroll(sr[1])
            sr_lwp_2d = sr[1]
        else:
            sr_lwp_1d = sr[1]
            sr_lwp_2d = rollback(sr[1])
    if cr.InheritsFrom("TH2"):
        cr_1d = unroll(cr)
        cr_2d = cr
    else:
        cr_1d = cr
        cr_2d = rollback(cr)
    if not ("Up" in name or "Down" in name):
        cr_1d.Write()
        cr_1d.Write()
    # Combine some bins (if needed)
    h_sr            = sr_1d.Clone()
    h_cr            = cr_1d.Clone()
    if combine_bins:
        h_sr            = combinebins(h_sr,            h_sr.GetName()+"_combined")
        h_cr            = combinebins(h_cr,            cr_1d.GetName()+"_combined")
        if binned_k:
            h_cr_data = combinebins(data, name+"_cr_data")
        else:
            h_cr_data = data.Clone(name+"_cr_data_nocomb")
        h_est     = combinebins(data, name)
    else:
        h_cr_data = data.Clone(name+"_cr_data")
        h_est     = data.Clone(name)
    h_sr_nobinfix = h_sr.Clone(h_sr.GetName()+"_nobinfix")
    
    # --------------------- Data Control Region --------------------
    # Subtract MC from data counts
    for hist in sub:
        if combine_bins and binned_k:
            h_cr_data.Add(combinebins(hist, hist.GetName()+"_combined"), -1)
        else:
            h_cr_data.Add(hist, -1)
    if not ("Up" in name or "Down" in name): h_cr_data.Write(name+"_data_nonzeroed")
    # Bin merging does not take into account data counts
    # --> MC subtraction may yield 0/negative yields
    # Zero counts for negative content, let the error stay the same
    # N.B: This only happens for the closure test
    for binx in range(1, h_cr_data.GetNbinsX()+1):
        if h_cr_data.GetBinContent(binx)<0:
            h_cr_data.SetBinContent(binx,0)
    if not ("Up" in name or "Down" in name): h_cr_data.Write(name+"_data_zeroed")
    ##    # Zero bins with negative counts
    ##    # Coming from NLO generator scale variations (flipped sign)
    ##    # See thread: https://hypernews.cern.ch/HyperNews/CMS/get/generators/3510.html
    ##    # TODO: Should we zero the event weight during filling of a histogram if original event weight was positive?
    ##    for binx in range(1, h_cr_data.GetNbinsX()+1):
    ##        if h_cr_data.GetBinContent(binx)<0:
    ##            h_cr_data.SetBinContent(binx,0)
    ##            h_cr_data.SetBinError(binx,0)
    ##    if not ("Up" in name or "Down" in name): h_cr_data.Write(name+"_data_zeroed")
    
    
    # ----------------- Data/MC Control Region---------------------
    # Check which empty (low, 0 or negative) bins needs to be combined for CR
    # Both data and MC --> Requested by convenors, that only check in MC
    threshold_cr = 1 # Specify this, so the bin combination will yield more than that much events
    if merge_CR_bins:
        if not ("Up" in name or "Down" in name):
            h_cr.Write(name+"_mc")
        combinations = [ [5,10], [10,15], [6,11], [11,12], [11,16], [5,10,15], [4,5], [9,10], [14,15], [9,10,14,15], [4,5,9,10], [8,9,10], [14,15,19], [13,14,15], [16,20], [17,21], [18,19], [13,14,15,18,19], [17,18,19], [12,13,14,15,17,18,19], [18,19,22], [21,22], [17,18,19,21,22], [20,21,22], [16,17,18,19,20,21,22] ]
        empty_bins = []
        for binx in range(1, h_cr.GetNbinsX()+1):
            c_data = h_cr_data. GetBinContent(binx)
            c_mc   = h_cr.GetBinContent(binx)
            if (c_data<=threshold_cr or c_mc<=threshold_cr):
            #if c_mc<=threshold_cr:
                empty_bins.append(binx)
        found_combos = []
        selected_combos = {}
        if len(empty_bins):
            # Check which combination is the best:
            for i in list(reversed(range(len(empty_bins)))):
                # Check if the bin is already in a previously combined bin
                selected = -1
                for j in range(len(found_combos)):
                    if empty_bins[i] in found_combos[j]:
                        best_combo = found_combos[j]
                        selected = j
                # If not, look for the best combination (lowest sum)
                if selected == -1:
                    bestsum_mc   = 9.99e6
                    bestsum_data = 9.99e6
                    best_combo = []
                    extend_previous = []
                    for combo in combinations:
                        if empty_bins[i] in combo:
                            # Check if the combination can be considered
                            # i.e. the bins in them are not already merged
                            # Check also if new bin(s) can be merged with previously merged ones
                            # The previously merged bins won't count in the sum when comparing
                            # with new combinations
                            can_consider = True
                            extension_candidates = []
                            for j in range(len(found_combos)):
                                for cand_bin in combo:
                                    if cand_bin in found_combos[j]:
                                        can_consider = False
                                extends_prev_found = True
                                for prev_found_bin in found_combos[j]:
                                    if not (prev_found_bin in combo):
                                        extends_prev_found = False
                                if extends_prev_found:
                                    extension_candidates.append(j)
                            if can_consider or len(extension_candidates):
                                combsum_data = [0] * (1+len(extension_candidates))
                                combsum_mc   = [0] * (1+len(extension_candidates))
                                for binx in combo:
                                    for j in range(len(combsum_data)):
                                        if j==0 or not (binx in found_combos[extension_candidates[j-1]]):
                                            combsum_data[j] += h_cr_data .GetBinContent(binx)
                                            combsum_mc[j]   += h_cr.GetBinContent(binx)
                                # Save the combination which yields the lowest sum
                                for k in range(len(combsum_data)):
                                    #if not ("Up" in name or "Down" in name):
                                    #    print ("k=%d/%d, data=%4.2f mc=%4.2f" % (k, len(combsum_data), combsum_data[k], combsum_mc[k]) )
                                    threshold = (threshold_cr if j==0 else 0) # Becasue previously already passed threshold
                                    # We favor a combination with previously merged bins
                                    if combsum_data[k]>=threshold and combsum_mc[k]>=threshold:
                                        if combsum_mc[k]<bestsum_mc and combsum_data[k]<bestsum_data:
                                            bestsum_mc   = combsum_mc[k]
                                            bestsum_data = combsum_data[k]
                                            best_combo = combo
                                            extend_previous = extension_candidates
                    # If new combination is found, append it
                    if len(extend_previous)==0:
                        selected = len(found_combos)
                        found_combos.append(best_combo)
                        #if not ("Up" in name or "Down" in name) and len(found_combos):
                        #    print "- New combination found for "+name+" bin "+str(empty_bins[i])+":",
                        #    for j in range(len(best_combo)):
                        #        print ("%d" % best_combo[j]),
                        #        if j!=(len(best_combo)-1): print ",",
                        #    print
                    # Otherwise extend the previous one
                    else:
                        for prev in extend_previous:
                            selected = prev
                            found_combos[prev] = best_combo
                        #if not ("Up" in name or "Down" in name) and len(found_combos):
                        #    print "- Extending previous combination for "+name+" bin "+str(empty_bins[i])+":",
                        #    for j in range(len(best_combo)):
                        #        print ("%d" % best_combo[j]),
                        #        if j!=(len(best_combo)-1): print ",",
                        #    print
                selected_combos[empty_bins[i]] = selected
        # Creating final best combinations
        final_combos = {}
        for i in list(reversed(range(len(empty_bins)))):
            final_combos[empty_bins[i]] = found_combos[selected_combos[empty_bins[i]]]
        if not ("Up" in name or "Down" in name) and len(found_combos):
            print "- Merged CR bins for "+name+":",
            print final_combos
        # Add rest of the unmerged bins as they are
        for binx in range(1, h_cr.GetNbinsX()+1):
            if not binx in empty_bins:
                for j in range(len(found_combos)):
                    if binx in found_combos[j]:
                        final_combos[binx] = found_combos[j]
                if not binx in final_combos:
                    final_combos[binx] = [binx]
    
    
    # ------------------------ Signal Region -----------------------
    # Fix low stat bins for Signal region
    # (if looser region is available to extrapolate from)
    if len(sr)>1:
        if combine_bins:
            h_sr_lwp = combinebins(sr_lwp_1d, sr_lwp_1d.GetName()+"_combined")
        else:
            h_sr_lwp = sr_lwp_1d.Clone()
        if not ("Up" in name or "Down" in name):
            h_sr.Write(name+"_sr")
        doExtrap = 1 if ("_extrapUp" in name) else (-1 if ("_extrapDown" in name) else 0)
        fix_low_stat_bins(fout, h_sr, h_sr_lwp, doExtrap)
        if not ("Up" in name or "Down" in name):
            h_sr           .Write(name+"_sr_fixed")
        #if name == "MultiJet" or name == "MultiJet_triggerUp" or name == "MultiJet_triggerDown":
        #    for binx in range(1, h_est.GetNbinsX()+1):
        #        if h_sr_nobinfix.GetBinContent(binx)==0:
        #            print ("Missing bin=%2d  SR(MC)=%4.2f+-%4.2f - %s" % (binx, h_sr.GetBinContent(binx), h_sr.GetBinError(binx), h_sr.GetName()) )
    
    
    # ----------------------- Final Estimate -----------------------
    # ----------------------- Previos Method -----------------------
    # Previous method with transfer factors (kappas)
    # Keep plots and histos for now
    nbin = h_sr.GetNbinsX()
    srname = sr_1d.GetName().split("_")[1].replace("s","S'").replace("q","Q'")
    crname = cr_1d.GetName().split("_")[1]
    h_kappa_binned        = ROOT.TH1D("kappa_binned_"+srname+"_"+crname, ";;"+srname+"/"+crname+" transfer factor", nbin,0,nbin)
    add_bin_labels(h_kappa_binned, combine_bins)
    h_kappa_binned.Divide(h_sr, h_cr)
    nbin_nocomb = sr_1d.GetNbinsX() # using unmerged bins
    h_kappa_binned_nocomb = ROOT.TH1D("kappa_binned_nocomb_"+srname+"_"+crname, ";;"+srname+"/"+crname+" transfer factor", nbin_nocomb,0,nbin_nocomb)
    add_bin_labels(h_kappa_binned_nocomb, 0)
    h_kappa_binned_nocomb.Divide(sr_1d, cr_1d)
    # Zero also negative kappas (for binned version only)
    for binx in range(1, h_sr.GetNbinsX()+1):
        if h_kappa_binned.GetBinContent(binx)<0:
            h_kappa_binned.SetBinContent(binx,0)
            h_kappa_binned.SetBinError(binx,0)
    for binx in range(1, sr_1d.GetNbinsX()+1):
        if h_kappa_binned_nocomb.GetBinContent(binx)<0:
            h_kappa_binned_nocomb.SetBinContent(binx,0)
            h_kappa_binned_nocomb.SetBinError(binx,0)
    # Factorized form
    h_kappa_fact = calc_factorized_kappa(sr_2d, cr_2d, "kappa_factorized_"+srname+"_"+crname, srname+"/"+crname)
    if not ("Up" in name or "Down" in name):
        h_kappa_binned.Write()
        h_kappa_binned_nocomb.Write()
        h_kappa_fact.Write()
    # Draw kappa comparison plot
    if not ("Up" in name or "Down" in name):
        processname = name.split("_")[0]
        can_tf = custom_can(h_kappa_binned, "kappa_"+processname+"_"+srname+"_"+crname+"_"+opt.box, 0,0)
        #can_tf.SetLogy(1)
        h_kappa_binned_nocomb.SetMarkerColor(633)
        h_kappa_binned_nocomb.SetLineColor  (633)
        #h_kappa_binned_nocomb.GetYaxis().SetRangeUser(1e-2,1e3)
        h_kappa_binned_nocomb.Draw("PE0")
        ###   h_kappa_fact.Draw("SAME PE0")
        leg = ROOT.TLegend(0.3,0.7,0.6,0.85, "No bin merging")
        ###   leg.AddEntry(h_kappa_fact,   "#color[1]{Factorized form}", "pe")
        leg.AddEntry(h_kappa_binned, "#color[633]{SR/CR transfer factors}",        "pe")
        leg.SetTextSize(0.04)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw("SAME")
        #save_plot(can_tf, "", "Plots/kappa/"+DATE+"/"+can_tf.GetName())
    # Final estimate
    if not merge_CR_bins:
        # with the kappa of choice
        if binned_k:
            h_est.Multiply(h_kappa_binned)
        else:
            h_est.Multiply(h_kappa_fact)
        # For the factorized kappa version
        # Merge only the bins after the estimate
        if (not binned_k) and combine_bins:
            h_est = combinebins(h_est, name)
    
    
    # ------------------------- New Method -------------------------
    # New, mathematically equivalent way for final estimate
    # It is the same as for the inclsuive analysis
    # This formulation allows the merging of bins in the control region
    if merge_CR_bins:
        prev_corr = 0
        prev_corr_err = 0
        for binx in range(1, h_est.GetNbinsX()+1):
            c_sr_mc   = h_sr.GetBinContent(binx)
            if opt.scenario>0:
                e_sr_mc   = 0
            else:
                e_sr_mc   = h_sr.GetBinError  (binx) * ((lumi_Run2/lumi) ** 0.5)
            # Use the previously found best combination for the merged CR bins
            c_cr_data = 0
            c_cr_mc   = 0
            e_cr_data = 0
            e_cr_mc   = 0
            for unmerged_bin in final_combos[binx]:
                c_cr_data += h_cr_data .GetBinContent(unmerged_bin)
                c_cr_mc   += h_cr.GetBinContent(unmerged_bin)
                e_cr_data  = (h_cr_data .GetBinError(unmerged_bin) ** 2 * lumi_Run2/lumi + e_cr_data ** 2 ) ** 0.5
                if opt.scenario>0:
                    e_cr_mc    = (h_cr.GetBinError(unmerged_bin) ** 2 * lumi_Run2/lumi + e_cr_mc ** 2   ) ** 0.5
            #print ("binx=%2d SR(MC)=%4.2f CR(Data)=%4.2f CR(MC)=%4.2f" % (binx, c_sr_mc, c_cr_data, c_cr_mc) )
            corr     = c_cr_data / c_cr_mc
            corr_err = div_err(c_cr_data, e_cr_data, c_cr_mc, e_cr_mc)
            # The bin merging did not merge sufficient bins (we can simply take the value of the previous correction)
            # This procedure gives the same results as that used for Run2 (where bin merging involved the previous bin)
            if corr>2.0:
                corr = 2.0
                #corr = prev_corr
                corr_err = prev_corr_err
            prev_corr = corr
            prev_corr_err = corr_err
            estimate = c_sr_mc * corr
            error    = mult_err(c_sr_mc, e_sr_mc, corr, corr_err)
            #if name == "MultiJet":
            #    print ("bin=%2d  SR(MC)=%4.2f+-%4.2f  CR(Data)=%4.2f+-%4.2f  CR(MC)=%4.2f+-%4.2f --> EST=%4.2f+-%4.2f - %s" % (binx, c_sr_mc, e_sr_mc, c_cr_data, e_cr_data, c_cr_mc, e_cr_mc, estimate, error, name) )
            #if name == "Top_s_T":
            #    print ("bin=%2d  SR(MC)=%4.2f+-%4.2f  CR(Data)=%4.2f+-%4.2f  CR(MC)=%4.2f+-%4.2f --> EST=%4.2f+-%4.2f - %s" % (binx, c_sr_mc, e_sr_mc, c_cr_data, e_cr_data, c_cr_mc, e_cr_mc, estimate, error, name) )
            h_est.SetBinContent(binx,estimate)
            h_est.SetBinError  (binx,error)
            #if not ("Up" in name or "Down" in name):
            #    print ("Bin=%2d -  Est = %4.2f +- %4.2f = (%4.2f +- %4.2f) * (%4.2f +- %4.2f) / (%4.2f +- %4.2f)" % (binx, estimate, error, c_sr_mc, e_sr_mc, c_cr_data, e_cr_data, c_cr_mc, e_cr_mc) )
    
    
    # If specified, scale the estimate by a factor, eg. with the double ratio
    if scale != 1.0: h_est.Scale(scale)
    if not ("Up" in name or "Down" in name):
        h_est.Write(name+"_est")
    
    
    return h_est

def zinv_est(name, vvh_sr, fout, vh_g_data_promptdirect, vh_g_cr, l_data, l_subtract, l_cr, combine_bins, binned_k, DR, eDR):
    ZInv_est = []
    # G region based estimate
    if len(vvh_sr)>1:
        ZInv_est.append(bg_est(name,                vh_g_data_promptdirect[0], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0], combine_bins, binned_k, DR))
    else:
        ZInv_est.append(bg_est(name,                vh_g_data_promptdirect[0], [], [vvh_sr[0][0]],               fout, vh_g_cr[0], combine_bins, binned_k, DR))
    # Systematics on double ratio, purity, direct photon fraction and extrapolation
    if len(vvh_sr)>1:
        ZInv_est.append(bg_est(name+"_doubleratioUp",   vh_g_data_promptdirect[0], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0],  combine_bins, binned_k, DR+eDR))
        ZInv_est.append(bg_est(name+"_doubleratioDown", vh_g_data_promptdirect[0], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0],  combine_bins, binned_k, DR-eDR))
        ZInv_est.append(bg_est(name+"_purityUp",        vh_g_data_promptdirect[1], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0],  combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_purityDown",      vh_g_data_promptdirect[2], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0],  combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_dirfracUp",       vh_g_data_promptdirect[3], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0],  combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_dirfracDown",     vh_g_data_promptdirect[4], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0],  combine_bins, binned_k, DR))
    else:
        ZInv_est.append(bg_est(name+"_doubleratioUp",   vh_g_data_promptdirect[0], [], [vvh_sr[0][0]],               fout, vh_g_cr[0],  combine_bins, binned_k, DR+eDR))
        ZInv_est.append(bg_est(name+"_doubleratioUp",   vh_g_data_promptdirect[0], [], [vvh_sr[0][0]],               fout, vh_g_cr[0],  combine_bins, binned_k, DR-eDR))
        ZInv_est.append(bg_est(name+"_purityUp",        vh_g_data_promptdirect[1], [], [vvh_sr[0][0]],               fout, vh_g_cr[0],  combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_purityDown",      vh_g_data_promptdirect[2], [], [vvh_sr[0][0]],               fout, vh_g_cr[0],  combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_dirfracUp",       vh_g_data_promptdirect[3], [], [vvh_sr[0][0]],               fout, vh_g_cr[0],  combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_dirfracDown",     vh_g_data_promptdirect[4], [], [vvh_sr[0][0]],               fout, vh_g_cr[0],  combine_bins, binned_k, DR))
    ztonunu_photon_est = ZInv_est[0].Clone(name+"_photon_est")
    # Standard systematics
    for i in range(1, len(vvh_sr[0])):
        systematic = vvh_sr[0][i].GetName().split("_")[-2]
        if len(vvh_sr)>1:
            ZInv_est.append(bg_est(name+"_"+systematic, vh_g_data_promptdirect[0], [], [vvh_sr[0][i], vvh_sr[1][i]], fout, vh_g_cr[i], combine_bins, binned_k, DR))
        else:
            ZInv_est.append(bg_est(name+"_"+systematic, vh_g_data_promptdirect[0], [], [vvh_sr[0][i]],               fout, vh_g_cr[i], combine_bins, binned_k, DR))
    # Extrapolation uncertainty (only if loose region is available, add string to name is enough)
    if len(vvh_sr)>1:
        ZInv_est.append(bg_est(name+"_extrapUp",        vh_g_data_promptdirect[0], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0], combine_bins, binned_k, DR))
        ZInv_est.append(bg_est(name+"_extrapDown",      vh_g_data_promptdirect[0], [], [vvh_sr[0][0], vvh_sr[1][0]], fout, vh_g_cr[0], combine_bins, binned_k, DR))
    # L region based estimate (as an additional systematic)
    if len(vvh_sr)>1:
        ztonunu_lepton_est = bg_est(name+"_lepest", l_data, l_subtract, [vvh_sr[0][0], vvh_sr[1][0]], fout, l_cr, combine_bins, binned_k)
    else:
        ztonunu_lepton_est = bg_est(name+"_lepest", l_data, l_subtract, [vvh_sr[0][0]],               fout, l_cr, combine_bins, binned_k)
    ztonunu_lepton_est.GetXaxis().SetTitle("Z(#nu#nu) 1-lepton estimate")
    ztonunu_leptonestUp   = ztonunu_lepton_est.Clone("h_"+name+"_leptonestUp")
    ztonunu_leptonestDown = ztonunu_lepton_est.Clone("h_"+name+"_leptonestDown")
    avg_closure_syst = 0
    for unrolled_bin in range(1, ztonunu_lepton_est.GetNbinsX()+1):
        zinv_est     = ztonunu_photon_est.GetBinContent(unrolled_bin)
        zinv_est_err = ztonunu_photon_est.GetBinError(unrolled_bin)
        lep_est      = ztonunu_lepton_est.GetBinContent(unrolled_bin)
        lep_est_syst = (max(zinv_est, lep_est) - min(zinv_est, lep_est)) / 2.0
        if lep_est>0:
            avg_closure_syst += lep_est_syst/lep_est*100
        ztonunu_leptonestUp  .SetBinContent(unrolled_bin, max(0, zinv_est + lep_est_syst))
        ztonunu_leptonestDown.SetBinContent(unrolled_bin, max(0, zinv_est - lep_est_syst))
        ztonunu_leptonestUp  .SetBinError  (unrolled_bin, zinv_est_err)
        ztonunu_leptonestDown.SetBinError  (unrolled_bin, zinv_est_err)
    avg_closure_syst /= 22
    if name == "ZInv": print "Average ZInv closure test systematics: "+str(avg_closure_syst)
    if combine_bins:
        ZInv_est.append(combinebins(ztonunu_leptonestUp,   name+"_leptonestUp"))
        ZInv_est.append(combinebins(ztonunu_leptonestDown, name+"_leptonestDown"))
    else:
        ZInv_est.append(ztonunu_leptonestUp  .Clone(name+"_leptonestUp"))
        ZInv_est.append(ztonunu_leptonestDown.Clone(name+"_leptonestDown"))
    for hist in ZInv_est: hist.SetDirectory(0)
    # Make G/L based estimate comparison plot
    can = ROOT.TCanvas(name)
    can.Divide(1,2)
    pad = can.cd(1)
    pad.SetPad(0,0.4,1,1)
    pad.SetBottomMargin(0.02)
    pad.SetTopMargin   (0.10)
    if combine_bins:
        ztonunu_mc = combinebins(vvh_sr[0][0], "ztonunu_mc_combined")
    else:
        ztonunu_mc = vvh_sr[0][0].Clone("ztonunu_mc")
    ztonunu_mc.GetXaxis().SetLabelSize(0)
    ztonunu_mc.GetXaxis().SetLabelColor(0)
    ymax = {"WAna_nj35": 24, "WAna_nj45": 12, "WAna_nj46": 24, "WAna_nj6": 10, "WAna_nj7": 10, "TopAna" : 5}
    if "q" in name:
        if "Top" in opt.box:
            ymax[opt.box] = 25 * ymax[opt.box]
        else:
            ymax[opt.box] = 2.5 * ymax[opt.box]
    ztonunu_mc.GetYaxis().SetRangeUser(0, ymax[opt.box])
    ztonunu_mc.GetYaxis().SetTitle("Events / bin")
    ztonunu_mc.GetYaxis().SetTitleSize(0.1)
    ztonunu_mc.GetYaxis().SetTitleOffset(0.8)
    ztonunu_mc.SetMarkerStyle(0)
    ztonunu_mc.SetLineWidth(2)
    ztonunu_mc.SetLineColor(633)
    ztonunu_mc.Draw("HISTE")
    ztonunu_photon_est.SetMarkerStyle(20)
    ztonunu_photon_est.Draw("SAME PE0")
    ztonunu_lepton_est.SetMarkerStyle(21)
    ztonunu_lepton_est.SetMarkerColor(418)
    ztonunu_lepton_est.SetLineColor(418)
    ztonunu_lepton_est.Draw("SAME LE0")
    if "WAna_nj45" in opt.box:
        draw_mr_bins([ztonunu_photon_est, ztonunu_lepton_est, ztonunu_mc], 0, ymax[opt.box], combine_bins, keep, mrbins_TeV, r2bins, {}, 0, 0.25, 0.06)
    elif "WAna_nj6" in opt.box:
        draw_mr_bins([ztonunu_lepton_est], 0, ymax[opt.box], combine_bins, keep, mrbins_TeV, r2bins, {}, 0, 1.0, 0.06)
    else:
        draw_mr_bins([ztonunu_photon_est], 0, ymax[opt.box], combine_bins, keep, mrbins_TeV, r2bins, {}, 0, 1.0, 0.06)
    legtitle = "Signal region"
    if "s_" in name:
        legtitle = "Signal-like validation region"
    elif "q_" in name:
        legtitle = "Multijet validation region"
    leg = ROOT.TLegend(0.5,0.57,0.7,0.87, legtitle+", "+BOX)
    #leg.SetNColumns(2)
    leg.SetTextSize(0.05)
    #leg.AddEntry(ztonunu_mc,         "ZToNuNu MC",       "LE")
    #leg.AddEntry(ztonunu_photon_est, "#gamma estimate",  "PE")
    #leg.AddEntry(ztonunu_lepton_est, "W(l#nu) estimate", "PE")
    leg.AddEntry(ztonunu_mc,         "Z(#nu#nu) MC",       "LE")
    #leg.AddEntry(0,                  ("#font[82]{%7.1f}" % ztonunu_mc.Integral()), "")
    leg.AddEntry(ztonunu_photon_est, "#gamma estimate",  "PE")
    #leg.AddEntry(0,                  ("#font[82]{%7.1f}" % ztonunu_photon_est.Integral()), "")
    leg.AddEntry(ztonunu_lepton_est, "W(l#nu) estimate", "PE")
    #leg.AddEntry(0,                  ("#font[82]{%7.1f}" % ztonunu_lepton_est.Integral()), "")
    leg.Draw("SAME")
    pad = can.cd(2)
    pad.SetPad(0,0,1,0.4)
    pad.SetTopMargin(0.04)
    pad.SetBottomMargin(0.55)
    pad.SetGridy(1)
    # Add ratio plots
    ztonunu_photon_ratio = ROOT.TH1D(name+"_ratio", ";;Estimate/MC", ztonunu_mc.GetNbinsX(),0,ztonunu_mc.GetNbinsX())
    ztonunu_photon_ratio.Divide(ztonunu_photon_est, ztonunu_mc)
    ztonunu_photon_ratio.GetXaxis().SetLabelOffset(0.02)
    ztonunu_photon_ratio.GetXaxis().SetLabelSize(0.125)
    ztonunu_photon_ratio.GetYaxis().SetLabelSize(0.125)
    ztonunu_photon_ratio.GetYaxis().SetNdivisions(502)
    ztonunu_photon_ratio.GetXaxis().SetTitleSize(0.15)
    #ztonunu_photon_ratio.GetXaxis().SetTitleOffset(1.0) # set ofset in add_bin_labels last option
    ztonunu_photon_ratio.GetYaxis().SetTitleSize(0.15)
    ztonunu_photon_ratio.GetYaxis().SetTitleOffset(0.5)
    ztonunu_photon_ratio.GetYaxis().SetRangeUser(0,2)
    ztonunu_photon_ratio.SetMarkerStyle(20)
    add_bin_labels(ztonunu_photon_ratio, combine_bins, mrbins, r2bins, 2.0)
    ztonunu_photon_ratio.GetXaxis().LabelsOption("v")
    ztonunu_photon_ratio.Draw("PE0")
    ztonunu_lepton_ratio = ROOT.TH1D("ztonunu_lepton_ratio", ";;Estmate/MC", ztonunu_mc.GetNbinsX(),0,ztonunu_mc.GetNbinsX())
    ztonunu_lepton_ratio.Divide(ztonunu_lepton_est, ztonunu_mc)
    ztonunu_lepton_ratio.SetMarkerStyle(21)
    ztonunu_lepton_ratio.SetMarkerColor(418)
    ztonunu_lepton_ratio.SetLineColor(418)
    ztonunu_lepton_ratio.Draw("SAME PE0")
    add_cms_era(can, ERA, keep, opt.energy, lumi, prefix)
    save_plot(can, "", plotdir+"/"+name+"_Estimate_"+opt.box)
    return ZInv_est

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

def fit_fraction(chiso, temp1, temp2, plotdir, bin, mcpurity=0.0):
    temp1.Write()
    temp2.Write()
    EBEE = ("EB" if "EB" in chiso.GetName() else "EE")
    data    = chiso.Clone("data_" +EBEE+"_"+bin)
    can = ROOT.TCanvas("CHiso_"+EBEE+"_fit_"+bin)
    can.SetTopMargin(0.08)
    data.SetMarkerStyle(20)
    data.Draw("PE0")
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
        fitter.GetResult(0, prompt_val, prompt_err)
        fitter.GetResult(1, fake_val, fake_err)
        #print "prompt: "+str(prompt_val)+" fake_value: "+str(fake_val)
        result = fitter.GetPlot();
        h_prompt = mctemp1.Clone("prompt_"+EBEE+"_"+bin)
        h_fake   = mctemp2.Clone("fake_"+EBEE+"_"+bin)
        h_prompt.Scale(prompt_val*result.Integral(1,20)/h_prompt.Integral(1,20))
        h_fake  .Scale(fake_val  *result.Integral(1,20)/h_fake  .Integral(1,20))
        h_prompt.SetLineColor(3)
        h_prompt.SetFillColor(3)
        h_fake  .SetLineColor(2)
        h_fake  .SetFillColor(2)
        #h_prompt.Draw("same hist")
        #h_fake  .Draw("same hist")
        stack = ROOT.THStack("stack_"+EBEE+"_"+bin,"")
        stack.Add(h_fake)
        stack.Add(h_prompt)
        stack.Draw("SAME HIST")
        result.SetLineColor(4)
        data.Draw("SAME PE0")
        #result.Draw("SAME")
        leg = ROOT.TLegend(0.2,0.70,0.9,0.9)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)
        leg.SetHeader("Bin: "+bin+", "+EBEE)
        leg.AddEntry(data,     "Data",           "PE")
        leg.AddEntry(h_prompt, "Prompt photons", "F")
        leg.AddEntry(h_fake,   "Fake photons",   "F")
        #leg.AddEntry(result,   "Total Fit",      "L")
        leg.Draw("SAME")
        lat1 = ROOT.TLatex(4, 0.5*data.GetMaximum(), ("Purity: %.2f #pm %.2f" % (prompt_val, prompt_err)))
        lat1.Draw("SAME")
        lat2 = ROOT.TLatex(4, 0.4*data.GetMaximum(), ("MC Purity: %.2f" % (mcpurity)))
        lat2.Draw("SAME")
        add_cms_era(can, ERA+1 if "WAna" in opt.box else ERA+2, keep, opt.energy, lumi, prefix)
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

def calc_syst_err(vh, name, additional_syst=0.0, add_stat=True):
    scen_1_drop = ["extrap", "trigger", "ak8scale", "doubleratio",
                   "elefastsim", "muonfastsim","btagfastsim", "wtagfastsim",
                   "wmistagfastsim", "toptagfastsim", "topmistagfastsim", "leptonest"]
    nominal = vh[0].Clone(name)
    # Add the statistical error
    syst_up = vh[0].Clone(name+"_up")
    syst_dn = vh[0].Clone(name+"_dn")
    for binx in range(1, nominal.GetNbinsX()+1):
        err_up  = nominal.GetBinError(binx) ** 2 if add_stat else 0
        err_dn  = nominal.GetBinError(binx) ** 2 if add_stat else 0
        err_up += (additional_syst*nominal.GetBinContent(binx)) ** 2
        err_dn += (additional_syst*nominal.GetBinContent(binx)) ** 2
        for i in range(1, len(vh)):
            consider_syst = True
            if opt.scenario==0:
                consider_syst = False
            if opt.scenario==1:
                for syst_to_drop in scen_1_drop:
                    if syst_to_drop in vh[i].GetName():
                        consider_syst = False
            if consider_syst:
                syst_err = vh[i].GetBinContent(binx)-nominal.GetBinContent(binx)
                if syst_err>0:
                    err_up += syst_err ** 2
                else:
                    err_dn += syst_err ** 2
        err_up = err_up ** 0.5
        err_dn = err_dn ** 0.5
        syst_up.SetBinError(binx, err_up)
        syst_dn.SetBinError(binx, err_dn)
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
    #h_new = ROOT.TH2D(h.GetName()+"_2D", ";M_{R} (GeV);R^{2}", len(mrbins)-1,0,len(mrbins)-1, len(r2bins)-1,0,len(r2bins)-1)
    h_new = ROOT.TH2D(h.GetName()+"_2D", ";M_{R} (GeV);R^{2}", len(mrbins)-1,array('d',mrbins), len(r2bins)-1,array('d',r2bins))
    h_new.SetDirectory(0)
    for i in range(h.GetNbinsX()):
        bin_mr = (i/(len(r2bins)-1))+1
        bin_r2 = (i%(len(mrbins)-1))+1
        h_new.SetBinContent(bin_mr, bin_r2, h.GetBinContent(i+1))
        h_new.SetBinError  (bin_mr, bin_r2, h.GetBinError  (i+1))
    return h_new

# Make transfer factor plots with systematics
def make_kappa_plot(process, sr, cr, combine_bins, additional_syst):
    # Unroll and merge bins for nominal distributions
    if sr[0].InheritsFrom("TH2"):
        sr_1d = unroll(sr[0])
    else:
        sr_1d = sr[0]
    if cr[0].InheritsFrom("TH2"):
        cr_1d = unroll(cr[0])
    else:
        cr_1d = cr[0]
    h_sr            = sr_1d.Clone()
    h_cr            = cr_1d.Clone()
    if combine_bins:
        h_sr = combinebins(h_sr, h_sr.GetName()+"_combined")
        h_cr = combinebins(h_cr, h_cr.GetName()+"_combined")
    nbin = h_sr.GetNbinsX()
    srname = sr_1d.GetName().split("_")[1].replace("s","S'").replace("q","Q'")
    crname = cr_1d.GetName().split("_")[1]
    h_kappa = ROOT.TH1D("kappa_"+srname+"_"+crname, ";;"+srname+"/"+crname+" transfer factor", nbin,0,nbin)
    add_bin_labels(h_kappa, combine_bins)
    h_kappa.Divide(h_sr, h_cr)
    # Calculate systematic error interval
    vh_kappa = []
    for i in range(1, len(sr)):
        # Unroll and merge bins similar to nominal
        if sr[i].InheritsFrom("TH2"):
            sr_syst_1d = unroll(sr[i])
        else:
            sr_syst_1d = sr[i]
        if cr[i].InheritsFrom("TH2"):
            cr_syst_1d = unroll(cr[i])
        else:
            cr_syst_1d = cr[i]
        h_sr_syst = sr_syst_1d.Clone()
        h_cr_syst = cr_syst_1d.Clone()
        if combine_bins:
            h_sr_syst = combinebins(h_sr_syst, h_sr_syst.GetName()+"_combined")
            h_cr_syst = combinebins(h_cr_syst, h_cr_syst.GetName()+"_combined")
        h_kappa_syst = ROOT.TH1D("kappa_"+srname+"_"+crname+systematics[i], ";;"+srname+"/"+crname+" transfer factor", nbin,0,nbin)
        h_kappa_syst.Divide(h_sr, h_cr)
        vh_kappa.append(h_kappa_syst)
    # Calculate total (systematic+stat) error for each bin
    for binx in range(1, h_kappa.GetNbinsX()+1):
        nom = h_kappa.GetBinContent(binx)
        # Stat error + optional additional flat systematic
        err_up = h_kappa.GetBinError(binx) ** 2 + ((h_kappa.GetBinContent(binx)*additional_syst) ** 2)
        err_dn = h_kappa.GetBinError(binx) ** 2 + ((h_kappa.GetBinContent(binx)*additional_syst) ** 2)
        for i in range(len(vh_kappa)):
            syst_err = nom-vh_kappa[i].GetBinContent(binx)
            if syst_err>0:
                err_up += syst_err ** 2
            else:
                err_dn += syst_err ** 2
        err_up = err_up ** 0.5
        err_dn = err_dn ** 0.5
        # Move central value to middle of the full interval of the total error
        h_kappa.SetBinContent(binx, nom+(err_up-err_dn)/2)
        h_kappa.SetBinError  (binx,     (err_up+err_dn)/2)
    # Zero negative kappas
    #for binx in range(1, h_sr.GetNbinsX()+1):
    #    if h_kappa.GetBinContent(binx)<0:
    #        h_kappa.SetBinContent(binx,0)
    #        h_kappa.SetBinError(binx,0)
    # Draw kappa comparison plot
    processname = process.split(",")[0].replace(" (G)","").replace(" (L)","")
    can_tf = custom_can(h_kappa, "kappa_"+processname+"_"+srname+"_"+crname+"_"+opt.box, 0,0, 500,500, 90,20,50)
    h_kappa.SetMarkerColor(633)
    h_kappa.SetLineColor  (633)
    h_kappa.Draw("PE0")
    h_kappa.GetXaxis().SetLabelSize(0.04)
    leg = ROOT.TLegend(0.3,0.7,0.6,0.85, "")
    leg.AddEntry(h_kappa, "#color[633]{"+process+"}", "pe")
    leg.SetTextSize(0.04)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw("SAME")
    ROOT.gPad.Update()
    draw_mr_bins([h_kappa], 0, h_kappa.GetMaximum(), combine_bins, keep, mrbins_TeV, r2bins)
    add_cms_era(can_tf, ERA, keep, opt.energy, lumi, prefix)
    save_plot(can_tf, "", "Plots/kappa/"+DATE+"/"+can_tf.GetName())

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
    if not os.path.exists("syst_"+opt.dir+"/hadd/gjets.root"):
        gjets_files = []
        for name in gjets: gjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/gjets.root"]+gjets_files,             "syst_"+opt.dir+"/hadd/log/gjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/dyjets.root"):
        dyjets_files = []
        for name in dyjets: dyjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/dyjets.root"]+dyjets_files,           "syst_"+opt.dir+"/hadd/log/dyjets.log")
    if not os.path.exists("syst_"+opt.dir+"/hadd/nondyjets.root"):
        nondyjets_files = []
        for name in nondyjets: nondyjets_files.append(opt.dir+"/hadd/"+name+".root")
        logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/nondyjets.root"]+nondyjets_files,     "syst_"+opt.dir+"/hadd/log/nondyjets.log")

    # also make sure example signals are hadded
    for signal_name, dirs in example_signals.iteritems():
        if not os.path.exists("syst_"+opt.dir+"/hadd/signal_"+signal_name+".root"):
            signal_files = []
            for name in dirs: signal_files.append(opt.dir+"/hadd/"+name+".root")
            if len(signal_files)>1:
                logged_call(["hadd", "-f", "-v", "syst_"+opt.dir+"/hadd/signal_"+signal_name+".root"]+signal_files, "syst_"+opt.dir+"/hadd/log/signal_"+signal_name+".log")
            else:
                logged_call(["cp", "-p"]+signal_files+["syst_"+opt.dir+"/hadd/signal_"+signal_name+".root"], "syst_"+opt.dir+"/hadd/log/signal_"+signal_name+".log")
            input_files = []
            for sig in dirs: input_files += glob.glob(ntuple+"/"+sig+"/*.root")
            for i in range(len(input_files)):
                f = ROOT.TFile.Open(input_files[i])
                if "T2tt" in signal_name:
                    h = f.Get("npvLowHigh_T2tt")
                else:
                    h = f.Get("npvLowHigh_T1tttt")
                if i==0:
                    npvLowHighHist_allevt = h.Clone(h.GetName()+"_allevt")
                    npvLowHighHist_allevt.SetDirectory(0)
                else:
                    npvLowHighHist_allevt.Add(h)
                f.Close()
            f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+signal_name+".root","UPDATE")
            npvLowHighHist_allevt.Write()
            f.Close()

# acceptance extrapolation
if "T2tt" in opt.model:
    highest_mass = 1200
    if opt.energy == 14:
        extension_points = range(1250, 2550, 50)
    elif opt.energy == 27:
        extension_points = range(1250, 2050, 50)
else:
    highest_mass = 2300
    extension_points = range(2350, 3550, 50)

if opt.box=="WAna_nj45":
    postfix = "_NJet4to5"
    facc = "results/Plotter_out_2018_08_08_syst_replot.root"
elif opt.box=="WAna_nj6":
    postfix = "_NJet6"
    facc = "results/Plotter_out_2018_08_08_syst_replot.root"
elif opt.box=="TopAna":
    postfix = ""
    facc = "results/Plotter_out_2018_06_08_syst_TopAna_replot.root"
if opt.model == "T2tt":
    subdir = "SignalSelectionEfficiency_vs_MLSP_vs_MStop"
else:
    subdir = "SignalSelectionEfficiency_vs_MLSP_vs_MGluino"

fin = ROOT.TFile.Open(facc)
h_acc = fin.Get(subdir+"/"+opt.model+postfix)
h_acc.SetDirectory(0)
keep.append(h_acc)
fin.Close()

def fit_turnon_improved_(h):
    if (h.GetEntries()>0):
        # Estimate Parameters
        p0 = 0.
        first_point = 0.
        for i in range(h.GetNbinsX(),0,-1):
            c = h.GetBinContent(i)
            if c>0:
                if first_point == 0:
                    first_point = h.GetBinCenter(i)
                    p0 = c
        p1 = h.GetBinContent(h.GetMaximumBin()) - p0
        bin1 = 0
        bin2=9999
        mid = (p0 + p1) / 2.0
        for i in range(h.GetNbinsX(),0,-1):
            cont = h.GetBinContent(i)
            if cont<mid and cont>0 and bin1==0: bin1 =i
            if cont>mid: bin2 = i
        x1 = h.GetBinCenter(bin1)
        x2 = h.GetBinCenter(bin2)
        y1 = h.GetBinContent(bin1)
        y2 = h.GetBinContent(bin2)
        p2 = x1 + (x2-x1)*(mid-y1)/(y2-y1)
        # Fit
        f_fit = ROOT.TF1(h.GetName()+"_turnon", "[0]+[1]/(1+exp(([2]-x)/[3]))")
        f_fit.SetParameter(0,p0)
        f_fit.SetParameter(1,p1)
        f_fit.SetParameter(2,p2)
        f_fit.SetParameter(3,400)
        f_fit.SetParLimits(0,0,0.2)
        f_fit.SetParLimits(1,0.9*p1,1.1*p1)
        #f_fit.SetParLimits(2,0,1)
        #f_fit.SetParLimits(3,0.25,15)
        f_fit.SetLineColor(h.GetMarkerColor())
        #std::cout<<h.GetName()<<" "<<p0<<" "<<p1<<" "<<p2<<std::endl
        h.Fit(h.GetName()+"_turnon", "MQ")
        #double Chi2NDoF = f_fit.GetChisquare()/f_fit.GetNDF()
        #double Const = f_fit.GetParameter(0)
        #double Rise = f_fit.GetParameter(1)
        #double Mid = f_fit.GetParameter(2)
        #double Width = f_fit.GetParameter(3)
        #std::cout<<printf( "Chi2/NDoF: %10.1f  Const: %1.3f  Rise: %1.3f  Mid: %3.1f  Width: %2.2lf", Chi2NDoF, Const, Rise, Mid, Width)<<std::endl
        # gPad.Update()
        # if (fabs(Mid-p2)>4) sleep(5)
        return h.GetFunction(h.GetName()+"_turnon")
    else: return 0

def turnon1_with_slope_(h, keep, xmin = 0, xmax = 3000, ymin=0, ymax=0.2):
    # Estimate Parameters
    p0 = 0.
    first_point = 0.
    for i in range(h.GetNbinsX(),0,-1):
        c = h.GetBinContent(i)
        if c>0:
            if first_point == 0:
                first_point = h.GetBinCenter(i)
                p0 = c
    p1 = h.GetBinContent(h.GetMaximumBin()) - p0
    bin1 = 0
    bin2=9999
    mid = 0.05
    for i in range(h.GetNbinsX(),0,-1):
        cont = h.GetBinContent(i)
        if cont<mid and cont>0 and bin1==0: bin1 =i
        if cont>mid: bin2 = i
    x1 = h.GetBinCenter(bin1)
    x2 = h.GetBinCenter(bin2)
    y1 = h.GetBinContent(bin1)
    y2 = h.GetBinContent(bin2)
    p2 = x1 + (x2-x1)*(mid-y1)/(y2-y1)
    maximum = h.GetMaximum()
    f_fit = ROOT.TF1(h.GetName()+"_turnon", "([0]+[1]*x)/(1+exp(([2]-x)/[3]))")
    f_fit.SetParameter(0,0.1)
    f_fit.SetParameter(1,1e-5)
    f_fit.SetParameter(2,p2)
    f_fit.SetParameter(3,100)
    f_fit.SetParLimits(0,0.9*maximum,1.1*maximum)
    #f_fit.SetParLimits(1,0,1e-4)
    f_fit.SetParLimits(1,5e-6,1e-5)
    #f_fit.SetParLimits(2,0,100)
    f_fit.SetParLimits(2,1.0*p2,1.3*p2)
    f_fit.SetParLimits(3,0,200)
    f_fit.SetLineWidth(2)
    h.Fit(h.GetName()+"_turnon","MQ")
    keep.append(f_fit)
    return f_fit

def turnon2_with_slope_(h, keep, mslp):
    # Estimate Parameters
    p0 = 0.
    first_point = 0.
    for i in range(h.GetNbinsX(),0,-1):
        c = h.GetBinContent(i)
        if c>0:
            if first_point == 0:
                first_point = h.GetBinCenter(i)
                p0 = c
    p1 = h.GetBinContent(h.GetMaximumBin()) - p0
    bin1 = 0
    bin2=9999
    mid = 0.05
    for i in range(h.GetNbinsX(),0,-1):
        cont = h.GetBinContent(i)
        if cont<mid and cont>0 and bin1==0: bin1 =i
        if cont>mid: bin2 = i
    x1 = h.GetBinCenter(bin1)
    x2 = h.GetBinCenter(bin2)
    y1 = h.GetBinContent(bin1)
    y2 = h.GetBinContent(bin2)
    p2 = x1 + (x2-x1)*(mid-y1)/(y2-y1)
    maximum = h.GetMaximum()
    f_fit = ROOT.TF1(h.GetName()+"_turnon", "([0]+[1]*x)/(1+exp(([2]-x)/[3]))")
    if mslp==0:
        f_fit.SetParameter(0,0.03)
        f_fit.SetParLimits(0,0.02,0.04)
        f_fit.SetParameter(1,2e-6)
        f_fit.SetParLimits(1,2e-6,2e-5)
        f_fit.SetParameter(2,700)
        f_fit.SetParLimits(2,600,800)
        f_fit.SetParameter(3,300)
        f_fit.SetParLimits(3,200,400)
    else:
        f_fit.SetParameter(0,0.1)
        if mlsp==400:
            f_fit.FixParameter(0,0.093)
        elif mlsp==500:
            f_fit.FixParameter(0,0.095)
        #f_fit.SetParLimits(0,0.85*maximum,1.1*maximum)
        if mlsp<300:
            f_fit.FixParameter(1,3e-6)
        else:
            f_fit.FixParameter(1,6e-6)
        f_fit.SetParameter(2,p2)
        f_fit.SetParLimits(2,0.9*p2,1.3*p2)
        if mlsp>=400 and mlsp<=500:
            f_fit.FixParameter(3,150)
        if mlsp>=600 and mlsp<=700:
            f_fit.FixParameter(3,125)
        else:
            f_fit.SetParameter(3,150)
            f_fit.SetParLimits(3,100,300)
    f_fit.SetLineWidth(2)
    h.Fit(h.GetName()+"_turnon","MQ")
    keep.append(f_fit)
    return f_fit

# Acceptance fits
# T2tt - Wn45 --> flat (use no fit)
# T2tt - Wn6  --> decreasing
# T2tt - Top  --> increasing, use turnon fit
# T1tttt - Wn45 --> flat (use no fit)
# T1tttt - Wn6  --> deacreasing (below lsp 1300), flat (above)
# T1tttt - Top  --> increasing
# T5ttcc - Wn45 --> decreasing
# T5ttcc - Wn6  --> decreasing
# T5ttcc - Top  --> increasing, use turnon fit
f_acc_fits = ROOT.TFile.Open("acceptance_"+opt.model+"_"+opt.box+".root","RECREATE")
binx_high = h_acc.GetXaxis().FindBin(highest_mass)
vf_acc = {}
for biny in range(1, h_acc.GetNbinsY()+1):
    if h_acc.GetBinContent(binx_high, biny)>0:
        mlsp = int(h_acc.GetYaxis().GetBinCenter(biny))
        h_slice = h_acc.ProjectionX("_"+str(mlsp),biny,biny)
        h_ext = ROOT.TH1D("acc_"+opt.model+"_lsp"+str(mlsp),";"+h_acc.GetXaxis().GetTitle()+";Acceptance",401,-12.5,10012.5);
        for binx in range(1, h_slice.GetNbinsX()+1):
            h_ext.SetBinContent(binx,h_slice.GetBinContent(binx))
            h_ext.SetBinError  (binx,h_slice.GetBinError  (binx))
        h_ext.GetXaxis().SetRangeUser(0,extension_points[-1])
        h_ext.SetMarkerStyle(20)
        can = custom_can(h_ext, "acc_"+str(mlsp))
        h_ext.Draw("PE1")
        turnon1 = False
        turnon2 = False
        func = "pol1"
        if opt.model == "T5ttcc":
            if "WAna_nj45" in opt.box:
                rangemin = min(2000,max(1200,mlsp+800))
            elif "WAna_nj6" in opt.box:
                rangemin = max(1000,mlsp+600)
            elif "TopAna" in opt.box:
                #func = "pol0"
                #if mlsp>1200:
                #    rangemin = highest_mass-100
                #elif mlsp>1000:
                #    rangemin = highest_mass-200
                #else:
                #    rangemin = highest_mass-300
                turnon2 = True
        elif opt.model == "T1tttt":
            if "WAna_nj45" in opt.box:
                func = "pol0"
                if mlsp>1500:
                    rangemin = highest_mass-150
                else:
                    rangemin = highest_mass-200
            elif "WAna_nj6" in opt.box:
                rangemin = min(2300,max(1300,mlsp+900))
            elif "TopAna" in opt.box:
                rangemin = min(1800, max(1300,mlsp+900))
        elif opt.model == "T2tt":
            if "WAna_nj45" in opt.box:
                func = "pol0"
                if mlsp>300:
                    rangemin = highest_mass-100
                else:
                    rangemin = highest_mass-200
            elif "WAna_nj6" in opt.box:
                rangemin = max(800,800+(mlsp-300)/2)
            elif "TopAna" in opt.box:
                turnon1 = True
                #func = "[0]+[1]/(1+exp(([2]-x)/[3]))"
                #rangemin = 0
                #func = "log(x*[0]+[1])+[2]"
                #rangemin = max(800,800+(mlsp-300)/2)
                #func = "pol0"
                #rangemin = highest_mass
        if turnon1:            
            #f_fit = fit_turnon_improved_(h_ext)
            f_fit = turnon1_with_slope_(h_ext, keep)
        elif turnon2:
            f_fit = turnon2_with_slope_(h_ext, keep, mlsp)
        else:
            f_fit = ROOT.TF1("fit_"+str(mlsp), func, rangemin,extension_points[-1])
            if opt.model == "T2tt" and "WAna_nj6" in opt.box and mlsp>400:
                f_fit.FixParameter(1,-4e-5)
            if opt.model == "T1tttt" and "WAna_nj6" in opt.box:
                f_fit.SetParLimits(1,-7e-5,-5e-6)
                if mlsp>=1000:
                    f_fit.FixParameter(1,-6.1e-6)
            h_ext.Fit("fit_"+str(mlsp), "RMQ")
            f_fit.SetLineWidth(2)
            f_fit.Draw("SAME")
            keep.append(f_fit)
        vf_acc[mlsp] = f_fit
        can.Write()
        keep.append(can)
        keep.append(h_ext)

f_acc_fits.Close()

#sys.exit()

# ----------------- Harvest histograms -------------------

# HL-LHC study
Run2_dir = "results/run_2018_06_08_syst_TopAna"
if not "TopAna" in opt.box: Run2_dir = "results/run_2018_08_08_syst"

tmp_dir = opt.dir # Somehow G region triggers were messed up recently, so revert back to previous good version
#if "WAna" in opt.box: tmp_dir = "results/run_2018_06_08_syst_NoTopVeto" 

# Get relative xsec and luminsoity ratios
rel_scales = {}
# Get the 13 TeV Cross sections first
for line in open('common/BackGroundXSec.txt','r'):
    line = line.replace('\n','')
    shortname = line.split()[0]
    xsec = float(line.split()[2])
    rel_scales[shortname] = xsec
# Then take the ratio with the 14/27 TeV one and multiply by relative lumi
if opt.energy != 14 and opt.energy != 27:
    print ("ERROR: Invalid energy, there are no cross-sections available for specified energy. Use -e 14 or -e 27")
    sys.exit()
for line in open('common/BackGroundXSec_'+str(opt.energy)+'TeV.txt','r'):
    line = line.replace('\n','')
    shortname = line.split()[0]
    xsec = float(line.split()[2])
    xsec_Run2 = rel_scales[shortname]
    rel_scales[shortname] = xsec / xsec_Run2 * lumi / lumi_Run2

# Do the same for signal
rel_scales_gluino = {}
rel_scales_stop   = {}
# Get the 13 TeV Cross sections first
for line in open('data/gluino13TeV.txt', 'r'):
    line = line.replace('\n','')
    mass = int(line.split(',')[0])
    xsec = float(line.split(',')[1])
    rel_scales_gluino[mass] = xsec
for line in open('data/stop13TeV.txt', 'r'):
    line = line.replace('\n','')
    mass = int(line.split(',')[0])
    xsec = float(line.split(',')[1])
    rel_scales_stop[mass] = xsec
# Then take the ratio with the 14/27 TeV one and multiply by relative lumi
for line in open('data/gluino'+str(opt.energy)+'TeV.txt', 'r'):
    line = line.replace('\n','')
    mass = int(line.split(',')[0])
    if mass<=3000:
        xsec = float(line.split(',')[1])
        xsec_Run2 = rel_scales_gluino[mass]
        rel_scales_gluino[mass] = xsec / xsec_Run2 * lumi / lumi_Run2
for line in open('data/stop'+str(opt.energy)+'TeV.txt', 'r'):
    line = line.replace('\n','')
    mass = int(line.split(',')[0])
    if mass<=2000:
        xsec = float(line.split(',')[1])
        xsec_Run2 = rel_scales_stop[mass]
        rel_scales_stop[mass] = xsec / xsec_Run2 * lumi / lumi_Run2

rel_scales_signal = rel_scales_gluino
if "T2" in opt.model: rel_scales_signal = rel_scales_stop

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
S_data = load(f,"MRR2_S_data"+BIN,"_data", combine_bins)
f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/data.root")
Q_data = load(f,"MRR2_Q_data"+BIN,"_data", 0)
W_data = load(f,"MRR2_W_data"+BIN,"_data", 0)
L_data = load(f,"MRR2_L_data"+BIN,"_data", 0)
T_data = load(f,"MRR2_T_data"+BIN,"_data", 0)
s_data = load(f,"MRR2_s_data"+BIN,"_data", combine_bins)
q_data = load(f,"MRR2_q_data"+BIN,"_data", combine_bins)
npvHist = load(f,"nvtx","_data")
npvHist.Scale(1/npvHist.Integral())

##PREV    # Signal
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root")
##PREV    S_signal = []
##PREV    counter = 0
##PREV    for ikey in range(0, f.GetListOfKeys().GetEntries()):
##PREV        name = f.GetListOfKeys().At(ikey).GetName()
##PREV        if name.startswith("MRR2_S_signal") and not "Up" in name and not "Down" in name:
##PREV            if not "_nj" in name:
##PREV                counter+=1
##PREV                S_syst = []
##PREV                for syst in systematics:
##PREV                    S_syst.append(load(f, name+BIN+syst, "_sig", combine_bins))
##PREV                S_signal.append(S_syst)
##PREV        if opt.TEST>0:
##PREV            if counter==opt.TEST:
##PREV                break
##PREV    # Histos for pileup acceptance systematic
##PREV    if "T2tt" in opt.model:
##PREV        npvLowHighHist        = loadclone(f,"npvLowHigh_T2tt","_sig")
##PREV        npvLowHighHist_allevt = loadclone(f,"npvLowHigh_T2tt_allevt","_sig")
##PREV    else:
##PREV        npvLowHighHist        = loadclone(f,"npvLowHigh_T1tttt","_sig")
##PREV        npvLowHighHist_allevt = loadclone(f,"npvLowHigh_T1tttt_allevt","_sig")
##PREV    # Merge statistics in Mglu/Mstop vs Mlsp
##PREV    npvLowHighHist       .Rebin3D(4,4,1)
##PREV    npvLowHighHist_allevt.Rebin3D(4,4,1)
##PREV    
##PREV    # Background
##PREV    # top + ttbar
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ttbar.root")
##PREV    S_TT = []
##PREV    S_LooseWP_TT = []
##PREV    for syst in systematics:
##PREV        S_TT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_TT", 0))
##PREV        S_LooseWP_TT.append(load(f, "MRR2_S_LooseWP_bkg"+BIN+syst,"_TT", 0))
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/top.root")
##PREV    for i in range(len(systematics)):
##PREV        # Fix problem with nonexistent scale weights for single top
##PREV        syst = systematics[i]
##PREV        if "scale" in syst: syst = ""
##PREV        S_TT[i].Add(load(f,"MRR2_S_bkg"+BIN+syst,"_T", 0))
##PREV        S_LooseWP_TT[i].Add(load(f, "MRR2_S_LooseWP_bkg"+BIN+syst,"_T", 0))
##PREV    f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/ttbar.root")
##PREV    Q_TT = []
##PREV    W_TT = []
##PREV    L_TT = []
##PREV    T_TT = []
##PREV    s_TT = []
##PREV    q_TT = []
##PREV    #s_LooseWP_TT = []
##PREV    for syst in systematics:
##PREV        Q_TT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_TT", 0))
##PREV        W_TT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_TT", 0))
##PREV        L_TT.append(load(f,"MRR2_L_bkg"+BIN+syst,"_TT", 0))
##PREV        T_TT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_TT", 0))
##PREV        s_TT.append(load(f,"MRR2_s_bkg"+BIN+syst,"_TT", 0))
##PREV        q_TT.append(load(f,"MRR2_q_bkg"+BIN+syst,"_TT", 0))
##PREV        #s_LooseWP_TT.append(load(f, "MRR2_s_LooseWP_bkg"+BIN+syst,"_TT", 0))
##PREV    f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/top.root")
##PREV    for i in range(len(systematics)):
##PREV        # Fix problem with nonexistent scale weights for single top
##PREV        syst = systematics[i]
##PREV        if "scale" in syst: syst = ""
##PREV        Q_TT[i].Add(load(f,"MRR2_Q_bkg"+BIN+syst,"_T", 0))
##PREV        W_TT[i].Add(load(f,"MRR2_W_bkg"+BIN+syst,"_T", 0))
##PREV        L_TT[i].Add(load(f,"MRR2_L_bkg"+BIN+syst,"_T", 0))
##PREV        T_TT[i].Add(load(f,"MRR2_T_bkg"+BIN+syst,"_T", 0))
##PREV        s_TT[i].Add(load(f,"MRR2_s_bkg"+BIN+syst,"_T", 0))
##PREV        q_TT[i].Add(load(f,"MRR2_q_bkg"+BIN+syst,"_T", 0))
##PREV        #s_LooseWP_TT[i].Add(load(f, "MRR2_s_LooseWP_bkg"+BIN+syst,"_T", 0))
##PREV    # multijet
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/multijet.root")
##PREV    S_MJ = []
##PREV    S_LooseWP_MJ = []
##PREV    for syst in systematics:
##PREV        S_MJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_MJ", 0))
##PREV        S_LooseWP_MJ.append(load(f, "MRR2_S_LooseWP_bkg"+BIN+syst,"_MJ", 0))
##PREV    f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/multijet.root")
##PREV    Q_MJ = []
##PREV    W_MJ = []
##PREV    L_MJ = []
##PREV    T_MJ = []
##PREV    s_MJ = []
##PREV    q_MJ = []
##PREV    #s_LooseWP_MJ = []
##PREV    for syst in systematics:
##PREV        Q_MJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_MJ", 0))
##PREV        W_MJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_MJ", 0))
##PREV        L_MJ.append(load(f,"MRR2_L_bkg"+BIN+syst,"_MJ", 0))
##PREV        T_MJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_MJ", 0))
##PREV        s_MJ.append(load(f,"MRR2_s_bkg"+BIN+syst,"_MJ", 0))
##PREV        q_MJ.append(load(f,"MRR2_q_bkg"+BIN+syst,"_MJ", 0))
##PREV        #s_LooseWP_MJ.append(load(f, "MRR2_s_LooseWP_bkg"+BIN+syst,"_MJ", 0))
##PREV    # wjets
##PREV    f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/wjets.root")
##PREV    S_WJ = []
##PREV    S_LooseWP_WJ = []
##PREV    for syst in systematics:
##PREV        S_WJ.append(load(f,"MRR2_S_bkg"+BIN+syst,"_WJ", 0))
##PREV        S_LooseWP_WJ.append(load(f, "MRR2_S_LooseWP_bkg"+BIN+syst,"_WJ", 0))
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/wjets.root")
##PREV    Q_WJ = []
##PREV    W_WJ = []
##PREV    L_WJ = []
##PREV    T_WJ = []
##PREV    L_WJ = []
##PREV    s_WJ = []
##PREV    q_WJ = []
##PREV    #s_LooseWP_WJ = []
##PREV    for syst in systematics:
##PREV        Q_WJ.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_WJ", 0))
##PREV        W_WJ.append(load(f,"MRR2_W_bkg"+BIN+syst,"_WJ", 0))
##PREV        L_WJ.append(load(f,"MRR2_L_bkg"+BIN+syst,"_WJ", 0))
##PREV        T_WJ.append(load(f,"MRR2_T_bkg"+BIN+syst,"_WJ", 0))
##PREV        s_WJ.append(load(f,"MRR2_s_bkg"+BIN+syst,"_WJ", 0))
##PREV        q_WJ.append(load(f,"MRR2_q_bkg"+BIN+syst,"_WJ", 0))
##PREV        #s_LooseWP_WJ.append(load(f, "MRR2_s_LooseWP_bkg"+BIN+syst,"_WJ", 0))
##PREV    # ztoinv
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/ztoinv.root")
##PREV    S_ZI = []
##PREV    S_LooseWP_ZI = []
##PREV    for syst in systematics:
##PREV        S_ZI.append(load(f,"MRR2_S_bkg"+BIN+syst,"_ZI", 0))
##PREV        S_LooseWP_ZI.append(load(f, "MRR2_S_LooseWP_bkg"+BIN+syst,"_ZI", 0))
##PREV    f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/ztoinv.root")
##PREV    Q_ZI = []
##PREV    W_ZI = []
##PREV    L_ZI = []
##PREV    T_ZI = []
##PREV    s_ZI = []
##PREV    q_ZI = []
##PREV    #s_LooseWP_ZI = []
##PREV    for syst in systematics:
##PREV        Q_ZI.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_ZI", 0))
##PREV        W_ZI.append(load(f,"MRR2_W_bkg"+BIN+syst,"_ZI", 0))
##PREV        L_ZI.append(load(f,"MRR2_L_bkg"+BIN+syst,"_ZI", 0))
##PREV        T_ZI.append(load(f,"MRR2_T_bkg"+BIN+syst,"_ZI", 0))
##PREV        s_ZI.append(load(f,"MRR2_s_bkg"+BIN+syst,"_ZI", 0))
##PREV        q_ZI.append(load(f,"MRR2_q_bkg"+BIN+syst,"_ZI", 0))
##PREV        #s_LooseWP_ZI.append(load(f, "MRR2_s_LooseWP_bkg"+BIN+syst,"_ZI", 0))
##PREV    # other
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/other.root")
##PREV    S_OT = []
##PREV    S_LooseWP_OT = []
##PREV    for syst in systematics:
##PREV        S_OT.append(load(f,"MRR2_S_bkg"+BIN+syst,"_OT", 0))
##PREV        S_LooseWP_OT.append(load(f, "MRR2_S_LooseWP_bkg"+BIN+syst,"_OT", 0))
##PREV    f = ROOT.TFile.Open("syst_"+tmp_dir+"/hadd/other.root")
##PREV    Q_OT = []
##PREV    W_OT = []
##PREV    L_OT = []
##PREV    T_OT = []
##PREV    s_OT = []
##PREV    q_OT = []
##PREV    #s_LooseWP_OT = []
##PREV    for syst in systematics:
##PREV        Q_OT.append(load(f,"MRR2_Q_bkg"+BIN+syst,"_OT", 0))
##PREV        W_OT.append(load(f,"MRR2_W_bkg"+BIN+syst,"_OT", 0))
##PREV        L_OT.append(load(f,"MRR2_L_bkg"+BIN+syst,"_OT", 0))
##PREV        T_OT.append(load(f,"MRR2_T_bkg"+BIN+syst,"_OT", 0))
##PREV        s_OT.append(load(f,"MRR2_s_bkg"+BIN+syst,"_OT", 0))
##PREV        q_OT.append(load(f,"MRR2_q_bkg"+BIN+syst,"_OT", 0))
##PREV        #s_LooseWP_OT.append(load(f, "MRR2_s_LooseWP_bkg"+BIN+syst,"_OT", 0))
##PREV    
##PREV    # Selected signals for Results plot
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_T1tttt.root")
##PREV    S_T1tttt = load(f, "MRR2_S_signal_2000_300"+BIN, "_sig", 0)
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_T2tt.root")
##PREV    S_T2tt = load(f, "MRR2_S_signal_1200_100"+BIN, "_sig", 0)
##PREV    f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_T5ttcc.root")
##PREV    S_T5ttcc = load(f, "MRR2_S_signal_2000_300"+BIN, "_sig", 0)




# Signal
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root")
S_signal = []
counter = 0
for ikey in range(0, f.GetListOfKeys().GetEntries()):
    name = f.GetListOfKeys().At(ikey).GetName()
    if name.startswith("MRR2_S_signal") and not "Up" in name and not "Down" in name:
        if not "_nj" in name:
            counter+=1
            mass = int(name.split("_")[-2])
            S_syst = []
            for syst in systematics:
                S_syst.append(load_and_scale_signal_1d(f, rel_scales_signal[mass], keep, name+BIN+syst, "_sig", combine_bins))
            S_signal.append(S_syst)
    if opt.TEST>0:
        if counter==opt.TEST:
            break
# Histos for pileup acceptance systematic
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_"+opt.model+".root")
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
S_TT = []
S_LooseWP_TT = []
for syst in systematics:
    S_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_S_bkg"+BIN+syst,"_TT", 0))
    S_LooseWP_TT.append(load_and_scale(ttbar, rel_scales, keep, "MRR2_S_LooseWP_bkg"+BIN+syst,"_TT", 0))
for i in range(len(systematics)):
    # Fix problem with nonexistent scale weights for single top
    syst = systematics[i]
    if "scale" in syst: syst = ""
    S_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_S_bkg"+BIN+syst,"_T", 0))
    S_LooseWP_TT[i].Add(load_and_scale(top, rel_scales, keep, "MRR2_S_LooseWP_bkg"+BIN+syst,"_T", 0))
Q_TT = []
W_TT = []
L_TT = []
T_TT = []
s_TT = []
q_TT = []
#s_LooseWP_TT = []
for syst in systematics:
    Q_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_Q_bkg"+BIN+syst,"_TT", 0))
    W_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_W_bkg"+BIN+syst,"_TT", 0))
    L_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_L_bkg"+BIN+syst,"_TT", 0))
    T_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_T_bkg"+BIN+syst,"_TT", 0))
    s_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_s_bkg"+BIN+syst,"_TT", 0))
    q_TT.append(load_and_scale(ttbar, rel_scales, keep,"MRR2_q_bkg"+BIN+syst,"_TT", 0))
    #s_LooseWP_TT.append(load_and_scale(ttbar, rel_scales, keep, "MRR2_s_LooseWP_bkg"+BIN+syst,"_TT", 0))
for i in range(len(systematics)):
    # Fix problem with nonexistent scale weights for single top
    syst = systematics[i]
    if "scale" in syst: syst = ""
    Q_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_Q_bkg"+BIN+syst,"_T", 0))
    W_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_W_bkg"+BIN+syst,"_T", 0))
    L_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_L_bkg"+BIN+syst,"_T", 0))
    T_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_T_bkg"+BIN+syst,"_T", 0))
    s_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_s_bkg"+BIN+syst,"_T", 0))
    q_TT[i].Add(load_and_scale(top, rel_scales, keep,"MRR2_q_bkg"+BIN+syst,"_T", 0))
    #s_LooseWP_TT[i].Add(load_and_scale(top, rel_scales, keep, "MRR2_s_LooseWP_bkg"+BIN+syst,"_T", 0))
# multijet
S_MJ = []
S_LooseWP_MJ = []
for syst in systematics:
    S_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_S_bkg"+BIN+syst,"_MJ", 0))
    S_LooseWP_MJ.append(load_and_scale(multijet, rel_scales, keep, "MRR2_S_LooseWP_bkg"+BIN+syst,"_MJ", 0))
Q_MJ = []
W_MJ = []
L_MJ = []
T_MJ = []
s_MJ = []
q_MJ = []
#s_LooseWP_MJ = []
for syst in systematics:
    Q_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_Q_bkg"+BIN+syst,"_MJ", 0))
    W_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_W_bkg"+BIN+syst,"_MJ", 0))
    L_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_L_bkg"+BIN+syst,"_MJ", 0))
    T_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_T_bkg"+BIN+syst,"_MJ", 0))
    s_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_s_bkg"+BIN+syst,"_MJ", 0))
    q_MJ.append(load_and_scale(multijet, rel_scales, keep,"MRR2_q_bkg"+BIN+syst,"_MJ", 0))
    #s_LooseWP_MJ.append(load_and_scale(multijet, rel_scales, keep, "MRR2_s_LooseWP_bkg"+BIN+syst,"_MJ", 0))
# wjets
S_WJ = []
S_LooseWP_WJ = []
for syst in systematics:
    S_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_S_bkg"+BIN+syst,"_WJ", 0))
    S_LooseWP_WJ.append(load_and_scale(wjets, rel_scales, keep, "MRR2_S_LooseWP_bkg"+BIN+syst,"_WJ", 0))
Q_WJ = []
W_WJ = []
L_WJ = []
T_WJ = []
L_WJ = []
s_WJ = []
q_WJ = []
#s_LooseWP_WJ = []
for syst in systematics:
    Q_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_Q_bkg"+BIN+syst,"_WJ", 0))
    W_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_W_bkg"+BIN+syst,"_WJ", 0))
    L_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_L_bkg"+BIN+syst,"_WJ", 0))
    T_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_T_bkg"+BIN+syst,"_WJ", 0))
    s_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_s_bkg"+BIN+syst,"_WJ", 0))
    q_WJ.append(load_and_scale(wjets, rel_scales, keep,"MRR2_q_bkg"+BIN+syst,"_WJ", 0))
    #s_LooseWP_WJ.append(load_and_scale(wjets, rel_scales, keep, "MRR2_s_LooseWP_bkg"+BIN+syst,"_WJ", 0))
# ztoinv
S_ZI = []
S_LooseWP_ZI = []
for syst in systematics:
    S_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_S_bkg"+BIN+syst,"_ZI", 0))
    S_LooseWP_ZI.append(load_and_scale(ztoinv, rel_scales, keep, "MRR2_S_LooseWP_bkg"+BIN+syst,"_ZI", 0))
Q_ZI = []
W_ZI = []
L_ZI = []
T_ZI = []
s_ZI = []
q_ZI = []
#s_LooseWP_ZI = []
for syst in systematics:
    Q_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_Q_bkg"+BIN+syst,"_ZI", 0))
    W_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_W_bkg"+BIN+syst,"_ZI", 0))
    L_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_L_bkg"+BIN+syst,"_ZI", 0))
    T_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_T_bkg"+BIN+syst,"_ZI", 0))
    s_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_s_bkg"+BIN+syst,"_ZI", 0))
    q_ZI.append(load_and_scale(ztoinv, rel_scales, keep,"MRR2_q_bkg"+BIN+syst,"_ZI", 0))
    #s_LooseWP_ZI.append(load_and_scale(ztoinv, rel_scales, keep, "MRR2_s_LooseWP_bkg"+BIN+syst,"_ZI", 0))
# other
S_OT = []
S_LooseWP_OT = []
for syst in systematics:
    S_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_S_bkg"+BIN+syst,"_OT", 0))
    S_LooseWP_OT.append(load_and_scale(other, rel_scales, keep, "MRR2_S_LooseWP_bkg"+BIN+syst,"_OT", 0))
Q_OT = []
W_OT = []
L_OT = []
T_OT = []
s_OT = []
q_OT = []
#s_LooseWP_OT = []
for syst in systematics:
    Q_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_Q_bkg"+BIN+syst,"_OT", 0))
    W_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_W_bkg"+BIN+syst,"_OT", 0))
    L_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_L_bkg"+BIN+syst,"_OT", 0))
    T_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_T_bkg"+BIN+syst,"_OT", 0))
    s_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_s_bkg"+BIN+syst,"_OT", 0))
    q_OT.append(load_and_scale(other, rel_scales, keep,"MRR2_q_bkg"+BIN+syst,"_OT", 0))
    #s_LooseWP_OT.append(load_and_scale(other, rel_scales, keep, "MRR2_s_LooseWP_bkg"+BIN+syst,"_OT", 0))

# Selected signals for Results plot
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_T1tttt.root")
S_T1tttt = load_and_scale_signal_1d(f, rel_scales_gluino[2000], keep, "MRR2_S_signal_2000_300"+BIN, "_sig", 0)
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_T2tt.root")
S_T2tt   = load_and_scale_signal_1d(f, rel_scales_stop  [1200], keep, "MRR2_S_signal_1200_100"+BIN, "_sig", 0)
f = ROOT.TFile.Open("syst_"+opt.dir+"/hadd/signal_T5ttcc.root")
S_T5ttcc = load_and_scale_signal_1d(f, rel_scales_gluino[2000], keep, "MRR2_S_signal_2000_300"+BIN, "_sig", 0)

# ---------------- Z(nunu) estimate ---------------------

# ----------------- Photon Estimate ---------------------

print "Calculate photon based Z(nunu) estimate"

# Loading plots
# Templates for purity
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/gjets.root")
CHIsoTemplate_Prompt_EB = loadclone(f, "CHIsoTemplate_Prompt_g_EB", "_MC")
CHIsoTemplate_Prompt_EE = loadclone(f, "CHIsoTemplate_Prompt_g_EE", "_MC")
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/data.root")
CHIsoTemplate_Fake_EB = loadclone(f, "CHIsoTemplate_Fake_g_EB", "_data")
CHIsoTemplate_Fake_EE = loadclone(f, "CHIsoTemplate_Fake_g_EE", "_data")
# Distributions to fit
purity_in_Gm1 = False
if purity_in_Gm1:
    CHIso_GNoIso_EB = loadclone(f, "MR_R2_CHIso_gNoIso_EB"+BIN, "_data")
    CHIso_GNoIso_EE = loadclone(f, "MR_R2_CHIso_gNoIso_EE"+BIN, "_data")
else:
    CHIso_GNoIso_EB = loadclone(f, "MR_R2_CHIso_GNoIso_EB"+BIN, "_data")
    CHIso_GNoIso_EE = loadclone(f, "MR_R2_CHIso_GNoIso_EE"+BIN, "_data")
# Total photon counts
G_data_EB = loadclone(f, "MR_R2_G_EB"+BIN, "_data")
G_data_EE = loadclone(f, "MR_R2_G_EE"+BIN, "_data")

# Direct photon fraction
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/bkg.root")
IsDirect_G_EB = loadclone(f, "MR_R2_IsDirect_G_EB", "_MC")
IsDirect_G_EE = loadclone(f, "MR_R2_IsDirect_G_EE", "_MC")

# Double ratio
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/data.root")
G_data    = load(f, "MRR2_G_data"+BIN,"_data", 0)
Z_data    = load(f, "MRR2_Z_data","_data", 0)
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/nondyjets.root")
Z_NONDY   = load(f, "MRR2_Z_bkg", "_NONDY")
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/dyjets.root")
Z_DY      = load(f, "MRR2_Z_bkg", "_DY")

# Transfer factors
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/gjets.root")
GDirectPrompt_GJ = []
for syst in systematics:
    GDirectPrompt_GJ.append(load(f, "MRR2_G_DirectPrompt_bkg"+BIN+syst,"_MJ", 0))

# MC truth purity
f = ROOT.TFile.Open("results/Plotter_out_2018_05_29.root" if "WAna" in opt.box else "results/Plotter_out_2018_05_29_TopAna.root")
MCPurity_EB = loadclone(f, "PhotonPurity_vs_R2NoPho_vs_MRNoPho/IsPormpt_vs_R2NoPho_vs_MRNoPho/Background_G_Barrel", "")
MCPurity_EE = loadclone(f, "PhotonPurity_vs_R2NoPho_vs_MRNoPho/IsPormpt_vs_R2NoPho_vs_MRNoPho/Background_G_Endcap", "")

# --------------------------------
#       HL-LHC Data scaling

#load for scaling data(Run2 background)
f = ROOT.TFile.Open("syst_"+Run2_dir+"/hadd/bkg.root")
T_MC_Run2 = load(f,"MRR2_T_bkg"+BIN,"_Run2", 0)
Q_MC_Run2 = load(f,"MRR2_Q_bkg"+BIN,"_Run2", 0)
W_MC_Run2 = load(f,"MRR2_W_bkg"+BIN,"_Run2", 0)
G_MC_Run2 = load(f,"MRR2_G_bkg"+BIN,"_Run2", 0)
Z_MC_Run2 = load(f,"MRR2_Z_bkg",    "_Run2", 0)
L_MC_Run2 = load(f,"MRR2_L_bkg"+BIN,"_Run2", 0)
S_MC_Run2 = load(f,"MRR2_S_bkg"+BIN,"_Run2", 0)
q_MC_Run2 = load(f,"MRR2_q_bkg"+BIN,"_Run2", 0)
s_MC_Run2 = load(f,"MRR2_s_bkg"+BIN,"_Run2", 0)

T_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_T_bkg"+BIN,"_HL", 0)
Q_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_Q_bkg"+BIN,"_HL", 0)
W_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_W_bkg"+BIN,"_HL", 0)
G_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_G_bkg"+BIN,"_HL", 0)
Z_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_Z_bkg",    "_HL", 0)
L_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_L_bkg"+BIN,"_HL", 0)
S_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_S_bkg"+BIN,"_HL", 0)
q_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_q_bkg"+BIN,"_HL", 0)
s_MC = load_and_scale(bkg, rel_scales, keep,"MRR2_s_bkg"+BIN,"_HL", 0)

print "- Scaling Run2 data"

#data_HL = data_2016 * (MC_HL/MC_2016)
def scale_data(h_data, h_HL, h_Run2):
    #print h_data.GetName()
    h_tmp = h_data.Clone(h_data.GetName()+"_template")
    avgw_HL   = 0
    avgw_Run2 = 0
    for binx in range(1,h_tmp.GetNbinsX()+1):
        n_data = h_data.GetBinContent(binx)
        e_data = h_data.GetBinError  (binx)
        n_HL   = h_HL  .GetBinContent(binx)
        e_HL   = h_HL  .GetBinError  (binx)
        n_Run2 = h_Run2.GetBinContent(binx)
        e_Run2 = h_Run2.GetBinError  (binx)
        #if n_data == 0: n_data = 1.83 / 2.0
        if n_HL   == 0: n_HL   = 1.83 / 2.0 * avgw_HL
        if n_Run2 == 0: n_Run2 = 1.83 / 2.0 * avgw_Run2
        avgw_HL   = e_HL   ** 2 / n_HL
        avgw_Run2 = e_Run2 ** 2 / n_Run2
        h_tmp.SetBinContent(binx, n_data*(n_HL/n_Run2))
        #print("%2d  -  %.2f * (%.2f / %.2f) = %.2f" % (binx, n_data, n_HL, n_Run2, n_data*(n_HL/n_Run2)))
    #print ("Integral (HL-LHC) = %f" % (h_tmp .Integral()))
    #print ("Integral (Run2)   = %f" % (h_data.Integral()))
    #print ("Ratio             = %f" % (h_tmp.Integral()/h_data.Integral()))
    #print ("Lumi Ratio        = %f" % (lumi/lumi_Run2))
    h_data.Reset("ICESM")
    h_data.FillRandom(h_tmp,int(h_tmp.Integral()))

scale_data(T_data, T_MC, T_MC_Run2)
scale_data(Q_data, Q_MC, Q_MC_Run2)
scale_data(W_data, W_MC, W_MC_Run2)
scale_data(G_data, G_MC, G_MC_Run2)
#scale_data(Z_data, Z_MC, Z_MC_Run2)
scale_data(L_data, L_MC, L_MC_Run2)
scale_data(S_data, S_MC, S_MC_Run2)
scale_data(q_data, q_MC, q_MC_Run2)
scale_data(s_data, s_MC, s_MC_Run2)

print "- Scaling the Data to HL-LHC values done."

#sys.exit()

scale_down_factors = {
    "toppt"       : 0.3333,
    "isr"         : 0.5, 
    "alphas"      : 0.5,
    "facrenscale" : 0.5,
    "lostlep"     : 0.1,
    "jes"         : 0.5,
    "jer"         : 0.5,
    "met"         : 0.5,
    "btag"        : 0.5
}

# Scale down systematics by the specified factors
if opt.scenario == 1:
    for iSignal in range(len(S_signal)):
        scale_down_syst(S_signal[iSignal], scale_down_factors, True)
    
    scale_down_syst(S_TT, scale_down_factors)
    scale_down_syst(Q_TT, scale_down_factors)
    scale_down_syst(W_TT, scale_down_factors)
    scale_down_syst(L_TT, scale_down_factors)
    scale_down_syst(T_TT, scale_down_factors)
    scale_down_syst(s_TT, scale_down_factors)
    scale_down_syst(q_TT, scale_down_factors)
    scale_down_syst(S_MJ, scale_down_factors)
    scale_down_syst(Q_MJ, scale_down_factors)
    scale_down_syst(W_MJ, scale_down_factors)
    scale_down_syst(L_MJ, scale_down_factors)
    scale_down_syst(T_MJ, scale_down_factors)
    scale_down_syst(s_MJ, scale_down_factors)
    scale_down_syst(q_MJ, scale_down_factors)
    scale_down_syst(S_WJ, scale_down_factors)
    scale_down_syst(Q_WJ, scale_down_factors)
    scale_down_syst(W_WJ, scale_down_factors)
    scale_down_syst(L_WJ, scale_down_factors)
    scale_down_syst(T_WJ, scale_down_factors)
    scale_down_syst(L_WJ, scale_down_factors)
    scale_down_syst(s_WJ, scale_down_factors)
    scale_down_syst(q_WJ, scale_down_factors)
    scale_down_syst(S_ZI, scale_down_factors)
    scale_down_syst(Q_ZI, scale_down_factors)
    scale_down_syst(W_ZI, scale_down_factors)
    scale_down_syst(L_ZI, scale_down_factors)
    scale_down_syst(T_ZI, scale_down_factors)
    scale_down_syst(s_ZI, scale_down_factors)
    scale_down_syst(q_ZI, scale_down_factors)
    scale_down_syst(S_OT, scale_down_factors)
    scale_down_syst(Q_OT, scale_down_factors)
    scale_down_syst(W_OT, scale_down_factors)
    scale_down_syst(L_OT, scale_down_factors)
    scale_down_syst(T_OT, scale_down_factors)
    scale_down_syst(s_OT, scale_down_factors)
    scale_down_syst(q_OT, scale_down_factors)
    
    scale_down_syst(S_LooseWP_TT, scale_down_factors)
    scale_down_syst(S_LooseWP_MJ, scale_down_factors)
    scale_down_syst(S_LooseWP_WJ, scale_down_factors)
    scale_down_syst(S_LooseWP_ZI, scale_down_factors)
    scale_down_syst(S_LooseWP_OT, scale_down_factors)
    
    scale_down_syst(GDirectPrompt_GJ, scale_down_factors)

# --------------------------------
#             Purity

if purity_in_Gm1:
    f_zinv_est = ROOT.TFile.Open("zinv_est_"+opt.box+"_Gm1.root","RECREATE")
    plotdir = "Plots/z_inv_est/"+DATE+"/Gm1"
else:
    f_zinv_est = ROOT.TFile.Open("zinv_est_"+opt.box+".root","RECREATE")
    plotdir = "Plots/z_inv_est/"+DATE

print "- Fitting purity"

# In bins of MR
first = True
purity_MR_EB = ROOT.TH1D("purity_MR_EB"," ;M_{R} (GeV);Photon purity",5,0.5,5.5)
purity_MR_EE = ROOT.TH1D("purity_MR_EE"," ;M_{R} (GeV);Photon purity",5,0.5,5.5)
purity_MR_EB_MC = ROOT.TH1D("purity_MR_EB_MC","MC;MR bin;Photon purity",5,0.5,5.5)
purity_MR_EE_MC = ROOT.TH1D("purity_MR_EE_MC","MC;MR bin;Photon purity",5,0.5,5.5)
mcpurity_MR_EB = MCPurity_EB.Project3D("zx").ProfileX()
mcpurity_MR_EE = MCPurity_EE.Project3D("zx").ProfileX()
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
    mcpur_EB = mcpurity_MR_EB.GetBinContent(binx+2)
    mcpur_EE = mcpurity_MR_EE.GetBinContent(binx+2)

    pur_EB, pur_EB_err = fit_fraction(chiso_EB, temp1_EB, temp2_EB, plotdir, binname, mcpur_EB)
    pur_EE, pur_EE_err = fit_fraction(chiso_EE, temp1_EE, temp2_EE, plotdir, binname, mcpur_EE)
    purity_MR_EB.SetBinContent(binx, pur_EB)
    purity_MR_EB.SetBinError  (binx, pur_EB_err)
    purity_MR_EE.SetBinContent(binx, pur_EE)
    purity_MR_EE.SetBinError  (binx, pur_EE_err)
    purity_MR_EB_MC.SetBinContent(binx, mcpur_EB)
    purity_MR_EE_MC.SetBinContent(binx, mcpur_EE)
    purity_MR_EB_MC.SetBinError  (binx, 0)
    purity_MR_EE_MC.SetBinError  (binx, 0)

pur_MR = custom_can(purity_MR_EB, "Purity_vs_MR_"+opt.box, 0,0, 500,500, 90,20,20,75)
purity_MR_EB.GetXaxis().SetNdivisions(505)
purity_MR_EB.GetYaxis().SetRangeUser(0,2)
for i in range(len(mrbins_TeV)-1): purity_MR_EB.GetXaxis().SetBinLabel(i+1, "["+str(mrbins_TeV[i])+","+str(mrbins_TeV[i+1])+"]")
purity_MR_EB.GetXaxis().SetLabelSize(0.05)
purity_MR_EB.GetXaxis().SetTitleOffset(1)
purity_MR_EB.Draw("PE")
purity_MR_EE.SetLineColor(2)
purity_MR_EE.SetMarkerColor(2)
purity_MR_EB_MC.SetMarkerStyle(25)
purity_MR_EB_MC.SetMarkerSize (2)
purity_MR_EE_MC.SetLineColor(2)
purity_MR_EE_MC.SetMarkerColor(2)
purity_MR_EE_MC.SetMarkerStyle(25)
purity_MR_EE_MC.SetMarkerSize (2)
purity_MR_EE.Draw("SAME PE")
purity_MR_EB_MC.Draw("SAME PE")
purity_MR_EE_MC.Draw("SAME PE")
leg = ROOT.TLegend(0.28,0.75,0.9,0.9, BOX)
leg.SetNColumns(2)
leg.SetTextSize(0.04)
leg.AddEntry(purity_MR_EB,    "Barrel - Data", "LPE")
leg.AddEntry(purity_MR_EB_MC, "MC",            "PE")
leg.AddEntry(purity_MR_EE,    "Endcap - Data", "LPE")
leg.AddEntry(purity_MR_EE_MC, "MC",            "PE")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
add_cms_era(pur_MR, ERA, keep, opt.energy, lumi, prefix)
save_plot(pur_MR, "", plotdir+"/"+pur_MR.GetName())

# In bins of R2
first = True
purity_R2_EB = ROOT.TH1D("purity_R2_EB"," ;R^{2};Photon purity",5,0.5,5.5)
purity_R2_EE = ROOT.TH1D("purity_R2_EE"," ;R^{2};Photon purity",5,0.5,5.5)
purity_R2_EB_MC = ROOT.TH1D("purity_R2_EB_MC","MC;R2 bin;Photon purity",5,0.5,5.5)
purity_R2_EE_MC = ROOT.TH1D("purity_R2_EE_MC","MC;R2 bin;Photon purity",5,0.5,5.5)
mcpurity_R2_EB = MCPurity_EB.Project3D("zy").ProfileX()
mcpurity_R2_EE = MCPurity_EE.Project3D("zy").ProfileX()
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
    mcpur_EB = mcpurity_R2_EB.GetBinContent(biny+2)
    mcpur_EE = mcpurity_R2_EE.GetBinContent(biny+2)
    pur_EB, pur_EB_err = fit_fraction(chiso_EB, temp1_EB, temp2_EB, plotdir, binname, mcpur_EB)
    pur_EE, pur_EE_err = fit_fraction(chiso_EE, temp1_EE, temp2_EE, plotdir, binname, mcpur_EE)
    purity_R2_EB.SetBinContent(biny, pur_EB)
    purity_R2_EB.SetBinError  (biny, pur_EB_err)
    purity_R2_EE.SetBinContent(biny, pur_EE)
    purity_R2_EE.SetBinError  (biny, pur_EE_err)
    purity_R2_EB_MC.SetBinContent(biny, mcpur_EB)
    purity_R2_EE_MC.SetBinContent(biny, mcpur_EE)
    purity_R2_EB_MC.SetBinError  (biny, 0)
    purity_R2_EE_MC.SetBinError  (biny, 0)

pur_R2 = custom_can(purity_R2_EB, "Purity_vs_R2_"+opt.box, 0,0, 500,500, 90,20,20,75)
purity_R2_EB.GetXaxis().SetNdivisions(505)
purity_R2_EB.GetYaxis().SetRangeUser(0,2)
for i in range(len(r2bins)-1): purity_R2_EB.GetXaxis().SetBinLabel(i+1, "["+str(r2bins[i])+","+str(r2bins[i+1])+"]")
purity_R2_EB.GetXaxis().SetLabelSize(0.05)
purity_R2_EB.GetXaxis().SetTitleOffset(1)
purity_R2_EB.Draw("PE")
purity_R2_EE.SetLineColor(2)
purity_R2_EE.SetMarkerColor(2)
purity_R2_EB_MC.SetMarkerStyle(25)
purity_R2_EB_MC.SetMarkerSize (2)
purity_R2_EE_MC.SetLineColor(2)
purity_R2_EE_MC.SetMarkerColor(2)
purity_R2_EE_MC.SetMarkerStyle(25)
purity_R2_EE_MC.SetMarkerSize (2)
purity_R2_EE.Draw("SAME PE")
purity_R2_EB_MC.Draw("SAME PE")
purity_R2_EE_MC.Draw("SAME PE")
leg = ROOT.TLegend(0.28,0.75,0.9,0.9, BOX)
leg.SetNColumns(2)
leg.SetTextSize(0.04)
leg.AddEntry(purity_R2_EB,    "Barrel - Data", "LPE")
leg.AddEntry(purity_R2_EB_MC, "MC",            "PE")
leg.AddEntry(purity_R2_EE,    "Endcap - Data", "LPE")
leg.AddEntry(purity_R2_EE_MC, "MC",            "PE")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
add_cms_era(pur_R2, ERA, keep, opt.energy, lumi, prefix)
save_plot(pur_R2, "", plotdir+"/"+pur_R2.GetName())

# Measure average for EB/EE
first = True
purity_EB = ROOT.TH1D("purity_EB"," ;;Average photon purity",1,0,1)
purity_EE = ROOT.TH1D("purity_EE"," ;;Average photon purity",1,0,1)
purity_EB_MC = ROOT.TH1D("purity_EB_MC","MC;;Average photon purity",1,0,1)
purity_EE_MC = ROOT.TH1D("purity_EE_MC","MC;;Average photon purity",1,0,1)
mcpurity_EB = MCPurity_EB.Project3D("z")
mcpurity_EE = MCPurity_EE.Project3D("z")
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
mcpur_EB = mcpurity_EB.GetBinContent(2)/mcpurity_EB.Integral()
mcpur_EE = mcpurity_EE.GetBinContent(2)/mcpurity_EE.Integral()
pur_EB, pur_EB_err = fit_fraction(chiso_EB, temp1_EB, temp2_EB, plotdir, binname, mcpur_EB)
pur_EE, pur_EE_err = fit_fraction(chiso_EE, temp1_EE, temp2_EE, plotdir, binname, mcpur_EE)
purity_EB.SetBinContent(1, pur_EB)
purity_EB.SetBinError  (1, pur_EB_err)
purity_EE.SetBinContent(1, pur_EE)
purity_EE.SetBinError  (1, pur_EE_err)
purity_EB_MC.SetBinContent(1, mcpur_EB)
purity_EE_MC.SetBinContent(1, mcpur_EE)
purity_EB_MC.SetBinError  (1, 0)
purity_EE_MC.SetBinError  (1, 0)

pur = custom_can(purity_EB, "Purity_Average_"+opt.box)
purity_EB.GetXaxis().SetNdivisions(505)
purity_EB.GetYaxis().SetRangeUser(0,2)
purity_EB.Draw("PE")
purity_EE.SetLineColor(2)
purity_EE.SetMarkerColor(2)
purity_EB_MC.SetMarkerStyle(25)
purity_EB_MC.SetMarkerSize (2)
purity_EE_MC.SetLineColor(2)
purity_EE_MC.SetMarkerColor(2)
purity_EE_MC.SetMarkerStyle(25)
purity_EE_MC.SetMarkerSize (2)
purity_EE.Draw("SAME PE")
purity_EB_MC.Draw("SAME PE")
purity_EE_MC.Draw("SAME PE")
leg = ROOT.TLegend(0.28,0.75,0.9,0.9, "")
leg.SetNColumns(2)
leg.SetTextSize(0.04)
leg.AddEntry(purity_EB,    "Barrel - Data", "LPE")
leg.AddEntry(purity_EB_MC, "MC",            "PE")
leg.AddEntry(purity_EE,    "Endcap - Data", "LPE")
leg.AddEntry(purity_EE_MC, "MC",            "PE")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
save_plot(pur, "", plotdir+"/"+pur.GetName())

# Measure average purity (to subtract QCD etc. in data)
nevt_G_prompt = 0
nevt_G = 0
for binx in range(1, G_data_EB.GetNbinsX()+1):
    for biny in range(1, G_data_EB.GetNbinsY()+1):
        nevt_G        += G_data_EB.GetBinContent(binx,biny) + G_data_EE.GetBinContent(binx,biny)
        nevt_G_prompt += G_data_EB.GetBinContent(binx,biny) * pur_EB + G_data_EE.GetBinContent(binx,biny) * pur_EE
avg_purity_data = nevt_G_prompt / nevt_G
avg_purity_mc   = (mcpurity_EB.GetBinContent(2)+mcpurity_EE.GetBinContent(2)) / (mcpurity_EB.Integral()+mcpurity_EE.Integral())
print "Purity - EB  (MC):  "+str(pur_EB)+" ("+str(mcpur_EB)+")"
print "Purity - EE  (MC):  "+str(pur_EE)+" ("+str(mcpur_EE)+")"
print "Purity - Avg (MC):  "+str(avg_purity_data)+" ("+str(avg_purity_mc)+")"

# --------------------------------
#      Direct photon fraction
print "- Calculate direct photon fraction"

# Use an average direct photon fraction for EB/EE (It is ~0.9)
direct_frac_EB = IsDirect_G_EB.Project3D("z").GetMean()
direct_frac_EE = IsDirect_G_EE.Project3D("z").GetMean()

print "Direct fraction - EB: "+str(direct_frac_EB)
print "Direct fraction - EE: "+str(direct_frac_EE)

# --------------------------------
#  Direct & Prompt photon counts

# Estimate Prompt direct photons in the G control region
print "- Estimate prompt direct photon counts in G region"
h_npromptdirect               = ROOT.TH1D("npromptdirect",               ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
h_npromptdirect_purityUp      = ROOT.TH1D("npromptdirect_purityUp",      ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
h_npromptdirect_purityDown    = ROOT.TH1D("npromptdirect_purityDown",    ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
h_npromptdirect_dirfracUp     = ROOT.TH1D("npromptdirect_dirfracUp",     ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
h_npromptdirect_dirfracDown   = ROOT.TH1D("npromptdirect_dirfracDown",   ";Bin;Z(#nu#nu) photon estimate",nrazorbin,0,nrazorbin)
sum_EB = 0
sum_EE = 0
for binx in range(1, G_data_EB.GetNbinsX()+1):
    for biny in range(1, G_data_EB.GetNbinsY()+1):
        # Use average purity (earlier version)
        #purity_EB = pur_EB
        #purity_EE = pur_EE
        # switch to use purity in bins of R2 (new version)
        purity_EB = purity_R2_EB.GetBinContent(biny)
        purity_EE = purity_R2_EE.GetBinContent(biny)
        # Add prompt direct photons in EB and EE
        npromptdirect  = G_data_EB.GetBinContent(binx,biny) * purity_EB * direct_frac_EB
        npromptdirect += G_data_EE.GetBinContent(binx,biny) * purity_EE * direct_frac_EE
        sum_EB += G_data_EB.GetBinContent(binx,biny) * purity_EB * direct_frac_EB
        sum_EE += G_data_EE.GetBinContent(binx,biny) * purity_EE * direct_frac_EE
        npromptdirect_purityUp     = G_data_EB.GetBinContent(binx,biny) * (purity_EB + purity_err) * direct_frac_EB
        npromptdirect_purityUp    += G_data_EE.GetBinContent(binx,biny) * (purity_EE + purity_err) * direct_frac_EE
        npromptdirect_purityDown   = G_data_EB.GetBinContent(binx,biny) * (purity_EB - purity_err) * direct_frac_EB
        npromptdirect_purityDown  += G_data_EE.GetBinContent(binx,biny) * (purity_EE - purity_err) * direct_frac_EE
        npromptdirect_dirfracUp    = G_data_EB.GetBinContent(binx,biny) * purity_EB * (direct_frac_EB + dirfrac_err)
        npromptdirect_dirfracUp   += G_data_EE.GetBinContent(binx,biny) * purity_EE * (direct_frac_EE + dirfrac_err)
        npromptdirect_dirfracDown  = G_data_EB.GetBinContent(binx,biny) * purity_EB * (direct_frac_EB - dirfrac_err)
        npromptdirect_dirfracDown += G_data_EE.GetBinContent(binx,biny) * purity_EE * (direct_frac_EE - dirfrac_err)
        # Calculate errors
        # add statistical error
        npromptdirect_err              = (G_data_EB.GetBinError(binx,biny) * purity_EB * direct_frac_EB) ** 2
        npromptdirect_err             += (G_data_EE.GetBinError(binx,biny) * purity_EE * direct_frac_EE) ** 2
        npromptdirect_err_purityUp     = (G_data_EB.GetBinError(binx,biny) * (purity_EB + purity_err) * direct_frac_EB) ** 2
        npromptdirect_err_purityUp    += (G_data_EE.GetBinError(binx,biny) * (purity_EE + purity_err) * direct_frac_EE) ** 2
        npromptdirect_err_purityDown   = (G_data_EB.GetBinError(binx,biny) * (purity_EB - purity_err) * direct_frac_EB) ** 2
        npromptdirect_err_purityDown  += (G_data_EE.GetBinError(binx,biny) * (purity_EE - purity_err) * direct_frac_EE) ** 2
        npromptdirect_err_dirfracUp    = (G_data_EB.GetBinError(binx,biny) * purity_EB * (direct_frac_EB + dirfrac_err)) ** 2
        npromptdirect_err_dirfracUp   += (G_data_EE.GetBinError(binx,biny) * purity_EE * (direct_frac_EE + dirfrac_err)) ** 2
        npromptdirect_err_dirfracDown  = (G_data_EB.GetBinError(binx,biny) * purity_EB * (direct_frac_EB - dirfrac_err)) ** 2
        npromptdirect_err_dirfracDown += (G_data_EE.GetBinError(binx,biny) * purity_EE * (direct_frac_EE - dirfrac_err)) ** 2
        npromptdirect_err              = npromptdirect_err             ** 0.5
        npromptdirect_err_purityUp     = npromptdirect_err_purityUp    ** 0.5
        npromptdirect_err_purityDown   = npromptdirect_err_purityDown  ** 0.5
        npromptdirect_err_dirfracUp    = npromptdirect_err_dirfracUp   ** 0.5
        npromptdirect_err_dirfracDown  = npromptdirect_err_dirfracDown ** 0.5
        # fill measured direct prompt photon counts
        unrolled_bin = (binx-1)*G_data_EB.GetNbinsY()+biny
        h_npromptdirect              .SetBinContent(unrolled_bin, npromptdirect)
        h_npromptdirect_purityUp     .SetBinContent(unrolled_bin, npromptdirect_purityUp)
        h_npromptdirect_purityDown   .SetBinContent(unrolled_bin, npromptdirect_purityDown)
        h_npromptdirect_dirfracUp    .SetBinContent(unrolled_bin, npromptdirect_dirfracUp)
        h_npromptdirect_dirfracDown  .SetBinContent(unrolled_bin, npromptdirect_dirfracDown)
        h_npromptdirect              .SetBinError  (unrolled_bin, npromptdirect_err)
        h_npromptdirect_purityUp     .SetBinError  (unrolled_bin, npromptdirect_err_purityUp)
        h_npromptdirect_purityDown   .SetBinError  (unrolled_bin, npromptdirect_err_purityDown)
        h_npromptdirect_dirfracUp    .SetBinError  (unrolled_bin, npromptdirect_err_dirfracUp)
        h_npromptdirect_dirfracDown  .SetBinError  (unrolled_bin, npromptdirect_err_dirfracDown)
print "EB:  "+str(sum_EB)
print "EE:  "+str(sum_EE)
print "SUM: "+str(h_npromptdirect.Integral())

vh_npromptdirect = []
vh_npromptdirect.append(h_npromptdirect)
vh_npromptdirect.append(h_npromptdirect_purityUp)
vh_npromptdirect.append(h_npromptdirect_purityDown)
vh_npromptdirect.append(h_npromptdirect_dirfracUp)
vh_npromptdirect.append(h_npromptdirect_dirfracDown)

# --------------------------------
#     Double ratio = k_Z / k_G

print "- Calculate double ratio"
# Get integrals and errors
eZ     = ROOT.Double(0)
eZ_NDY = ROOT.Double(0)
eZ_DY  = ROOT.Double(0)
eG     = ROOT.Double(0)
eG_GJ  = ROOT.Double(0)
nZ     = Z_data             .IntegralAndError(0,-1,eZ)
nZ_NDY = Z_NONDY            .IntegralAndError(0,-1,eZ_NDY)
nZ_DY  = Z_DY               .IntegralAndError(0,-1,eZ_DY)
nG     = h_npromptdirect    .IntegralAndError(0,-1,eG)
nG_GJ  = GDirectPrompt_GJ[0].IntegralAndError(0,-1,eG_GJ)
# For DY in Z region, use MC
k_Z = (nZ - nZ_NDY)/nZ_DY
ek_Z = div_err(nZ-nZ_NDY, (eZ**2 + eZ_NDY**2)**0.5, nZ_DY, eZ_DY) 
# For GJets in G, instead use data measurement
k_G = nG/nG_GJ
ek_G = div_err(nG, eG, nG_GJ, eG_GJ) 
DR  = k_Z / k_G
eDR = div_err(k_Z, ek_Z, k_G, ek_G)
print "Z_data:  "+str(nZ)    +" +- "+str(eZ)
print "Z_nonDY: "+str(nZ_NDY)+" +- "+str(eZ_NDY)
print "Z_DY:    "+str(nZ_DY) +" +- "+str(eZ_DY)
print "G_data (prompt, direct): "+str(nG)   +" +- "+str(eG)
print "G_GJ   (prompt, direct): "+str(nG_GJ)+" +- "+str(eG_GJ)
print ("k_Z = %4.3f +- %4.3f (%4.3f / %4.3f), k_G = %4.3f +- %4.3f (%4.3f / %4.3f), double ratio = %3.2f +- %3.2f" %
       (k_Z, ek_Z, nZ - nZ_NDY, nZ_DY, k_G, ek_G, nG, nG_GJ, DR, eDR))
#sys.exit()

# ------------------ Estimate ---------------------------

print "Perform Z(nunu) estimate"

L_subtract = [L_TT[0],  L_MJ[0],           L_ZI[0], L_OT[0]]
ZInv_est   = zinv_est("ZInv",   [S_ZI, S_LooseWP_ZI], f_zinv_est, vh_npromptdirect, GDirectPrompt_GJ, L_data, L_subtract, L_WJ[0], combine_bins, binned_k, DR, eDR)

# Do the same for closure tests
s_ZInv_est = zinv_est("s_ZInv", [s_ZI],               f_zinv_est, vh_npromptdirect, GDirectPrompt_GJ, L_data, L_subtract, L_WJ[0], combine_bins, binned_k, DR, eDR)
q_ZInv_est = zinv_est("q_ZInv", [q_ZI],               f_zinv_est, vh_npromptdirect, GDirectPrompt_GJ, L_data, L_subtract, L_WJ[0], combine_bins, binned_k, DR, eDR)

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
    evt_weight = (lumi * 1000 * xsecs[float(scan_point.split("_")[0])] / N)
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
        #print ("npv=%.1f fitPred=%.3f npvWeight=%f avgAcc+=%f accErr+=%f" %(npv, fitPred, npvWeight, npvWeight * fitPred, npvWeight * fitError)) 
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

# --------------- Signal extrapolation ------------------
#                  to higher masses


xsecs_upgrade = {}
for line in open('./data/stop'+str(opt.energy)+'TeV.txt' if ('T2' in opt.model) else 'data/gluino'+str(opt.energy)+'TeV.txt','r'):
    line = line.replace('\n','')
    xsecs_upgrade[float(line.split(',')[0])]=float(line.split(',')[1]) #pb

S_signal_extension = []
for i in range(len(S_signal)):
    if "MRR2_S_signal_"+str(highest_mass) in S_signal[i][0].GetName():
        mlsp = int(S_signal[i][0].GetName().split("_")[4])
        extend = []
        for j in range(len(S_signal[i])):
            for k in range(len(extension_points)):
                h_ext = S_signal[i][j].Clone(S_signal[i][j].GetName().replace("MRR2_S_signal_"+str(highest_mass), "MRR2_S_signal_"+str(extension_points[k])))
                xsec_ratio = xsecs_upgrade[extension_points[k]]/xsecs_upgrade[highest_mass]
                acceptance_ratio = vf_acc[mlsp].Eval(extension_points[k])/vf_acc[mlsp].Eval(highest_mass)
                h_ext.Scale(xsec_ratio*acceptance_ratio)
                if len(extend)<len(extension_points): extend.append([])
                extend[k].append(h_ext)
        S_signal_extension = S_signal_extension + extend

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
##  hist = ROOT.TH1D("h",";Data/Prediction;A. U.",20,-10,10)
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
    fname_closure = "closure_test_"+opt.box+"_binned_k.root"
else:
    fname_closure = "closure_test_"+opt.box+"_factorized_k.root"

fout = ROOT.TFile(fname_closure,"RECREATE")

# S'
# s_MJ_est = (Q_data - Q_notMJ) * s_MJ / Q_MJ
# s_WJ_est = (W_data - W_notWJ) * s_WJ / W_WJ
# s_TT_est = (T_data - T_notTT) * s_TT / T_TT
# s_ZI_est = (L_data - L_notWJ) * s_ZI / L_ZI
print "- in S'"

s_Top_est      = []
s_MultiJet_est = []
s_WJets_est    = []
if not use_G: s_ZInv_est     = []
s_Other_est    = []
for i in range(0, len(systematics)):
    s_TT_est = bg_est("Top_s_T"       +systematics[i], T_data, [          T_MJ[0], T_WJ[0], T_ZI[0], T_OT[0]], [s_TT[i]], fout, T_TT[i], combine_bins, binned_k)
    s_MJ_est = bg_est("MultiJet_s_Q"  +systematics[i], Q_data, [Q_TT[0],           Q_WJ[0], Q_ZI[0], Q_OT[0]], [s_MJ[i]], fout, Q_MJ[i], combine_bins, binned_k)
    s_WJ_est = bg_est("WJets_s_W"     +systematics[i], W_data, [W_TT[0],  W_MJ[0],          W_ZI[0], W_OT[0]], [s_WJ[i]], fout, W_WJ[i], combine_bins, binned_k)
    if not use_G:
        s_ZI_est = bg_est("ZInv_s_L"      +systematics[i], L_data, [L_TT[0],  L_MJ[0],          L_ZI[0], L_OT[0]], [s_WJ[i]], fout, L_WJ[i], combine_bins, binned_k)
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
    if not use_G:
        s_ZInv_est    .append(s_ZI_est)
    s_Other_est   .append(s_OT_est)

s_TT_nom, s_TT_up, s_TT_dn = calc_syst_err(s_Top_est,     "s_Top_syst")
s_MJ_nom, s_MJ_up, s_MJ_dn = calc_syst_err(s_MultiJet_est,"s_MultiJet_syst", qcd_syst)
s_WJ_nom, s_WJ_up, s_WJ_dn = calc_syst_err(s_WJets_est,   "s_WJets_syst")
s_ZI_nom, s_ZI_up, s_ZI_dn = calc_syst_err(s_ZInv_est,    "s_ZInv_syst", dy_syst)
s_OT_nom, s_OT_up, s_OT_dn = calc_syst_err(s_Other_est,   "s_Other_syst")
s_syst_err = s_TT_nom.Clone("s_TotalErr")
##s_syst_err.SetFillColor(1)
##s_syst_err.SetFillStyle(3002)
#s_syst_err.SetFillColor(ROOT.kGray) # 920
#s_syst_err.SetFillStyle(1001)
s_syst_err.SetFillColor(13)
s_syst_err.SetFillStyle(3001)
s_syst_err.SetMarkerStyle(0)
s_stat_err = s_TT_nom.Clone("s_StatErr")
s_stat_err.SetFillColor(1)
s_stat_err.SetFillStyle(3004)
s_stat_err.SetMarkerStyle(0)
s_stat_err.SetMarkerColor(0)
for binx in range(1, s_syst_err.GetNbinsX()+1):
    # Sum total stat+syst error for all bkg components
    err_up = (s_TT_up.GetBinError(binx)**2 + s_MJ_up.GetBinError(binx)**2 + s_WJ_up.GetBinError(binx)**2 + s_ZI_up.GetBinError(binx)**2 + s_OT_up.GetBinError(binx)**2) ** 0.5
    err_dn = (s_TT_dn.GetBinError(binx)**2 + s_MJ_dn.GetBinError(binx)**2 + s_WJ_dn.GetBinError(binx)**2 + s_ZI_dn.GetBinError(binx)**2 + s_OT_dn.GetBinError(binx)**2) ** 0.5
    # Sum nominal counts
    nom = s_TT_nom.GetBinContent(binx) + s_MJ_nom.GetBinContent(binx) + s_WJ_nom.GetBinContent(binx) + s_ZI_nom.GetBinContent(binx) + s_OT_nom.GetBinContent(binx)
    # Make asymmetric interval by simply shifting the mean
    s_syst_err.SetBinContent(binx, nom + (err_up - err_dn)/2)
    s_syst_err.SetBinError(binx, (err_up+err_dn)/2)
    # Save stat error separately
    totstat  = s_TT_nom.GetBinError(binx) ** 2
    totstat += s_MJ_nom.GetBinError(binx) ** 2
    totstat += s_WJ_nom.GetBinError(binx) ** 2
    totstat += s_ZI_nom.GetBinError(binx) ** 2
    totstat += s_OT_nom.GetBinError(binx) ** 2
    s_stat_err.SetBinContent(binx, nom)
    s_stat_err.SetBinError  (binx, totstat ** 0.5)

can = custom_can(s_data, "Closure_test_s_"+opt.box)
can.SetLogy(1)
s_data.GetYaxis().SetRangeUser(1.01e-1,1e6)
s_data.GetYaxis().SetTitle("Events / bin")
s_data.Draw("PE0")
s_TT_nom.SetLineColor(600)
s_MJ_nom.SetLineColor(600)
s_WJ_nom.SetLineColor(600)
s_ZI_nom.SetLineColor(600)
s_OT_nom.SetLineColor(600)
s_TT_nom.SetFillColor(418)
s_MJ_nom.SetFillColor(618)
s_WJ_nom.SetFillColor(633)
s_ZI_nom.SetFillColor(433)
s_OT_nom.SetFillColor(865)
s_stack = ROOT.THStack("s_TotBkgEst","")
s_stack.Add(s_OT_nom)
s_stack.Add(s_ZI_nom)
s_stack.Add(s_WJ_nom)
if "TopAna" in opt.box:
    s_stack.Add(s_MJ_nom)
    s_stack.Add(s_TT_nom)
else:
    s_stack.Add(s_TT_nom)
    s_stack.Add(s_MJ_nom)
s_stack.Draw("SAME HIST")
s_syst_err.Draw("SAME E2")
s_stat_err.Draw("SAME E2")
s_data.Draw("SAMEPE0")
draw_mr_bins([s_data, s_stack], 1.01e-1,1e4, combine_bins, keep, mrbins_TeV, r2bins)
leg = ROOT.TLegend(0.68,0.51,0.98,0.86, BOX)
leg.AddEntry(s_data,   "#color[1]{Data}",                       "pe")
if "TopAna" in opt.box:
    leg.AddEntry(s_TT_nom, "#color[418]{t#bar{t} or single t}", "f")
    leg.AddEntry(s_MJ_nom, "#color[618]{Multijet}",             "f")
else:
    leg.AddEntry(s_MJ_nom, "#color[618]{Multijet}",             "f")
    leg.AddEntry(s_TT_nom, "#color[418]{t#bar{t} or single t}", "f")
leg.AddEntry(s_WJ_nom, "#color[633]{W(l#nu)+jets}}", "f")
leg.AddEntry(s_ZI_nom, "#color[433]{Z(#nu#nu)+jets}","f")
leg.AddEntry(s_OT_nom, "#color[865]{Other}",                    "f")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
ROOT.gPad.Update()
#can.Write()
can = add_stack_ratio_plot(can, 0,s_data.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins)
add_cms_era(can, ERA, keep, opt.energy, lumi, prefix)
lat_s = ROOT.TLatex(0.75,2.8e5, "Signal-like validation region")
lat_s.Draw("SAME")
zinv_region = "G_region" if use_G else "L_region"
if binned_k:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/binned_k/"+zinv_region+"/"+can.GetName(), 1)
else:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/factorized_k/"+zinv_region+"/"+can.GetName(), 1)

# Q'
# q_MJ_est = (Q_data - Q_notMJ) * q_MJ / Q_MJ
# q_WJ_est = (W_data - W_notWJ) * q_WJ / W_WJ
# q_TT_est = (T_data - T_notTT) * q_TT / T_TT
# q_ZI_est = (L_data - L_notWJ) * q_ZI / L_ZI
print "- in Q'"

q_Top_est      = []
q_MultiJet_est = []
q_WJets_est    = []
if not use_G: q_ZInv_est     = []
q_Other_est    = []
for i in range(0, len(systematics)):
    q_TT_est = bg_est("Top_q_T"         +systematics[i], T_data, [          T_MJ[0], T_WJ[0], T_ZI[0], T_OT[0]], [q_TT[i]], fout, T_TT[i], combine_bins, binned_k)
    q_MJ_est = bg_est("MultiJet_q_Q"    +systematics[i], Q_data, [Q_TT[0],           Q_WJ[0], Q_ZI[0], Q_OT[0]], [q_MJ[i]], fout, Q_MJ[i], combine_bins, binned_k)
    q_WJ_est = bg_est("WJets_q_W"       +systematics[i], W_data, [W_TT[0],  W_MJ[0],          W_ZI[0], W_OT[0]], [q_WJ[i]], fout, W_WJ[i], combine_bins, binned_k)
    if not use_G: q_ZI_est = bg_est("ZInv_q_L"        +systematics[i], L_data, [L_TT[0],  L_MJ[0],          L_ZI[0], L_OT[0]], [q_WJ[i]], fout, L_WJ[i], combine_bins, binned_k)
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
    if not use_G: q_ZInv_est    .append(q_ZI_est)
    q_Other_est   .append(q_OT_est)

q_TT_nom, q_TT_up, q_TT_dn = calc_syst_err(q_Top_est,     "q_Top_syst")
q_MJ_nom, q_MJ_up, q_MJ_dn = calc_syst_err(q_MultiJet_est,"q_MultiJet_syst", qcd_syst)
q_WJ_nom, q_WJ_up, q_WJ_dn = calc_syst_err(q_WJets_est,   "q_WJets_syst")
q_ZI_nom, q_ZI_up, q_ZI_dn = calc_syst_err(q_ZInv_est,    "q_ZInv_syst", dy_syst)
q_OT_nom, q_OT_up, q_OT_dn = calc_syst_err(q_Other_est,   "q_Other_syst")
q_syst_err = q_TT_nom.Clone("q_TotalErr")
##q_syst_err.SetFillColor(1)
##q_syst_err.SetFillStyle(3002)
#q_syst_err.SetFillColor(ROOT.kGray) # 920
#q_syst_err.SetFillStyle(1001)
q_syst_err.SetFillColor(13)
q_syst_err.SetFillStyle(3001)
q_syst_err.SetMarkerStyle(0)
q_stat_err = q_TT_nom.Clone("q_StatErr")
q_stat_err.SetFillColor(1)
q_stat_err.SetFillStyle(3004)
q_stat_err.SetMarkerStyle(0)
q_stat_err.SetMarkerColor(0)
for binx in range(1, q_syst_err.GetNbinsX()+1):
    # Sum total stat+syst error for all bkg components
    err_up = (q_TT_up.GetBinError(binx)**2 + q_MJ_up.GetBinError(binx)**2 + q_WJ_up.GetBinError(binx)**2 + q_ZI_up.GetBinError(binx)**2 + q_OT_up.GetBinError(binx)**2) ** 0.5
    err_dn = (q_TT_dn.GetBinError(binx)**2 + q_MJ_dn.GetBinError(binx)**2 + q_WJ_dn.GetBinError(binx)**2 + q_ZI_dn.GetBinError(binx)**2 + q_OT_dn.GetBinError(binx)**2) ** 0.5
    # Sum nominal counts
    nom = q_TT_nom.GetBinContent(binx) + q_MJ_nom.GetBinContent(binx) + q_WJ_nom.GetBinContent(binx) + q_ZI_nom.GetBinContent(binx) + q_OT_nom.GetBinContent(binx)
    # Make asymmetric interval by simply shifting the mean
    q_syst_err.SetBinContent(binx, nom + (err_up - err_dn)/2)
    q_syst_err.SetBinError(binx, (err_up+err_dn)/2)
    # Save stat error separately
    totstat  = q_TT_nom.GetBinError(binx) ** 2
    totstat += q_MJ_nom.GetBinError(binx) ** 2
    totstat += q_WJ_nom.GetBinError(binx) ** 2
    totstat += q_ZI_nom.GetBinError(binx) ** 2
    totstat += q_OT_nom.GetBinError(binx) ** 2
    q_stat_err.SetBinContent(binx, nom)
    q_stat_err.SetBinError  (binx, totstat ** 0.5)  

can = custom_can(q_data, "Closure_test_q_"+opt.box)
can.SetLogy(1)
q_data.GetYaxis().SetRangeUser(1.01e-1,1e5)
q_data.GetYaxis().SetTitle("Events / bin")
q_data.Draw("PE0")
q_TT_nom.SetLineColor(600)
q_MJ_nom.SetLineColor(600)
q_WJ_nom.SetLineColor(600)
q_ZI_nom.SetLineColor(600)
q_OT_nom.SetLineColor(600)
q_TT_nom.SetFillColor(418)
q_MJ_nom.SetFillColor(618)
q_WJ_nom.SetFillColor(633)
q_ZI_nom.SetFillColor(433)
q_OT_nom.SetFillColor(865)
q_stack = ROOT.THStack("q_TotBkgEst","")
q_stack.Add(q_OT_nom)
q_stack.Add(q_TT_nom)
q_stack.Add(q_WJ_nom)
if opt.box=="WAna_nj45":
    q_stack.Add(q_MJ_nom)
    q_stack.Add(q_ZI_nom)
else:
    q_stack.Add(q_ZI_nom)
    q_stack.Add(q_MJ_nom)
q_stack.Draw("SAME HIST")
q_syst_err.Draw("SAME E2")
q_stat_err.Draw("SAME E2")
q_data.Draw("SAMEPE0")
draw_mr_bins([q_data, q_stack], 1.01e-1,1e4, combine_bins, keep, mrbins_TeV, r2bins)
leg = ROOT.TLegend(0.68,0.51,0.98,0.86, BOX)
leg.AddEntry(q_data,   "#color[1]{Data}",                            "pe")
if opt.box=="WAna_nj45":
    leg.AddEntry(q_MJ_nom, "#color[618]{Multijet}",                  "f")
    leg.AddEntry(q_ZI_nom, "#color[433]{Z(#nu#nu)+jets}", "f")
else:
    leg.AddEntry(q_ZI_nom, "#color[433]{Z(#nu#nu)+jets}", "f")
    leg.AddEntry(q_MJ_nom, "#color[618]{Multijet}",                  "f")
leg.AddEntry(q_WJ_nom, "#color[633]{W(l#nu)+jets}",       "f")
leg.AddEntry(q_TT_nom, "#color[418]{t#bar{t} or single t}",          "f")
leg.AddEntry(q_OT_nom, "#color[865]{Other}",                         "f")
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw("SAME")
ROOT.gPad.Update()
#can.Write()
can = add_stack_ratio_plot(can, 0,q_data.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins, {}, False, 0.26 if opt.box=="WAna_nj45" else 0.16)
add_cms_era(can, ERA, keep, opt.energy, lumi, prefix)
lat_q = ROOT.TLatex(0.75,3.5e4, "Multijet validation region")
lat_q.Draw("SAME")
if binned_k:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/binned_k/"+zinv_region+"/"+can.GetName(), 1)
else:
    save_plot(can, "", "Plots/closure_test/"+DATE+"/factorized_k/"+zinv_region+"/"+can.GetName(), 1)
f.Close()
#sys.exit()

# --------------- Background Estimation ------------------

# Formulas for bkg estimate:
# S_MJ_est = (Q_data - Q_notMJ) * S_MJ / Q_MJ
# T_MJ_est = (Q_data - Q_notMJ) * T_MJ / Q_MJ --> Not used, because not much MJ in T
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

systematics.append("_extrapUp")
systematics.append("_extrapDown")

fout = ROOT.TFile.Open("bkg_estimate_"+opt.box+".root","RECREATE")
Top_est      = []
MultiJet_est = []
WJets_est    = []
if not use_G: ZInv_est     = []
Other_est    = []
for i in range(0, len(systematics)):
    j = i if i<len(systematics)-2 else 0
    S_TT_est = bg_est("Top"     +systematics[i], T_data, [          T_MJ[0], T_WJ[0], T_ZI[0], T_OT[0]], [S_TT[j], S_LooseWP_TT[j]], fout, T_TT[j], combine_bins, binned_k)
    S_MJ_est = bg_est("MultiJet"+systematics[i], Q_data, [Q_TT[0],           Q_WJ[0], Q_ZI[0], Q_OT[0]], [S_MJ[j], S_LooseWP_MJ[j]], fout, Q_MJ[j], combine_bins, binned_k)
    S_WJ_est = bg_est("WJets"   +systematics[i], W_data, [W_TT[0],  W_MJ[0],          W_ZI[0], W_OT[0]], [S_WJ[j], S_LooseWP_WJ[j]], fout, W_WJ[j], combine_bins, binned_k)
    if not use_G: S_ZI_est = bg_est("ZInv"+systematics[i], L_data, [L_TT[0],  L_MJ[0],          L_ZI[0], L_OT[0]], [S_ZI[j], S_LooseWP_ZI[j]], fout, L_WJ[j], combine_bins, binned_k)
    doExtrap = 1 if systematics[i]=="_extrapUp" else (-1 if systematics[i]=="_extrapDown" else 0)
    if combine_bins:
        S_OT_est            = combinebins(S_OT[j], "Other"+systematics[i])
        LooseRegionName = S_LooseWP_OT[j].GetName()+"_combined"
        if i>=len(systematics)-2:
            LooseRegionName = LooseRegionName.replace("_OT",systematics[i]+"_OT")
        fix_low_stat_bins(fout, S_OT_est, combinebins(S_LooseWP_OT[j], LooseRegionName), doExtrap)
    else:
        S_OT_est            = S_OT[j].Clone("Other"+systematics[i])
        fix_low_stat_bins(fout, S_OT_est, S_LooseWP_OT[j], doExtrap)
    # Sometimes MC sum has negative counts (due to NLO MCs)
    for binx in range(1,S_OT_est.GetNbinsX()+1):
        if S_OT_est.GetBinContent(binx)<0:
            S_OT_est.SetBinContent(binx,0)
            S_OT_est.SetBinError  (binx,0)
        # remove MC stat error for the YR18 and stat only scenarios
        if opt.scenario>0: S_OT_est.SetBinError(binx,0)            
    Top_est     .append(S_TT_est)
    MultiJet_est.append(S_MJ_est)
    WJets_est   .append(S_WJ_est)
    if not use_G: ZInv_est.append(S_ZI_est)
    Other_est   .append(S_OT_est)

# Transfer factor plots
make_kappa_plot("Top, "     +BOX, S_TT, T_TT,             combine_bins, 0)
make_kappa_plot("Multijet, "+BOX, S_MJ, Q_MJ,             combine_bins, qcd_syst)
make_kappa_plot("WJets, "   +BOX, S_WJ, W_WJ,             combine_bins, 0)
make_kappa_plot("ZInv (G), "+BOX, S_ZI, GDirectPrompt_GJ, combine_bins, dy_syst)
make_kappa_plot("ZInv (L), "+BOX, S_ZI, L_WJ,             combine_bins, dy_syst)

#sys.exit()

# Unblinded plot
S_TT_nom, S_TT_up, S_TT_dn = calc_syst_err(Top_est,     "S_Top_syst", 0)
S_MJ_nom, S_MJ_up, S_MJ_dn = calc_syst_err(MultiJet_est,"S_MultiJet_syst", qcd_syst)
S_WJ_nom, S_WJ_up, S_WJ_dn = calc_syst_err(WJets_est,   "S_WJets_syst", 0)
S_ZI_nom, S_ZI_up, S_ZI_dn = calc_syst_err(ZInv_est,    "S_ZInv_syst", dy_syst)
S_OT_nom, S_OT_up, S_OT_dn = calc_syst_err(Other_est,   "S_Other_syst", 0)
S_syst_err = S_TT_nom.Clone("S_TotalErr")
S_stat_err = S_TT_nom.Clone("S_StatErr")
S_stat_err.SetLineColor(1)
S_stat_err.SetMarkerStyle(0)
S_stat_err.SetMarkerColor(0)
S_stat_err.SetFillColor(1)
S_stat_err.SetFillStyle(3004)
S_syst_err.SetLineColor(1)
S_syst_err.SetMarkerStyle(0)
S_syst_err.SetMarkerColor(0)
S_syst_err.SetFillColor(ROOT.kGray)
S_syst_err.SetFillStyle(1001)
##  ##S_syst_err.SetFillColor(1)
##  ##S_syst_err.SetFillStyle(3002)
##  #S_syst_err.SetFillColor(ROOT.kGray) # 920
##  #S_syst_err.SetFillStyle(1001)
##  S_syst_err.SetFillColor(13)
##  S_syst_err.SetFillStyle(3001)
##  S_syst_err.SetMarkerStyle(0)
##  S_stat_err.SetFillColor(1)
##  S_stat_err.SetFillStyle(3004)
##  S_stat_err.SetMarkerStyle(0)
##  S_stat_err.SetMarkerColor(0)
for binx in range(1, S_syst_err.GetNbinsX()+1):
    # Sum total stat+syst error for all bkg components
    err_up = (S_TT_up.GetBinError(binx)**2 + S_MJ_up.GetBinError(binx)**2 + S_WJ_up.GetBinError(binx)**2 + S_ZI_up.GetBinError(binx)**2 + S_OT_up.GetBinError(binx)**2) ** 0.5
    err_dn = (S_TT_dn.GetBinError(binx)**2 + S_MJ_dn.GetBinError(binx)**2 + S_WJ_dn.GetBinError(binx)**2 + S_ZI_dn.GetBinError(binx)**2 + S_OT_dn.GetBinError(binx)**2) ** 0.5
    # Sum nominal counts
    nom = S_TT_nom.GetBinContent(binx) + S_MJ_nom.GetBinContent(binx) + S_WJ_nom.GetBinContent(binx) + S_ZI_nom.GetBinContent(binx) + S_OT_nom.GetBinContent(binx)
    # Make asymmetric interval by simply shifting the mean
    S_syst_err.SetBinContent(binx, nom + (err_up - err_dn)/2)
    S_syst_err.SetBinError(binx, (err_up+err_dn)/2)
    # Save stat error separately
    totstat  = S_TT_nom.GetBinError(binx) ** 2
    totstat += S_MJ_nom.GetBinError(binx) ** 2
    totstat += S_WJ_nom.GetBinError(binx) ** 2
    totstat += S_ZI_nom.GetBinError(binx) ** 2
    totstat += S_OT_nom.GetBinError(binx) ** 2
    S_stat_err.SetBinContent(binx, nom)
    S_stat_err.SetBinError  (binx, totstat ** 0.5)

ymin = 1.01e0
if lumi<100:
    ymax = 1e6
elif lumi<=3000:
    ymax = 1e6
else:
    ymin = 1.01e1
    ymax = 1e8

S_data.SetTitle(" ")
can = custom_can(S_data, "BkgEstimate_"+opt.box, 0, 0)
can.SetLogy(1)
S_data.GetXaxis().SetLabelSize(0.04)
S_data.GetYaxis().SetRangeUser(ymin,ymax)
S_data.GetYaxis().SetTitle("Events / bin")
S_data.SetMarkerColor(0)
S_data.Draw("AXIS")
S_TT_nom.SetLineColor(600)
S_MJ_nom.SetLineColor(600)
S_WJ_nom.SetLineColor(600)
S_ZI_nom.SetLineColor(600)
S_OT_nom.SetLineColor(600)
S_TT_nom.SetFillColor(418)
S_MJ_nom.SetFillColor(618)
S_WJ_nom.SetFillColor(633)
S_ZI_nom.SetFillColor(433)
S_OT_nom.SetFillColor(865)
S_stack = ROOT.THStack("S_TotBkgEst","")
S_stack.Add(S_OT_nom)
S_stack.Add(S_ZI_nom)
S_stack.Add(S_WJ_nom)
S_stack.Add(S_MJ_nom)
S_stack.Add(S_TT_nom)
S_stack.Draw("SAME HIST")
S_syst_err.Draw("SAME E2")
S_stat_err.Draw("SAME E2")
#S_data.Draw("SAMEPE0")
S_T1tttt.SetLineColor(619)
S_T2tt  .SetLineColor(401)
S_T5ttcc.SetLineColor(601)
S_T1tttt.SetLineStyle(7)
S_T2tt  .SetLineStyle(7)
S_T5ttcc.SetLineStyle(7)
S_T1tttt.SetLineWidth(3)
S_T2tt  .SetLineWidth(3)
S_T5ttcc.SetLineWidth(3)
S_T1tttt.Draw("SAME HIST")
S_T2tt  .Draw("SAME HIST")
S_T5ttcc.Draw("SAME HIST")
#leg = ROOT.TLegend(0.55,0.67,0.95,0.92, BOX) # in case no ratio
leg = ROOT.TLegend(0.55,0.56,0.95,0.86, BOX)
leg.SetNColumns(2)
#leg.AddEntry(S_data,     "#color[1]{Data}",                      "pe")
leg.AddEntry(S_T1tttt,   "#color[619]{T1tttt}",                   "l")
leg.AddEntry(S_TT_nom,   "#color[418]{t#bar{t} or single t}",     "f")
leg.AddEntry(S_T2tt,     "#color[401]{T2tt}",                     "l")
leg.AddEntry(S_MJ_nom,   "#color[618]{Multijet}",                 "f")
leg.AddEntry(S_T5ttcc,   "#color[601]{T5ttcc}",                   "l")
leg.AddEntry(S_WJ_nom,   "#color[633]{W(l#nu)+jets}",  "f")
leg.AddEntry(0,          "",                                      "")
leg.AddEntry(S_ZI_nom,   "#color[433]{Z(#nu#nu)+jets}","f")
leg.AddEntry(0,          "",                                      "")
leg.AddEntry(S_OT_nom,   "#color[865]{Other}",                    "f")
#leg.AddEntry(0,          "",                                      "")
#leg.AddEntry(S_syst_err, "Stat. + syst. unc.",                    "f");
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
#leg.SetTextSize(0.03)
leg.Draw("SAME")
draw_mr_bins([S_data, S_stack], ymin,ymax, combine_bins, keep, mrbins_TeV, r2bins)
ROOT.gPad.Update()
ROOT.gPad.RedrawAxis()
if opt.energy==14 and "WAna_nj45" in opt.box:
    can = add_stack_ratio_plot(can, 0,S_data.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins, { 3 : 4e3})
elif opt.energy==14 and "WAna_nj6" in opt.box:
    can = add_stack_ratio_plot(can, 0,S_data.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins, { 2 : 1.5e4, 3: 3e3 })
elif opt.energy==27 and "TopAna" in opt.box:
    can = add_stack_ratio_plot(can, 0,S_data.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins, { 4 : 4e4 })
else:
    can = add_stack_ratio_plot(can, 0,S_data.GetNbinsX(), keep, 1, combine_bins, mrbins_TeV, r2bins)
add_cms_era(can, ERA, keep, opt.energy, lumi, prefix)
save_plot(can, "", "Plots/result/"+DATE+"/"+can.GetName().replace("_Ratio",""), 1)
f.Close()
#sys.exit()

# Now save a different root file for each signal point
if opt.nocards or opt.TEST!=0: sys.exit()

print "Looping on Signal points and creating data cards"
if not os.path.exists("syst_"+opt.dir+"/cards"):
    special_call(["mkdir", "-p", "syst_"+opt.dir+"/cards"], 0)

cards = []
for signal_syst in S_signal_extension:
    scan_point = signal_syst[0].GetName()[:-4].replace("MRR2_S_signal_","").replace(BIN,"")
    #root_filename = "syst_"+opt.dir+"/cards/RazorBoost_"+opt.box+"_"+opt.model+"_"+scan_point+".root"
    root_filename = "syst_"+opt.dir+"/cards/RazorBoost_SMS-"+opt.model+"_"+scan_point+"_"+opt.box+".root"
    if not opt.nocards:
        fout = ROOT.TFile.Open(root_filename,"RECREATE")
        #print "  Creating root file: "+root_filename
        # Add signal systematics
        for i in range(0, len(systematics)-2):
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
        if opt.scenario==0:
            card.write(
'''
------------------------------------------------------------
'''
            )
        elif opt.scenario==1:
            card.write(
'''
------------------------------------------------------------
lumi			lnN	1.01	1.01	1.01	1.01	1.01	1.01
toppt			shape	-	1.0	1.0	1.0	1.0	1.0
isr			shape	1.0	-	-	-	-	-
jes			shape	1.0	1.0	1.0	1.0	1.0	1.0
jer 			shape	1.0	1.0	1.0	1.0	1.0	1.0
met			shape	1.0	1.0	1.0	1.0	1.0	1.0
facrenscale		shape	1.0	1.0	1.0	1.0	1.0	1.0
alphas			shape	1.0	1.0	1.0	1.0	1.0	1.0
elereco			shape	1.0	1.0	1.0	1.0	1.0	1.0
eleid			shape	1.0	1.0	1.0	1.0	1.0	1.0
eleiso			shape	1.0	1.0	1.0	1.0	1.0	1.0
muontrk			shape	1.0	1.0	1.0	1.0	1.0	1.0
muonidiso		shape	1.0	1.0	1.0	1.0	1.0	1.0
lostlep			shape	1.0	1.0	1.0	1.0	1.0	1.0
btag			shape	1.0	1.0	1.0	1.0	1.0	1.0
wtag			shape	1.0	1.0	1.0	1.0	1.0	1.0
wmistag			shape	1.0	1.0	1.0	1.0	1.0	1.0
wmasstag		shape	-	-	-	1.0	1.0	-
wantitag		shape	-	-	1.0	-	-	-
toptag			shape	1.0	1.0	1.0	1.0	1.0	1.0
topmistag		shape	1.0	1.0	1.0	1.0	1.0	1.0
top0bmasstag		shape	-	-	-	1.0	-	-
topmasstag		shape	-	-	-	-	1.0	-
topantitag		shape	-	-	1.0	-	-	-
purity			shape	-	-	-	-	1.0	-
dirfrac			shape	-	-	-	-	1.0	-
'''
#'''
#extrap			shape	-	1.0	1.0	1.0	1.0	1.0
#trigger			shape	1.0	1.0	1.0	1.0	1.0	1.0
#ak8scale		shape	-	1.0	1.0	1.0	1.0	1.0
#doubleratio		shape	-	-	-	-	1.0	-
#elefastsim		shape	1.0	-	-	-	-	-
#muonfastsim		shape	1.0	-	-	-	-	-
#btagfastsim		shape	1.0	-	-	-	-	-
#wtagfastsim		shape	1.0	-	-	-	-	-
#wmistagfastsim		shape	1.0	-	-	-	-	-
#toptagfastsim		shape	1.0	-	-	-	-	-
#topmistagfastsim	shape	1.0	-	-	-	-	-
#leptonest		shape	-	-	-	-	1.0	-
#            card.write("dytoll\t\t\tlnN\t-\t-\t-\t-\t%2.2f\t-\n" % (1.0+dy_syst))
#'''
            )
            if "T1tttt" in opt.model:
                card.write("pileup\t\t\tlnN\t0.995/1.005\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\n")
            elif "T5ttcc" in opt.model:
                card.write("pileup\t\t\tlnN\t\t0.987/1.013\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\n")
            elif "T2tt" in opt.model:
                card.write("pileup\t\t\tlnN\t0.99/1.01\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\t0.978/1.026\n")
            card.write("qcd\t\t\tlnN\t-\t-\t%2.3f\t-\t-\t-\n" % (1.0+qcd_syst))
            card.write("dytoll\t\t\tlnN\t-\t-\t-\t-\t%2.3f\t-\n" % (1.0+dy_syst))
        elif opt.scenario==2:
            card.write(
'''
------------------------------------------------------------
lumi			lnN	1.025	1.025	1.025	1.025	1.025	1.025
pileup			shape	1.0	1.0	1.0	1.0	1.0	1.0
toppt			shape	-	1.0	1.0	1.0	1.0	1.0
isr			shape	1.0	-	-	-	-	-
jes			shape	1.0	1.0	1.0	1.0	1.0	1.0
jer 			shape	1.0	1.0	1.0	1.0	1.0	1.0
met			shape	1.0	1.0	1.0	1.0	1.0	1.0
trigger			shape	1.0	1.0	1.0	1.0	1.0	1.0
facrenscale		shape	1.0	1.0	1.0	1.0	1.0	1.0
alphas			shape	1.0	1.0	1.0	1.0	1.0	1.0
ak8scale		shape	-	1.0	1.0	1.0	1.0	1.0
elereco			shape	1.0	1.0	1.0	1.0	1.0	1.0
eleid			shape	1.0	1.0	1.0	1.0	1.0	1.0
eleiso			shape	1.0	1.0	1.0	1.0	1.0	1.0
elefastsim		shape	1.0	-	-	-	-	-
muontrk			shape	1.0	1.0	1.0	1.0	1.0	1.0
muonidiso		shape	1.0	1.0	1.0	1.0	1.0	1.0
muonfastsim		shape	1.0	-	-	-	-	-
lostlep			shape	1.0	1.0	1.0	1.0	1.0	1.0
btag			shape	1.0	1.0	1.0	1.0	1.0	1.0
btagfastsim		shape	1.0	-	-	-	-	-
wtag			shape	1.0	1.0	1.0	1.0	1.0	1.0
wtagfastsim		shape	1.0	-	-	-	-	-
wmistag			shape	1.0	1.0	1.0	1.0	1.0	1.0
wmistagfastsim		shape	1.0	-	-	-	-	-
wmasstag		shape	-	-	-	1.0	1.0	-
wantitag		shape	-	-	1.0	-	-	-
toptag			shape	1.0	1.0	1.0	1.0	1.0	1.0
toptagfastsim		shape	1.0	-	-	-	-	-
topmistag		shape	1.0	1.0	1.0	1.0	1.0	1.0
topmistagfastsim	shape	1.0	-	-	-	-	-
top0bmasstag		shape	-	-	-	1.0	-	-
topmasstag		shape	-	-	-	-	1.0	-
topantitag		shape	-	-	1.0	-	-	-
extrap			shape	-	1.0	1.0	1.0	1.0	1.0
doubleratio		shape	-	-	-	-	1.0	-
purity			shape	-	-	-	-	1.0	-
dirfrac			shape	-	-	-	-	1.0	-
leptonest		shape	-	-	-	-	1.0	-
'''
            )
            card.write("qcd\t\t\tlnN\t-\t-\t%2.2f\t-\t-\t-\n" % (1.0+qcd_syst))
            card.write("dytoll\t\t\tlnN\t-\t-\t-\t-\t%2.2f\t-\n" % (1.0+dy_syst))
        
        card.close()

print "All data cards ready"
print "Done."
