import FWCore.ParameterSet.Config as cms
import ROOT

numPileupBins = 100
maxPileupBin  = 100

# 2015 - 76X
#from SimGeneral.MixingModule.mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi import mix
#outputfile = "pileup/Mar02_Golden_JSON/mc_pileup.root"

# 2016 - 80X
from SimGeneral.MixingModule.mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi import mix
outputfile = "pileup/ICHEP16_Golden_JSON/mc_pileup.root"

# Write histogram to specified file
f = ROOT.TFile.Open(outputfile, "recreate")
h = ROOT.TH1D("pileup","pileup", numPileupBins, 0, maxPileupBin)

for bin in range (1, len(mix.input.nbPileupEvents.probValue)+1):
    h.SetBinContent(bin, mix.input.nbPileupEvents.probValue[bin-1])

h.Write()
f.Close()


#  Usage:
#      cd /data/jkarancs/CMSSW/ntuple
#      cmsenv
#      cd /data/jkarancs/CMSSW/SusyAnalysis/Analyzer/CMSSW_7_6_3_patch2/src/BoostedRazorAnalysis/Analyzer
#      python /data/jkarancs/CMSSW/SusyAnalysis/Analyzer/CMSSW_7_6_3_patch2/src/BoostedRazorAnalysis/Analyzer/scripts/make_pileup_dist_mc.py
