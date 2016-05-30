import FWCore.ParameterSet.Config as cms
import ROOT

from SimGeneral.MixingModule.mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi import mix
numPileupBins = 100
maxPileupBin  = 100
outputfile = "pileup/Mar02_Silver_JSON/mc_pileup.root"

# Write histogram to specified file
f = ROOT.TFile.Open(outputfile, "recreate")
h = ROOT.TH1D("pileup","pileup", numPileupBins, 0, maxPileupBin)

for bin in range (1, len(mix.input.nbPileupEvents.probValue)+1):
    h.SetBinContent(bin, mix.input.nbPileupEvents.probValue[bin-1])

h.Write()
f.Close()

