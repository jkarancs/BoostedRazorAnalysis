# BoostedRazorAnalysis
Run II Search for Supersymmetry using Razor variables in the boosted regime

## Checkout recipe

```Shell
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_15_patch1
cd CMSSW_7_4_15_patch1/src
cmsenv
git clone https://github.com/jkarancs/BoostedRazorAnalysis
scram b -j 20
```

## Analyzer example

```Shell
cd BoostedRazorAnalysis
eosmount data/eos_mount
cd Analyzer
make Analyzer
./Analyzer ../data/filelists/backgrounds/TTJets_madgraph.txt 0 0 2500 Bkg_TTJets_madgraph.root Bkg_TTJets_madgraph Pileup_True
```

Use this script to print command to run on all samples

```Shell
source scripts/run_all.csh
```

Remarks:
    - The cross-sections and total weight is taken straight from ntuple files (set 2nd and 3rd argument to 0)
    - The Analyzer program is used to count events that pass the Analysis selection
    - Additionally it also saves a skimmed TTree with same content for the selected events

## Contact Information

Janos Karancsi (janos.karancsi@cern.ch)

