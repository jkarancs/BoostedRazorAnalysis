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
cd BoostedRazorAnalysis/Analyzer
source scripts/setup.csh
make
./Analyzer Bkg_TTJets_madgraph.root filelists/backgrounds/TTJets_madgraph.txt
```
scripts/setup.csh is used only the first time to mount eos on lxplus, create softlinks, update filelists

```Shell
source scripts/run_all.csh
```

scripts/run_all.csh shows commands to run on everything

Remarks:
    - The Analyzer program is basically an event looper, histograms/methods/cuts etc. are defined in [Name]_Analysis.h
    - All settings (printed also while running) are defined in settings.h it also includes the Analysis code you want
    - common methods (for all analyses/studies) are defined in common/AnalysisBase.h
    - The cross-sections and total weight is taken straight from ntuple files
    - Additionally there's an option to save a skimmed TTree with same content for the selected events
    - counts are saved for all common and specific analysis cuts in the order they are defined

## Contact Information

Janos Karancsi (janos.karancsi@cern.ch)

