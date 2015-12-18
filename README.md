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

## Run the Analyzer

When running for the first time, run the setup script in order  
to mount eos on lxplus, create softlinks and create/update filelists
```Shell
cd BoostedRazorAnalysis/Analyzer
source scripts/setup.csh
make
```

Before running check/set the options for the Analyzer in settings.h

More info about Analyzer:
   * The Analyzer program is basically an event looper, histograms/methods/cuts etc. are defined in [Name]_Analysis.h
   * All settings (printed also while running) are defined in settings.h it also includes the Analysis code you want
   * common methods (for all analyses/studies) are defined in common/AnalysisBase.h
   * The cross-sections and total weight is taken straight from ntuple files
   * Reweighting and systematics weight methods are also given in AnalysisBase
   * Additionally there's an option (in settings.h) to save a skimmed TTree with same content for the selected events
   * counts are saved for all common and specific analysis cuts in the order they are defined

Run your anaylsis code with
```Shell
./Analyzer <output filename> <filelist.txt or root file(s) ... >
eg:
./Analyzer Bkg_TTJets_madgraph.root filelists/backgrounds/TTJets_madgraph.txt
```


Run systematics (first, set doSystematics = true in settings.h)  
when running systematics, additionally must specify options:  
systematicsFileName and numSyst (i.e. the line to read in syst file)
```Shell
./Analyzer systematicsFileName=systematics/2015_12_18.txt numSyst=1   Bkg_TTJets_madgraph.root filelists/backgrounds/TTJets_madgraph.txt
./Analyzer systematicsFileName=systematics/2015_12_18.txt numSyst=2   Bkg_TTJets_madgraph.root filelists/backgrounds/TTJets_madgraph.txt
...
./Analyzer systematicsFileName=systematics/2015_12_18.txt numSyst=100 Bkg_TTJets_madgraph.root filelists/backgrounds/TTJets_madgraph.txt
```

You can remake the systematics file with:
```Shell
python scripts/write_syst.py <nline>
```

## Contact Information

Janos Karancsi (janos.karancsi@cern.ch)

