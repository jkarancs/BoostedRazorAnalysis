# BoostedRazorAnalysis
Run II Search for Supersymmetry using Razor variables in the boosted regime

## Checkout recipe

```Shell
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv
git cms-addpkg RecoLuminosity/LumiDB
git clone https://github.com/jkarancs/BoostedRazorAnalysis
scram b -j 20
```

## Run the Analyzer

When running for the first time, run the setup script in order  
to mount eos on lxplus, create softlinks and create/update filelists
```Shell
cd BoostedRazorAnalysis/Analyzer
python scripts/setup.py
```

Before running check/set the options for the Analyzer in settings.h

More info about Analyzer:
   * The Analyzer program is basically an event looper, histograms/methods/cuts etc. are defined in Analysis_[Name].h
   * All settings (printed also while running) are defined in settings_[Name].h it also includes the Analysis_[Name].h
   * common methods (for all analyses/studies) are defined in common/AnalysisBase.h
   * The cross-sections and total weight is taken straight from ntuple files
   * Reweighting and systematics weight methods are also given in AnalysisBase
   * Additionally there's an option (in settings_[Name].h) to save a skimmed TTree with same content for the selected events
   * counts are saved for all common and specific analysis cuts in the order they are defined

Run your anaylsis code with
```Shell
make
./Analyzer <output filename> <filelist.txt or root file(s) ... >
eg:
./Analyzer Bkg_TTJets_madgraph.root filelists/backgrounds/TTJets_madgraph.txt
```

There is a py script to run the Analyzer to run on filelists of datasets
with lot of options to use
```Shell
python scripts/run_all.py --help
```

Run a full analysis on all samples (not using --run option prints commands)
```Shell
python scripts/run_all.py --full --run
python scripts/run_all.py --full --run --batch --nevt=2000000
```

Run systematics (first, set varySystematics = true in settings.h)  
when running systematics, make sure you define one histogram for each variation
syst.index in the Analyzer.cc files can be used to select a histogram from a vector

You can remake the systematics file with:
```Shell
python scripts/write_syst.py <nline>
```

## Contact Information

Janos Karancsi (janos.karancsi@cern.ch)

