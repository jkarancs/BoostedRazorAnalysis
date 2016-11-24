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

There is a py script to run the Analyzer and use the results to create plots all in once
on a few datasets (--test option) and fewer events (--quick option)
```Shell
python scripts/run_all.py --quick --test --run
```

Or run a full analysis on all samples (--run option) and similarly create plots in the end
```Shell
python scripts/run_all.py --full --plot --run
```

All options:
   * --run: without specifying this iption, the script does nothing (dry run, prints commands to be run only)
   * --test: run a test with a few datasets: TTbar, QCD, small JetHT dataset
   * --quick: run a quick test with 1/100 events
   * --full: run on all datasets
   * --plot: after the Analyses, use the Plotter to produce plots in the output root files (also uses the Analysis methods)
   * --skim: Make a skimmed ntuple in the directory specified inside the script

Run systematics (first, set varySystematics = true in settings.h)  
when running systematics, make sure you define one histogram for each variation
syst.index in the Analyzer.cc files can be used to select a histogram from a vector

You can remake the systematics file with:
```Shell
python scripts/write_syst.py <nline>
```

## Contact Information

Janos Karancsi (janos.karancsi@cern.ch)

