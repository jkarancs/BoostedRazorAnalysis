# BoostedRazorAnalysis
Run II Search for Supersymmetry using Razor variables in the boosted regime

## Checkout recipe

```Shell
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_0_28
cd CMSSW_8_0_28src
cmsenv
git cms-addpkg RecoLuminosity/LumiDB
git clone https://github.com/jkarancs/BoostedRazorAnalysis
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-https.sh)
scram b -j 20
cd BoostedRazorAnalysis/Analyzer
```

## Run the Analyzer

When running for the first time, run the setup script in order  
to mount eos on lxplus, create softlinks and create/update filelists
```Shell
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
make clean; make Analyzer
./Analyzer <output filename> <filelist.txt or root file(s) ... >
eg:
./Analyzer test.root ntuple/Latest/TT_powheg-pythia8/Skim_1.root
./Analyzer ttbar.root filelists/backgrounds/TT_powheg-pythia8.txt
```

There is a py script to run the Analyzer to run on filelists of datasets
with lot of options to use (use option --help)
```Shell
python scripts/run_all.py --help
```

Run a quick (1/Nth events) interactive test, using 4 cpus on selected/all samples
```Shell
python scripts/run_all.py --full --nproc=4 --quick=1000 --outdir=test/quicktest --run
python scripts/run_all.py --plot --nproc=4 --quick=100 --outdir=test/quicktest_data_bkg --run filelists/data/JetHT_Run2016* filelists/backgrounds/*
```

If you are on lxplus, it is highly recommended to run full analysis on the batch
First a quick test on a faster queue
N.B: you need to generate a temp filelist for split jobs, this is a bit slow, so next time you can reuse previous list
For this example, please make sure the varySystematics is set to False in the settings file you include in the Analyzer.cc
```Shell
python scripts/run_all.py --full --plot --batch --run --outdir=results/run_2018_03_23 --queue=8nh --optim --nevt=2000000
```
Run systematics (first, set varySystematics = true in settings.h)  
when running systematics, make sure you define one histogram for each variation
syst.index in the Analyzer.cc files can be used to select a histogram from a vector

```Shell
python scripts/run_all.py --full --plot --batch --nohadd --run --outdir=results/run_2018_03_24_syst --queue=2nd --optim --nevt=300000
```

## Contact Information

Janos Karancsi (janos.karancsi@cern.ch)

