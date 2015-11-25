# BoostedRazorAnalysis
Run II Search for Supersymmetry using Razor variables in the boosted regime

# Checkout recipe
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_15_patch1
cd CMSSW_7_4_15_patch1/src
cmsenv
git clone https://github.com/jkarancs/BoostedRazorAnalysis
scram b -j 20

# Analyzer example
cd BoostedRazorAnalysis
eosmount data/eos_mount
cd Analyzer
make Analyzer
./Analyzer ../data/filelists/backgrounds/TTJets_HT.txt 1 1 1 out.root TTjets_HT Pileup_True
