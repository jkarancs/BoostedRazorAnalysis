setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_15_patch1
cd CMSSW_7_4_15_patch1/src
cmsenv
git cms-addpkg RecoLuminosity/LumiDB
scram b -j 20
ln -s /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV CERT15

mkdir Nov13_Silver_JSON
RecoLuminosity/LumiDB/scripts/pileupCalc.py -i CERT15/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt --inputLumiJSON CERT15/PileUp/pileup_JSON_11-19-2015.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 100 --numPileupBins 100 Nov13_Silver_JSON/data_pileup.root
# increase cross-section by 5% (pileup goes down)
RecoLuminosity/LumiDB/scripts/pileupCalc.py -i CERT15/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt --inputLumiJSON CERT15/PileUp/pileup_JSON_11-19-2015.txt --calcMode true --minBiasXsec 72450 --maxPileupBin 100 --numPileupBins 100 Nov13_Silver_JSON/data_pileup_down.root
# decrease cross-section by 5% (pileup goes up)
RecoLuminosity/LumiDB/scripts/pileupCalc.py -i CERT15/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt --inputLumiJSON CERT15/PileUp/pileup_JSON_11-19-2015.txt --calcMode true --minBiasXsec 65550 --maxPileupBin 100 --numPileupBins 100 Nov13_Silver_JSON/data_pileup_up.root

# then just copy this to BoostedRazorAnalysis/Analyzer/pileup
