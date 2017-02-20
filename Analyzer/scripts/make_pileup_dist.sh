#!/bin/bash

# 2015
#CERT_dir=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV
#PU_file=$CERT_dir/PileUp/pileup_JSON_11-19-2015.txt # = pileup_latest.txt 13 May 2016

# Mar 02 Silver JSON
#JSON=$CERT_dir/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt
#OUT_dir=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/pileup/Mar02_Silver_JSON

# Mar 02 Golden JSON
#JSON=$CERT_dir/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt
#OUT_dir=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/pileup/Mar02_Golden_JSON

# 2015 MB XSec
#XSec_MB_down=65550.0  # -5% xsec
#XSec_MB_nom=69000.0   # nominal xsec
#XSec_MB_up=72450.0    # +5% xsec


# 2016
CERT_dir=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV
PU_file=$CERT_dir/PileUp/pileup_latest.txt # Dec06 is latest (size 9942473)

# ICHEP16 Golden JSON
#JSON=$CERT_dir/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt
#OUT_dir=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/pileup/ICHEP16_Golden_JSON
# Oct21 Golden JSON
#JSON=$CERT_dir/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON.txt
#OUT_dir=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/pileup/Oct21_Golden_JSON
# Dec02 Final ReReco Golden JSON
JSON=$CERT_dir/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
OUT_dir=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/pileup/Dec02_Golden_JSON

# 2016 MB XSec
# https://hypernews.cern.ch/HyperNews/CMS/get/b2g/883/1.html
# --> https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/613/2/1/1/1.html
XSec_MB_down=66016.8  # -4.6% xsec
XSec_MB_nom=69200.0   # nominal xsec
XSec_MB_up=72383.2    # +4.6% xsec


cmsenv
cd $CMSSW_BASE/src

if [ ! -d RecoLuminosity ]; then
    mkdir -p ../tmp
    mv * ../tmp/
    git cms-addpkg RecoLuminosity/LumiDB
    mv ../tmp/* .
    rm -r ../tmp
    scram b -j 20
fi

mkdir -p $OUT_dir
RecoLuminosity/LumiDB/scripts/pileupCalc.py --inputLumiJSON $PU_file -i $JSON --calcMode true --minBiasXsec $XSec_MB_up --maxPileupBin 100 --numPileupBins 100 $OUT_dir/data_pileup_up.root
RecoLuminosity/LumiDB/scripts/pileupCalc.py --inputLumiJSON $PU_file -i $JSON --calcMode true --minBiasXsec $XSec_MB_nom --maxPileupBin 100 --numPileupBins 100 $OUT_dir/data_pileup.root
RecoLuminosity/LumiDB/scripts/pileupCalc.py --inputLumiJSON $PU_file -i $JSON --calcMode true --minBiasXsec $XSec_MB_down --maxPileupBin 100 --numPileupBins 100 $OUT_dir/data_pileup_down.root

cd -
