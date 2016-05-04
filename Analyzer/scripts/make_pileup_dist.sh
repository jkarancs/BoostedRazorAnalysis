CERT_dir=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV
PU_file=$CERT_dir/PileUp/pileup_JSON_11-19-2015.txt # = pileup_latest.txt 13 May 2016
JSON=$CERT_dir/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt
OUT_dir=$CMSSW_BASE/src/BoostedRazorAnalysis/Analyzer/pileup/Mar02_Silver_JSON
XSec_MB_up=65550.0      # 5% more pile-up
XSec_MB_nom=69000.0     # nominal
XSec_MB_down=72450.0    # 5% less pile-up

cmsenv
cd $CMSSW_BASE/src

if [ ! -d RecoLuminosity ]; then
    mkdir ../tmp
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
