#include "common/DataStruct.h"
#include "common/treestream.h"

#include "Analysis_T.h" // Specify here the implementations for your Analysis

struct settings {

#include "common/selectVariables_skim.h"
//#include "common/selectVariables_fast.h"

  //-----------------------------------------------------------------------------
  // -- Constants
  //-----------------------------------------------------------------------------
  settings() :
    runOnSkim                ( true  ),
    saveSkimmedNtuple        ( false ),
    doPileupReweighting      ( true  ),
    applyWTagSF              ( false ),
    applyBTagSF              ( false ),
    //applyHadTopTagSF         ( false ),
    scaleQCD                 ( false ),
    doHTReweighting          ( false ),
    varySystematics          ( false ),
    systematicsFileName      ( "systematics/2016_10_31_1SigmaUpDown_NoPdf.txt" ),
    treeName                 ( runOnSkim ? "B2GTree"   : "B2GTTreeMaker/B2GTree" ),
    totWeightHistoName       ( runOnSkim ? "totweight" : "EventCounter/totweight" ), // saved in ntuple
    mcPileupHistoName        ( runOnSkim ? "pileup_mc" : "EventCounter/pileup" ),    // saved in ntuple
    useJSON                  ( false ), // Default is not appliying any JSON, but Golden JSON (tighter selection) can be applied on top of the Default Silver JSON
    jsonFileName             ( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON.txt" ),
    pileupDir                ( "pileup/Oct21_Golden_JSON/" ),
    intLumi                  ( 27655.802 /* brilcalc - Oct21 Golden JSON */ ), // Tot int lumi in (pb^-1),
    lumiUncertainty          ( 0.062 ),
    useXSecFileForBkg        ( true  ), // true: use file below, false: use value in the ntuple (evt_XSec)
    xSecFileName             ( "common/BackGroundXSec.txt" ),
    triggerEffScaleFactor    (  1.00 ),
    triggerEffUncertainty    (  0.01 )
  {
    if (runOnSkim) {
      totWeightHistoNamesSignal.push_back("totweight_T1tttt"); // lsp mass vs gluino mass scan, also used for T5ttcc and T5tttt
      totWeightHistoNamesSignal.push_back("totweight_T2tt");   // T2tt
    } else {
      totWeightHistoNamesSignal.push_back("EventCounter/h_totweight_T1tttt");
      totWeightHistoNamesSignal.push_back("EventCounter/h_totweight_T2tt");
    }
  };
  ~settings(){};

  const bool runOnSkim;
  const bool saveSkimmedNtuple;
  const bool doPileupReweighting;
  const bool applyWTagSF;
  const bool applyBTagSF;
  //const bool applyHadTopTagSF;
  const bool scaleQCD;
  const bool doHTReweighting;
  const bool varySystematics;
  const std::string systematicsFileName;
  const std::string treeName;
  const std::string totWeightHistoName;
  const std::string mcPileupHistoName;
  const bool useJSON;
  const std::string jsonFileName;
  const std::string pileupDir;
  const double intLumi;
  const double lumiUncertainty;
  const bool useXSecFileForBkg;
  const std::string xSecFileName;
  const double triggerEffScaleFactor;
  const double triggerEffUncertainty;
  std::vector<std::string> totWeightHistoNamesSignal;

} settings;
