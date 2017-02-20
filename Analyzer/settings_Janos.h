// VER 0 - Spring15 MC, ICHEP Dataset
// VER 1 - Moriond17 datasets
#define VER 0

#include "common/DataStruct_Oct31.h"
#include "common/treestream.h"
#include "Analysis_Janos.h" // Specify here the implementations for your Analysis

struct settings {

#include "common/selectVariables_fast_Oct31.h"

  //-----------------------------------------------------------------------------
  // -- Constants
  //-----------------------------------------------------------------------------
  settings() :
    runOnSkim                ( true  ),
    saveSkimmedNtuple        ( false ),
    doPileupReweighting      ( true  ),
    scaleQCD                 ( false ),
    doHTReweighting          ( false ),
    applySmearing            ( true  ),
    applyScaleFactors        ( true  ),
    nSigmaScaleFactors       ( 13    ), // Count the number of sigmas you use in Analysis_*.h - 4 ele, 7 mu, 1 W, 1 b
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
    useXSecFileForBkg        ( false ), // true: use file below, false: use value in the ntuple (evt_XSec)
    xSecFileName             ( "common/BackGroundXSec.txt" )
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
  const bool scaleQCD;
  const bool doHTReweighting;
  const bool applySmearing;
  const bool applyScaleFactors;
  const int  nSigmaScaleFactors;
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
  std::vector<std::string> totWeightHistoNamesSignal;

} settings;
