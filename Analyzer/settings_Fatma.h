// VER 0 - Spring15 MC, ICHEP Dataset
// VER 1 - Moriond17 datasets
// SKIM - 1: save skimmed ntuple, 0: run on already skimmed ntuple
#define VER     1
#define SKIM    0

#if VER == 0
#include "common/DataStruct_Oct31.h"
#else
#include "common/DataStruct_Jan12.h"
#endif
#include "common/treestream.h"
#include "Analysis_T.h" // Specify here the implementations for your Analysis

struct settings {
#if VER == 0
#if SKIM == 1
#include "common/selectVariables_skim_Oct31.h"
#else
#include "common/selectVariables_fast_Oct31.h"
#endif

#else
#if SKIM == 1
#include "common/selectVariables_skim_Jan12.h"
#else
#include "common/selectVariables_fast_Jan12.h"
#endif
#endif

  //-----------------------------------------------------------------------------
  // -- Constants
  //-----------------------------------------------------------------------------
  settings() :
    runOnSkim                ( 1-SKIM),
    saveSkimmedNtuple        ( SKIM  ),
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
    useJSON                  ( false ), // by default: no need to apply, but can be useful if some lumisections need to be excluded additionally
#if VER == 0
    jsonFileName             ( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/"
			       "Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON.txt" ),
    pileupDir                ( "pileup/Oct21_Golden_JSON/" ),
    intLumi                  ( 27655.802 /* brilcalc - Oct21 Golden JSON */ ), // Tot int lumi in (pb^-1),
    lumiUncertainty          ( 0.062 ),
#else
    jsonFileName             ( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/"
			       "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt" ),
    pileupDir                ( "pileup/Dec02_Golden_JSON/" ),
    intLumi                  ( 35867 /* brilcalc - Dec02 Golden JSON */ ), // Tot int lumi in (pb^-1),
    lumiUncertainty          ( 0.026 ),
#endif
    useXSecFileForBkg        ( true  ), // true: use file below, false: use value in the ntuple (evt_XSec)
    xSecFileName             ( "common/BackGroundXSec.txt" )
  {
    totWeightHistoNamesSignal.push_back(runOnSkim ? "totweight_T1tttt" : "EventCounter/h_totweight_T1tttt"); // lsp mass vs gluino mass scan, also used for T5ttcc and T5tttt
    totWeightHistoNamesSignal.push_back(runOnSkim ? "totweight_T2tt"   : "EventCounter/h_totweight_T2tt");   // T2tt
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
