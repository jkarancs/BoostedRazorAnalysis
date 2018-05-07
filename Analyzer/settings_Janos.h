// VER 0 - Spring16 MC, ICHEP Dataset
// VER 1 - Moriond17 datasets
// VER 2 - Moriond17 + 03Feb2017 ReMiniAOD datasets
// VER 3 - VER2 with updates (added photons, corridor event veto, ISR, gen MET)
// VER 4 - VER3 with updates (added taus, updated isolations)
// SKIM - 1: save skimmed ntuple, 0: run on already skimmed ntuple
#define VER     4
#define SKIM    0
#define SYST    1
#define TOP     1

#if VER == 1
#include "common/DataStruct_Jan12.h"
#elif VER == 2
#include "common/DataStruct_May10.h"
#elif VER == 3
#include "common/DataStruct_Sep26.h"
#elif VER == 4
#include "common/DataStruct_Nov30.h"
#endif
#include "common/treestream.h"
#include "Analysis_Janos.h" // Specify here the implementations for your Analysis
//#include "Analysis_T.h" // Specify here the implementations for your Analysis

struct settings {
#if VER == 1
#if SKIM == 1
#include "common/selectVariables_skim_Jan12.h"
#else
#include "common/selectVariables_fast_Jan12.h"
#endif

#elif VER == 2
#if SKIM == 1
//#include "common/selectVariables_skim_May10.h"
#include "common/selectVariables_skim_May10_photon.h"
#else
#include "common/selectVariables_fast_May10.h"
//#include "common/selectVariables_fast_May10_photon.h"
#endif

#elif VER == 3
#if SKIM == 1
#include "common/selectVariables_skim_Sep26.h"
#else
#include "common/selectVariables_fast_Sep26.h"
#endif

#elif VER == 4
#if SKIM == 1
#include "common/selectVariables_skim_Nov30.h"
#else
#include "common/selectVariables_fast_Nov30.h"
#endif

#endif

  //-----------------------------------------------------------------------------
  // -- Constants
  //-----------------------------------------------------------------------------
  settings() :
    runOnSkim                ( 1-SKIM),
    saveSkimmedNtuple        ( SKIM  ),
    doTopPtReweighting       ( true  ),
    doISRReweighting         ( true  ),
    doPileupReweighting      ( true  ),
    doAK8JetPtRescaling      ( true  ),
    doNJetReweighting        ( true  ),
    applySmearing            ( true  ),
    applyScaleFactors        ( true  ),
    nSigmaScaleFactors       ( 22    ), // Count the number of sigmas you use in Analysis_*.h - 4 ele, 3 mu, 6 W, 2b, 7 top
    varySystematics          ( SYST  ),
    systematicsFileName      ( "systematics/2018_05_02_1SigmaUpDown_NoPdf.txt" ),
//  systematicsFileName      ( "systematics/2017_12_26_AllUpDown_NoPdf.txt" ),
//  systematicsFileName      ( "systematics/2018_03_19_JESOnly.txt" ),
    useJSON                  ( false ), // by default: no need to apply, but can be useful if some lumisections need to be excluded additionally
#if VER != 0
    jsonFileName             ( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/"
			       "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt" ),
    pileupDir                ( "pileup/Dec02_Golden_JSON/" ),
    intLumi                  ( 35867 /* brilcalc - Dec02 Golden JSON */ ), // Tot int lumi in (pb^-1),
    lumiUncertainty          ( 0.025  ),
#endif
    useXSecFileForBkg        ( true   ), // true: use file below, false: use value in the ntuple (evt_XSec)
    xSecFileName             ( "common/BackGroundXSec.txt" ) {};
  ~settings(){};

  const bool runOnSkim;
  const bool saveSkimmedNtuple;
  const bool doTopPtReweighting;
  const bool doISRReweighting;
  const bool doPileupReweighting;
  const bool doAK8JetPtRescaling;
  const bool doNJetReweighting;
  const bool applySmearing;
  const bool applyScaleFactors;
  const int  nSigmaScaleFactors;
  const bool varySystematics;
  const std::string systematicsFileName;
  const bool useJSON;
  const std::string jsonFileName;
  const std::string pileupDir;
  const double intLumi;
  const double lumiUncertainty;
  const bool useXSecFileForBkg;
  const std::string xSecFileName;

} settings;
