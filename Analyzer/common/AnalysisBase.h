#ifndef VER
#define VER 0
#endif

#include <iostream>
#include <functional>
#include <map>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TStopwatch.h"
#include "TRandom3.h"
#include <unistd.h>

#include "utils.h"
#include "GluinoXSec.h"
#include "StopXSec.h"
#include "Razor.h"

#include "BTagCalibrationStandalone.cpp"

// _____________________________________________________________
//        AnalysisBase: Methods common in all analysis

class AnalysisBase
{
public:
  AnalysisBase(const bool& isData, const bool& isSignal, const std::string& dirname) :
    isData(isData),
    isSignal(isSignal),
    isBackground(!isData&&!isSignal),
    sample(dirname)
  {
    sw_1_  = new TStopwatch;
    sw_1k_  = new TStopwatch;
    sw_10k_ = new TStopwatch;
    sw_job_ = new TStopwatch;

    syst = {
      "lumi", "toppt", "isr", "pileup",
      "alphas", "facscale", "renscale", "facrenscale", "lostlep",
      "trigger", "jes", "jer", "met", "ak8scale",
      "elereco", "eleid", "eleiso", "elefastsim",
      "muontrk", "muonidiso", "muonfastsim",
      "btag", "btagfastsim",
      "wtag", "wtagfastsim", "wmistag", "wmistagfastsim", "wmasstag","wantitag",
      "toptag", "toptagfastsim", "topmistag", "topmistagfastsim",
      "top0bmasstag","topmasstag","topantitag"
    };
    //for (int i=1; i<=100; ++i) {
    //  std::stringsteam ss;
    //  ss<<"pdf"<<i;
    //  bkg_syst.push_bask(ss.str());
    //}

    other_trigger_eff = 1.0;
  }
  ~AnalysisBase() {
    delete sw_1_;
    delete sw_1k_;
    delete sw_10k_;
    delete sw_job_;
  }

  typedef struct Cut { std::string name; std::function<bool()> func; } Cut;
  std::vector<Cut> baseline_cuts;

  // Functions used by the Analyzer
  void define_preselections(const DataStruct&);

  void calculate_common_variables(DataStruct&, const unsigned int&, int);

  void init_common_histos(const bool&);

  void fill_common_histos(DataStruct&, const bool&, const bool&, const unsigned int&, const double&);

  double get_xsec_from_ntuple(const std::vector<std::string>&, const bool&);

  std::pair<double, double> get_xsec_totweight_from_txt_file(const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const bool&);

  void calc_weightnorm_histo_from_ntuple(const std::vector<std::string>&, const double&, const std::vector<std::string>&,
					 const bool&, const bool&, TDirectory*, bool);

  void init_pileup_reweighting(const std::string&, const bool&, const std::vector<std::string>&);

  double get_toppt_weight(DataStruct&, const double&, const unsigned int&, const bool&);

  double get_isr_weight(DataStruct&, const double&, const unsigned int&, const bool&);

  double get_pileup_weight(DataStruct&, const double&, const unsigned int&, const bool&);

  void rescale_smear_jet_met(DataStruct&, const bool&, const unsigned int&, const double&, 
			     const double&, const double&, const bool&, const double&);

  double get_ht_weight(DataStruct&, const double&);

  double get_alphas_weight(const std::vector<float>&, const double&, const int&);

  double get_scale_weight(const std::vector<float>&, const std::vector<double>&, const double&, const unsigned int&);

  double get_syst_weight(const double&, const double&, const double&, const double&);

  double get_syst_weight(const double&, const double&, const double&);

  void job_monitoring(const int&, const int&, const std::string&, const float);

  void init_syst_input();

  double calc_top_tagging_sf(DataStruct&, const double&, const double&, const double&, const double&, const bool&);
  double calc_fake_top_0b_mass_tagging_sf(DataStruct&, const double&);
  double calc_fake_top_mass_tagging_sf(DataStruct&, const double&);
  double calc_fake_top_anti_tagging_sf(DataStruct&, const double&);

  double calc_w_tagging_sf(DataStruct&, const double&, const double&, const double&, const double&, const bool&);
  double calc_fake_w_mass_tagging_sf(DataStruct&, const double&);
  double calc_fake_w_anti_tagging_sf(DataStruct&, const double&);

  std::pair<double, double> calc_b_tagging_sf(DataStruct&, const double&, const double&, const bool&);

  std::tuple<double, double, double> calc_ele_sf(DataStruct&, const double&, const double&, const double&, const double&, const bool&);

  std::tuple<double, double, double> calc_muon_sf(DataStruct&, const double&, const double&, const double&, const bool&);
  
  double calc_lostlep_weight(DataStruct&, const double&);

  double calc_trigger_efficiency(DataStruct&, const double&);

  int calc_mrr2_bin(DataStruct&, const char&);

  std::map<char, std::vector<double> > scale_factors;
  std::map<char, double> sf_weight;
  std::vector<double> all_weights;

  double other_trigger_eff;

  const bool isData;
  const bool isSignal;
  const bool isBackground;
  const std::string sample;
  std::vector<std::string> syst;

  std::map<char, std::vector<Cut> > analysis_cuts;

  bool apply_cut(char, std::string);
  bool apply_cut(char, unsigned int);
  bool apply_ncut(char, std::string);
  bool apply_ncut(char, unsigned int);
  bool apply_cuts(char, std::vector<std::string>);
  bool apply_cuts(char, std::vector<unsigned int>);
  bool apply_all_cuts(char);
  bool apply_all_cuts_except(char, std::string);
  bool apply_all_cuts_except(char, unsigned int);
  bool apply_all_cuts_except(char, std::vector<std::string>);
  bool apply_all_cuts_except(char, std::vector<unsigned int>);

private:

  TStopwatch *sw_1_, *sw_1k_, *sw_10k_, *sw_job_;
  TRandom3 rnd_;
  std::map<std::string, int> bad_files;

  BTagCalibration* btag_calib_full_;
  BTagCalibration* btag_calib_fast_;
  BTagCalibrationReader* btag_sf_full_loose_;
  BTagCalibrationReader* btag_sf_fast_loose_;
  BTagCalibrationReader* btag_sf_full_medium_;
  BTagCalibrationReader* btag_sf_fast_medium_;

  TF1* puppisd_corrGEN_      = 0;
  TF1* puppisd_corrRECO_cen_ = 0;
  TF1* puppisd_corrRECO_for_ = 0;
};

// _____________________________________________________________
//         Analysis: Analysis specific methods/histos

class Analysis : public AnalysisBase
{
public:
  Analysis(const bool isData, const bool& isSignal, const std::string& dirname) : 
    AnalysisBase(isData, isSignal, dirname)
  {}
  ~Analysis() {}

  void calculate_variables(DataStruct&, const unsigned int&);

  bool pass_skimming(DataStruct&);

  void define_selections(const DataStruct&);

  virtual bool signal_selection(const DataStruct&);

  void apply_scale_factors(DataStruct&, const unsigned int&, const std::vector<std::vector<double> >&);

  void define_histo_options(const double&, const DataStruct&, const unsigned int&, const unsigned int&, bool);

  void init_analysis_histos(const unsigned int&, const unsigned int&);

  void fill_analysis_histos(DataStruct&, const unsigned int&, const double&);

  void load_analysis_histos(std::string);

  void save_analysis_histos(bool);

private:

  typedef struct Sample { std::string postfix; std::string legend; std::string color; std::vector<std::string> dirs; } Sample;
  typedef struct PostfixOptions { size_t index; std::string postfixes; std::string legends; std::string colors; } PostfixOptions;
  PostfixOptions get_pf_opts_(std::vector<std::vector<Sample> > lists, std::string);
};


//_______________________________________________________
//                 Define baseline cuts
void
AnalysisBase::define_preselections(const DataStruct& data)
{
  baseline_cuts.clear();

  // Apply the same cuts as it is in the ntuple - Only for check
  // cut is an std::function, which we can define easily with a lambda function

  //baseline_cuts.push_back({ .name="ntuple_filter", .func = [&data]() { 
  //      		      // Define cut function here:
  //      		      if ( !(data.jetsAK8.size>=2) ) return 0;
  //      		      if ( !(data.jetsAK8.Pt[0]>350) ) return 0;
  //      		      if ( !(data.jetsAK8.Pt[1]>350) ) return 0;
  //      		      return 1;
  //      		    } });
  
  // Recommended event filters by MET group - Updated to 80X Recommendations
  // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=101#Analysis_Recommendations_for_ana
  // 
  // Select at least one good vertex (|z|<24, |rho|<2, ndof>=4)
  // NGoodVtx defined in:
  // https://github.com/jkarancs/B2GTTrees/blob/v8.0.x_v2.1_Oct24/plugins/B2GEdmExtraVarProducer.cc#L528-L531
  // baseline_cuts.push_back({ .name="met_filter_NGoodVtx",          .func = [&data] { return data.evt.NGoodVtx>0; } });
  baseline_cuts.push_back({ .name="Clean_goodVertices",      .func = [&data] { return data.filter.goodVertices; } });
  
  // Other filters (in 80X MiniAODv2)
  // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=101#What_is_available_in_MiniAOD
  baseline_cuts.push_back({ .name="Clean_CSC_Halo_Tight",    .func = [&data,this] { return isSignal ? 1 : data.filter.globalTightHalo2016Filter; } });
  baseline_cuts.push_back({ .name="Clean_HBHE_Noise",        .func = [&data] { return data.filter.HBHENoiseFilter; } });
  baseline_cuts.push_back({ .name="Clean_HBHE_IsoNoise",     .func = [&data] { return data.filter.HBHENoiseIsoFilter; } });
  baseline_cuts.push_back({ .name="Clean_Ecal_Dead_Cell_TP", .func = [&data] { return data.filter.EcalDeadCellTriggerPrimitiveFilter; } });
  baseline_cuts.push_back({ .name="Clean_EE_Bad_Sc",         .func = [&data,this] { return isData ? data.filter.eeBadScFilter : 1; } });
  // Not in MiniAODv2 (producer added)
  baseline_cuts.push_back({ .name="Clean_Bad_Muon",          .func = [&data] { return data.filter.BadPFMuonFilter; } });
  baseline_cuts.push_back({ .name="Clean_Bad_Charged",       .func = [&data] { return data.filter.BadChargedCandidateFilter; } });
#if VER > 2
  baseline_cuts.push_back({ .name="Clean_Corridor_Veto",     .func = [&data,this] { return isSignal ? !data.filter.isPathologicalFastsimEvent : 1; } });
#endif
}


//_______________________________________________________
//                 Object selections


/*
  Jet ID (Oct31/Jan12 ntuple):
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
  
  Latest Recommendation (Exactly the same for |eta| <2.4):
  https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016?rev=4
  
  For AK4 Jet Selection Choose:
  - Loose jet ID
  - pt > 30
  - |eta| < 2.4

  For AK8 Jet Selection Choose:
  - Loose jet ID
  - pt > 200 (this cut was lowered to 170 for skimming)
  - |eta| < 2.4

*/
#define JET_AK4_PT_CUT  30
#define JET_AK4_ETA_CUT 2.4
#define JET_AK8_PT_CUT  200
#define JET_AK8_ETA_CUT 2.4

/*
  MET:
  https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes

  MET definition:
  - Use Type1 corrected Particle-flow MET
  + Data: slimmedMETsMuEGClean
  + MC:   slimmedMETsMuClean
*/

// 0: Muon + E/Gamma cleaned met (new, to be used)
// 1: Muon cleaned only (old default)
// 2: Puppi MET
#define MET_CHOICE 0

/*
  Latest b-tagging WPs/SFs:
  https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco?rev=5#Supported_Algorithms_and_Operati

  Choose:
  - CombinedSecondaryVertex v2
  - CSVv2 >= 0.5426 (Loose - for Veto)
  - CSVv2 >= 0.8484 (Medium - for Tag)

*/
#define B_SUBJET_CSV_LOOSE_CUT 0.5426
#define B_CSV_LOOSE_CUT        0.5426
#define B_CSV_MEDIUM_CUT       0.8484
#define B_CSV_TIGHT_CUT        0.9535

/* 
   Boosted W-tagging:

   Latest WPs/SFs:
   https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=54#Working_points_and_scale_factors

   W-Tag selection:
   - AK8 CHS jets
   - pt > 200
   - |eta| < 2.4
   - 65 <= Puppi SD Mass (corrected) < 105
   - Medium WP: Puppi tau21 < 0.4

*/

#define W_PT_CUT            200
#define W_ETA_CUT           2.4
#define W_SD_MASS_CUT_LOW   65
#define W_SD_MASS_CUT_HIGH  105
#define W_TAU21_LOOSE_CUT   0.55
#define W_TAU21_TIGHT_CUT   0.4 // There is a Tighter WP: 0.35

#define W_TAG_HP_SF       1.00
#define W_TAG_HP_SF_ERR   0.06
#define W_TAG_LP_SF       0.96
#define W_TAG_LP_SF_ERR   0.11
#define W_TAG_JMS_SF      1.00
#define W_TAG_JMS_SF_ERR  0.0094
#define W_TAG_JMR_SF      1.00
#define W_TAG_JMR_SF_ERR  0.20
#define W_TAG_SIGMA_MC    10.1

/*

Top Tagging working points (No subjet B tag):
  
Latest WPs/SFs:
https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging?rev=14#13_TeV_working_points_CMSSW_8_0
  
Choose:
- Tightest selection
- AK8 CHS jets
- 105 < Puppi SD Mass < 210
- Puppi tau32 < 0.46


Top Tagging working point (with subjet b tag):
  
Latest WPs/SFs:
https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging?rev=14#13_TeV_working_points_CMSSW_8_0
https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco?rev=14#Boosted_event_topologies
  
Choose:
- Use Tightest selection
- AK8 CHS jets
- 105 <= Puppi SD Mass < 210
- Puppi tau32 < 0.46
- max subjet BTag CSV > 0.5426
*/

#define USE_BTAG 1

#if USE_BTAG == 0

#define TOP_PT_CUT            400
#define TOP_SD_MASS_CUT_LOW   105
#define TOP_SD_MASS_CUT_HIGH  210
#define TOP_TAU32_CUT        0.46

#define TOP_TAG_SF           1.07
#define TOP_TAG_SF_ERR_UP    0.08
#define TOP_TAG_SF_ERR_DOWN  0.04

#else

#define TOP_PT_CUT            400
#define TOP_SD_MASS_CUT_LOW   105
#define TOP_SD_MASS_CUT_HIGH  210
#define TOP_TAU32_CUT        0.46
#define TOP_BTAG_CSV         0.5426

#define TOP_TAG_SF           1.05
#define TOP_TAG_SF_ERR_UP    0.07
#define TOP_TAG_SF_ERR_DOWN  0.04

#endif

/*
  Latest Electron IDs:
  [1] Cut Based  - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Working_points_for_2016_data_for
  [2] MVA        - https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2?rev=30#Recommended_MVA_recipes_for_2016
  [3] SUSY (Use) - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#Electrons

  Latest Isolation WPs:
  [4] SUSY MiniIso Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO_AN1
  [5] RazorInclusive - https://github.com/RazorCMS/RazorAnalyzer/blob/master/src/RazorAnalyzer.cc#L1210

  Latest Impact Point Cut:
  [6] SUSY Loose/Tight IP2D (Use) - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO_AN1
  [7] POG  Tight - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Offline_selection_criteria

  For Veto (Regions except Z) Choose:
  - Spring15 Cut based Veto ID without relIso (EA) cut
  - for pt<20  absolute iso03 (EA) < 5 
  - Mini-Isolation (EA)/pt < 0.1 (Medium WP [4])
  - pt >= 5
  - |eta| < 2.5
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [6])
  OR when using MVA:
  - Spring16 MVA Loose ID
  - for pt>=20 Mini-Isolation (EA)/pt < 0.2 (RazorInclusive cut [5])
  - pt >= 5
  - |eta| < 2.5
  - 3D IP sig < 4

  For Selection (Z) Choose:
  - Spring15 Cut based Medium ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [6])
  OR when using MVA:
  - Spring16 MVA Loose ID
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [6])


  For Tight Selection (TriggerEff Only) Choose:
  - Spring15 Cut based Tight ID (including relIso (EA) cut)
  - pt >= 30
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D and IP3D [6])

  For Loose selection (deprecated) Choose:
  - Spring15 Cut based Loose ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [6])

*/

#define USE_MVA_ID             1
#define USE_POG_ID             0

#if USE_POG_ID > 0
#define ELE_VETO_PT_CUT        5
#define ELE_VETO_ETA_CUT       2.5
#define ELE_VETO_MINIISO_CUT   0.1
#define ELE_VETO_ABSISO_CUT    5   // For skim only
#if USE_MVA_ID == 1
#define ELE_VETO_IP_D0_CUT     0.05
#define ELE_VETO_IP_DZ_CUT     0.1
#define ELE_VETO_IP_3D_CUT     4   // For skim only
#else
#define ELE_VETO_IP_D0_CUT     0.2
#define ELE_VETO_IP_DZ_CUT     0.5
#endif
#else
#define ELE_VETO_PT_CUT        5
#define ELE_VETO_ETA_CUT       2.5
#define ELE_VETO_MINIISO_CUT   0.2
#define ELE_VETO_ABSISO_CUT    5
#if USE_MVA_ID == 1
#define ELE_VETO_IP_D0_CUT     0.05
#define ELE_VETO_IP_DZ_CUT     0.1
#define ELE_VETO_IP_3D_CUT     4
#else
#define ELE_VETO_IP_D0_CUT     0.2
#define ELE_VETO_IP_DZ_CUT     0.5
#endif
#endif

#define ELE_LOOSE_PT_CUT       10
#define ELE_LOOSE_ETA_CUT      2.5
#define ELE_LOOSE_MINIISO_CUT  0.1
#define ELE_LOOSE_IP_D0_CUT    0.2
#define ELE_LOOSE_IP_DZ_CUT    0.5

#define ELE_SELECT_PT_CUT      10
#define ELE_SELECT_ETA_CUT     2.5
#define ELE_SELECT_MINIISO_CUT 0.1
#define ELE_SELECT_IP_D0_CUT   0.05
#define ELE_SELECT_IP_DZ_CUT   0.1

#define ELE_TIGHT_PT_CUT       30
#define ELE_TIGHT_ETA_CUT      2.5
#define ELE_TIGHT_IP_D0_CUT    0.05
#define ELE_TIGHT_IP_DZ_CUT    0.1
#define ELE_TIGHT_IP_SIG_CUT   4

/*
  Latest Muon IDs (Loose/Medium):
  [1] POG Loose/Medium - https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=28#Short_Term_Instructions_for_Mori

  Latest Isolation WPs:
  [2] SUSY MiniISo Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO
  [3] POG Tight RelIso - https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=28#Muon_Isolation
  [4] RazorInclusive - https://github.com/RazorCMS/RazorAnalyzer/blob/master/src/RazorAnalyzer.cc#L1602

  Latest Impact Point Cut (Loose/Tight):
  [5] SUSY Loose/Tight IP2D - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO
  
  For Veto (All except Z) Choose:
  - POG recommended Loose ID (No Iso/IP)
  - for pt>=20 Mini-Isolation (EA)/pt < 0.2 (RazorInclusive cut [4])
  - for pt<20  absolute iso04 (EA) < 10
  - pt >= 5
  - |eta| < 2.4
  - 3D IP sig < 4

  For Loose Choose:
  - POG recommended Loose ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.2 (Tight WP [2])
  - pt >= 10
  - |eta| < 2.4
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [5])

  For Selection (Z) Choose:
  - POG recommended Medium ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.2 (tight WP [2])
  - pt >= 5
  - |eta| < 2.4
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [5])

  For Tight Selection (TriggerEff Only) Choose:
  - POG recommended Tight ID (No Iso/IP)
  - comb. rel. Isolation (R=0.4) < 0.15 (tight WP [3])
  - pt >= 30
  - |eta| < 2.4
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D and IP3D [5])

*/

#if USE_POG_ID > 0
#define MU_VETO_PT_CUT         5
#define MU_VETO_ETA_CUT        2.4
#define MU_VETO_MINIISO_CUT    0.4
#define MU_VETO_ABSISO_CUT     10  // For skim only
#define MU_VETO_IP_D0_CUT      0.2
#define MU_VETO_IP_DZ_CUT      0.5
#define MU_VETO_IP_3D_CUT      4   // For skim only
#else
#define MU_VETO_PT_CUT         5
#define MU_VETO_ETA_CUT        2.4
#define MU_VETO_MINIISO_CUT    0.2
#define MU_VETO_ABSISO_CUT     10
#define MU_VETO_IP_D0_CUT      0.2
#define MU_VETO_IP_DZ_CUT      0.5
#define MU_VETO_IP_3D_CUT      4
#endif

#define MU_LOOSE_PT_CUT        10
#define MU_LOOSE_ETA_CUT       2.4
#define MU_LOOSE_MINIISO_CUT   0.2
#define MU_LOOSE_IP_D0_CUT     0.2
#define MU_LOOSE_IP_DZ_CUT     0.5

#define MU_SELECT_PT_CUT       10
#define MU_SELECT_ETA_CUT      2.4
#define MU_SELECT_MINIISO_CUT  0.2
#define MU_SELECT_IP_D0_CUT    0.05
#define MU_SELECT_IP_DZ_CUT    0.1

#define MU_TIGHT_PT_CUT        30
#define MU_TIGHT_ETA_CUT       2.4
#define MU_TIGHT_RELISO_CUT    0.15
#define MU_TIGHT_IP_D0_CUT     0.05
#define MU_TIGHT_IP_DZ_CUT     0.1
#define MU_TIGHT_IP_SIG_CUT    4

/*
  Latest Tau IDs:
  https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV?rev=35

  For Veto Choose:
  - POG Loose ID:
  + byLooseCombinedIsolationDeltaBetaCorr3Hits"
  - pt >= 18 (Same as MINIAOD)

  OR use:
  - Isolated charged trk = 0
  + pt>=5 (ele/mu), 10 (hadrons)
  + isolation (dR=0.3) / pt < 0.2 (ele/mu), 0.1 (hadrons)
  + |dz| < 0.1

*/

#define USE_ISO_TRK_VETO  0
#define COMBINE_MRR2_BINS 1

/*
  Latest Photon IDs:
  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=36

  For Selection Choose:
  - Spring16 Cut based Medium photon ID
  - Pass electron veto
  - pt >= 80
  - |eta| < 2.5

*/


#define PHOTON_SELECT_PT_CUT        80
#define PHOTON_SELECT_ETA_CUT       2.5


//_______________________________________________________
//              Rescale jet 4-momenta

// Variables to save the original values before applying any systematics on them
std::vector<float> AK4_E, AK4_Pt;
std::vector<float> AK8_E, AK8_Pt, AK8_softDropMass;//, AK8_trimmedMass, AK8_prunedMass, AK8_filteredMass;
std::vector<float> AK8_softDropMass0;
std::vector<float> AK8_softDropMass1;
std::vector<float> AK8_softDropMass2;
std::vector<float> AK8_softDropMassCorr; // Correction for W tagging
std::vector<float> softDropMassCorr;     // POG Correction + uncertainties for W tagging
std::vector<float> softDropMass;         // JES Correction + uncertainties for top tagging

std::vector<float> AK4_JERSmearFactor,     AK8_JERSmearFactor;
std::vector<float> AK4_JERSmearFactorUp,   AK8_JERSmearFactorUp;
std::vector<float> AK4_JERSmearFactorDown, AK8_JERSmearFactorDown;
std::vector<float> AK8_JMR_random;

TVector3 met;
TVector3 dmet_JESUp,  dmet_JESDown;
TVector3 dmet_JERUp,  dmet_JERDown;
TVector3 dmet_RestUp, dmet_RestDown;

bool recalc_megajets;

void
AnalysisBase::rescale_smear_jet_met(DataStruct& data, const bool& applySmearing, const unsigned int& syst_index,
				    const double& nSigmaJES, const double& nSigmaJER, const double& nSigmaRestMET,
				    const bool& rescaleAK8, const double& nSigmaRescaleAK8)
{
  // Apply Jet Energy Scale (JES) and Jet Energy Resolution (JER) corrections
  // For AK8 jets which are used for W tagging (only):
  // - Additionally or instead apply jet mass scale (JMS) and jet mass resolutin (JMR) corrections
  recalc_megajets = (nSigmaJES!=0) || (nSigmaJER!=0);

  // Initialization (needed for later variations
  if (syst_index==0) {
    // Save the original values for later (before applying any systematics)
    AK4_E            = data.jetsAK4.E;
    AK4_Pt           = data.jetsAK4.Pt;
    AK8_E            = data.jetsAK8.E;
    AK8_Pt           = data.jetsAK8.Pt;
#if VER == 0
    AK8_softDropMass = data.jetsAK8.softDropMass;
#elif VER == 1
    AK8_softDropMass = data.jetsAK8.softDropMassPuppi;
#else
    //AK8_softDropMass = data.jetsAK8.softDropMassPuppi; 
    AK8_softDropMass = data.jetsAK8.corrSDMassPuppi; // Make sure to change also sd_mass_top if needed!!!
    AK8_softDropMass0 = data.jetsAK8.uncorrSDMassPuppi; // No correction (needed for W)
    AK8_softDropMass1 = data.jetsAK8.corrSDMassPuppi;   // L1L2L3 on subjet
    AK8_softDropMass2 = data.jetsAK8.softDropMassPuppi; // L2L3 on fatjet
#endif
    //AK8_trimmedMass  = data.jetsAK8.trimmedMass;
    //AK8_prunedMass   = data.jetsAK8.prunedMass;
    //AK8_filteredMass = data.jetsAK8.filteredMass;

    // Correction for Puppi SoftDrop Mass
    // (Needed for W tagging)
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=43#Working_points_and_scale_factors
    if (puppisd_corrGEN_==0) {
#if VER < 2
      AK8_softDropMassCorr = AK8_softDropMass;
#else
      AK8_softDropMassCorr = data.jetsAK8.uncorrSDMassPuppi;
#endif
    } else {
      AK8_softDropMassCorr.clear();
      for (size_t i=0; i<data.jetsAK8.size; ++i) {
#if VER == 0
	double puppi_pt  = data.jetsAK8.Pt[i];
	double puppi_eta = data.jetsAK8.Eta[i];
	double puppi_sd_mass_w = data.jetsAK8.softDropMass[i];
#else
	double puppi_pt  = data.jetsAK8.PtPuppi[i];
	double puppi_eta = data.jetsAK8.EtaPuppi[i];
#if VER == 1
	double puppi_sd_mass_w = data.jetsAK8.softDropMassPuppi[i];
#else
	double puppi_sd_mass_w = data.jetsAK8.uncorrSDMassPuppi[i];
#endif
#endif
	double corr = puppisd_corrGEN_->Eval(puppi_pt);
	if(std::abs(puppi_eta)<=1.3) corr *= puppisd_corrRECO_cen_->Eval(puppi_pt);
	else corr *= puppisd_corrRECO_for_->Eval(puppi_pt);
	
	AK8_softDropMassCorr.push_back(puppi_sd_mass_w * corr);
      }
    } // end sys_index==0

    // Calculate the JER/JMR smear factors
    if (applySmearing) {
      AK4_JERSmearFactor    .clear();
      AK4_JERSmearFactorUp  .clear();
      AK4_JERSmearFactorDown.clear();
      for (size_t i=0; i<data.jetsAK4.size; ++i) {
        double JERSmear     = data.jetsAK4.SmearedPt[i]/data.jetsAK4.Pt[i];
        double JERSmearUp   = 1 + (JERSmear-1) * (data.jetsAK4.JERSFUp[i]  -1) / (data.jetsAK4.JERSF[i]-1);
        double JERSmearDown = 1 + (JERSmear-1) * (data.jetsAK4.JERSFDown[i]-1) / (data.jetsAK4.JERSF[i]-1);
        AK4_JERSmearFactor    .push_back(JERSmear);
        AK4_JERSmearFactorUp  .push_back(JERSmearUp);
        AK4_JERSmearFactorDown.push_back(JERSmearDown);
      }
      AK8_JERSmearFactor    .clear();
      AK8_JERSmearFactorUp  .clear();
      AK8_JERSmearFactorDown.clear();
      AK8_JMR_random.clear();
      for (size_t i=0; i<data.jetsAK8.size; ++i) {
        double JERSmear     = data.jetsAK8.SmearedPt[i]/data.jetsAK8.Pt[i];
        double JERSmearUp   = 1 + (JERSmear-1) * (data.jetsAK8.JERSFUp[i]  -1) / (data.jetsAK8.JERSF[i]-1);
        double JERSmearDown = 1 + (JERSmear-1) * (data.jetsAK8.JERSFDown[i]-1) / (data.jetsAK8.JERSF[i]-1);
        AK8_JERSmearFactor    .push_back(JERSmear);
        AK8_JERSmearFactorUp  .push_back(JERSmearUp);
        AK8_JERSmearFactorDown.push_back(JERSmearDown);
	// Apply random gaussian smearing to worsen mass resolution (It cannot improve with this method)
	// Recipe is the same as the stochastic smearing explained here:
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
	// Generate the random number once, and vary the systematics only on the SF
	double sigma_rel = W_TAG_SIGMA_MC / (AK8_softDropMassCorr[i] * W_TAG_JMS_SF);
	double random = rnd_.Gaus(0,sigma_rel);
	AK8_JMR_random.push_back(random);
      }
    }

    //_________________________________________________
    //                      MET

    // Choose the met to use
#if MET_CHOICE == 1
    // Muon cleaned only PF MET
    data.met.Pt       = data.met.MuCleanOnly_Pt;
    data.met.Phi      = data.met.MuCleanOnly_Phi;
    data.syst_met.Pt  = data.syst_met.MuCleanOnly_Pt;
    data.syst_met.Phi = data.syst_met.MuCleanOnly_Phi;
#elif MET_CHOICE == 2
    // Puppi MET
    data.met.Pt       = data.puppimet.Pt;
    data.met.Phi      = data.puppimet.Phi;
    data.syst_met.Pt  = data.syst_puppimet.Pt;
    data.syst_met.Phi = data.syst_puppimet.Phi;
#endif

    // Save the original MET
    met.SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);

#if VER > 0
    // MET Uncertainties
    /*
      Met uncertainty vector indices:

      enum METUncertainty {
      JetResUp =0, JetResDown =1,
      JetEnUp =2, JetEnDown =3,
      MuonEnUp =4, MuonEnDown =5,
      ElectronEnUp =6, ElectronEnDown =7,
      TauEnUp =8, TauEnDown =9,
      UnclusteredEnUp =10, UnclusteredEnDown =11,
      PhotonEnUp =12, PhotonEnDown =13,
      }
    */
    float maxdpt_up = 0, maxdpt_down = 0;
    float dphi_up   = 0, dphi_down   = 0;
    float ptsum_up  = 0, ptsum_down  = 0;
    // Consider JES/JER modulation separately
    // Add the rest of the systematic pt modulations in quadrature
    // Use the phi direction of the largest remaining systematic
    for (size_t i=0; i<data.syst_met.size; ++i) {
      TVector3 met_syst;
      met_syst.SetPtEtaPhi(data.syst_met.Pt[i], 0, data.syst_met.Phi[i]);
      TVector3 dmet;
      dmet = met_syst - met;
      if (i==0) {
	dmet_JERUp = dmet;
      } else if (i==1) {
	dmet_JERDown = dmet;
      } else if (i==2) {
	dmet_JESUp = dmet;
      } else if (i==3) {
	dmet_JESDown = dmet;
      } else if (i%2==0) {
	// Rest Up
	if (dmet.Pt()>maxdpt_up) {
	  maxdpt_up = dmet.Pt();
	  dphi_up   = dmet.Phi();
	  ptsum_up  = std::sqrt(ptsum_up*ptsum_up + dmet.Perp2());
	}
      } else {
	// Rest Down
	if (dmet.Pt()>maxdpt_down) {
	  maxdpt_down = dmet.Pt();
	  dphi_down   = dmet.Phi();
	  ptsum_down  = std::sqrt(ptsum_down*ptsum_down + dmet.Perp2());
	}
      }
    }
    dmet_RestUp.  SetPtEtaPhi(ptsum_up,   0, dphi_up);
    dmet_RestDown.SetPtEtaPhi(ptsum_down, 0, dphi_down);
#endif
  }


  // Apply systematic variations
  // Even if Sigmas=0, we still smear jets!
  // AK4 jets
  while(data.jetsAK4.Loop()) {
    size_t i = data.jetsAK4.it;
    double scaleJES = get_syst_weight(1.0, data.jetsAK4.jecUncertainty[i], nSigmaJES);
    data.jetsAK4.Pt[i] = AK4_Pt[i] * scaleJES;
    data.jetsAK4.E[i]  = AK4_E[i]  * scaleJES;
    if (applySmearing) {
      double scaleJER = get_syst_weight(AK4_JERSmearFactor[i], AK4_JERSmearFactorUp[i], AK4_JERSmearFactorDown[i], nSigmaJER);
      data.jetsAK4.Pt[i] *= scaleJER;
      data.jetsAK4.E[i]  *= scaleJER;
    }
  }
  // AK8 jets
  softDropMassCorr.clear();
  softDropMass    .clear();
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    double scaleJES = get_syst_weight(1.0, data.jetsAK8.jecUncertainty[i], nSigmaJES);
    data.jetsAK8.Pt[i] = AK8_Pt[i] * scaleJES;
    data.jetsAK8.E[i]  = AK8_E[i]  * scaleJES;
    double scaleJER = 1;
    if (applySmearing) {
      scaleJER = get_syst_weight(AK8_JERSmearFactor[i], AK8_JERSmearFactorUp[i], AK8_JERSmearFactorDown[i], nSigmaJER);
      data.jetsAK8.Pt[i] *= scaleJER;
      data.jetsAK8.E[i]  *= scaleJER;
    }

    double rescale_v4 = 1.0;
    if (rescaleAK8) {
      // Apply AK8 jet rescaling to match data distribution for Z/G/DY MCs
      // Get the scale/offset from: scripts/calc_pt_scaling.C
      // Rescale the whole 4-momentum, because the Z-mass peak (in AK8 jet) also looks shifted
      //double scale = 0.867 + 0.039*nSigmaRescaleAK8, offset = 33.2 + 36.4 *nSigmaRescaleAK8;
#if USE_ISO_TRK_VETO > 0
      double scale = 0.858 + 0.047*nSigmaRescaleAK8, offset = 30.8 + 38.6 *nSigmaRescaleAK8;
#else
      //double scale = 0.840 + 0.037*nSigmaRescaleAK8, offset = 50.5 + 35.0 *nSigmaRescaleAK8;
      double scale = 0.848 + 0.047*nSigmaRescaleAK8, offset = 39.6 + 38.7 *nSigmaRescaleAK8;
#endif
      // Use original pt without any scale/smear to determine the right scaling
      double orig_pt = AK8_Pt[i];
      if (applySmearing) orig_pt *= AK8_JERSmearFactor[i];
      rescale_v4 = scale + offset/orig_pt;
      data.jetsAK8.Pt[i] *= rescale_v4;
      data.jetsAK8.E[i]  *= rescale_v4;
    }
    //AK8_Ht += data.jetsAK8.Pt[i];

    // For W jet mass apply combination of both JES+JMS and JMR (JER not needed)
    // JES uncertainty is added on top of the JMS one (in quadrature)
    double comb_unc = std::sqrt(W_TAG_JMS_SF_ERR*W_TAG_JMS_SF_ERR + data.jetsAK8.jecUncertainty[i]*data.jetsAK8.jecUncertainty[i]);
    double scaleJMS = get_syst_weight(W_TAG_JMS_SF, comb_unc, nSigmaJES);
    //double scaled_corrected_mass = AK8_softDropMassCorr[i] * scaleJMS * rescale_v4;
    double scaled_corrected_mass = AK8_softDropMassCorr[i] * scaleJMS;
    // Apply random gaussian smearing to worsen mass resolution (It cannot improve with this method)
    if (applySmearing) {
      double scale_factor = get_syst_weight(W_TAG_JMR_SF, W_TAG_JMR_SF_ERR, nSigmaJER);
      if (scale_factor > 1) {
	double scaleJMR = 1 + AK8_JMR_random[i] * std::sqrt( scale_factor*scale_factor - 1 );
	scaled_corrected_mass *= scaleJMR;
      }
    }
    softDropMassCorr.push_back(scaled_corrected_mass);

    // For Top jet mass, apply only JES + JER for now
    // (Since there's no other recommendation)
#if VER == 0
    //double scaled_top_mass = AK8_softDropMass[i] * scaleJES * rescale_v4;
    double scaled_top_mass = AK8_softDropMass[i] * scaleJES;
    if (applySmearing) scaled_top_mass *= scaleJER;
#else
    //double scaled_top_mass = AK8_softDropMass[i] * scaleJES * rescale_v4;
    double scaled_top_mass = AK8_softDropMass[i] * scaleJES;
    if (applySmearing) scaled_top_mass *= scaleJER;
#endif
    softDropMass.push_back(scaled_top_mass);

  }

  TVector3 dmet(0,0,0);
#if VER > 0
  // MET Uncertainties
  if (!isSignal) {
    // Background
    if (nSigmaJES   >=0) dmet += std::abs(nSigmaJES) * dmet_JESUp;
    else                 dmet += std::abs(nSigmaJES) * dmet_JESDown;
    if (applySmearing) {
      if (nSigmaJER   >=0) dmet += std::abs(nSigmaJER) * dmet_JERUp;
      else                 dmet += std::abs(nSigmaJER) * dmet_JERDown;
    }
    if (nSigmaRestMET>=0) dmet += std::abs(nSigmaRestMET) * dmet_RestUp;
    else                  dmet += std::abs(nSigmaRestMET) * dmet_RestDown;
  }
#if VER > 2
  else {
    // Signal - Use only a one sided variation
    // Take the gen MET as the Up variation
    // and use the PF MET as the nominal/down variation -> i.e do nothing for Sigma < 0
    // --> In the end use the average of the acceptances as the nominal value
    // But do not use the average PF/gen met as the nominal value
    if (nSigmaRestMET>0) {
      TVector3 genmet;
      genmet.SetPtEtaPhi(data.met.genPt, 0, data.met.genPhi);
      dmet = genmet - met;
      dmet = std::abs(nSigmaRestMET) * dmet;
    }
  }
#endif
#endif
  TVector3 shifted_met = met + dmet;
  data.met.Pt[0]  = shifted_met.Pt();
  data.met.Phi[0] = shifted_met.Phi();
}


//_______________________________________________________
//                 Define common variables


// AK4 jets
/*
  convention:

  iObject  -  gives the index of the nth selected object in the original collection

  example:
  for (size_t i=0; i<nJet; ++i) h_pt->Fill( data.jetsAK4.Pt[iJet[i]] );  
  or
  if (nJet>0) vh_pt[0] -> Fill( data.jetsAK4.Pt[iJet[0]] );
  if (nJet>1) vh_pt[1] -> Fill( data.jetsAK4.Pt[iJet[1]] );


  itObject  -  gives the index in the selected collection

  example:
  for (size_t it=0; it<data.jetsAK4.size; ++it)
  if (passLooseJet[it]) vh_pt[itJet[it]]->Fill( data.jetsAK4.Pt[it] );

*/

std::vector<size_t > iJet;
std::vector<size_t > iLooseBTag;
std::vector<size_t > iMediumBTag;
std::vector<size_t > iTightBTag;
std::vector<size_t > itJet;
std::vector<size_t > itLooseBTag;
std::vector<size_t > itMediumBTag;
std::vector<size_t > itTightBTag;
std::vector<bool> passLooseJet;
std::vector<bool> passLooseJetNoLep;
std::vector<bool> passLooseJetNoPho;
std::vector<bool> passLooseBTag;
std::vector<bool> passMediumBTag;
std::vector<bool> passTightBTag;
std::vector<double> AK4_photonDR;
std::vector<double> AK4_photonPtRatio;
std::vector<double> AK4_eleDR;
std::vector<double> AK4_elePtRatio;
std::vector<double> AK4_muonDR;
std::vector<double> AK4_muonPtRatio;
unsigned int nJet;
unsigned int nJetNoLep;
unsigned int nJetNoPho;
unsigned int nLooseBTag;
unsigned int nMediumBTag;
unsigned int nMediumBTagNoPho;
unsigned int nTightBTag;
double AK4_Ht, AK4_HtOnline, AK4_HtNoLep;
double minDeltaPhi; // Min(DeltaPhi(Jet_i, MET)), i=1,2,3,4
double dPhiRazor, dPhiRazorNoPho;

// AK8 jets
std::vector<size_t > iJetAK8;
std::vector<size_t > iWMassTag;
std::vector<size_t > iLooseWTag;
std::vector<size_t > iTightWTag;
std::vector<size_t > iTightWAntiTag;
std::vector<size_t > iHadTopMassTag;
std::vector<size_t > iHadTop1BMassTag;
std::vector<size_t > iHadTopTag;
std::vector<size_t > iHadTop0BMassTag;
std::vector<size_t > iHadTop0BAntiTag;
std::vector<size_t > itJetAK8;
std::vector<size_t > itWMassTag;
std::vector<size_t > itLooseWTag;
std::vector<size_t > itTightWTag;
std::vector<size_t > itTightWAntiTag;
std::vector<size_t > itHadTopMassTag;
std::vector<size_t > itHadTop1BMassTag;
std::vector<size_t > itHadTopTag;
std::vector<size_t > itHadTop0BMassTag;
std::vector<size_t > itHadTop0BAntiTag;
std::vector<double> tau21;
std::vector<double> tau31;
std::vector<double> tau32;
std::vector<double> AK8_photonDR;
std::vector<double> AK8_photonPtRatio;
std::vector<double> AK8_eleDR;
std::vector<double> AK8_elePtRatio;
std::vector<double> AK8_muonDR;
std::vector<double> AK8_muonPtRatio;
std::vector<float> softDropMassW;
std::vector<float> softDropMassTop;
std::vector<float> softDropMassTop1;
std::vector<float> softDropMassTop2;
#if VER == 0
std::vector<double> maxSubjetCSV;
#endif
std::vector<bool> passSubjetBTag;
std::vector<bool> passLooseJetAK8;
std::vector<bool> passWMassTag;
std::vector<bool> passLooseWTag;
std::vector<bool> passTightWTag;
std::vector<bool> passTightWAntiTag;
std::vector<bool> passHadTopTag;
std::vector<bool> passHadTopMassTag;
std::vector<bool> passHadTop1BMassTag;
std::vector<bool> passHadTop0BMassTag;
std::vector<bool> passHadTop0BAntiTag;
std::vector<bool> passOtherWMassTag;
std::vector<bool> passOtherTightWTag;
std::vector<bool> passOtherTightWAntiTag;
std::vector<bool> passOtherHadTopTag;
std::vector<bool> passOtherHadTopMassTag;
std::vector<bool> passOtherHadTop0BMassTag;
std::vector<bool> passOtherHadTop0BAntiTag;
std::vector<bool> hasGenW;
std::vector<bool> hasGenTop;
unsigned int nJetAK8;
unsigned int nWMassTag;
unsigned int nLooseWTag;
unsigned int nTightWTag;
unsigned int nHadWTag1;
unsigned int nHadWTag3;
unsigned int nTightWAntiTag;
unsigned int nSubjetBTag;
unsigned int nHadTopTag;
unsigned int nHadTopTag2;
unsigned int nHadTopTag3;
unsigned int nHadTopTag4;
unsigned int nHadTopMassTag;
unsigned int nHadTop1BMassTag;
unsigned int nHadTop0BMassTag;
unsigned int nHadTop0BAntiTag;
double AK8_Ht;
double minDeltaR_W_b;

// Event Letpons
std::vector<size_t > iEleVeto;
std::vector<size_t > iMuVeto;
std::vector<size_t > iTauVeto;
std::vector<size_t > iEleLoose;
std::vector<size_t > iMuLoose;
std::vector<size_t > iEleSelect;
std::vector<size_t > iMuSelect;
std::vector<size_t > itEleVeto;
std::vector<size_t > itMuVeto;
std::vector<size_t > itTauVeto;
std::vector<size_t > itEleLoose;
std::vector<size_t > itMuLoose;
std::vector<size_t > itEleSelect;
std::vector<size_t > itMuSelect;
std::vector<size_t > iEleTight;
std::vector<size_t > iMuTight;
std::vector<size_t > itEleTight;
std::vector<size_t > itMuTight;
std::vector<bool> passEleVeto;
std::vector<bool> passMuVeto;
std::vector<bool> passTauVeto;
std::vector<bool> passEleLoose;
std::vector<bool> passMuLoose;
std::vector<bool> passEleSelect;
std::vector<bool> passMuSelect;
std::vector<bool> passEleTight;
std::vector<bool> passMuTight;
std::vector<double> eleJetDR;
std::vector<double> eleJetDPhi;
std::vector<double> eleJetPt;
std::vector<double> muJetDR;
std::vector<double> muJetDPhi;
std::vector<double> muJetPt;
unsigned int nEleVetoNoIso;
unsigned int nEleVeto;
unsigned int nEleLoose;
unsigned int nEleSelect;
unsigned int nEleTight;
unsigned int nMuVetoNoIso;
unsigned int nMuVeto;
unsigned int nMuLoose;
unsigned int nMuSelect;
unsigned int nMuTight;
unsigned int nTauVeto;
unsigned int nLepVetoNoIso;
unsigned int nLepVeto;
unsigned int nLepLoose;
unsigned int nLepSelect;
unsigned int nLepTight;
std::vector<size_t > iPhotonPreSelect;
std::vector<size_t > iPhotonSelect;
std::vector<size_t > iPhotonFake;
std::vector<size_t > itPhotonPreSelect;
std::vector<size_t > itPhotonSelect;
std::vector<size_t > itPhotonFake;
std::vector<bool> passPhotonPreSelect;
std::vector<bool> passPhotonSelect;
std::vector<bool> passPhotonPrompt;
std::vector<bool> passPhotonFake;
std::vector<double> ChargedHadronIsoEACorr;
unsigned int nPhotonPreSelect, nPhotonSelect, nPhotonFake;
double MT, MT_vetolep, MT_allvetolep;
double vetoLepPt, vetoLepEta, vetoLepPhi;
double MET_1l, MTR_1l, R_1l, R2_1l, minDeltaPhi_1l, M_1l;
double MET_1vl, MTR_1vl, R_1vl, R2_1vl, minDeltaPhi_1vl, M_1vl;
double MET_ll, MTR_ll, R_ll, R2_ll, minDeltaPhi_ll, M_ll;
double MET_pho, MR_pho, MTR_pho, R_pho, R2_pho, minDeltaPhi_pho; 
double dPhi_ll_met, dPhi_ll_jet;
std::vector<TLorentzVector> hemis_AK4;
std::vector<TLorentzVector> hemis_AK4_nophoton;
std::vector<TLorentzVector> saved_hemis_AK4;
std::vector<TLorentzVector> saved_hemis_AK4_nophoton;

// gen particles
std::vector<bool> passGenHadW;
std::vector<bool> passGenTop;
std::vector<bool> passGenLepton;
std::vector<bool> genLepPassLepID;
std::vector<bool> genHadWPassWTag;
std::vector<bool> genTopPassTopTag;
std::vector<bool> genHadWPassWMassTag;
std::vector<bool> genTopPassTopMassTag;
std::vector<bool> genHadWHasJetMatch;
std::vector<bool> genTopHasJetMatch;
std::vector<float> genHadWMatchedAK8JetPt;
std::vector<float> genHadWMatchedAK8JetEta;
std::vector<float> genTopMatchedAK8JetPt;
std::vector<float> genTopMatchedAK8JetEta;
std::vector<size_t > iGenHadW,   itGenHadW;
std::vector<size_t > iGenTop,    itGenTop;
std::vector<size_t > iGenZ,		   itGenZ;
std::vector<size_t > iGenWtagZ,   itGenWtagZ;
std::vector<size_t > iGenMassW,   itGenMassW;
std::vector<size_t > iGenMassTop, itGenMassTop;
unsigned int nGenHadW;
unsigned int nGenTop;
unsigned int nGenZ;
unsigned int nGenWtagZ;
unsigned int nGenMassW;
unsigned int nGenMassTop;
int nWTag;
int nmWTag;
int npreWTag;
int nTopTag;
int nmTopTag;
int npreTopTag;
int nFakePhoton;
int nDirectPromptPhoton;
int nFragmentationPromptPhoton;

std::vector<TLorentzVector> veto_electrons;
std::vector<TLorentzVector> veto_muons;
std::vector<TLorentzVector> veto_leptons_noiso, veto_leptons, selected_leptons, tight_leptons;
std::vector<float> r_iso_tight_leptons;
TLorentzVector lep_pair;
//std::vector<TLorentzVector> veto_muons_noiso, veto_muons, selected_muons;
//std::vector<bool> veto_lep_in_jet;
//std::vector<bool> veto_mu_in_jet, selected_mu_in_jet;

void
AnalysisBase::calculate_common_variables(DataStruct& data, const unsigned int& syst_index, int debug = 0)
{
  if (debug) std::cout<<"AnalysisBase::calc_common_var: start (new event)"<<std::endl;

  // It only makes sense to calculate certain variables only once 
  // if they don't depend on jet energy (eg. leptons)
  if (debug) std::cout<<"AnalysisBase::calc_common_var: init start"<<std::endl;
  if (syst_index == 0) {
    veto_electrons     .clear();
    veto_muons         .clear();
    veto_leptons_noiso .clear();
    veto_leptons       .clear();
    selected_leptons   .clear();
    tight_leptons      .clear();
    r_iso_tight_leptons.clear();
    lep_pair.SetPxPyPzE(0,0,0,0);

    // Loop on AK8 jets
    tau21         .assign(data.jetsAK8.size, 9999);
    tau31         .assign(data.jetsAK8.size, 9999);
    tau32         .assign(data.jetsAK8.size, 9999);
#if VER == 0
    maxSubjetCSV .assign(data.jetsAK8.size, 0);
#endif
    passSubjetBTag.assign(data.jetsAK8.size, 0);
    nSubjetBTag = 0;
    while(data.jetsAK8.Loop()) {
      size_t i = data.jetsAK8.it;
      // N-subjettiness
#if VER == 0
      if (data.jetsAK8.tau1[i]>0) tau21[i] = data.jetsAK8.tau2[i]/data.jetsAK8.tau1[i];
      if (data.jetsAK8.tau1[i]>0) tau31[i] = data.jetsAK8.tau3[i]/data.jetsAK8.tau1[i];
      if (data.jetsAK8.tau2[i]>0) tau32[i] = data.jetsAK8.tau3[i]/data.jetsAK8.tau2[i];
      // Maximum Subjet btag discriminator
      maxSubjetCSV[i] = -9999;
      int i_sj0 = data.jetsAK8.vSubjetIndex0[i], i_sj1 = data.jetsAK8.vSubjetIndex1[i];
      if (i_sj0 != -1) if (data.subjetsAK8.CSVv2[i_sj0] > maxSubjetCSV[i]) maxSubjetCSV[i] = data.subjetsAK8.CSVv2[i_sj0];
      if (i_sj1 != -1) if (data.subjetsAK8.CSVv2[i_sj1] > maxSubjetCSV[i]) maxSubjetCSV[i] = data.subjetsAK8.CSVv2[i_sj1];
#if USE_BTAG == 1
      if (passSubjetBTag[i] = (maxSubjetCSV[i] >= TOP_BTAG_CSV) ) nSubjetBTag++;
#else
      if (passSubjetBTag[i] = (maxSubjetCSV[i] >= B_SUBJET_CSV_LOOSE_CUT) ) nSubjetBTag++;
#endif
#else
      if (data.jetsAK8.tau1Puppi[i]>0) tau21[i] = data.jetsAK8.tau2Puppi[i]/data.jetsAK8.tau1Puppi[i];
      if (data.jetsAK8.tau1Puppi[i]>0) tau31[i] = data.jetsAK8.tau3Puppi[i]/data.jetsAK8.tau1Puppi[i];
      if (data.jetsAK8.tau2Puppi[i]>0) tau32[i] = data.jetsAK8.tau3Puppi[i]/data.jetsAK8.tau2Puppi[i];
      // Maximum Subjet btag discriminator
#if USE_BTAG == 1
      if (passSubjetBTag[i] = (data.jetsAK8.maxSubjetCSVv2[i] >= TOP_BTAG_CSV) ) nSubjetBTag++;
#else
      if (passSubjetBTag[i] = (data.jetsAK8.maxSubjetCSVv2[i] >= B_SUBJET_CSV_LOOSE_CUT) ) nSubjetBTag++;
#endif
#endif
    }
    if (debug) std::cout<<"AnalysisBase::calc_common_var: end init AK8"<<std::endl;

    // Event Letpons
    iEleVeto     .clear();
    iEleLoose    .clear();
    iEleSelect   .clear();
    iEleTight    .clear();
    itEleVeto    .assign(data.ele.size, (size_t)-1);
    itEleLoose   .assign(data.ele.size, (size_t)-1);
    itEleSelect  .assign(data.ele.size, (size_t)-1);
    itEleTight   .assign(data.ele.size, (size_t)-1);
    passEleVeto  .assign(data.ele.size, 0);
    passEleLoose .assign(data.ele.size, 0);
    passEleSelect.assign(data.ele.size, 0);
    passEleTight .assign(data.ele.size, 0);
    nEleVetoNoIso = nEleVeto = nEleLoose = nEleSelect = nEleTight = 0;
    while(data.ele.Loop()) {
      size_t i = data.ele.it;
      TLorentzVector ele_v4; ele_v4.SetPtEtaPhiE(data.ele.Pt[i], data.ele.Eta[i], data.ele.Phi[i], data.ele.E[i]);
      float pt = data.ele.Pt[i];
      float abseta = std::abs(data.ele.Eta[i]);
      float miniIso = data.ele.MiniIso[i];
#if USE_POG_ID == 0
      float absIso  = data.ele.Iso03[i]*pt;
#endif
      float absd0 = std::abs(data.ele.Dxy[i]);
      float absdz = std::abs(data.ele.Dz[i]);
      float ipsig = std::abs(data.ele.DB[i])/data.ele.DBerr[i];
#if USE_MVA_ID == 1
      // https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#Electrons
      //  VLoose WP, pT > 10 GeV. A/B values are -0.48/-0.85, -0.67/-0.91, -0.49/-0.83 for |eta| [0, 0.8], [0.8, 1.479], [1.479, 2.5]
      //  In addition, for pT 5-10, the following values are used on the HZZ MVA: 0.46, -0.03, 0.06 for |eta| [0, 0.8], [0.8, 1.479], [1.479, 2.5] 
      int categ = data.ele.vidMvaHZZcateg[i];
      double mva = categ<3 ? data.ele.vidMvaHZZvalue[i] : data.ele.vidMvaGPvalue[i];
      std::vector<double> cut = { 
	// pt<10
	0.46, -0.03, 0.06, 
	// pt>=10
	std::min(-0.48, std::max(-0.85 , -0.48 + (-0.85 - -0.48)/10.0*(pt-15))),
	std::min(-0.67, std::max(-0.91 , -0.67 + (-0.91 - -0.67)/10.0*(pt-15))),
	std::min(-0.49, std::max(-0.83 , -0.49 + (-0.83 - -0.49)/10.0*(pt-15)))
      };
      bool id_veto_noiso = (mva>cut[categ]);
#else
      bool id_veto_noiso = (data.ele.vidVetonoiso[i] == 1.0);
#endif
      bool id_loose_noiso = (data.ele.vidLoosenoiso[i] == 1.0);
      bool id_select_noiso = (data.ele.vidMediumnoiso[i] == 1.0);
      bool id_tight  = (data.ele.vidTight[i] == 1.0);
      //bool id_veto = (data.ele.vidVeto[i] == 1.0);
      //bool id_select = (data.ele.vidTight[i] == 1.0);
      // Veto
#if USE_POG_ID > 0
      if (passEleVeto[i] = 
	  ( id_veto_noiso &&
	    pt      >= ELE_VETO_PT_CUT &&
	    abseta  <  ELE_VETO_ETA_CUT && //!(abseta>=1.442 && abseta< 1.556) &&
	    absd0   <  ELE_VETO_IP_D0_CUT &&
	    absdz   <  ELE_VETO_IP_DZ_CUT) ) {
	veto_leptons_noiso.push_back(ele_v4);
	nEleVetoNoIso++;
	if (passEleVeto[i] = (miniIso <  ELE_VETO_MINIISO_CUT) ) {
	  iEleVeto.push_back(i);
	  itEleVeto[i] = nEleVeto++;
	  veto_electrons.push_back(ele_v4);
	  veto_leptons.push_back(ele_v4);
	  //veto_lep_in_jet.push_back(data.ele.IsPartOfNearAK4Jet[i]);
	}
      }
#else
      if (passEleVeto[i] = 
	  ( id_veto_noiso &&
	    pt      >= ELE_VETO_PT_CUT &&
	    abseta  <  ELE_VETO_ETA_CUT && //!(abseta>=1.442 && abseta< 1.556) &&
	    ipsig   <  ELE_VETO_IP_3D_CUT) ) {
	veto_leptons_noiso.push_back(ele_v4);
	nEleVetoNoIso++;
	if (passEleVeto[i] = 
	    (pt>20 ? miniIso <  ELE_VETO_MINIISO_CUT : absIso < ELE_VETO_ABSISO_CUT) ) {
	  iEleVeto.push_back(i);
	  itEleVeto[i] = nEleVeto++;
	  veto_electrons.push_back(ele_v4);
	  veto_leptons.push_back(ele_v4);
	  //veto_lep_in_jet.push_back(data.ele.IsPartOfNearAK4Jet[i]);
	}
      }
#endif
      // Loose
      if (passEleLoose[i] = 
	  ( id_loose_noiso &&
	    pt        >= ELE_LOOSE_PT_CUT &&
	    abseta    <  ELE_LOOSE_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
	    miniIso   <  ELE_LOOSE_MINIISO_CUT &&
	    absd0     <  ELE_LOOSE_IP_D0_CUT &&
	    absdz     <  ELE_LOOSE_IP_DZ_CUT) ) {
	iEleLoose.push_back(i);
	itEleLoose[i] = nEleLoose++;
      }
      // Select
      if (passEleSelect[i] = 
	  ( id_select_noiso &&
	    pt        >= ELE_SELECT_PT_CUT &&
	    abseta    <  ELE_SELECT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
	    miniIso   <  ELE_SELECT_MINIISO_CUT &&
	    absd0     <  ELE_SELECT_IP_D0_CUT &&
	    absdz     <  ELE_SELECT_IP_DZ_CUT) ) {
	selected_leptons.push_back(ele_v4);
	iEleSelect.push_back(i);
	itEleSelect[i] = nEleSelect++;
      }
      // Tight
      if (passEleTight[i] = 
	  ( id_tight &&
	    pt        >= ELE_TIGHT_PT_CUT &&
	    abseta    <  ELE_TIGHT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
	    absd0     <  ELE_TIGHT_IP_D0_CUT &&
	    absdz     <  ELE_TIGHT_IP_DZ_CUT &&
	    ipsig     <  ELE_TIGHT_IP_SIG_CUT) ) {
	tight_leptons.push_back(ele_v4);
	iEleTight.push_back(i);
	itEleTight[i] = nEleTight++;
	float r_iso = std::max(0.05, std::min(0.2, 10./pt));
	r_iso_tight_leptons.push_back(r_iso);
      }
    }
    if (debug) std::cout<<"AnalysisBase::calc_common_var: end init ele"<<std::endl;

    // Number of Veto/Select Muons
    iMuVeto      .clear();
    iMuLoose     .clear();
    iMuSelect    .clear();
    iMuTight     .clear();
    itMuVeto     .assign(data.mu.size,  (size_t)-1);
    itMuLoose    .assign(data.mu.size,  (size_t)-1);
    itMuSelect   .assign(data.mu.size,  (size_t)-1);
    itMuTight    .assign(data.mu.size,  (size_t)-1);
    passMuVeto   .assign(data.mu.size,  0);
    passMuLoose  .assign(data.mu.size,  0);
    passMuSelect .assign(data.mu.size,  0);
    passMuTight  .assign(data.mu.size,  0);
    nMuVetoNoIso = nMuVeto = nMuLoose = nMuSelect = nMuTight = 0;
    while(data.mu.Loop()) {
      size_t i = data.mu.it;
      TLorentzVector mu_v4; mu_v4.SetPtEtaPhiE(data.mu.Pt[i], data.mu.Eta[i], data.mu.Phi[i], data.mu.E[i]);
      float pt = data.mu.Pt[i];
      float abseta = std::abs(data.mu.Eta[i]);
      float miniIso = data.mu.MiniIso[i];
      float relIso = data.mu.Iso04[i];
#if USE_POG_ID == 0
      float absIso  = data.mu.Iso04[i]*pt;
#endif
      float absd0 = std::abs(data.mu.Dxy[i]);
      float absdz = std::abs(data.mu.Dz[i]);
      float ipsig = std::abs(data.mu.DB[i])/data.mu.DBerr[i];
      bool id_veto_noiso   = (data.mu.IsLooseMuon[i] == 1.0);
      bool id_loose_noiso  = (data.mu.IsLooseMuon[i] == 1.0);
      bool id_select_noiso = (data.mu.IsMediumMuon[i] == 1.0);
      bool id_tight_noiso  = (data.mu.IsTightMuon[i] == 1.0);
      // Veto
#if USE_POG_ID > 0
      if (passMuVeto[i] =
	  (id_veto_noiso &&
	   pt      >= MU_VETO_PT_CUT &&
	   abseta  <  MU_VETO_ETA_CUT &&
	   absd0   <  MU_VETO_IP_D0_CUT &&
	   absdz   <  MU_VETO_IP_DZ_CUT) ) {
	veto_leptons_noiso.push_back(mu_v4);
	//veto_muons_noiso.push_back(mu_v4);
	nMuVetoNoIso++;
	if (passMuVeto[i] = (miniIso <  MU_VETO_MINIISO_CUT) ) {
	  iMuVeto.push_back(i);
	  itMuVeto[i] = nMuVeto++;
	  veto_muons.push_back(mu_v4);
	  veto_leptons.push_back(mu_v4);
	  //veto_muons.push_back(mu_v4);
	  //veto_lep_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
	  //veto_mu_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
	}
      }
#else
      if (passMuVeto[i] =
	  (id_veto_noiso &&
	   pt      >= MU_VETO_PT_CUT &&
	   abseta  <  MU_VETO_ETA_CUT &&
	   ipsig   <  MU_VETO_IP_3D_CUT) ) {
	veto_leptons_noiso.push_back(mu_v4);
	//veto_muons_noiso.push_back(mu_v4);
	nMuVetoNoIso++;
	if (passMuVeto[i] =
	    ( pt>20 ? miniIso <  MU_VETO_MINIISO_CUT : absIso < MU_VETO_ABSISO_CUT ) ) {
	  iMuVeto.push_back(i);
	  itMuVeto[i] = nMuVeto++;
	  veto_muons.push_back(mu_v4);
	  veto_leptons.push_back(mu_v4);
	  //veto_muons.push_back(mu_v4);
	  //veto_lep_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
	  //veto_mu_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
	}
      }
#endif
      // Loose
      if (passMuLoose[i] =
	  ( id_loose_noiso &&
	    pt      >= MU_LOOSE_PT_CUT &&
	    abseta  <  MU_LOOSE_ETA_CUT &&
	    miniIso <  MU_LOOSE_MINIISO_CUT &&
	    absd0   <  MU_LOOSE_IP_D0_CUT &&
	    absdz   <  MU_LOOSE_IP_DZ_CUT) ) {
	iMuLoose.push_back(i);
	itMuLoose[i] = nMuLoose++;
      }
      // Select
      if (passMuSelect[i] =
	  ( id_select_noiso &&
	    pt      >= MU_SELECT_PT_CUT &&
	    abseta  <  MU_SELECT_ETA_CUT &&
	    miniIso <  MU_SELECT_MINIISO_CUT &&
	    absd0   <  MU_SELECT_IP_D0_CUT &&
	    absdz   <  MU_SELECT_IP_DZ_CUT) ) {
	selected_leptons.push_back(mu_v4);
	//selected_muons.push_back(mu_v4);
	iMuSelect.push_back(i);
	itMuSelect[i] = nMuSelect++;
	//selected_mu_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
      }
      // Tight
      if (passMuTight[i] =
	  ( id_tight_noiso &&
	    pt      >= MU_TIGHT_PT_CUT &&
	    abseta  <  MU_TIGHT_ETA_CUT &&
	    relIso  <  MU_TIGHT_RELISO_CUT &&
	    absd0   <  MU_TIGHT_IP_D0_CUT &&
	    absdz   <  MU_TIGHT_IP_DZ_CUT &&
	    ipsig   <  MU_TIGHT_IP_SIG_CUT) ) {
	tight_leptons.push_back(mu_v4);
	iMuTight.push_back(i);
	itMuTight[i] = nMuTight++;
	float r_iso = std::max(0.05, std::min(0.2, 10./pt));
	r_iso_tight_leptons.push_back(r_iso);
      }
    } // end of muon loop
    if (debug) std::cout<<"AnalysisBase::calc_common_var: end init mu"<<std::endl;

    nLepVetoNoIso = nEleVetoNoIso + nMuVetoNoIso;
    nLepVeto      = nEleVeto   + nMuVeto;
    nLepLoose     = nEleLoose  + nMuLoose;
    nLepSelect    = nEleSelect + nMuSelect;
    nLepTight     = nEleTight  + nMuTight;

    // Number of Veto Taus
    iTauVeto      .clear();
    itTauVeto     .assign(data.tau.size,  (size_t)-1);
    passTauVeto   .assign(data.tau.size,  0);
    nTauVeto = 0;
    while(data.tau.Loop()) {
      size_t i = data.tau.it;
      //TLorentzVector tau_v4; tau_v4.SetPtEtaPhiE(data.tau.Pt[i], data.tau.Eta[i], data.tau.Phi[i], data.tau.E[i]);
      //float pt = data.tau.Pt[i];
      //float abseta = std::abs(data.tau.Eta[i]);
      bool id_veto = (data.tau.IsLoose[i] == 1);
      // Veto
      if (passTauVeto[i] = 
	  ( id_veto //&&
	    //	    pt      >= TAU_VETO_PT_CUT &&
	    //	    abseta  <  TAU_VETO_ETA_CUT 
	    ) ) {
	iTauVeto.push_back(i);
	itTauVeto[i] = nTauVeto++;
      }
    }
    if (debug) std::cout<<"AnalysisBase::calc_common_var: end init tau"<<std::endl;

    // MT
    MT = 9999;
    if (nLepSelect==1) {
      if (nEleSelect==1) {
	MT = sqrt( 2*data.ele.Pt[iEleSelect[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.ele.Phi[iEleSelect[0]])) );
      } else if (nMuSelect==1) {
	MT = sqrt( 2*data.mu.Pt[iMuSelect[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.mu.Phi[iMuSelect[0]])) );
      }
    }
    MT_vetolep = 9999;
    if (nLepVeto==1) {
      if (nEleVeto==1) {
	MT_vetolep = sqrt( 2*data.ele.Pt[iEleVeto[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.ele.Phi[iEleVeto[0]])) );
      } else if (nMuVeto==1) {
	MT_vetolep = sqrt( 2*data.mu.Pt[iMuVeto[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.mu.Phi[iMuVeto[0]])) );
      }
    }
    MT_allvetolep = 9999;
    vetoLepPt  = -9999;
    vetoLepPhi = -9999;
    vetoLepEta = -9999;
    if (nEleVeto+nMuVeto+nTauVeto==1) {
      if (nEleVeto==1) {
	vetoLepPt  = data.ele.Pt[iEleVeto[0]];
	vetoLepPhi = data.ele.Phi[iEleVeto[0]];
	vetoLepEta = data.ele.Eta[iEleVeto[0]];
      } else if (nMuVeto==1) {
	vetoLepPt  = data.mu.Pt[iMuVeto[0]];
	vetoLepPhi = data.mu.Phi[iMuVeto[0]];
	vetoLepEta = data.mu.Eta[iMuVeto[0]];
      } else if (nTauVeto==1) {
	vetoLepPt  = data.tau.Pt[iTauVeto[0]];
	vetoLepPhi = data.tau.Phi[iTauVeto[0]];
	vetoLepEta = data.tau.Eta[iTauVeto[0]];
      }
      MT_allvetolep = sqrt( 2*vetoLepPt*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-vetoLepPhi)) );
    }
    
    // M_ll, dPhi_ll_met
    M_ll = -9999;
    dPhi_ll_met = 9999;
    if (nLepSelect==2) {
      TLorentzVector lep1, lep2;
      if (nEleSelect==2) {
	lep1.SetPtEtaPhiE(data.ele.Pt[iEleSelect[0]], data.ele.Eta[iEleSelect[0]], data.ele.Phi[iEleSelect[0]], data.ele.E[iEleSelect[0]]);
	lep2.SetPtEtaPhiE(data.ele.Pt[iEleSelect[1]], data.ele.Eta[iEleSelect[1]], data.ele.Phi[iEleSelect[1]], data.ele.E[iEleSelect[1]]);
	lep_pair = lep1+lep2;
	M_ll = lep_pair.M();
	dPhi_ll_met = std::abs(TVector2::Phi_mpi_pi(lep_pair.Phi() - data.met.Phi[0]));
      } else if (nMuSelect==2) {
	lep1.SetPtEtaPhiE(data.mu.Pt[iMuSelect[0]], data.mu.Eta[iMuSelect[0]], data.mu.Phi[iMuSelect[0]], data.mu.E[iMuSelect[0]]);
	lep2.SetPtEtaPhiE(data.mu.Pt[iMuSelect[1]], data.mu.Eta[iMuSelect[1]], data.mu.Phi[iMuSelect[1]], data.mu.E[iMuSelect[1]]);
	lep_pair = lep1+lep2;
	M_ll = lep_pair.M();
	dPhi_ll_met = std::abs(TVector2::Phi_mpi_pi(lep_pair.Phi() - data.met.Phi[0]));
      }
    }
    if (debug) std::cout<<"AnalysisBase::calc_common_var: end init razor"<<std::endl;
  } // end if (syst_index==0)
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end init"<<std::endl;
    
  // Event Letpons
#if VER < 4
  // Do a first loop to reverse engineer the rho value (which is used later)
  double rho_from_other_iso = 0;
  while(data.pho.Loop()) {
    size_t i = data.pho.it;
    float abseta = std::abs(data.pho.Eta[i]);
    double EA_neu = 0, EA_pho = 0;
    if (abseta>0.0 && abseta<=1.0) {
      EA_neu = 0.0597;
      EA_pho = 0.1210;
    } else if (abseta>1.0 && abseta<=1.479) {
      EA_neu = 0.0807;
      EA_pho = 0.1107;
    } else if (abseta>1.479 && abseta<=2.0) {
      EA_neu = 0.0629;
      EA_pho = 0.0699;
    } else if (abseta>2.0 && abseta<=2.2) {
      EA_neu = 0.0197;
      EA_pho = 0.1056;
    } else if (abseta>2.2 && abseta<=2.3) {
      EA_neu = 0.0184;
      EA_pho = 0.1457;
    } else if (abseta>2.3 && abseta<=2.4) {
      EA_neu = 0.0284;
      EA_pho = 0.1719;
    } else if (abseta>2.4 && abseta<=5.0) {
      EA_neu = 0.0591;
      EA_pho = 0.1998;
    }
    // rho was not saved either, but can be deduced from the corrected and uncorrected iso
    double rho = 0;
    if (data.pho.NeutralHadronIsoEAcorrectedsp15[i]>0) {
      rho = (data.pho.NeutralHadronIso[i] - data.pho.NeutralHadronIsoEAcorrectedsp15[i])/EA_neu;
    } else if (data.pho.PhotonIsoEAcorrectedsp15[i]) {
      rho = (data.pho.PhotonIso[i] - data.pho.PhotonIsoEAcorrectedsp15[i])/EA_pho;
    }
    if (rho_from_other_iso==0&&rho>0) rho_from_other_iso = rho;
  }
#endif
  // Count photons
  std::vector<TLorentzVector> photons_PreSelect, selected_photons;
  std::vector<TLorentzVector> faked_photons;
  iPhotonPreSelect   .clear();
  iPhotonSelect      .clear();
  iPhotonFake        .clear();
  itPhotonPreSelect  .assign(data.pho.size, (size_t)-1);
  itPhotonSelect     .assign(data.pho.size, (size_t)-1);
  itPhotonFake       .assign(data.pho.size, (size_t)-1);
  passPhotonPreSelect.assign(data.pho.size, 0);
  passPhotonSelect   .assign(data.pho.size, 0);
  passPhotonFake   .assign(data.pho.size, 0);
  ChargedHadronIsoEACorr.clear();
  nPhotonPreSelect = nPhotonSelect = nPhotonFake = 0;
  while(data.pho.Loop()) {
    size_t i = data.pho.it;
    TLorentzVector pho_v4; pho_v4.SetPtEtaPhiE(data.pho.Pt[i], data.pho.Eta[i], data.pho.Phi[i], data.pho.E[i]);
    float pt = data.pho.Pt[i];
    float abseta = std::abs(data.pho.Eta[i]);
    bool ele_veto = (data.pho.ElectronVeto[i]==1);
    bool id_select = data.pho.PassMediumID[i];
#if VER < 4
    // For cut based ID, need all cut variables at hand
    // But the EA corrected value was not saved for the Charged isolation (in VER 3)
    // --> so calculate it
    double EA_ch = 0;
    if      (abseta>0.0 && abseta<=1.0)   EA_ch  = 0.0360;
    else if (abseta>1.0 && abseta<=1.479) EA_ch  = 0.0377;
    else if (abseta>1.479 && abseta<=2.0) EA_ch  = 0.0306;
    else if (abseta>2.0 && abseta<=2.2)   EA_ch  = 0.0283;
    else if (abseta>2.2 && abseta<=2.3)   EA_ch  = 0.0254;
    else if (abseta>2.3 && abseta<=2.4)   EA_ch  = 0.0217;
    else if (abseta>2.4 && abseta<=5.0)   EA_ch  = 0.0167;
    double ChargedHadronIsoEACorrSpring16 = std::max(0., data.pho.ChargedHadronIso[i] - rho_from_other_iso*EA_ch);
    ChargedHadronIsoEACorr.push_back(ChargedHadronIsoEACorrSpring16);
#else 
    ChargedHadronIsoEACorr.push_back(data.pho.ChargedHadronIsoEAcorrectedsp15[i]);
#endif
    // Implementing cuts for Spring16 v2.2 Medium ID
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=36#Recommended_Working_points_for_2
    // Using the table, not the presentation, same as this VID cfg:
    // https://github.com/ikrav/cmssw/blob/65fa87654e6744efdff2e16d55a5b86bbdccd48d/RecoEgamma/PhotonIdentification/python/Identification/cutBasedPhotonID_Spring16_V2p2_cff.py#L50-L78
    // Exclude Sigma_ietaieta cut in order to be able to measure QCD background
    bool id_PreSel = false;
    if ( std::abs(data.pho.SCEta[i])<1.479 ) {
      // Barrel cuts (EB)
      id_PreSel = 
  	data.pho.HoverE[i]                          < 0.0396 &&
  	//data.pho.SigmaIEtaIEta[i]                   < 0.01022 &&
  	//ChargedHadronIsoEACorr[i]                   < 0.441 &&
  	data.pho.NeutralHadronIsoEAcorrectedsp15[i] < 2.725+0.0148*pt+0.000017*pt*pt &&
  	data.pho.PhotonIsoEAcorrectedsp15[i]        < 2.571+0.0047*pt;
    } else {
      // Encap cuts (EE)
      id_PreSel = 
  	data.pho.HoverE[i]                          < 0.0219 &&
  	//data.pho.SigmaIEtaIEta[i]                   < 0.03001 &&
  	//ChargedHadronIsoEACorr[i]                   < 0.442 &&
  	data.pho.NeutralHadronIsoEAcorrectedsp15[i] < 1.715+0.0163*pt+0.000014*pt*pt &&
  	data.pho.PhotonIsoEAcorrectedsp15[i]        < 3.863+0.0034*pt;
    }
    // Medium ID without Sigma_ietaieta cut
    if ( passPhotonPreSelect[i] = 
	 ( id_PreSel &&
	   ele_veto &&
	   pt        >= PHOTON_SELECT_PT_CUT &&
	   abseta    <  PHOTON_SELECT_ETA_CUT ) ) {
      iPhotonPreSelect.push_back(i);
      itPhotonPreSelect[i] = nPhotonPreSelect++;
      photons_PreSelect.push_back(pho_v4);
      // Fake photons (those that fail the SigmaIeteIeta cut)
      if ( passPhotonFake[i] = 
	   (data.pho.SigmaIEtaIEta[i] >= (std::abs(data.pho.SCEta[i])<1.479 ? 0.01022 : 0.03001) ) ) {
	faked_photons.push_back(pho_v4);
      	iPhotonFake.push_back(i);
      	itPhotonFake[i] = nPhotonFake++;
      }
    }
    // Photons passing full ID
    if ( passPhotonSelect[i] = 
	 ( id_select &&
	   ele_veto &&
	   pt        >= PHOTON_SELECT_PT_CUT &&
	   abseta    <  PHOTON_SELECT_ETA_CUT ) ) {
      selected_photons.push_back(pho_v4);
      iPhotonSelect.push_back(i);
      itPhotonSelect[i] = nPhotonSelect++;
    }
    // Sanity check for cut based implementation
    //  bool id_select_manual= false;
    //  if ( std::abs(data.pho.SCEta[i])<1.479 ) {
    //    // Barrel cuts (EB)
    //    id_select_manual = 
    //      data.pho.HoverE[i]                          < 0.0396 &&
    //      data.pho.SigmaIEtaIEta[i]                   < 0.01022 &&
    //      data.pho.ChargedHadronIsoEAcorrectedsp15[i] < 0.441 &&
    //      //data.pho.ChargedHadronIso[i]                < 0.441 &&
    //      data.pho.NeutralHadronIsoEAcorrectedsp15[i] < 2.725+0.0148*pt+0.000017*pt*pt &&
    //      data.pho.PhotonIsoEAcorrectedsp15[i]        < 2.571+0.0047*pt;
    //  } else {
    //    // Encap cuts (EE)
    //    id_select_manual = 
    //      data.pho.HoverE[i]                          < 0.0219 &&
    //      data.pho.SigmaIEtaIEta[i]                   < 0.03001 &&
    //      data.pho.ChargedHadronIsoEAcorrectedsp15[i] < 0.442 &&
    //      //data.pho.ChargedHadronIso[i]                < 0.442 &&
    //      data.pho.NeutralHadronIsoEAcorrectedsp15[i] < 1.715+0.0163*pt+0.000014*pt*pt&&
    //      data.pho.PhotonIsoEAcorrectedsp15[i]        < 3.863+0.0034*pt;
    //  }
    //  if (id_select_manual != id_select) {
    //    std::cout<<"POG ID - Pass VID="<<id_select<<std::endl;
    //    if (std::abs(data.pho.SCEta[i]) < 1.479) {
    //      std::cout<<(std::abs(data.pho.SCEta[i])<1.479                                             )<<" SCEta="<<std::abs(data.pho.SCEta[i])<<std::endl;
    //      std::cout<<(data.pho.HoverE[i] < 0.0396                                                   )<<" HoverE="<<data.pho.HoverE[i]<<std::endl;
    //      std::cout<<(data.pho.SigmaIEtaIEta[i] < 0.01022                                           )<<" SigmaIEtaIEta="<<data.pho.SigmaIEtaIEta[i]<<std::endl;   
    //      std::cout<<(data.pho.ChargedHadronIso[i] < 0.441                                          )<<" ChargedHadronIso="<<data.pho.ChargedHadronIso[i]<<std::endl;
    //      std::cout<<(data.pho.ChargedHadronIsoEAcorrectedsp15[i] < 0.441                           )<<" ChargedHadronIsoEACorr="<<data.pho.ChargedHadronIsoEAcorrectedsp15[i]<<std::endl;
    //      std::cout<<(data.pho.NeutralHadronIsoEAcorrectedsp15[i] < (2.725+0.0148*pt+0.000017*pt*pt))<<" NeutralHadronIsoEAcorrectedsp15="<<data.pho.NeutralHadronIsoEAcorrectedsp15[i]<<" cut="<<(2.725+0.0148*pt+0.000017*pt*pt)<<std::endl;
    //      std::cout<<(data.pho.PhotonIsoEAcorrectedsp15[i] < 2.571+0.0047*pt                        )<<" PhotonIsoEAcorrectedsp15="<<data.pho.PhotonIsoEAcorrectedsp15[i]<<" cut="<<(2.571+0.0047*pt)<<std::endl<<std::endl;
    //    } else {
    //      std::cout<<(std::abs(data.pho.SCEta[i])>=1.479                                            )<<" SCEta="<<std::abs(data.pho.SCEta[i])<<std::endl;
    //      std::cout<<(data.pho.HoverE[i] < 0.0219                                                   )<<" HoverE="<<data.pho.HoverE[i]<<std::endl;
    //      std::cout<<(data.pho.SigmaIEtaIEta[i] < 0.03001                                           )<<" SigmaIEtaIEta="<<data.pho.SigmaIEtaIEta[i]<<std::endl;   
    //      std::cout<<(data.pho.ChargedHadronIso[i] < 0.442                                          )<<" ChargedHadronIso="<<data.pho.ChargedHadronIso[i]<<std::endl;
    //      std::cout<<(data.pho.ChargedHadronIsoEAcorrectedsp15[i] < 0.442                           )<<" ChargedHadronIsoEACorr="<<data.pho.ChargedHadronIsoEAcorrectedsp15[i]<<std::endl;
    //      std::cout<<(data.pho.NeutralHadronIsoEAcorrectedsp15[i] < (1.715+0.0163*pt+0.000014*pt*pt))<<" NeutralHadronIsoEAcorrectedsp15="<<data.pho.NeutralHadronIsoEAcorrectedsp15[i]<<" cut="<<(2.725+0.0148*pt+0.000017*pt*pt)<<std::endl;
    //      std::cout<<(data.pho.PhotonIsoEAcorrectedsp15[i] < 3.863+0.0034*pt                        )<<" PhotonIsoEAcorrectedsp15="<<data.pho.PhotonIsoEAcorrectedsp15[i]<<" cut="<<(2.571+0.0047*pt)<<std::endl<<std::endl;
    //    }
    //  }
  }
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end pho"<<std::endl;

  // Add the lepton (pair) to MET
  TVector3 met_1l, met_1vl, met_ll;
  met_1l. SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);
  met_1vl.SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);
  met_ll. SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);
  if (nLepSelect==1) {
    TVector3 lep_met;
    lep_met.SetPtEtaPhi(selected_leptons[0].Pt(), 0, selected_leptons[0].Phi());
    met_1l += lep_met;
  }
  if (nLepVeto==1) {
    TVector3 lep_met;
    lep_met.SetPtEtaPhi(veto_leptons[0].Pt(), 0, veto_leptons[0].Phi());
    met_1vl += lep_met;
  }
  if (M_ll!=-9999) {
    TVector3 lep_pair_met;
    lep_pair_met.SetPtEtaPhi(lep_pair.Pt(), 0, lep_pair.Phi());
    met_ll += lep_pair_met;
  }

  // Add the photon to MET
  TVector3 met_pho;
  met_pho.SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);
  if (nPhotonSelect==1) {
    TVector3 pho_met;
    pho_met.SetPtEtaPhi(selected_photons[0].Pt(), 0, selected_photons[0].Phi());
    met_pho += pho_met;
  } else if (nPhotonPreSelect==1) {
    // Add also preselected photons
    TVector3 pho_met;
    pho_met.SetPtEtaPhi(photons_PreSelect[0].Pt(), 0, photons_PreSelect[0].Phi());
    met_pho += pho_met;
  }

  // Add the fake photon to MET
  TVector3 met_fakepho;
  met_fakepho.SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);
  if (nPhotonFake==1) {
    TVector3 pho_met;
    pho_met.SetPtEtaPhi(faked_photons[0].Pt(), 0, faked_photons[0].Phi());
    met_fakepho += pho_met;
  }
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end met+lep/pho"<<std::endl;

  // Rest of the vairables need to be recalculated each time the jet energy is changed
  // eg. Jet selection, W/top tags, HT (obviously), etc. that depends on jet pt

  // AK4 jets
  iJet       .clear();
  iLooseBTag .clear();
  iMediumBTag.clear();
  iTightBTag .clear();
  itJet              .assign(data.jetsAK4.size, (size_t)-1);
  itLooseBTag        .assign(data.jetsAK4.size, (size_t)-1);
  itMediumBTag       .assign(data.jetsAK4.size, (size_t)-1);
  itTightBTag        .assign(data.jetsAK4.size, (size_t)-1);
  passLooseJet       .assign(data.jetsAK4.size, 0);
  passLooseJetNoLep  .assign(data.jetsAK4.size, 0);
  passLooseJetNoPho  .assign(data.jetsAK4.size, 0);
  passLooseBTag      .assign(data.jetsAK4.size, 0);
  passMediumBTag     .assign(data.jetsAK4.size, 0);
  passTightBTag      .assign(data.jetsAK4.size, 0);
  AK4_photonDR       .assign(data.jetsAK4.size, 9999);
  AK4_photonPtRatio  .assign(data.jetsAK4.size, 0);
  AK4_eleDR          .assign(data.jetsAK4.size, 9999);
  AK4_elePtRatio     .assign(data.jetsAK4.size, 0);
  AK4_muonDR         .assign(data.jetsAK4.size, 9999);
  AK4_muonPtRatio    .assign(data.jetsAK4.size, 0);
  nJet = nJetNoLep = nJetNoPho = 0;
  nLooseBTag  = 0;
  nMediumBTag = nMediumBTagNoPho = 0;
  nTightBTag  = 0;
  AK4_Ht = AK4_HtOnline = AK4_HtNoLep = 0;
  minDeltaPhi = minDeltaPhi_1l = minDeltaPhi_1vl = minDeltaPhi_ll = minDeltaPhi_pho = dPhi_ll_jet = 9999;
  //std::vector<bool> add_lepton_to_ht(veto_leptons.size(),1);
  //std::vector<bool> remove_muon_from_ht(selected_muons.size(),0);
  while(data.jetsAK4.Loop()) {
    size_t i = data.jetsAK4.it;
    TLorentzVector jet_v4; jet_v4.SetPtEtaPhiE(data.jetsAK4.Pt[i], data.jetsAK4.Eta[i], data.jetsAK4.Phi[i], data.jetsAK4.E[i]);
    // Jet ID
    if ( passLooseJet[i] = 
	 ( data.jetsAK4.looseJetID[i] == 1 &&
	   data.jetsAK4.Pt[i]         >= JET_AK4_PT_CUT &&
	   std::abs(data.jetsAK4.Eta[i])  <  JET_AK4_ETA_CUT ) ) {

      // Photon jet veto (not applied)
      // - Calculate DR and Pt ratio of photons overlapping with the jet
      // - Remove jet if DR<0.1, otherwise
      //   if the photon is within the cone, subtract it's v4
      // - Finally, recheck that the remnant jet passes the ID cuts
      //size_t jMatch = -1;
      for (size_t j=0; j<selected_photons.size(); ++j) {
	double dR = selected_photons[j].DeltaR(jet_v4);
	if (dR<AK4_photonDR[i]) {
	  AK4_photonDR[i] = dR;
	  AK4_photonPtRatio[i] = selected_photons[j].Pt()/jet_v4.Pt();
	  //jMatch = j;
	}
      }
      //if (AK4_photonDR[i]<0.4) {
      //  jet_v4 -= selected_photons[jMatch];
      //  data.jetsAK4.Pt[i]  = jet_v4.Pt();
      //  data.jetsAK4.Eta[i] = jet_v4.Eta();
      //  data.jetsAK4.Phi[i] = jet_v4.Phi();
      //}
      //if ( passLooseJet[i] = 
      //     ( AK4_photonDR[i]>0.1 &&
      //       data.jetsAK4.Pt[i]         >= JET_AK4_PT_CUT &&
      //       std::abs(data.jetsAK4.Eta[i])  <  JET_AK4_ETA_CUT ) ) {

      // Lepton-jet overlap
      for (size_t j=0; j<veto_electrons.size(); ++j) {
	double dR = veto_electrons[j].DeltaR(jet_v4);
	if (dR<AK4_eleDR[i]) {
	  AK4_eleDR[i] = dR;
	  AK4_elePtRatio[i] = veto_electrons[j].Pt()/jet_v4.Pt();
	}
      }
      for (size_t j=0; j<veto_muons.size(); ++j) {
	double dR = veto_muons[j].DeltaR(jet_v4);
	if (dR<AK4_muonDR[i]) {
	  AK4_muonDR[i] = dR;
	  AK4_muonPtRatio[i] = veto_muons[j].Pt()/jet_v4.Pt();
	}
      }

      nJet++;
      AK4_Ht += data.jetsAK4.Pt[data.jetsAK4.it];

      // B tagging
      if (passLooseBTag[i]  = (data.jetsAK4.CSVv2[i] >= B_CSV_LOOSE_CUT ) ) {
	iLooseBTag.push_back(i);
	itLooseBTag[i] = nLooseBTag++;
      }
      if (passMediumBTag[i] = (data.jetsAK4.CSVv2[i] >= B_CSV_MEDIUM_CUT) ) {
	iMediumBTag.push_back(i);
	itMediumBTag[i] = nMediumBTag++;
      }
      if (passTightBTag[i]  = (data.jetsAK4.CSVv2[i] >= B_CSV_TIGHT_CUT ) ) {
	iTightBTag.push_back(i);
	itTightBTag[i] = nTightBTag++;
      }

	
      // minDeltaPhi
      if (nJet<=4) {
	double dphi = std::abs(TVector2::Phi_mpi_pi(data.met.Phi[0] - data.jetsAK4.Phi[i]));
	if (dphi<minDeltaPhi) minDeltaPhi = dphi;
	// with added (veto) lepton
	double dphi_met1l = std::abs(TVector2::Phi_mpi_pi(met_1l.Phi() - data.jetsAK4.Phi[i]));
	if (dphi_met1l<minDeltaPhi_1l) minDeltaPhi_1l = dphi_met1l;
	double dphi_met1vl = std::abs(TVector2::Phi_mpi_pi(met_1vl.Phi() - data.jetsAK4.Phi[i]));
	if (dphi_met1vl<minDeltaPhi_1vl) minDeltaPhi_1vl = dphi_met1vl;
	// with added lepton pair
	double dphi_metll = std::abs(TVector2::Phi_mpi_pi(met_ll.Phi() - data.jetsAK4.Phi[i]));
	if (dphi_metll<minDeltaPhi_ll) minDeltaPhi_ll = dphi_metll;
	// with added photon
	double dphi_metpho = std::abs(TVector2::Phi_mpi_pi(met_pho.Phi() - data.jetsAK4.Phi[i]));
	if (dphi_metpho<minDeltaPhi_pho) minDeltaPhi_pho = dphi_metpho;
	// jet lep-pair angle
	if (M_ll!=-9999) {
	  double dphi_ll = std::abs(TVector2::Phi_mpi_pi(lep_pair.Phi() - data.jetsAK4.Phi[i]));
	  if (dphi_ll<dPhi_ll_jet) dPhi_ll_jet = dphi_ll;
	}
      }

      // Exclude jets that have tight leptons in the isolation cone for the DeltaPhi calculation
      float minDR = 9999;
      float r_iso = -9999;
      for (size_t i=0, n=tight_leptons.size(); i<n; ++i) {
	double DR = tight_leptons[i].DeltaR(jet_v4);
	if (DR<minDR) {
	  minDR=DR;
	  r_iso = r_iso_tight_leptons[i];
	}
      }

      if (passLooseJetNoLep[i] = (minDR>=r_iso)) {
	iJet.push_back(i);
	itJet[i] = nJetNoLep++;
	AK4_HtNoLep += data.jetsAK4.Pt[data.jetsAK4.it];
      }

      // Exclude jets that have overlapping photon
      if (nPhotonSelect==1) {
	double dR = selected_photons[0].DeltaR(jet_v4);
	double pt_ratio = selected_photons[0].Pt()/jet_v4.Pt();
	//if (passLooseJetNoPho[i] = (dR>=0.4)) {
	//if (passLooseJetNoPho[i] = (dR>=0.4 && pt_ratio>=0.5&&pt_ratio<2.0)) {
	if (passLooseJetNoPho[i] = !(dR<0.4 && pt_ratio>=0.5&&pt_ratio<2.0)) {
	  nJetNoPho++;
	  if (data.jetsAK4.CSVv2[i] >= B_CSV_MEDIUM_CUT) nMediumBTagNoPho++;
	} //else std::cout<<"dR="<<dR<<" pt (pho, jet)="<<selected_photons[0].Pt()<<" "<<jet_v4.Pt()<<std::endl;	  
      } else if (nPhotonPreSelect==1) {
	// Also do this for the preselected photon
	double dR = photons_PreSelect[0].DeltaR(jet_v4);
	double pt_ratio = photons_PreSelect[0].Pt()/jet_v4.Pt();
	//if (passLooseJetNoPho[i] = (dR>=0.4)) {
	//if (passLooseJetNoPho[i] = (dR>=0.4 && pt_ratio>=0.5&&pt_ratio<2.0)) {
	if (passLooseJetNoPho[i] = !(dR<0.4 && pt_ratio>=0.5&&pt_ratio<2.0)) {
	  nJetNoPho++;
	  if (data.jetsAK4.CSVv2[i] >= B_CSV_MEDIUM_CUT) nMediumBTagNoPho++;
	}
      }

      //} // End photon jet veto

    } // End Jet Selection
    
    // Online jet selection for HT (+ testing Additional Loose Jet ID)
    if ( //data.jetsAK4.looseJetID[i] == 1 &&
	data.jetsAK4.Pt[i]         >  30 &&
	std::abs(data.jetsAK4.Eta[i])  <  3.0 ) {
      AK4_HtOnline += data.jetsAK4.Pt[data.jetsAK4.it];

      // Lepton (complete) isolation from jets
      // float minDR_lep = 9999; 
      // int ilep_match = -1;
      // for (size_t ilep=0, nlep=veto_leptons.size(); ilep<nlep; ++ilep) {
      //   float dR_lep = jet_v4.DeltaR(veto_leptons[ilep]);
      //   if (dR_lep < minDR_lep) {
      //     minDR_lep = dR_lep;
      //     ilep_match = ilep;
      //   }
      // }
      // if (minDR_lep<0.4 && veto_lep_in_jet[ilep_match]) add_lepton_to_ht[ilep_match] = 0;

      // Muons inside the jet
      //float minDR_mu = 9999; 
      //int imu_match = -1;
      //for (size_t imu=0, nmu=selected_muons.size(); imu<nmu; ++imu) {
      //  float dR_mu = jet_v4.DeltaR(selected_muons[imu]);
      //  if (dR_mu < minDR_mu) {
      //    minDR_mu = dR_mu;
      //    imu_match = imu;
      //  }
      //}
      //if (minDR_mu<0.4 && selected_mu_in_jet[imu_match]) remove_muon_from_ht[imu_match] = 1;

    }
  } // End AK4 Jet Loop

  // Nearest jet to lepton
  eleJetDR     .assign(data.ele.size,  9999);
  eleJetDPhi   .assign(data.ele.size,  9999);
  eleJetPt     .assign(data.ele.size, -9999);
  while(data.ele.Loop()) {
    size_t i = data.ele.it;
    TLorentzVector ele_v4; ele_v4.SetPtEtaPhiE(data.ele.Pt[i], data.ele.Eta[i], data.ele.Phi[i], data.ele.E[i]);
    // Nearest jet info
    while(data.jetsAK4.Loop()) {
      size_t j = data.jetsAK4.it;
      if (passLooseJet[j]) {
	TLorentzVector jet_v4; jet_v4.SetPtEtaPhiE(data.jetsAK4.Pt[j], data.jetsAK4.Eta[j], data.jetsAK4.Phi[j], data.jetsAK4.E[j]);
	double dR = ele_v4.DeltaR(jet_v4);
	if (dR<eleJetDR[i]) {
	  eleJetDR[i] = dR;
	  eleJetPt[i] = data.jetsAK4.Pt[j];
	}
	double dPhi = std::abs(ele_v4.DeltaPhi(jet_v4));
	if (dPhi<eleJetDPhi[i]) eleJetDPhi[i] = dPhi;
      }
    }
  }
  muJetDR     .assign(data.mu.size,  9999);
  muJetDPhi   .assign(data.mu.size,  9999);
  muJetPt     .assign(data.mu.size, -9999);
  while(data.mu.Loop()) {
    size_t i = data.mu.it;
    TLorentzVector mu_v4; mu_v4.SetPtEtaPhiE(data.mu.Pt[i], data.mu.Eta[i], data.mu.Phi[i], data.mu.E[i]);
    // Nearest jet info
    while(data.jetsAK4.Loop()) {
      size_t j = data.jetsAK4.it;
      if (passLooseJet[j]) {
	TLorentzVector jet_v4; jet_v4.SetPtEtaPhiE(data.jetsAK4.Pt[j], data.jetsAK4.Eta[j], data.jetsAK4.Phi[j], data.jetsAK4.E[j]);
	double dR = mu_v4.DeltaR(jet_v4);
	if (dR<muJetDR[i]) {
	  muJetDR[i] = dR;
	  muJetPt[i] = data.jetsAK4.Pt[j];
	}
	double dPhi = std::abs(mu_v4.DeltaPhi(jet_v4));
	if (dPhi<muJetDPhi[i]) muJetDPhi[i] = dPhi;
      }
    }
  }
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end AK4"<<std::endl;


  // Add isolated leptons to HT computation
  //for (size_t ilep=0, nlep=veto_leptons.size(); ilep<nlep; ++ilep)
  //  if (add_lepton_to_ht[ilep]) AK4_Ht += veto_leptons[ilep].Pt();

  //for (size_t imu=0, nmu=selected_muons.size(); imu<nmu; ++imu)
  //  if (remove_muon_from_ht[imu]) AK4_Ht -= selected_muons[imu].Pt();
  
  // AK8 jets
  iJetAK8         .clear();
  iWMassTag       .clear();
  iLooseWTag      .clear();
  iTightWTag      .clear();
  iTightWAntiTag  .clear();
  iHadTopMassTag  .clear();
  iHadTop1BMassTag  .clear();
  iHadTopTag      .clear();
  iHadTop0BMassTag.clear();
  iHadTop0BAntiTag.clear();
  softDropMassW   .clear();
  softDropMassTop .clear();
  softDropMassTop1.clear();
  softDropMassTop2.clear();
  itJetAK8             .assign(data.jetsAK8.size, (size_t)-1);
  itWMassTag           .assign(data.jetsAK8.size, (size_t)-1);
  itLooseWTag          .assign(data.jetsAK8.size, (size_t)-1);
  itTightWTag          .assign(data.jetsAK8.size, (size_t)-1);
  itTightWAntiTag      .assign(data.jetsAK8.size, (size_t)-1);
  itHadTopMassTag      .assign(data.jetsAK8.size, (size_t)-1);
  itHadTop1BMassTag    .assign(data.jetsAK8.size, (size_t)-1);
  itHadTopTag          .assign(data.jetsAK8.size, (size_t)-1);
  itHadTop0BMassTag    .assign(data.jetsAK8.size, (size_t)-1);
  itHadTop0BAntiTag    .assign(data.jetsAK8.size, (size_t)-1);
  passLooseJetAK8      .assign(data.jetsAK8.size, 0);
  passWMassTag         .assign(data.jetsAK8.size, 0);
  passLooseWTag        .assign(data.jetsAK8.size, 0);
  passTightWTag        .assign(data.jetsAK8.size, 0);
  passTightWAntiTag    .assign(data.jetsAK8.size, 0);
  passHadTopTag        .assign(data.jetsAK8.size, 0);
  passHadTopMassTag    .assign(data.jetsAK8.size, 0);
  passHadTop1BMassTag  .assign(data.jetsAK8.size, 0);
  passHadTop0BMassTag  .assign(data.jetsAK8.size, 0);
  passHadTop0BAntiTag  .assign(data.jetsAK8.size, 0);
  passOtherWMassTag       .assign(data.jetsAK8.size, 0);
  passOtherTightWTag      .assign(data.jetsAK8.size, 0);
  passOtherTightWAntiTag  .assign(data.jetsAK8.size, 0);
  passOtherHadTopTag      .assign(data.jetsAK8.size, 0);
  passOtherHadTopMassTag  .assign(data.jetsAK8.size, 0);
  passOtherHadTop0BMassTag.assign(data.jetsAK8.size, 0);
  passOtherHadTop0BAntiTag.assign(data.jetsAK8.size, 0);
  AK8_photonDR         .assign(data.jetsAK8.size, 9999);
  AK8_photonPtRatio    .assign(data.jetsAK8.size, 0);
  AK8_eleDR            .assign(data.jetsAK8.size, 9999);
  AK8_elePtRatio       .assign(data.jetsAK8.size, 0);
  AK8_muonDR           .assign(data.jetsAK8.size, 9999);
  AK8_muonPtRatio      .assign(data.jetsAK8.size, 0);
  nJetAK8          = 0;
  nWMassTag        = 0;
  nLooseWTag       = 0;
  nTightWTag       = 0;
  nHadWTag1        = 0;
  nHadWTag3        = 0;
  nTightWAntiTag   = 0;
  nSubjetBTag      = 0;
  nHadTopMassTag   = 0;
  nHadTop1BMassTag = 0;
  nHadTopTag       = 0;
  nHadTopTag2      = 0;
  nHadTopTag3      = 0;
  nHadTopTag4      = 0;
  nHadTop0BMassTag = 0;
  nHadTop0BAntiTag = 0;
  AK8_Ht   = 0;
  minDeltaR_W_b = 9999;
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" start"<<std::endl;
    TLorentzVector AK8_v4; AK8_v4.SetPtEtaPhiE(data.jetsAK8.Pt[i], data.jetsAK8.Eta[i], data.jetsAK8.Phi[i], data.jetsAK8.E[i]);
#if VER == 0
    double sd_mass_w = isData ? data.jetsAK8.softDropMass[i] : softDropMassCorr[i];
    double sd_mass_top = data.jetsAK8.softDropMass[i];
#elif VER == 1
    double sd_mass_w = isData ? data.jetsAK8.softDropMassPuppi[i] : softDropMassCorr[i];
    double sd_mass_top = data.jetsAK8.softDropMassPuppi[i];
#else
    // For W   tagging in MC we use: GEN/RECO corrected +scaled+smeared softdrop mass
    // For top tagging in MC we use: L1L2L3 subjet corrected +scaled+smeared softdrop mass
    double puppi_pt  = data.jetsAK8.PtPuppi[i];
    double puppi_eta = data.jetsAK8.EtaPuppi[i];
    double corr = puppisd_corrGEN_->Eval(puppi_pt);
    if(std::abs(puppi_eta)<=1.3) corr *= puppisd_corrRECO_cen_->Eval(puppi_pt);
    else corr *= puppisd_corrRECO_for_->Eval(puppi_pt);
    double sd_mass_w_uncorr = data.jetsAK8.uncorrSDMassPuppi[i];
    double sd_mass_w        = isData ? data.jetsAK8.uncorrSDMassPuppi[i] * corr : softDropMassCorr[i];
    double sd_mass_top      = isData ? data.jetsAK8.corrSDMassPuppi[i]          : softDropMass[i];
    double sd_mass_top1     = isData ? data.jetsAK8.corrSDMassPuppi[i]          : AK8_softDropMass1[i];
    double sd_mass_top2     = isData ? data.jetsAK8.softDropMassPuppi[i]        : AK8_softDropMass2[i];
#endif
    softDropMassW   .push_back(sd_mass_w);
    softDropMassTop .push_back(sd_mass_top);
    softDropMassTop1.push_back(sd_mass_top1);
    softDropMassTop2.push_back(sd_mass_top2);
    if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" mass ok"<<std::endl;

    // Jet ID
    if ( passLooseJetAK8[i] = 
	 ( data.jetsAK8.looseJetID[i] == 1 &&
	   data.jetsAK8.Pt[i]         >= JET_AK8_PT_CUT &&
	   std::abs(data.jetsAK8.Eta[i])  <  JET_AK8_ETA_CUT ) ) {

      if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" id ok"<<std::endl;
      // Photon jet veto (not applied):
      // - Calculate DR and Pt ratio of photons overlapping with the jet
      // - Remove jet if DR<0.8
      for (size_t j=0; j<selected_photons.size(); ++j) {
	double dR = selected_photons[j].DeltaR(AK8_v4);
	if (dR<AK8_photonDR[i]) {
	  AK8_photonDR[i] = dR;
	  AK8_photonPtRatio[i] = selected_photons[j].Pt()/AK8_v4.Pt();
	}
      }
      //if ( passLooseJetAK8[i] = (AK8_photonDR[i]>=0.8) ) {
      iJetAK8.push_back(i);
      itJetAK8[i] = nJetAK8++;

      // Lepton-jet overlap
      for (size_t j=0; j<veto_electrons.size(); ++j) {
	double dR = veto_electrons[j].DeltaR(AK8_v4);
	if (dR<AK8_eleDR[i]) {
	  AK8_eleDR[i] = dR;
	  AK8_elePtRatio[i] = veto_electrons[j].Pt()/AK8_v4.Pt();
	}
      }
      for (size_t j=0; j<veto_muons.size(); ++j) {
	double dR = veto_muons[j].DeltaR(AK8_v4);
	if (dR<AK8_muonDR[i]) {
	  AK8_muonDR[i] = dR;
	  AK8_muonPtRatio[i] = veto_muons[j].Pt()/AK8_v4.Pt();
	}
      }

      // Tagging Variables
      double pt      = data.jetsAK8.Pt[i];
      double abseta  = data.jetsAK8.Eta[i];
      double tau_21 = tau21[i];
      double tau_32 = tau32[i];
      if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" overlap ok"<<std::endl;

      // _______________________________________________________
      //                   Hadronic W Tag definition

      if (passWMassTag[i] = 
	  ( pt        >= W_PT_CUT &&
	    abseta    <  W_ETA_CUT &&
	    sd_mass_w >= W_SD_MASS_CUT_LOW && 
	    sd_mass_w <  W_SD_MASS_CUT_HIGH) ) {
	iWMassTag.push_back(i);
	itWMassTag[i] = nWMassTag++;
	// Loose/Tight W Tag Working points
	if (passLooseWTag[i] = (tau_21 < W_TAU21_LOOSE_CUT) ) {
	  iLooseWTag.push_back(i);
	  itLooseWTag[i] = nLooseWTag++;
	}
	if (tau_21 < 0.35) nHadWTag1++;
	if (tau_21 < 0.55) nHadWTag3++;
	if (passTightWTag[i] = (tau_21 < W_TAU21_TIGHT_CUT) ) {
	  iTightWTag.push_back(i);
	  itTightWTag[i] = nTightWTag++;
	  // DR between W and b
	  while(data.jetsAK4.Loop()) {
	    size_t j = data.jetsAK4.it;
	    TLorentzVector AK4_v4; AK4_v4.SetPtEtaPhiE(data.jetsAK4.Pt[j], data.jetsAK4.Eta[j], data.jetsAK4.Phi[j], data.jetsAK4.E[j]);
	    if (passMediumBTag[j]) {
	      double dR = AK4_v4.DeltaR(AK8_v4);
	      if (dR<minDeltaR_W_b) minDeltaR_W_b = dR;
	    }
	  }
	} else {
	  passTightWAntiTag[i] = 1;
	  iTightWAntiTag.push_back(i);
	  itTightWAntiTag[i] = nTightWAntiTag++;
	}
      }
      // Same without mass smearing
      if (passOtherWMassTag[i] = 
          ( pt        >= W_PT_CUT &&
            abseta    <  W_ETA_CUT &&
            sd_mass_w_uncorr >= W_SD_MASS_CUT_LOW && 
            sd_mass_w_uncorr <  W_SD_MASS_CUT_HIGH) ) {
        passOtherTightWTag[i] = (tau_21 < W_TAU21_TIGHT_CUT);
        passOtherTightWAntiTag[i] = !passOtherTightWTag[i];
      }
      if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" wtag ok"<<std::endl;

      // _______________________________________________________
      //                  Hadronic Top Tag definition

      minDeltaR_W_b = 9999;
      if (passHadTopMassTag[i] = 
	  ( pt          >= TOP_PT_CUT && 
	    sd_mass_top >= TOP_SD_MASS_CUT_LOW &&
	    sd_mass_top <  TOP_SD_MASS_CUT_HIGH) ) {
	itHadTopMassTag[i] = nHadTopMassTag++;
	iHadTopMassTag.push_back(i);
#if USE_BTAG == 1
	if (passHadTop1BMassTag[i] = passSubjetBTag[i]) {
#endif
	  itHadTop1BMassTag[i] = nHadTop1BMassTag++;
	  iHadTop1BMassTag.push_back(i);
	  if (passHadTopTag[i] = (tau_32 < TOP_TAU32_CUT) ) {
	    itHadTopTag[i] = nHadTopTag++;
	    iHadTopTag.push_back(i);
	  }
	  if (tau_32 < 0.54) nHadTopTag2++;
	  if (tau_32 < 0.65) nHadTopTag3++;
	  if (tau_32 < 0.80) nHadTopTag4++;
#if USE_BTAG == 1
	} else {
	  while(data.jetsAK4.Loop()) {
	    size_t j = data.jetsAK4.it;
	    TLorentzVector AK4_v4; AK4_v4.SetPtEtaPhiE(data.jetsAK4.Pt[j], data.jetsAK4.Eta[j], data.jetsAK4.Phi[j], data.jetsAK4.E[j]);
	    if (passMediumBTag[j]) {
	      double dR = AK4_v4.DeltaR(AK8_v4);
	      if (dR<minDeltaR_W_b) minDeltaR_W_b = dR;
	    }
	  }
	  if(minDeltaR_W_b > 0.8) {
	    passHadTop0BMassTag[i] = 1; 
	    itHadTop0BMassTag[i] = nHadTop0BMassTag++;
	    iHadTop0BMassTag.push_back(i);
	  }
	  if (passHadTop0BAntiTag[i] = (tau_32 >= TOP_TAU32_CUT) ) {
	    itHadTop0BAntiTag[i] = nHadTop0BAntiTag++;
	    iHadTop0BAntiTag.push_back(i);
	  }
	}
#endif
      }
      // Other top mass correction
      minDeltaR_W_b = 9999;
      if (passOtherHadTopMassTag[i] = 
          ( pt           >= TOP_PT_CUT && 
            sd_mass_top2 >= TOP_SD_MASS_CUT_LOW &&
            sd_mass_top2 <  TOP_SD_MASS_CUT_HIGH) ) {
        if (passSubjetBTag[i]) {
          passOtherHadTopTag[i] = (tau_32 < TOP_TAU32_CUT);
        } else {
          while(data.jetsAK4.Loop()) {
            size_t j = data.jetsAK4.it;
            TLorentzVector AK4_v4; AK4_v4.SetPtEtaPhiE(data.jetsAK4.Pt[j], data.jetsAK4.Eta[j], data.jetsAK4.Phi[j], data.jetsAK4.E[j]);
            if (passMediumBTag[j]) {
              double dR = AK4_v4.DeltaR(AK8_v4);
              if (dR<minDeltaR_W_b) minDeltaR_W_b = dR;
            }
          }
          if(minDeltaR_W_b > 0.8) {
            passOtherHadTop0BMassTag[i] = minDeltaR_W_b > 0.8;
            passOtherHadTop0BAntiTag[i] = (tau_32 >= TOP_TAU32_CUT);
          }
        }
      }
      if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" top tag ok"<<std::endl;

      //} // End photon jet veto

    } // End Jet Selection

    // Online jet selection for AK8 HT
    if ( data.jetsAK8.Pt[i]         > 150 &&
	 std::abs(data.jetsAK8.Eta[i])  <  2.5 ) {
      // Ht
      AK8_Ht += data.jetsAK8.Pt[i];
    }

    if (debug>1) std::cout<<"AnalysisBase::calc_common_var: AK8 "<<i<<" end"<<std::endl;
  } // End AK8 Jet Loop
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end AK8"<<std::endl;
  
  // Loop on generator particles
  iGenHadW   .clear();
  iGenTop    .clear();
  iGenZ    .clear();
  iGenWtagZ  .clear();
  iGenMassW  .clear();
  iGenMassTop.clear();
  itGenHadW              .assign(data.gen.size, (size_t)-1);
  itGenTop               .assign(data.gen.size, (size_t)-1);
  itGenZ                 .assign(data.gen.size, (size_t)-1);
  itGenWtagZ             .assign(data.gen.size, (size_t)-1);
  itGenMassW             .assign(data.gen.size, (size_t)-1);
  itGenMassTop           .assign(data.gen.size, (size_t)-1);
  passGenHadW            .assign(data.gen.size, 0);
  passGenTop             .assign(data.gen.size, 0);
  passGenLepton          .assign(data.gen.size, 0);
  genLepPassLepID        .assign(data.gen.size, 0);
  genHadWHasJetMatch     .assign(data.gen.size, 0);
  genHadWPassWTag        .assign(data.gen.size, 0);
  genTopHasJetMatch      .assign(data.gen.size, 0);
  genTopPassTopTag       .assign(data.gen.size, 0);
  genHadWPassWMassTag    .assign(data.gen.size, 0);
  genTopPassTopMassTag   .assign(data.gen.size, 0);
  genHadWMatchedAK8JetPt .assign(data.gen.size, -9999);
  genHadWMatchedAK8JetEta.assign(data.gen.size, -9999);
  genTopMatchedAK8JetPt  .assign(data.gen.size, -9999);
  genTopMatchedAK8JetEta .assign(data.gen.size, -9999);
  passPhotonPrompt.assign(data.gen.size, 0);
  nGenHadW = nGenTop = 0;
  nGenMassW = nGenMassTop = 0;
  nGenZ = nGenWtagZ = 0;
  npreWTag = nWTag = nmWTag = npreTopTag = nTopTag = nmTopTag = 0;
  nFakePhoton = nDirectPromptPhoton = nFragmentationPromptPhoton = 0;
  float dR;
  float dR1;
  //bool passWTag = 0;
  //bool passpreWTag = 0;
  TLorentzVector photag;
  TLorentzVector genpho;
  TLorentzVector ztag_v4;
  TLorentzVector genz_v4;
  TLorentzVector wtag_v4;
  TLorentzVector genw_v4;
  TLorentzVector genb_v4; 
  std::vector<TLorentzVector> selected_genphoton;
  std::vector<TLorentzVector> selected_genw_v4;
  std::vector<TLorentzVector> selected_genb_v4;
  //bool passTopTag = 0;
  //bool passpreTopTag = 0;
  TLorentzVector gentop_v4;
  TLorentzVector jet_v4;
  hasGenW            .assign(data.jetsAK8.size, 0);
  hasGenTop          .assign(data.jetsAK8.size, 0);
  while(data.gen.Loop()) {
    size_t i = data.gen.it;
    TLorentzVector gen_v4; gen_v4.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);

    while(data.pho.Loop()){
      size_t j = data.pho.it;
      if( passPhotonSelect[j]){
	if( std::abs(data.gen.ID[i]) == 22 && (std::abs(data.gen.Mom0ID[i])==1||std::abs(data.gen.Mom0ID[i])==2||std::abs(data.gen.Mom0ID[i])==3||std::abs(data.gen.Mom0ID[i])==4||std::abs(data.gen.Mom0ID[i])==5||std::abs(data.gen.Mom0ID[i])==6||std::abs(data.gen.Mom0ID[i])==21)){
	  photag.SetPtEtaPhiE(data.pho.Pt[j], data.pho.Eta[j], data.pho.Phi[j], data.pho.E[j]);
	  genpho.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);
	  if(genpho.DeltaR(photag) < 0.1 && data.gen.Pt[i]/data.pho.Pt[j] > 0.5 && data.gen.Pt[i]/data.pho.Pt[j] < 2 && data.gen.Status[i] == 1) selected_genphoton.push_back(photag);
	}
	else nFakePhoton++;
      }
    }

    // Select only final version of the particles (their daughters have different IDs)
    if (data.gen.Dau0ID[i]!=data.gen.ID[i]&&data.gen.Dau1ID[i]!=data.gen.ID[i]) {
      
      // Generator leptons with W mother
      if (passGenLepton[i] = 
	  ( std::abs(data.gen.ID[i])==11 ||
	    std::abs(data.gen.ID[i])==13 ||
	    std::abs(data.gen.ID[i])==15 ) &&
	  std::abs(data.gen.Mom0ID[i])==24) {
	
	// Match to veto electrons
        if (abs(data.gen.ID[i])==11) {
          bool match = 0;
          while(data.ele.Loop()) {
            size_t j = data.ele.it;
            if (passEleVeto[j]) {
              TLorentzVector ele_v4; ele_v4.SetPtEtaPhiE(data.ele.Pt[j], data.ele.Eta[j], data.ele.Phi[j], data.ele.E[j]);
              if (gen_v4.DeltaR(ele_v4)<0.1) {
                match = 1;
                //std::cout<<"gen/ele - pt = "<<gen_v4.Pt()<<"/"<<ele_v4.Pt()<<" - eta = "<<gen_v4.Eta()<<"/"<<ele_v4.Eta()<<"- DR = "<<gen_v4.DeltaR(ele_v4)<<std::endl;
              }
            }
          }
          if (match) genLepPassLepID[i] = 1;
        }
        // Match to veto muons
        if (abs(data.gen.ID[i])==13) {
          bool match = 0;
          while(data.mu.Loop()) {
            size_t j = data.mu.it;
            if (passMuVeto[j]) {
              TLorentzVector mu_v4; mu_v4.SetPtEtaPhiE(data.mu.Pt[j], data.mu.Eta[j], data.mu.Phi[j], data.mu.E[j]);
              if (gen_v4.DeltaR(mu_v4)<0.1) {
                match = 1;
                //std::cout<<"gen/mu - pt = "<<gen_v4.Pt()<<"/"<<mu_v4.Pt()<<" - eta = "<<gen_v4.Eta()<<"/"<<mu_v4.Eta()<<"- DR = "<<gen_v4.DeltaR(mu_v4)<<std::endl;
              }
            }
          }
          if (match) genLepPassLepID[i] = 1;
        }
        // Match to veto taus
        if (abs(data.gen.ID[i])==15) {
          bool match = 0;
          while(data.tau.Loop()) {
            size_t j = data.tau.it;
            if (passTauVeto[j]) {
              TLorentzVector tau_v4; tau_v4.SetPtEtaPhiE(data.tau.Pt[j], data.tau.Eta[j], data.tau.Phi[j], data.tau.E[j]);
              if (gen_v4.DeltaR(tau_v4)<0.1) {
                match = 1;
                //std::cout<<"gen/tau - pt = "<<gen_v4.Pt()<<"/"<<tau_v4.Pt()<<" - eta = "<<gen_v4.Eta()<<"/"<<tau_v4.Eta()<<"- DR = "<<gen_v4.DeltaR(tau_v4)<<std::endl;
              }
            }
          }
          if (match) genLepPassLepID[i] = 1;
        }
      }
      
      // Apply cut |eta| < 2.4
      if (std::abs(data.gen.Eta[i])<2.4) {
	// gen bs
	if(abs(data.gen.ID[i])==5&&data.gen.Pt[i]>0) {
	  genb_v4.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);
	  selected_genb_v4.push_back(genb_v4);
	}
  
	// gen Zs
	if ( abs(data.gen.ID[i])==23 && data.gen.Pt[i] > 0) {
	  genz_v4.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);
	  while(data.jetsAK8.Loop()) {
	    size_t j = data.jetsAK8.it;
	    ztag_v4.SetPtEtaPhiE(data.jetsAK8.Pt[j], data.jetsAK8.Eta[j], data.jetsAK8.Phi[j], data.jetsAK8.E[j]);
	    dR = genz_v4.DeltaR(ztag_v4);
	    if (dR<0.8) {
	      iGenZ.push_back(j);
	      itGenZ[j] = nGenZ++;
	      if (passTightWTag[j]) {
		iGenWtagZ.push_back(j);
		itGenWtagZ[j] = nGenWtagZ++;
	      }
	    }
	  }
	}
        
	// gen Ws
	// Consider only hadronically decaying Ws
	if ( passGenHadW[i] =
	     ( abs(data.gen.ID[i])==24 &&
	       ! (abs(data.gen.Dau0ID[i])>=11&&abs(data.gen.Dau0ID[i])<=16) ) ) {
	  iGenHadW.push_back(i);
	  itGenHadW[i] = nGenHadW++;
	  //passWTag = 0;
	  //passpreWTag = 0;
	  genw_v4.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);
	  selected_genw_v4.push_back(genw_v4);
	  double minDR = 9999;
	  while(data.jetsAK8.Loop()) {
	    size_t j = data.jetsAK8.it;
	    wtag_v4.SetPtEtaPhiE(data.jetsAK8.Pt[j], data.jetsAK8.Eta[j], data.jetsAK8.Phi[j], data.jetsAK8.E[j]);
	    dR = genw_v4.DeltaR(wtag_v4);
	    if (dR<0.8) {
	      hasGenW[j] = true;
	      genHadWHasJetMatch[i] = true;
	      //passpreWTag = 1;
	      npreWTag++;
	      if (dR<minDR) {
		minDR = dR;
		genHadWMatchedAK8JetPt[i]  = data.jetsAK8.Pt[j];
		genHadWMatchedAK8JetEta[i] = data.jetsAK8.Eta[j];
		genHadWPassWTag[i] = passTightWTag[j];
		genHadWPassWMassTag[i] = passWMassTag[j];
		if (passTightWTag[j]) {
		  nWTag=1;
		  //iGenHadW.push_back(i);
		  //passWTag = 1;
		  /*while(data.gen.Loop()) {
		    k = data.gen.it;
		    if((abs(data.gen.ID[k])==5)){
		    genb_v4.SetPtEtaPhiE(data.gen.Pt[k], data.gen.Eta[k], data.gen.Phi[k], data.gen.E[k]);
		    dR1 = genb_v4.DeltaR(wtag_v4);
		    if (dR1<0.8) nWTag=-1;
		    }
		    }
		  */
		}
	      }
	    }
	  }
	}
      
	// gen tops
	// Consider also leptonic W decays, because lepton is usually energetic
	if(passGenTop[i] = (abs(data.gen.ID[i])==6)) {
	  iGenTop.push_back(i);
	  itGenTop[i] = nGenTop++;

	  //passTopTag = 0;
	  //passpreTopTag = 0;
	  gentop_v4.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);
	  double minDR = 9999;
	  while(data.jetsAK8.Loop()) {
	    size_t j = data.jetsAK8.it;
	    jet_v4.SetPtEtaPhiE(data.jetsAK8.Pt[j], data.jetsAK8.Eta[j], data.jetsAK8.Phi[j], data.jetsAK8.E[j]);
	    dR = gentop_v4.DeltaR(jet_v4);
	    if (dR<0.8) {
	      hasGenTop[j] = true;
	      genTopHasJetMatch[i] = true;
	      npreTopTag++;
	      if (dR<minDR) {
		minDR = dR;
		genTopMatchedAK8JetPt[i]  = data.jetsAK8.Pt[j];
		genTopMatchedAK8JetEta[i] = data.jetsAK8.Eta[j];
		genTopPassTopTag[i] = passHadTopTag[j];
		genTopPassTopMassTag[i] = passHadTopMassTag[j];
		//if (passHadTopTag[j]) {
		//  iGenTop.push_back(i);
		//  itGenTop[j] = nTopTag++;
		//}
		if (passHadTop1BMassTag[j]) {
		  iGenMassTop.push_back(j);
		  itGenMassTop[j] = nmTopTag++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  for(size_t i=0;i<selected_genw_v4.size();++i){
    for(size_t k=0;k<selected_genb_v4.size();++k){
      while(data.jetsAK8.Loop()) {
	size_t j = data.jetsAK8.it;
	wtag_v4.SetPtEtaPhiE(data.jetsAK8.Pt[j], data.jetsAK8.Eta[j], data.jetsAK8.Phi[j], data.jetsAK8.E[j]);
	//dR = wtag_v4.DeltaR(selected_genw_v4[i]);
	dR = selected_genw_v4[i].DeltaR(wtag_v4);
	if (dR<0.8) {
	  if (passWMassTag[j]) {
	    nmWTag=1;
	    dR1 = selected_genb_v4[k].DeltaR(wtag_v4);
	    if (dR1<0.8) nmWTag=-1;
	    else{iGenMassW.push_back(j); itGenMassW[j] = nGenMassW++;}
	  }
	  if (passTightWTag[j]) {
	    //passWTag = 1;
	    nWTag=1;
	    //dR1 = wtag_v4.DeltaR(selected_genb_v4[k]);
	    dR1 = selected_genb_v4[k].DeltaR(wtag_v4);
	    if (dR1<0.8) nWTag=-1;
	    else{ 
	      //iGenHadW.push_back(j); 
	      //itGenHadW[j] = nGenHadW++;
	    }
	  }
	  //passpreWTag = 1;
	  npreWTag++;
	}
      }
    }
  }
  for(size_t j=0;j<selected_genphoton.size();++j){
    while(data.gen.Loop()){
      size_t i = data.gen.it;
      photag.SetPtEtaPhiE(data.pho.Pt[j], data.pho.Eta[j], data.pho.Phi[j], data.pho.E[j]);
      genpho.SetPtEtaPhiE(data.gen.Pt[i], data.gen.Eta[i], data.gen.Phi[i], data.gen.E[i]);
      dR = genpho.DeltaR(photag);
      if(dR > 0.4 && (std::abs(data.gen.ID[i])==1||std::abs(data.gen.ID[i])==2||std::abs(data.gen.ID[i])==3||std::abs(data.gen.ID[i])==4||std::abs(data.gen.ID[i])==5||std::abs(data.gen.ID[i])==6||std::abs(data.gen.ID[i])==21)) nDirectPromptPhoton++;
      else nFragmentationPromptPhoton++;
    }
  }
  //if (syst_index==0&&applySmearing) {
  //  data.evt.MR  = data.evt.MR_Smear;
  //  data.evt.MTR = data.evt.MTR_Smear;
  //  data.evt.R   = data.evt.MTR/data.evt.MR;
  //  data.evt.R2  = data.evt.R*data.evt.R;    
  //} else if (syst_index!=0) {
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end gen"<<std::endl;

  // Recalculation of Razor variables
  // Has to be done after jet uncertainties applied
  // Get selected AK4 jets (input for megajets)

  // Calculating hemispheres is a CPU intensive task
  // So we only do it when really necessary (no pt variations)
  // And save them if no systematics are applied on jets
  if (syst_index==0 || recalc_megajets) {
    // Get input jets
    std::vector<TLorentzVector> selected_jets_AK4, nophoton_jets_AK4;
    while(data.jetsAK4.Loop()) {
      size_t i = data.jetsAK4.it;
      TLorentzVector jet_v4; jet_v4.SetPtEtaPhiE(data.jetsAK4.Pt[i], data.jetsAK4.Eta[i], data.jetsAK4.Phi[i], data.jetsAK4.E[i]);
      // Pass jet selection criteria
      if (passLooseJet[i]) selected_jets_AK4.push_back(jet_v4);
      if (passLooseJetNoPho[i]&&(nPhotonSelect==1||nPhotonPreSelect==1)) nophoton_jets_AK4.push_back(jet_v4);
    }
    // Calculate hemispheres
    if (selected_jets_AK4.size()>=2) hemis_AK4 = Razor::CombineJets(selected_jets_AK4);
    else hemis_AK4.clear();
    if (nPhotonSelect==1||nPhotonPreSelect==1) {
      if (nophoton_jets_AK4.size()>=2) hemis_AK4_nophoton = Razor::CombineJets(nophoton_jets_AK4);
      else hemis_AK4_nophoton.clear();
    } else hemis_AK4_nophoton = hemis_AK4;
    // Save result
    if (syst_index==0) {
      saved_hemis_AK4 = hemis_AK4;
      saved_hemis_AK4_nophoton = hemis_AK4_nophoton;
    }
  } else {
    // Load previously calculated result
    hemis_AK4 = saved_hemis_AK4;
    hemis_AK4_nophoton = saved_hemis_AK4_nophoton;
  }

  // Recalculate Razor (also with 1/2lep or pho added to MET)
  MET_1l = MET_1vl = MET_ll = MET_pho = -9999;
  data.evt.MR  = MR_pho = -9999;
  data.evt.MTR = MTR_1l = MTR_1vl = MTR_ll = MTR_pho = -9999;
  data.evt.R   = R_1l   = R_1vl   = R_ll   = R_pho   = -9999;
  data.evt.R2  = R2_1l  = R2_1vl  = R2_ll  = R2_pho  = -9999;
  dPhiRazor = dPhiRazorNoPho = 9999;
  if (hemis_AK4.size()==2) {
    // Normal Razor
    TVector3 shifted_met;
    shifted_met.SetPtEtaPhi(data.met.Pt[0], 0, data.met.Phi[0]);
    data.evt.MR  = Razor::CalcMR(hemis_AK4[0], hemis_AK4[1]);
    data.evt.MTR = Razor::CalcMTR(hemis_AK4[0], hemis_AK4[1], shifted_met);
    data.evt.R   = data.evt.MTR/data.evt.MR;
    data.evt.R2  = data.evt.R*data.evt.R;
    dPhiRazor = std::abs(TVector2::Phi_mpi_pi(hemis_AK4[0].Phi() - hemis_AK4[1].Phi()));
    // 1 (selected) lepton added
    if (nLepSelect==1) {
      MET_1l = met_1l.Pt();
      MTR_1l = Razor::CalcMTR(hemis_AK4[0], hemis_AK4[1], met_1l);
      R_1l   = MTR_1l/data.evt.MR;
      R2_1l  = R_1l*R_1l;
    }
    // 1 veto lepton added (default)
    if (nLepVeto==1) {
      MET_1vl = met_1vl.Pt();
      MTR_1vl = Razor::CalcMTR(hemis_AK4[0], hemis_AK4[1], met_1vl);
      R_1vl   = MTR_1vl/data.evt.MR;
      R2_1vl  = R_1vl*R_1vl;
    }
    // 2 leptons added
    if (M_ll!=-9999) {
      MET_ll = met_ll.Pt();
      MTR_ll = Razor::CalcMTR(hemis_AK4[0], hemis_AK4[1], met_ll);
      R_ll   = MTR_ll/data.evt.MR;
      R2_ll  = R_ll*R_ll;
    }
  }
  // Remove photon from both jet collections and add to MET
  if (hemis_AK4_nophoton.size()==2) {
    MET_pho = met_pho.Pt();
    MR_pho  = Razor::CalcMR(hemis_AK4_nophoton[0], hemis_AK4_nophoton[1]);
    MTR_pho = Razor::CalcMTR(hemis_AK4_nophoton[0], hemis_AK4_nophoton[1], met_pho);
    R_pho   = MTR_pho/MR_pho;
    R2_pho  = R_pho*R_pho;
    dPhiRazorNoPho = std::abs(TVector2::Phi_mpi_pi(hemis_AK4_nophoton[0].Phi() - hemis_AK4_nophoton[1].Phi()));
  }
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end razor"<<std::endl;
  if (debug) std::cout<<"AnalysisBase::calc_common_var: end"<<std::endl;
}

//_______________________________________________________
//  Apply analysis cuts in the specified search region

bool
AnalysisBase::apply_all_cuts(char region) {
  return apply_ncut(region, analysis_cuts[region].size());
}

bool
AnalysisBase::apply_ncut(char region, unsigned int ncut) {
  if (ncut>analysis_cuts[region].size()) return 0;
  for (unsigned int i=0; i<ncut; ++i) if ( ! analysis_cuts[region][i].func() ) return 0;
  return 1;
}

// Cuts to apply/exclude by cut name
bool
AnalysisBase::apply_cut(char region, std::string cut_name) {
  for (const auto& cut : analysis_cuts[region]) if (cut_name == cut.name) return cut.func();
  return 0;
}

bool
AnalysisBase::apply_cuts(char region, std::vector<std::string> cuts) {
  for (const auto& cut_in_region : analysis_cuts[region]) for (const auto& cut : cuts) 
    if (cut == cut_in_region.name) if (!cut_in_region.func()) return 0;
  return 1;
}

bool
AnalysisBase::apply_all_cuts_except(char region, std::string cut_to_skip) {
  bool result = true, found = false;
  for (const auto& cut : analysis_cuts[region]) {
    if (cut.name == cut_to_skip) { 
      found = true;
      continue;
    }
    if (!cut.func()) result = false;
  }
  // If a certain cut meant to be skipped (N-1) is not found for some reason
  // eg. mistyped, then end the job with ar error
  // This is for safety: We do not want to fill histograms wrongly by mistake
  if (!found) {
    std::cout<<"No cut to be skipped exsists in search region \""<<region<<"\" with name: \""<<cut_to_skip<<"\""<<std::endl;
    utils::error("AnalysisBase - the second argument for apply_all_cuts_except() is a non-sensical cut");
  }
  return result;
}

bool
AnalysisBase::apply_all_cuts_except(char region, std::vector<std::string> cuts_to_skip) {
  bool result = true;
  unsigned int found = 0;
  for (const auto& cut : analysis_cuts[region]) {
    for (const auto& cut_to_skip : cuts_to_skip) if (cut.name==cut_to_skip) { 
      ++found;
      continue;
    }
    if (!cut.func()) result = false;
  }
  // If a certain cut meant to be skipped is not found for some reason
  // eg. mistyped, then end the job with ar error
  // This is for safety: We do not want to fill histograms wrongly by mistake
  if (found!=cuts_to_skip.size()) {
    std::cout<<"A cut to be skipped does not exsist in seaerch region \""<<region<<"\" with names: ";
    for (const auto& cut : cuts_to_skip) std::cout<<cut<<", "; std::cout<<std::endl;
    utils::error("AnalysisBase - the second argument for apply_all_cuts_except() contains at least one non-sensical cut");
  }
  return result;
}


// Same functions but with cut index which is faster (can use an enum, to make it nicer)
bool
AnalysisBase::apply_cut(char region, unsigned int cut_index) { return analysis_cuts[region][cut_index].func(); }

bool
AnalysisBase::apply_cuts(char region, std::vector<unsigned int> cuts) {
  for (const unsigned int& cut : cuts) if ( ! analysis_cuts[region][cut].func() ) return 0;
  return 1;
}

bool
AnalysisBase::apply_all_cuts_except(char region, unsigned int cut_to_skip) {
  if (cut_to_skip>=analysis_cuts[region].size()) {
    std::cout<<"Index ("<<cut_to_skip<<") is too high for the cut to be skipped in search region '"<<region<<"'"<<std::endl;
    utils::error("AnalysisBase::apply_all_cuts_except(char region, unsigned int cut_to_skip)");
  }
  for (unsigned int i=0, n=analysis_cuts[region].size(); i<n; ++i) {
    if (i==cut_to_skip) continue;
    if ( ! analysis_cuts[region][i].func() ) return 0;
  }
  return 1;
}

bool
AnalysisBase::apply_all_cuts_except(char region, std::vector<unsigned int> cuts_to_skip) {
  for (unsigned int i=0, n=analysis_cuts[region].size(); i<n; ++i) {
    for (const unsigned int& cut_to_skip : cuts_to_skip) if (i!=cut_to_skip) 
      if ( ! analysis_cuts[region][i].func() ) return 0;
  }
  return 1;
}

//_______________________________________________________
//                 List of Histograms

TH1D* h_totweight;
TH1D* h_totweight_toppt;
TH1D* h_totweight_pileup;
std::vector<TH2D*> vh_totweight_signal;
std::vector<TH2D*> vh_xsec_signal;
std::vector<TH2D*> vh_weightnorm_signal;
TH1D* h_nisrjets;
TH1D* h_totweight_isr;
std::vector<TH3D*> vh_nisrjets_signal;
std::vector<TH3D*> vh_totweight_signal_isr;
TH1D* h_npvLowHigh_data;
std::vector<TH3D*> vh_npvLowHigh_signal;
TH1D* h_pileup_data;
TH1D* h_pileup_data_down;
TH1D* h_pileup_data_up;
TH1D* h_pileup_mc;
TH1D* h_pileup_weight;
TH1D* h_pileup_weight_down;
TH1D* h_pileup_weight_up;
TH1D* h_nvtx;
TH1D* h_nvtx_rw;
TH1D* h_read_speed_1k;
TH1D* h_read_speed_10k;
TH1D* h_read_speed_job;
TH2D* h_read_speed_vs_nevt_10k;
TH2D* h_read_speed_vs_nevt_job;
TH1D* h_runtime_job;
TH2D* h_runtime_vs_nevt_10k;
TH2D* h_runtime_vs_nevt_job;
TH2D* h_btag_eff_b_loose;
TH2D* h_btag_eff_c_loose;
TH2D* h_btag_eff_l_loose;
TH2D* h_btag_eff_b_medium;
TH2D* h_btag_eff_c_medium;
TH2D* h_btag_eff_l_medium;

TH1D* h_trigger_pass;
TH1D* h_trigger_total;
TH2D* h_trigger2d_pass;
TH2D* h_trigger2d_total;
TH2D* h_trigger2d_nolep_pass;
TH2D* h_trigger2d_nolep_total;

std::vector<TH1D*> vh_MRR2_data;
std::vector<TH1D*> vh_MRR2_data_nj45;
std::vector<TH1D*> vh_MRR2_data_nj6;
std::vector<std::vector<TH1D*> > vvh_MRR2_bkg;
std::vector<std::vector<TH1D*> > vvh_MRR2_bkg_nj45;
std::vector<std::vector<TH1D*> > vvh_MRR2_bkg_nj6;
std::map<uint32_t, std::vector<TH1D*> > m_vh_MRR2_sig;
std::map<uint32_t, std::vector<TH1D*> > m_vh_MRR2_sig_nj45;
std::map<uint32_t, std::vector<TH1D*> > m_vh_MRR2_sig_nj6;

// Histos for the Z(nunu) estimate
TH2D* h_MR_R2_G_EB; // For counts in data
TH2D* h_MR_R2_G_EE;
TH2D* h_MR_R2_G_EB_nj45;
TH2D* h_MR_R2_G_EB_nj6;
TH2D* h_MR_R2_G_EE_nj45;
TH2D* h_MR_R2_G_EE_nj6;
TH3D* h_CHIsoTemplate_Fake_g_EB; // Fake photon templates (data)
TH3D* h_CHIsoTemplate_Fake_g_EE;
TH3D* h_CHIsoTemplate_Fake_g_EB_MC; // Fake photon templates (mc)
TH3D* h_CHIsoTemplate_Fake_g_EE_MC;
TH3D* h_CHIsoTemplate_Prompt_g_EB; // Prompt photon template (GJets)
TH3D* h_CHIsoTemplate_Prompt_g_EE; 
TH3D* h_MR_R2_CHIso_GNoIso_EB; // For purity in data
TH3D* h_MR_R2_CHIso_GNoIso_EE;
TH3D* h_MR_R2_CHIso_GNoIso_EB_nj45;
TH3D* h_MR_R2_CHIso_GNoIso_EB_nj6;
TH3D* h_MR_R2_CHIso_GNoIso_EE_nj45;
TH3D* h_MR_R2_CHIso_GNoIso_EE_nj6;
TH3D* h_MR_R2_CHIso_gNoIso_EB; // For purity in data (G-1)
TH3D* h_MR_R2_CHIso_gNoIso_EE;
TH3D* h_MR_R2_CHIso_gNoIso_EB_nj45;
TH3D* h_MR_R2_CHIso_gNoIso_EB_nj6;
TH3D* h_MR_R2_CHIso_gNoIso_EE_nj45;
TH3D* h_MR_R2_CHIso_gNoIso_EE_nj6;
TH3D* h_MR_R2_IsPrompt_G_EB; // MC truth for photon purity
TH3D* h_MR_R2_IsPrompt_G_EE;
TH3D* h_MR_R2_IsDirect_G_EB; // MC truth for direct fraction
TH3D* h_MR_R2_IsDirect_G_EE;
TH2D* h_MR_R2_S; // Counts for Transfer factor: S(mc)/G(mc)
TH2D* h_MR_R2_S_nj45;
TH2D* h_MR_R2_S_nj6;
TH2D* h_MR_R2_G_DirectPrompt;
TH2D* h_MR_R2_G_DirectPrompt_nj45;
TH2D* h_MR_R2_G_DirectPrompt_nj6;
TH2D* h_MR_R2_Z; // Double ratio: [G(data)/G(mc)] / [Z(ll,data)/Z(ll,mc)]
TH2D* h_MR_R2_Z_nj45;
TH2D* h_MR_R2_Z_nj6;
TH2D* h_MR_R2_G;
TH2D* h_MR_R2_G_nj45;
TH2D* h_MR_R2_G_nj6;
// 1 Lepton estimate
TH2D* h_MR_R2_L;
TH2D* h_MR_R2_L_nj45;
TH2D* h_MR_R2_L_nj6;
// Q' closure
TH2D* h_MR_R2_q;
TH2D* h_MR_R2_q_nj45;
TH2D* h_MR_R2_q_nj6;

//_______________________________________________________
//              Define Histograms here
void
AnalysisBase::init_common_histos(const bool& varySystematics)
{
  // bins
  Double_t gluinoBins[202]; for (int i=0; i<202; ++i) gluinoBins[i] = (i-0.5)*25; 
  Double_t stopBins[402];   for (int i=0; i<402; ++i) stopBins[i] = (i-0.5)*5; 
  Double_t npvLowHighBins[3] = { 0,20,100 };
  Double_t isrWeightBins[3] = {0,1,2};
  Double_t isrJetBins[17]; for (int i=0; i<17; ++i) isrJetBins[i] = i-0.5; 
  // total weight
  h_totweight                     = new TH1D("totweight",           "MC;;Total (generator) event weight", 1,0,1);
  h_totweight_toppt               = new TH1D("totweight_toppt",     "MC;;Total toppt weight",             2,0,2);
  h_totweight_pileup              = new TH1D("totweight_pileup",    "MC;;Total pileup weight",            2,0,2);
  vh_totweight_signal    .push_back(new TH2D("totweight_T1tttt",    "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",        201,gluinoBins, 201,gluinoBins));
  vh_totweight_signal    .push_back(new TH2D("totweight_T2tt",      "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",        401,stopBins, 401,stopBins));
  // ISR reweighting
  h_nisrjets                      = new TH1D("nisrjets",            ";N_{ISR jets}", 16,isrJetBins);
  h_totweight_isr                 = new TH1D("totweight_isr",       "MC;;Total (generator) event weight", 2,isrWeightBins);
  vh_nisrjets_signal     .push_back(new TH3D("nisrjets_T1tttt",     ";M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight;N_{ISR jets}", 201,gluinoBins, 201,gluinoBins, 16,isrJetBins));
  vh_nisrjets_signal     .push_back(new TH3D("nisrjets_T2tt",       ";M_{#tilde{t}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight;N_{ISR jets}", 401,stopBins,  401,stopBins,    16,isrJetBins));
  vh_totweight_signal_isr.push_back(new TH3D("totweight_isr_T1tttt","T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",        201,gluinoBins, 201,gluinoBins, 2,isrWeightBins));
  vh_totweight_signal_isr.push_back(new TH3D("totweight_isr_T2tt",  "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",                              401,stopBins,   401,stopBins,   2,isrWeightBins));
  // signal weight
  vh_xsec_signal         .push_back(new TH2D("xsec_T1tttt",         "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Cross-section (pb)",  201,gluinoBins, 201,gluinoBins));
  vh_weightnorm_signal   .push_back(new TH2D("weightnorm_T1tttt",   "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);weight norm. factor", 201,gluinoBins, 201,gluinoBins));
  vh_xsec_signal         .push_back(new TH2D("xsec_T2tt",           "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Cross-section (pb)",  401,stopBins, 401,stopBins));
  vh_weightnorm_signal   .push_back(new TH2D("weightnorm_T2tt",     "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);weight norm. factor", 401,stopBins, 401,stopBins));
  // npv for extrapolations
  h_npvLowHigh_data               = new TH1D("npvLowHigh_data",     "Number of vertices - Data;N_{Vertices}", 2,npvLowHighBins);
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_T1tttt",   "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);N_{vertex}", 201,gluinoBins, 201,gluinoBins, 2,npvLowHighBins));
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_T2tt",     "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);N_{vertex}",                       401,stopBins, 401,stopBins,     2,npvLowHighBins));
  // pileup
  h_pileup_data                = new TH1D("pileup_data",        "Pile-up distribution - Data (Nominal);Pile-up", 100,0,100);
  h_pileup_data_down           = new TH1D("pileup_data_down",   "Pile-up distribution - Data (down);Pile-up",    100,0,100);
  h_pileup_data_up             = new TH1D("pileup_data_up",     "Pile-up distribution - Data (up);Pile-up",      100,0,100);
  h_pileup_mc                  = new TH1D("pileup_mc",          "Pile-up distribution - MC;Pile-up",             100,0,100);
  h_pileup_weight              = new TH1D("pileup_weight",      "Pile-up weights - Nominal MB X-sec (69 mb);Pile-up;Weight",    100,0,100);
  h_pileup_weight_down         = new TH1D("pileup_weight_down", "Pile-up weights - MB X-sec up 5% (72.45 mb);Pile-up;Weight",   100,0,100);
  h_pileup_weight_up           = new TH1D("pileup_weight_up",   "Pile-up weights - MB X-sec down 5% (65.55 mb);Pile-up;Weight", 100,0,100);
  h_nvtx                       = new TH1D("nvtx",               "Number of vertices - Nominal;N_{Vertices}",                      101,-0.5,100.5);
  h_nvtx_rw                    = new TH1D("nvtx_rw",            "Number of vertices - Pile-up reweighted (MC only);N_{Vertices}", 101,-0.5,100.5);
  // job_monitoring histos
  h_read_speed_1k              = new TH1D("read_speed_1k",          ";Read speed (Events/s);Measurement/1k Event",  1000,0,10000);
  h_read_speed_10k             = new TH1D("read_speed_10k",         ";Read speed (Events/s);Measurement/10k Event", 1000,0,10000);
  h_read_speed_job             = new TH1D("read_speed_job",         ";Read speed (Events/s);Measurement/Job",       1000,0,10000);
  h_read_speed_vs_nevt_10k     = new TH2D("read_speed_vs_nevt_10k", ";Entry;Read speed (Events/s)/10k Event",       100,0,10000000, 200,0,10000);
  h_read_speed_vs_nevt_job     = new TH2D("read_speed_vs_nevt_job", ";Total Entries;Read speed (Events/s)/Job",     100,0,10000000, 200,0,10000);
  h_runtime_job                = new TH1D("runtime_job",            ";Total job run-time (min)",                    600,0,600);
  h_runtime_vs_nevt_job        = new TH2D("runtime_vs_nevt_job",    ";Total Entries;Total job run-time (min)",      100,0,10000000, 600,0,600);

  // btagging efficiency
  double ptbins[11]  = { 20,30,50,70,100,140,200,300,600,1000,4000 };
  double effbins[3] = { -0.5,0.5,1.5 };
  h_btag_eff_b_loose            = new TH2D("btag_eff_b_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_c_loose            = new TH2D("btag_eff_c_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_l_loose            = new TH2D("btag_eff_l_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_b_medium           = new TH2D("btag_eff_b_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_c_medium           = new TH2D("btag_eff_c_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_l_medium           = new TH2D("btag_eff_l_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);

  // trigger efficiency
  double htbins[19]  = { 0, 200, 300, 400, 500, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1200, 1500, 2000, 4000, 10000 };
  double HTB[12] = {400, 500, 600, 700, 750, 800, 850, 900, 950, 1000, 1500, 10000};
  double PtB[9]  = {200, 300, 400, 450, 500, 550, 600, 1000, 10000};
  h_trigger_pass                = new TH1D("trigger_pass",    "Pass trigger;H_{T} (GeV)", 18,htbins);
  h_trigger_total               = new TH1D("trigger_total",          "Total;H_{T} (GeV)", 18,htbins);
  h_trigger2d_pass              = new TH2D("trigger2d_pass",  "Pass trigger;H_{T} (GeV);Leading AK8 jet p_{T} (GeV)", 11,HTB, 8,PtB);
  h_trigger2d_total             = new TH2D("trigger2d_total",        "Total;H_{T} (GeV);Leading AK8 jet p_{T} (GeV)", 11,HTB, 8,PtB);
  h_trigger2d_nolep_pass        = new TH2D("trigger2d_nolep_pass",  "Pass trigger;H_{T} (GeV);Leading AK8 jet p_{T} (GeV)", 11,HTB, 8,PtB);
  h_trigger2d_nolep_total       = new TH2D("trigger2d_nolep_total",        "Total;H_{T} (GeV);Leading AK8 jet p_{T} (GeV)", 11,HTB, 8,PtB);

  std::vector<std::string> regions = {"S", "s", "T", "W", "Q", "q", "Z", "L", "G", "G_DirectPrompt", "S_LooseWP"};

  if (varySystematics) {
    for (size_t i=0; i<regions.size(); ++i) {
      if (isData) {
	// Data
	vh_MRR2_data     .push_back(new TH1D((std::string("MRR2_")+regions[i]+"_data").c_str(),      ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	vh_MRR2_data_nj45.push_back(new TH1D((std::string("MRR2_")+regions[i]+"_data_nj45").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	vh_MRR2_data_nj6 .push_back(new TH1D((std::string("MRR2_")+regions[i]+"_data_nj6").c_str(),  ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
      } else if (!isSignal) {
	// Background
	vvh_MRR2_bkg     .push_back(std::vector<TH1D*>());
	vvh_MRR2_bkg_nj45.push_back(std::vector<TH1D*>());
	vvh_MRR2_bkg_nj6 .push_back(std::vector<TH1D*>());
	vvh_MRR2_bkg[i]     .push_back(new TH1D((std::string("MRR2_")+regions[i]+"_bkg").c_str(),      ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	vvh_MRR2_bkg_nj45[i].push_back(new TH1D((std::string("MRR2_")+regions[i]+"_bkg_nj45").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	vvh_MRR2_bkg_nj6[i] .push_back(new TH1D((std::string("MRR2_")+regions[i]+"_bkg_nj6").c_str(),  ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	for (size_t j=0; j<syst.size(); ++j) {
	  std::stringstream ss;
	  ss<<"MRR2_"<<regions[i]<<"_bkg_"<<syst[j];
	  vvh_MRR2_bkg[i].push_back(new TH1D((ss.str()+"Up").c_str(),   ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	  vvh_MRR2_bkg[i].push_back(new TH1D((ss.str()+"Down").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	  std::stringstream ss45;
	  ss45<<"MRR2_"<<regions[i]<<"_bkg_nj45_"<<syst[j];
	  vvh_MRR2_bkg_nj45[i].push_back(new TH1D((ss45.str()+"Up").c_str(),   ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	  vvh_MRR2_bkg_nj45[i].push_back(new TH1D((ss45.str()+"Down").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	  std::stringstream ss6;
	  ss6<<"MRR2_"<<regions[i]<<"_bkg_nj6_"<<syst[j];
	  vvh_MRR2_bkg_nj6[i].push_back(new TH1D((ss6.str()+"Up").c_str(),   ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	  vvh_MRR2_bkg_nj6[i].push_back(new TH1D((ss6.str()+"Down").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	}
      }
    }
  }
  // Signals --> Declare them later after the signal weight calculation
  // in calc_weightnorm_histo_from_ntuple()
  
  double mrbins[6]  = { 800, 1000, 1200, 1600, 2000, 4000 };
  double r2bins[6]  = { 0.08, 0.12, 0.16, 0.24, 0.4, 1.5 }; 
  double isdirectbins[3] = {-0.5,0.5,1.5};
  double chisobins[21] = { 0 };
  for (int i=0; i<21; ++i) chisobins[i] = float(i)/2.;
  
  // Histos for the Z(nunu) estimate
  h_MR_R2_G_EB                 = new TH2D("MR_R2_G_EB",                 "G region, EB;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_G_EE                 = new TH2D("MR_R2_G_EE",                 "G region, EE;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_G_EB_nj45            = new TH2D("MR_R2_G_EB_nj45",            "G region, EB, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G_EB_nj6             = new TH2D("MR_R2_G_EB_nj6",             "G region, EB, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G_EE_nj45            = new TH2D("MR_R2_G_EE_nj45 ",           "G region, EE, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G_EE_nj6             = new TH2D("MR_R2_G_EE_nj6",             "G region, EE, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_CHIsoTemplate_Fake_g_EB    = new TH3D("CHIsoTemplate_Fake_g_EB",    "G-1 region, EB;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_CHIsoTemplate_Fake_g_EE    = new TH3D("CHIsoTemplate_Fake_g_EE",    "G-1 region, EE;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_CHIsoTemplate_Fake_g_EB_MC = new TH3D("CHIsoTemplate_Fake_g_EB_MC", "G-1 region, EB (MC);M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_CHIsoTemplate_Fake_g_EE_MC = new TH3D("CHIsoTemplate_Fake_g_EE_MC", "G-1 region, EE (MC);M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_CHIsoTemplate_Prompt_g_EB  = new TH3D("CHIsoTemplate_Prompt_g_EB",  "G-1 region, EB;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_CHIsoTemplate_Prompt_g_EE  = new TH3D("CHIsoTemplate_Prompt_g_EE",  "G-1 region, EE;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_GNoIso_EB      = new TH3D("MR_R2_CHIso_GNoIso_EB",      "G region (w/o CH iso.), EB;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)",      5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_GNoIso_EE      = new TH3D("MR_R2_CHIso_GNoIso_EE",      "G region (w/o CH iso.), EE;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)",      5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_GNoIso_EB_nj45 = new TH3D("MR_R2_CHIso_GNoIso_EB_nj45", "G region (w/o CH iso.), EB nj45;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_GNoIso_EB_nj6  = new TH3D("MR_R2_CHIso_GNoIso_EB_nj6",  "G region (w/o CH iso.), EB nj6-;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_GNoIso_EE_nj45 = new TH3D("MR_R2_CHIso_GNoIso_EE_nj45 ","G region (w/o CH iso.), EE nj45;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_GNoIso_EE_nj6  = new TH3D("MR_R2_CHIso_GNoIso_EE_nj6",  "G region (w/o CH iso.), EE nj6-;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_gNoIso_EB      = new TH3D("MR_R2_CHIso_gNoIso_EB",      "G-1 region (w/o CH iso.), EB;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)",      5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_gNoIso_EE      = new TH3D("MR_R2_CHIso_gNoIso_EE",      "G-1 region (w/o CH iso.), EE;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)",      5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_gNoIso_EB_nj45 = new TH3D("MR_R2_CHIso_gNoIso_EB_nj45", "G-1 region (w/o CH iso.), EB nj45;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_gNoIso_EB_nj6  = new TH3D("MR_R2_CHIso_gNoIso_EB_nj6",  "G-1 region (w/o CH iso.), EB nj6-;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_gNoIso_EE_nj45 = new TH3D("MR_R2_CHIso_gNoIso_EE_nj45 ","G-1 region (w/o CH iso.), EE nj45;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_CHIso_gNoIso_EE_nj6  = new TH3D("MR_R2_CHIso_gNoIso_EE_nj6",  "G-1 region (w/o CH iso.), EE nj6-;M_{R} (GeV);R^{2};Photon Charged Isolation (GeV)", 5,mrbins, 5,r2bins, 20,chisobins);
  h_MR_R2_IsPrompt_G_EB        = new TH3D("MR_R2_IsPrompt_G_EB",        "G region, EB;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins, 2,isdirectbins);
  h_MR_R2_IsPrompt_G_EE        = new TH3D("MR_R2_IsPrompt_G_EE",        "G region, EE;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins, 2,isdirectbins);
  h_MR_R2_IsDirect_G_EB        = new TH3D("MR_R2_IsDirect_G_EB",        "G region, EB;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins, 2,isdirectbins);
  h_MR_R2_IsDirect_G_EE        = new TH3D("MR_R2_IsDirect_G_EE",        "G region, EE;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins, 2,isdirectbins);
  h_MR_R2_S                    = new TH2D("MR_R2_S",                    "S region;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_S_nj45               = new TH2D("MR_R2_S_nj45",               "S region, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_S_nj6                = new TH2D("MR_R2_S_nj6",                "S region, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G_DirectPrompt       = new TH2D("MR_R2_G_DirectPrompt",       "G region (direct,prompt #gamma);M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_G_DirectPrompt_nj45  = new TH2D("MR_R2_G_DirectPrompt_nj45",  "G region (direct,prompt #gamma), nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G_DirectPrompt_nj6   = new TH2D("MR_R2_G_DirectPrompt_nj6",   "G region (direct,prompt #gamma), nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_Z                    = new TH2D("MR_R2_Z",                    "Z region;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_Z_nj45               = new TH2D("MR_R2_Z_nj45",               "Z region, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_Z_nj6                = new TH2D("MR_R2_Z_nj6",                "Z region, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G                    = new TH2D("MR_R2_G",                    "G region;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_G_nj45               = new TH2D("MR_R2_G_nj45",               "G region, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_G_nj6                = new TH2D("MR_R2_G_nj6",                "G region, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_L                    = new TH2D("MR_R2_L",                    "L region;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_L_nj45               = new TH2D("MR_R2_L_nj45",               "L region, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_L_nj6                = new TH2D("MR_R2_L_nj6",                "L region, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_q                    = new TH2D("MR_R2_q",                    "q region;M_{R} (GeV);R^{2}",       5,mrbins, 5,r2bins);
  h_MR_R2_q_nj45               = new TH2D("MR_R2_q_nj45",               "q region, nj45;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
  h_MR_R2_q_nj6                = new TH2D("MR_R2_q_nj6",                "q region, nj6-;M_{R} (GeV);R^{2}", 5,mrbins, 5,r2bins);
}

//_______________________________________________________
//               Fill Histograms here
void
AnalysisBase::fill_common_histos(DataStruct& d, const bool& varySystematics, const bool& runOnSkim, const unsigned int& syst_index, const double& weight)
{
  if (runOnSkim) {
    if (syst_index == 0) {
      // btag efficiency - No event selection cuts to be applied
      // When making this plot, should remove all baseline cuts
      while(d.jetsAK4.Loop()) {
        size_t i = d.jetsAK4.it;
        if (passLooseJet[i]) {
          if (d.jetsAK4.HadronFlavour[i]==5) {
            h_btag_eff_b_loose ->Fill(d.jetsAK4.Pt[i], passLooseBTag[i]);
            h_btag_eff_b_medium->Fill(d.jetsAK4.Pt[i], passMediumBTag[i]);
          } else if (d.jetsAK4.HadronFlavour[i]==4) {
            h_btag_eff_c_loose ->Fill(d.jetsAK4.Pt[i], passLooseBTag[i]);
            h_btag_eff_c_medium->Fill(d.jetsAK4.Pt[i], passMediumBTag[i]);
          } else {
            h_btag_eff_l_loose ->Fill(d.jetsAK4.Pt[i], passLooseBTag[i]);
            h_btag_eff_l_medium->Fill(d.jetsAK4.Pt[i], passMediumBTag[i]);
          }
        }
      }
      // trigger efficiency, measured in single lepton datasets
      // SingleElectron dataset: Pass HLT_Ele27_WPTight_Gsf && 1 Electron
      // SingleMuon     dataset: Pass HLT_IsoMu24 && 1 Muon
      // SinglePhoton   dataset: Pass all single photon triggers && 1 Photon
      // Baseline cuts to be applied: 3 jets, 1 AK8 jet, MR & R^2
      bool pass_aux_trigger = 0;
      double R2 = d.evt.R2;
      double MR = d.evt.MR;
      int njet = nJet;
      if (TString(sample).Contains("MET")) {
        //if (d.hlt.PFMET120_PFMHT120_IDTight==1&&d.met.Pt[0]>200&&nLepVeto==0&&d.evt.NIsoTrk==0) pass_aux_trigger = 1;      
        if (d.hlt.PFMET120_PFMHT120_IDTight==1&&nLepVeto==0&&nTauVeto==0) pass_aux_trigger = 1;      
      } else if (TString(sample).Contains("SingleElectron")) {
        //if ((d.hlt.Ele23_WPLoose_Gsf==1||d.hlt.Ele27_WPTight_Gsf==1)&&nEleTight>=1&&nMuVeto==0) pass_aux_trigger = 1;
        if ((d.hlt.Ele23_WPLoose_Gsf==1||d.hlt.Ele27_WPTight_Gsf==1)&&nEleVeto==1&&nMuVeto==0) pass_aux_trigger = 1;
      } else if (TString(sample).Contains("SingleMuon")) {
        //if (d.hlt.IsoMu24==1&&nMuTight>=1&&nEleVeto==0) pass_aux_trigger = 1;
        if (d.hlt.IsoMu24==1&&nMuVeto==1&&nEleVeto==0) pass_aux_trigger = 1;
      } else if (TString(sample).Contains("SinglePhoton")) {
	bool OR_HLT_Photon = 
	  d.hlt.Photon22==1 ||
	  d.hlt.Photon30==1 ||
	  d.hlt.Photon36==1 ||
	  d.hlt.Photon50==1 ||
	  d.hlt.Photon75==1 ||
	  d.hlt.Photon90==1 ||
	  d.hlt.Photon120==1 ||
	  d.hlt.Photon165_HE10==1 ||
	  d.hlt.Photon175==1 ||
	  d.hlt.Photon250_NoHE==1 ||
	  d.hlt.Photon300_NoHE==1 ||
	  d.hlt.Photon500==1 ||
	  d.hlt.Photon600==1;
        //if (OR_HLT_Photon==1&&nPhotonSelect==1&&nLepVeto==0&&d.evt.NIsoTrk==0) pass_aux_trigger = 1;
        if (OR_HLT_Photon==1&&nPhotonSelect==1&&nLepVeto==0&&nTauVeto==0) pass_aux_trigger = 1;
	njet = nJetNoPho;
	R2 = R2_pho;
	MR = MR_pho;
      }
      if (pass_aux_trigger) {
        if (nJetAK8>=1 && njet>=3 && MR>=800 && R2>=0.08) {
          if (d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1) {
            h_trigger_pass  ->Fill(AK4_Ht);
            h_trigger2d_pass->Fill(AK4_Ht, d.jetsAK8.Pt[iJetAK8[0]]);
            h_trigger2d_nolep_pass->Fill(AK4_HtNoLep, d.jetsAK8.Pt[iJetAK8[0]]);
          }
          h_trigger_total   ->Fill(AK4_Ht);
          h_trigger2d_total ->Fill(AK4_Ht, d.jetsAK8.Pt[iJetAK8[0]]);
          h_trigger2d_nolep_total ->Fill(AK4_HtNoLep, d.jetsAK8.Pt[iJetAK8[0]]);
        }
      }
      
      // Fill Histos for the Z(nunu) estimate
      if (isData) {
	// All photon counts in G (data)
	if (apply_all_cuts('G')) {
	  bool EB = std::abs(d.pho.SCEta[iPhotonSelect[0]])<1.479;
	  if (EB) h_MR_R2_G_EB->Fill(MR_pho, R2_pho);
	  else    h_MR_R2_G_EE->Fill(MR_pho, R2_pho);
	  if (nJetNoPho>=4) {
	    if (nJetNoPho<6) {
	      if (EB) h_MR_R2_G_EB_nj45->Fill(MR_pho, R2_pho);
	      else    h_MR_R2_G_EE_nj45->Fill(MR_pho, R2_pho);
	    } else {
	      if (EB) h_MR_R2_G_EB_nj6->Fill(MR_pho, R2_pho);
	      else    h_MR_R2_G_EE_nj6->Fill(MR_pho, R2_pho);
	    }
	  }
	}
	// Fake photon templates
	if (apply_all_cuts_except('g',"1Pho")) {
	  if (nPhotonFake==1) {
	    bool EB = std::abs(d.pho.SCEta[iPhotonFake[0]])<1.479;
	    //if (EB) h_CHIsoTemplate_Fake_g_EB->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonFake[0]], sf_weight['g']);
	    //else    h_CHIsoTemplate_Fake_g_EE->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonFake[0]], sf_weight['g']);
	    if (EB) h_CHIsoTemplate_Fake_g_EB->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonFake[0]], sf_weight['g']);
	    else    h_CHIsoTemplate_Fake_g_EE->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonFake[0]], sf_weight['g']);
	  }
	}
      } else if (isBackground) {
	// Prompt photon templates
	if (apply_all_cuts_except('g',"1Pho")) {
	  if (nPhotonPreSelect==1) {
	    if (d.pho.SigmaIEtaIEta[iPhotonPreSelect[0]] < (std::abs(d.pho.SCEta[iPhotonPreSelect[0]])<1.479 ? 0.01022 : 0.03001)) {
	      bool EB = std::abs(d.pho.SCEta[iPhotonPreSelect[0]])<1.479;
	      if (d.pho.isPromptDirect[iPhotonPreSelect[0]]==1||d.pho.isPromptFrag[iPhotonPreSelect[0]]==1) {
		//if (EB) h_CHIsoTemplate_Prompt_g_EB->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]], sf_weight['g']);
		//else    h_CHIsoTemplate_Prompt_g_EE->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]], sf_weight['g']);
		if (EB) h_CHIsoTemplate_Prompt_g_EB->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]], sf_weight['g']);
		else    h_CHIsoTemplate_Prompt_g_EE->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]], sf_weight['g']);
	      } else {
		//if (EB) h_CHIsoTemplate_Fake_g_EB_MC->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]], sf_weight['g']);
		//else    h_CHIsoTemplate_Fake_g_EE_MC->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]], sf_weight['g']);
		if (EB) h_CHIsoTemplate_Fake_g_EB_MC->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]], sf_weight['g']);
		else    h_CHIsoTemplate_Fake_g_EE_MC->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]], sf_weight['g']);
	      }
	    }
	  }
	}
      }
      // Charged Hadron Isolation distributions for Purity measurement
      if (isData) {
	if (apply_all_cuts_except('G',"1Pho")) {
	  if (nPhotonPreSelect==1) {
	    if (d.pho.SigmaIEtaIEta[iPhotonPreSelect[0]] < (std::abs(d.pho.SCEta[iPhotonPreSelect[0]])<1.479 ? 0.01022 : 0.03001)) {
	      // Select 1 photon with CHIso N-1 cuts
	      bool EB = std::abs(d.pho.SCEta[iPhotonPreSelect[0]])<1.479;
	      //if (EB) h_MR_R2_CHIso_GNoIso_EB->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
	      //else    h_MR_R2_CHIso_GNoIso_EE->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
	      if (EB) h_MR_R2_CHIso_GNoIso_EB->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
	      else    h_MR_R2_CHIso_GNoIso_EE->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
	      if (nJetNoPho>=4) {
		if (nJetNoPho<6) {
		  //if (EB) h_MR_R2_CHIso_GNoIso_EB_nj45->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  //else    h_MR_R2_CHIso_GNoIso_EE_nj45->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  if (EB) h_MR_R2_CHIso_GNoIso_EB_nj45->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		  else    h_MR_R2_CHIso_GNoIso_EE_nj45->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		} else {
		  //if (EB) h_MR_R2_CHIso_GNoIso_EB_nj6->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  //else    h_MR_R2_CHIso_GNoIso_EE_nj6->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  if (EB) h_MR_R2_CHIso_GNoIso_EB_nj6->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		  else    h_MR_R2_CHIso_GNoIso_EE_nj6->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		}
	      }
	    }
	  }
	}
	// G-1 region
	if (apply_all_cuts_except('g',"1Pho")) {
	  if (nPhotonPreSelect==1) {
	    if (d.pho.SigmaIEtaIEta[iPhotonPreSelect[0]] < (std::abs(d.pho.SCEta[iPhotonPreSelect[0]])<1.479 ? 0.01022 : 0.03001)) {
	      // Select 1 photon with CHIso N-1 cuts
	      bool EB = std::abs(d.pho.SCEta[iPhotonPreSelect[0]])<1.479;
	      //if (EB) h_MR_R2_CHIso_gNoIso_EB->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
	      //else    h_MR_R2_CHIso_gNoIso_EE->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
	      if (EB) h_MR_R2_CHIso_gNoIso_EB->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
	      else    h_MR_R2_CHIso_gNoIso_EE->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
	      if (nJetNoPho>=4) {
		if (nJetNoPho<6) {
		  //if (EB) h_MR_R2_CHIso_gNoIso_EB_nj45->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  //else    h_MR_R2_CHIso_gNoIso_EE_nj45->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  if (EB) h_MR_R2_CHIso_gNoIso_EB_nj45->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		  else    h_MR_R2_CHIso_gNoIso_EE_nj45->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		} else {
		  //if (EB) h_MR_R2_CHIso_gNoIso_EB_nj6->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  //else    h_MR_R2_CHIso_gNoIso_EE_nj6->Fill(MR_pho, R2_pho, ChargedHadronIsoEACorr[iPhotonPreSelect[0]]);
		  if (EB) h_MR_R2_CHIso_gNoIso_EB_nj6->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		  else    h_MR_R2_CHIso_gNoIso_EE_nj6->Fill(MR_pho, R2_pho, d.pho.ChargedHadronIso[iPhotonPreSelect[0]]);
		}
	      }
	    }
	  }
	}
      } else if (isBackground) {
	// Direct photon ratio from MC
	if (apply_all_cuts('G')) {
	  //std::cout<<"Pass G: "<<MR_pho<<" "<<R2_pho<<" "<<d.evt.RunNumber<<" "<<d.evt.LumiBlock<<" "<<d.evt.EventNumber<<std::endl;
	  bool EB = std::abs(d.pho.SCEta[iPhotonSelect[0]])<1.479;
	  bool isPrompt = d.pho.isPromptDirect[iPhotonSelect[0]]==1||d.pho.isPromptFrag[iPhotonSelect[0]]==1;
	  // Prompt photon fraction
	  if (EB) h_MR_R2_IsPrompt_G_EB->Fill(MR_pho, R2_pho, isPrompt, sf_weight['G']);
	  else    h_MR_R2_IsPrompt_G_EE->Fill(MR_pho, R2_pho, isPrompt, sf_weight['G']);
	  // Direct photon fraction of prompt photons
	  if (isPrompt) {
	    if (EB) h_MR_R2_IsDirect_G_EB->Fill(MR_pho, R2_pho, d.pho.isPromptDirect[iPhotonSelect[0]], sf_weight['G']);
	    else    h_MR_R2_IsDirect_G_EE->Fill(MR_pho, R2_pho, d.pho.isPromptDirect[iPhotonSelect[0]], sf_weight['G']);
	  }
	  // Prompt direct photons in G control region
	  if (d.pho.isPromptDirect[iPhotonSelect[0]]==1) {
	    h_MR_R2_G_DirectPrompt->Fill(MR_pho, R2_pho, sf_weight['G']);
	    if (nJetNoPho>=4) {
	      if (nJetNoPho<6) h_MR_R2_G_DirectPrompt_nj45->Fill(MR_pho, R2_pho, sf_weight['G']);
	      else             h_MR_R2_G_DirectPrompt_nj6 ->Fill(MR_pho, R2_pho, sf_weight['G']);
	    }
	  }
	}
      }
      // Counts for transfer factors/ratios
      if (apply_all_cuts('S')) {
	h_MR_R2_S->Fill(d.evt.MR, d.evt.R2, sf_weight['S']);
	if (nJet>=4) {
	  if (nJet<6) h_MR_R2_S_nj45->Fill(d.evt.MR, d.evt.R2, sf_weight['S']);
	  else        h_MR_R2_S_nj6 ->Fill(d.evt.MR, d.evt.R2, sf_weight['S']);
	}
      }
      // photon/1-lepton/2-lepton estimates
      if (apply_all_cuts('Z')) {
	h_MR_R2_Z->Fill(d.evt.MR, R2_ll, sf_weight['Z']);
	if (nJet>=4) {
	  if (nJet<6) h_MR_R2_Z_nj45->Fill(d.evt.MR, R2_ll, sf_weight['Z']);
	  else        h_MR_R2_Z_nj6 ->Fill(d.evt.MR, R2_ll, sf_weight['Z']);
	}
      }
      if (apply_all_cuts('L')) {
	h_MR_R2_L->Fill(d.evt.MR, R2_1l, sf_weight['L']);
	if (nJet>=4) {
	  if (nJet<6) h_MR_R2_L_nj45->Fill(d.evt.MR, R2_1l, sf_weight['L']);
	  else        h_MR_R2_L_nj6 ->Fill(d.evt.MR, R2_1l, sf_weight['L']);
	}
      }
      if (apply_all_cuts('G')) {
	h_MR_R2_G->Fill(MR_pho, R2_pho, sf_weight['G']);
	if (nJetNoPho>=4) {
	  if (nJetNoPho<6) h_MR_R2_G_nj45->Fill(MR_pho, R2_pho, sf_weight['G']);
	  else             h_MR_R2_G_nj6 ->Fill(MR_pho, R2_pho, sf_weight['G']);
	}
      }
      // closure test in Q'
      if (apply_all_cuts('q')) {
	h_MR_R2_q->Fill(d.evt.MR, d.evt.R2, sf_weight['q']);
	if (nJet>=4) {
	  if (nJet<6) h_MR_R2_q_nj45->Fill(d.evt.MR, d.evt.R2, sf_weight['q']);
	  else        h_MR_R2_q_nj6 ->Fill(d.evt.MR, d.evt.R2, sf_weight['q']);
	}
      }


    } // end syst_index == 0

    if (varySystematics) {
      // Fill plots for systematics
      const std::vector<char> regions = {'S', 's', 'T', 'W', 'Q', 'q', 'Z', 'L', 'G'};
      if (isData) {
        if (TString(sample).Contains("JetHT")) {
          for (size_t i=0; i<regions.size(); ++i) {
	    int MRR2_bin = calc_mrr2_bin(d, regions[i]);
            if (this->apply_all_cuts(regions[i])) {
              vh_MRR2_data[i]->Fill(MRR2_bin);
	      if (nJet>=4) {
		if (nJet<6) vh_MRR2_data_nj45[i]->Fill(MRR2_bin);
		else        vh_MRR2_data_nj6[i] ->Fill(MRR2_bin);
	      }
            }
          }
	  // Save also vertex count (for signal systematics)
	  if (nJetAK8>=1 && nJet>=3 && d.evt.MR>=800 && d.evt.R2>=0.08) {
	    if (d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1) {
	      h_nvtx->Fill(d.evt.NGoodVtx);
	      h_npvLowHigh_data->Fill(d.evt.NGoodVtx);
	    }
	  }
        }
      } else if (isSignal) {
        if (apply_all_cuts('S')) {
          uint32_t mMother = TString(sample).Contains("T2tt") ? std::round(d.evt.SUSY_Stop_Mass/5.0)*5 : std::round(d.evt.SUSY_Gluino_Mass/25.0)*25;
          uint32_t mLSP    = TString(sample).Contains("T2tt") ? std::round(d.evt.SUSY_LSP_Mass /5.0)*5 : std::round(d.evt.SUSY_LSP_Mass   /25.0)*25;
          uint32_t signal_bin = mMother * 10000 + mLSP;
	  int MRR2_bin = calc_mrr2_bin(d, 'S');
          m_vh_MRR2_sig[signal_bin][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	  if (nJet>=4) {
	    if (nJet<6) m_vh_MRR2_sig_nj45[signal_bin][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	    else        m_vh_MRR2_sig_nj6 [signal_bin][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	  }
	  // npv for Signal
	  if (TString(sample).Contains("T2tt")) {
	    vh_npvLowHigh_signal[1]->Fill(d.evt.SUSY_Stop_Mass,   d.evt.SUSY_LSP_Mass, d.evt.NGoodVtx);
	  } else {
	    vh_npvLowHigh_signal[0]->Fill(d.evt.SUSY_Gluino_Mass, d.evt.SUSY_LSP_Mass, d.evt.NGoodVtx);
	  }
        }
      } else {
        // Backgrounds
        for (size_t i=0; i<regions.size(); ++i) {
          if (apply_all_cuts(regions[i])) {
	    int MRR2_bin = calc_mrr2_bin(d, regions[i]);
            vvh_MRR2_bkg[i][syst_index]->Fill(MRR2_bin, sf_weight[regions[i]]);
	    if (nJet>=4) {
	      if (nJet<6) vvh_MRR2_bkg_nj45[i][syst_index]->Fill(MRR2_bin, sf_weight[regions[i]]);
	      else        vvh_MRR2_bkg_nj6 [i][syst_index]->Fill(MRR2_bin, sf_weight[regions[i]]);
	    }
	    if (regions[i]=='G') {
	      // Save also the direct prompt photons for G region
	      // Needed for transfer factors
	      if (d.pho.isPromptDirect[iPhotonSelect[0]]==1) {
		vvh_MRR2_bkg[i+1][syst_index]->Fill(MRR2_bin, sf_weight[regions[i]]);
		if (nJet>=4) {
		  if (nJet<6) vvh_MRR2_bkg_nj45[i+1][syst_index]->Fill(MRR2_bin, sf_weight[regions[i]]);
		  else        vvh_MRR2_bkg_nj6 [i+1][syst_index]->Fill(MRR2_bin, sf_weight[regions[i]]);
		}
	      }
	    }
          }
        }
	// Signal regions with loosened tau32/21 requirements
#if TOP == 0
	if (apply_all_cuts_except('S',"1W")&&nWMassTag>=1) {
	  int MRR2_bin = calc_mrr2_bin(d, 'S');
	  vvh_MRR2_bkg[regions.size()+1][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	  if (nJet>=4) {
	    if (nJet<6) vvh_MRR2_bkg_nj45[regions.size()+1][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	    else        vvh_MRR2_bkg_nj6 [regions.size()+1][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	  }
	}
#else
	if (apply_all_cuts_except('S',"1Top")&&nHadTopTag4>=1) {
	  int MRR2_bin = calc_mrr2_bin(d, 'S');
	  vvh_MRR2_bkg[regions.size()+1][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	  if (nJet>=4) {
	    if (nJet<6) vvh_MRR2_bkg_nj45[regions.size()+1][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	    else        vvh_MRR2_bkg_nj6 [regions.size()+1][syst_index]->Fill(MRR2_bin, sf_weight['S']);
	  }
	}
#endif
      }
    } // end of systematics
  } else {
    // Skimming only histos
    nJet = nJetAK8 = 0;
    while(d.jetsAK4.Loop()) {
      size_t i = d.jetsAK4.it;
      // Jet ID
      if ( d.jetsAK4.looseJetID[i] == 1 &&
	   d.jetsAK4.Pt[i]         >= JET_AK4_PT_CUT &&
	   std::abs(d.jetsAK4.Eta[i])  <  JET_AK4_ETA_CUT ) {
	nJet++;
      }
    }
    while(d.jetsAK8.Loop()) {
      size_t i = d.jetsAK8.it;
      // Jet ID
      if ( d.jetsAK8.looseJetID[i] == 1 &&
	   d.jetsAK8.Pt[i]         >= JET_AK8_PT_CUT &&
	   std::abs(d.jetsAK8.Eta[i])  <  JET_AK8_ETA_CUT ) {
	nJetAK8++;
      }
    }
    // Save the number of vertices in Data and Signal
    if (isData) {
      // Use baseline selection
      if (nJetAK8>=1 && nJet>=3 && d.evt.MR>=800 && d.evt.R2>=0.08) {
	if (d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1) {
	  h_npvLowHigh_data->Fill(d.evt.NGoodVtx);
	}
      }
    } else if (isSignal) {
      if (TString(sample).Contains("T2tt")) {
	vh_npvLowHigh_signal[1]->Fill(d.evt.SUSY_Stop_Mass,   d.evt.SUSY_LSP_Mass, d.evt.NGoodVtx);
      } else {
	vh_npvLowHigh_signal[0]->Fill(d.evt.SUSY_Gluino_Mass, d.evt.SUSY_LSP_Mass, d.evt.NGoodVtx);
      }
    }
  }
}

//_______________________________________________________
//           Read cross-section from ntuple
double
AnalysisBase::get_xsec_from_ntuple(const std::vector<std::string>& filenames, const bool& runOnSkim)
{
  float evt_XSec=0, prev_XSec=0;
  for (const auto& filename : filenames) {
    TFile *f = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)f->Get(runOnSkim ? "B2GTree" : "B2GTTreeMaker/B2GTree");
    tree->GetBranch("evt_XSec")->SetAddress(&evt_XSec);
    tree->GetEntry(0);
    f->Close();
    if (prev_XSec!=0&&prev_XSec!=evt_XSec) {
      utils::error("AnalysisBase - Files added with different cross-sections. Please, add them separately!");
      return 0;
    }
    prev_XSec = evt_XSec;
  }
  return evt_XSec;
}

//_______________________________________________________
//           Read cross-section from txt file
std::pair<double, double>
AnalysisBase::get_xsec_totweight_from_txt_file(const std::string& txt_file)
{
  double XSec = 0, Totweight = 0;
  std::ifstream xsecFile(txt_file.c_str());
  if ( !xsecFile.good() ) {
    return std::make_pair(0,0);
    std::cout<<"Unable to open cross-section file: "<<txt_file<<std::endl;
    utils::error("Please provide the correct txt file for Cross-sections in settings.h!");
  } else {

    std::string line;
    std::string shortname, primary_dataset;
    double xsec, totweight;
    while ( std::getline(xsecFile, line) ) {
      std::stringstream nth_line;
      nth_line<<line;
      nth_line>>shortname;
      nth_line>>primary_dataset;
      nth_line>>xsec;
      nth_line>>totweight;
      // For skimmed samples, remove certain postfixes
      // Please, synchronize with setup.py script
      std::string dirname = sample;
      for (std::string pf : { "_2", "_ext1", "_ext2", "_ext3", "_ext4", "_backup", "_unskimmed" })
	if (TString(dirname).EndsWith(pf.c_str())) dirname.erase(dirname.size()-pf.size(), pf.size());
      if (dirname==shortname) {
	XSec = xsec;
	Totweight = totweight;
      }

    }
  }
  if (XSec == 0) {
    std::cout<<"No crossection found for "<<sample<<" in cross section file: "<<txt_file<<std::endl;
    utils::error("Please fix the cross-section file in settings.h!");
  }

  return std::make_pair(XSec, Totweight);
}

//_______________________________________________________
//          Read total weight from ntuple histos
double
AnalysisBase::get_totweight_from_ntuple(const std::vector<std::string>& filenames, const bool& runOnSkim)
{
  // Merging totweight histos
  for (const auto& filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    h_totweight->Add((TH1D*)f->Get(runOnSkim ? "totweight" : "EventCounter/totweight"));
    //std::cout<<f<<" "<<filename<<" "<<histoname.c_str()<<" "<<f->Get(histoname.c_str())<<" "<<h_totweight->GetBinContent(1)<<std::endl;
#if VER > 2
    // Get also total ISR jet info
    if (runOnSkim) {
      h_totweight_isr->Add((TH1D*)f->Get("totweight_isr"));
    }
#endif
    f->Close();
  }
  return h_totweight->GetBinContent(1);
}

//_______________________________________________________
//       Calculate weight normalization for signal
void
AnalysisBase::calc_weightnorm_histo_from_ntuple(const std::vector<std::string>& filenames, const double& intLumi, const std::vector<std::string>& vname_signal,
						const bool& runOnSkim, const bool& varySystematics, TDirectory* dir, bool verbose=1)
{
  // Find the index of the current signal
  int signal_index = TString(filenames[0]).Contains("T2tt"); // 0: Mlsp vs Mgluino - T1tttt, T1ttbb, T5ttcc, T5tttt; 1: Mlsp vs Mstop - T2tt

  // Merge totweight histos
  std::map<int, double> xsec_mother;
  for (const auto& filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    // Get total weight
    if (signal_index==0) {
      vh_totweight_signal[signal_index]->Add((TH2D*)f->Get(runOnSkim ? "totweight_T1tttt" : "EventCounter/h_totweight_T1tttt"));
    } else {
      vh_totweight_signal[signal_index]->Add((TH2D*)f->Get(runOnSkim ? "totweight_T2tt" :   "EventCounter/h_totweight_T2tt"));
    }
#if VER > 2
    // Get also total ISR jet weight
    if (runOnSkim) {
      if (signal_index==0) {
	vh_totweight_signal_isr[signal_index]->Add((TH3D*)f->Get("totweight_isr_T1tttt"));
      } else {
	vh_totweight_signal_isr[signal_index]->Add((TH3D*)f->Get("totweight_isr_T2tt"));
      }
    }
#endif
    f->Close();
  }

  // Set xsec for each gluino/stop mass bin
  // Read gluino/stop xsec from same file used in TTree step
  for (int binx=1, nbinx=vh_xsec_signal[signal_index]->GetNbinsX(); binx<=nbinx; ++binx) {
    double mMother = vh_xsec_signal[signal_index]->GetXaxis()->GetBinCenter(binx);
    xsec_mother[binx] = signal_index ? GetStopXSec(mMother).first : GetGluinoXSec(mMother).first; // first: mean xsec (pb), second: error (%)
    for (int biny=1, nbiny=vh_xsec_signal[signal_index]->GetNbinsY(); biny<=nbiny; ++biny)
      vh_xsec_signal[signal_index]->SetBinContent(binx, biny, xsec_mother[binx]);
  }
  // Calculate weight normalization
  // weightnorm = (settings.intLumi*xsec)/totweight;
  // Divide(h1,h2,c1,c2) --> c1*h1/(c2*h2)
  vh_weightnorm_signal[signal_index]->Divide(vh_xsec_signal[signal_index], vh_totweight_signal[signal_index], intLumi);
  std::map<uint32_t, std::string> signal_bins;
  if (verbose) {
    std::cout<<"- Signal:"<<std::endl;
    for (int binx=1, nbinx=vh_xsec_signal[signal_index]->GetNbinsX(); binx<=nbinx; ++binx) 
      for (int biny=1, nbiny=vh_xsec_signal[signal_index]->GetNbinsY(); biny<=nbiny; ++biny) {
        double mMother = vh_xsec_signal[signal_index]->GetXaxis()->GetBinCenter(binx);
        double mLSP = vh_xsec_signal[signal_index]->GetYaxis()->GetBinCenter(biny);
        double xsec  = vh_xsec_signal[signal_index]      ->GetBinContent(binx, biny);
        double totw  = vh_totweight_signal[signal_index] ->GetBinContent(binx, biny);
        double wnorm = vh_weightnorm_signal[signal_index]->GetBinContent(binx, biny);
        if (totw>0) {
	  std::cout<<(signal_index?"  Bin: M(s~)=":"  Bin: M(g~)=")<<mMother<<" M(LSP)="<<mLSP<<":   xsec="<<xsec<<" totweight="<<totw<<" weightnorm="<<wnorm<<std::endl;
	  uint32_t bin = mMother * 10000 + mLSP;
	  std::stringstream ss;
	  ss<<"_"<<mMother<<"_"<<mLSP;
	  signal_bins[bin] = ss.str();
	}
      }
    std::cout<<std::endl;
  }

  if (varySystematics&&isSignal) {
    dir->cd();
    // Declare systematics histos for Signals
    for (const auto& bin : signal_bins) {
      m_vh_MRR2_sig     [bin.first] = std::vector<TH1D*>();
      m_vh_MRR2_sig_nj45[bin.first] = std::vector<TH1D*>();
      m_vh_MRR2_sig_nj6 [bin.first] = std::vector<TH1D*>();
      m_vh_MRR2_sig     [bin.first].push_back(new TH1D((std::string("MRR2_S_signal")+bin.second).c_str(),         ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
      m_vh_MRR2_sig_nj45[bin.first].push_back(new TH1D((std::string("MRR2_S_signal")+bin.second+"_nj45").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
      m_vh_MRR2_sig_nj6 [bin.first].push_back(new TH1D((std::string("MRR2_S_signal")+bin.second+"_nj6").c_str(),  ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
      for (size_t j=0; j<syst.size(); ++j) {
	std::stringstream ss, ss45, ss6;
	ss  <<"MRR2_S_signal"<<bin.second<<"_"<<syst[j];
	ss45<<"MRR2_S_signal"<<bin.second<<"_nj45_"<<syst[j];
	ss6 <<"MRR2_S_signal"<<bin.second<<"_nj6_"<<syst[j];
	m_vh_MRR2_sig     [bin.first].push_back(new TH1D((ss.str()+"Up").c_str(),     ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	m_vh_MRR2_sig     [bin.first].push_back(new TH1D((ss.str()+"Down").c_str(),   ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	m_vh_MRR2_sig_nj45[bin.first].push_back(new TH1D((ss45.str()+"Up").c_str(),   ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	m_vh_MRR2_sig_nj45[bin.first].push_back(new TH1D((ss45.str()+"Down").c_str(), ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	m_vh_MRR2_sig_nj6 [bin.first].push_back(new TH1D((ss6.str()+"Up").c_str(),    ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
	m_vh_MRR2_sig_nj6 [bin.first].push_back(new TH1D((ss6.str()+"Down").c_str(),  ";MR/R^{2} bins (unrolled);Counts", 25,0,25));
      }
    }
  }
}


//_______________________________________________________
//             Load pile-up reweighting infos
void
AnalysisBase::init_pileup_reweighting(const std::string& pileupDir, const bool& runOnSkim, const std::vector<std::string>& filenames)
{
  // Get data histogram (generated by pileupCalc.py script)
  TFile* f_pileup_data = TFile::Open((pileupDir+"data_pileup.root").c_str());
  h_pileup_data->Add((TH1D*)f_pileup_data->Get("pileup"));
  f_pileup_data->Close();
  // Also get up/down variations
  TFile* f_pileup_data_down = TFile::Open((pileupDir+"data_pileup_down.root").c_str());
  h_pileup_data_down->Add((TH1D*)f_pileup_data_down->Get("pileup"));
  f_pileup_data_down->Close();
  TFile* f_pileup_data_up = TFile::Open((pileupDir+"data_pileup_up.root").c_str());
  h_pileup_data_up->Add((TH1D*)f_pileup_data_up->Get("pileup"));
  f_pileup_data_up->Close();
  // get mc histogram (used to generate mc pile-up)
  TFile* f_pileup_mc = TFile::Open((pileupDir+"mc_pileup.root").c_str());
  h_pileup_mc->Add((TH1D*)f_pileup_mc->Get("pileup"));
  f_pileup_mc->Close();
  // // Get mc histogram saved inside the ntuple (unfiltered pileup distribution)
  // std::cout<<h_pileup_mc->GetEntries()<<std::endl;
  // for (const auto& filename : filenames) {
  //   TFile* f_pileup_mc = TFile::Open(filename.c_str());
  //   h_pileup_mc->Add((TH1D*)f_pileup_mc->Get(runOnSkim ? "pileup_mc" : "EventCounter/pileup"));
  //   f_pileup_mc->Close();
  //   std::cout<<h_pileup_mc->GetEntries()<<std::endl;
  // }
  // Divide normalized data histo by normalized mc histo to get pileup weights for each bin
  h_pileup_weight     ->Divide(h_pileup_data,      h_pileup_mc, 1/h_pileup_data->Integral(),      1/h_pileup_mc->Integral());
  h_pileup_weight_down->Divide(h_pileup_data_down, h_pileup_mc, 1/h_pileup_data_down->Integral(), 1/h_pileup_mc->Integral());    
  h_pileup_weight_up  ->Divide(h_pileup_data_up,   h_pileup_mc, 1/h_pileup_data_up->Integral(),   1/h_pileup_mc->Integral());    
}

//_______________________________________________________
//              function to get scaled weight
double
AnalysisBase::get_syst_weight(const double& weight_nominal, const double& weight_up, const double& weight_down, const double& nSigma)
{
  double w = weight_nominal;
  if (nSigma == 0) {
    return w;
  } else {
    // Compute the weight according to the systematic variation considered
    // Use difference between nominal and up/down as 1 sigma variation 
    double dw_up = weight_up - weight_nominal;
    double dw_down = weight_nominal - weight_down;
    if (nSigma >= 0.) {
      w += nSigma*dw_up; 
    } else {
      w += nSigma*dw_down;
    }
    return w; 
  }
}

double
AnalysisBase::get_syst_weight(const double& weight_nominal, const double& uncertainty, const double& nSigma)
{
  double w = weight_nominal;
  // Use symmetrical difference for up/down variation
  if (nSigma!=0.) w *= 1.0 + nSigma * uncertainty;
  return w;
}

//_______________________________________________________
//                  Top pt reweighting
double
AnalysisBase::get_toppt_weight(DataStruct& data, const double& nSigmaTopPt, const unsigned int& syst_index, const bool& runOnSkim)
{
  double w_nom = 1;//, n=0;
  while(data.gen.Loop()) {
    size_t i = data.gen.it;
    // Select last copy of the particles only (i.e. their daughters are different)
    if (data.gen.Dau0ID[i]!=data.gen.ID[i]&&data.gen.Dau1ID[i]!=data.gen.ID[i]) {
      if(abs(data.gen.ID[i])==6) {
	double a = 0.0615, b = -0.0005;
	w_nom *= std::exp(a + b * data.gen.Pt[i]);
	//std::cout<<"evt="<<data.evt.EventNumber<<" i="<<i<<" top id="<<data.gen.ID[i]<<" dau0 id="<<data.gen.Dau0ID[i]<<" dau1 id="<<data.gen.Dau1ID[i]<<" pt="<<data.gen.Pt[i]<<" w="<<w_nom<<std::endl;
	//n+=1;
      }
    }
  }
  w_nom = std::sqrt(w_nom);
  //std::cout<<"N top = "<<n<<" w_nom = "<<w_nom<<std::endl<<std::endl;
  double w_toppt_up = 1;
  double w_toppt = std::sqrt(w_nom);
  double w_toppt_down = w_nom; // Nominal weight is typically below 1 (below top pt>123) --> use as down variation
  double w = get_syst_weight(w_toppt, w_toppt_up, w_toppt_down, nSigmaTopPt);
  if (syst_index==0&&!runOnSkim) {
    h_totweight_toppt->Fill(0);
    h_totweight_toppt->Fill(1, w_toppt);
  }
  return w;
}

//_______________________________________________________
//                    ISR reweighting
double
AnalysisBase::get_isr_weight(DataStruct& data, const double& nSigmaISR, const unsigned int& syst_index, const bool& runOnSkim)
{
  // Implementing the reweighting in this presentation:
  // https://indico.cern.ch/event/592621/contributions/2398559/attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf
  // Using the values found on slide 8 (T2tt and T1tttt)
  double w = 1;
  // ttbar ISR reweighting not needed, we do top pt reweighting!
  if (isSignal) {
    double d = TString(sample).Contains("T2tt") ? 1.121 : 1.143;
    int n = data.evt.NISRJets;
    double w_nom = 0;
    if      (n==0) w_nom = d;
    else if (n==1) w_nom = d * 0.920;
    else if (n==2) w_nom = d * 0.821;
    else if (n==3) w_nom = d * 0.715;
    else if (n==4) w_nom = d * 0.662;
    else if (n==5) w_nom = d * 0.561;
    else           w_nom = d * 0.511;
    double err = (1-w_nom)/2;
    double w_isr_up   = w_nom + err;
    double w_isr      = w_nom;
    double w_isr_down = w_nom - err;
    w = get_syst_weight(w_isr, w_isr_up, w_isr_down, nSigmaISR);
    if (syst_index==0&&!runOnSkim) {
      if (TString(sample).Contains("T2tt")) {
	vh_totweight_signal_isr[1]->Fill(data.evt.SUSY_Stop_Mass,   data.evt.SUSY_LSP_Mass, 0);
	vh_totweight_signal_isr[1]->Fill(data.evt.SUSY_Stop_Mass,   data.evt.SUSY_LSP_Mass, 1, w_isr);
      } else {
	vh_totweight_signal_isr[0]->Fill(data.evt.SUSY_Gluino_Mass, data.evt.SUSY_LSP_Mass, 0);
	vh_totweight_signal_isr[0]->Fill(data.evt.SUSY_Gluino_Mass, data.evt.SUSY_LSP_Mass, 1, w_isr);
      }
    }
  }
  return w;
}

//_______________________________________________________
//                  Get pile-up weight
double
AnalysisBase::get_pileup_weight(DataStruct& data, const double& nSigmaPU, const unsigned int& syst_index, const bool& runOnSkim)
{
  if (!isSignal) {
    // Background
    int pu_bin = data.pu.NtrueInt+1; // eg. pileup 0, is filled in bin 1
    double w_pileup = h_pileup_weight->GetBinContent(pu_bin);
    double w_pileup_up = h_pileup_weight_up->GetBinContent(pu_bin);
    double w_pileup_down = h_pileup_weight_down->GetBinContent(pu_bin);
    double w = get_syst_weight(w_pileup, w_pileup_up, w_pileup_down, nSigmaPU);
    if (syst_index==0&&!runOnSkim) {
      h_totweight_pileup->Fill(0);
      h_totweight_pileup->Fill(1, w_pileup);
    }
    return w;
  } else {
    // Signal
    // Do not reweight, but split to a low/high pileup region
    if (nSigmaPU==1) return data.evt.NGoodVtx>=h_npvLowHigh_data->GetBinLowEdge(2);
    else if (nSigmaPU==-1) return data.evt.NGoodVtx<h_npvLowHigh_data->GetBinLowEdge(2);
    return 1;
  }
}


//____________________________________________________
//                  HT reweighting

// Silver JSON
/*
  const double p0[2]     = { 1.16434, 1.00188 };
  const double p0_err[2] = { 0.00459931, 0.0266651 };
  const double p1[2]     = { -0.000142391, -7.80628e-05 };
  const double p1_err[2] = { 3.62929e-06, 1.11035e-05 };
*/

// Golden JSON
const double p0[2]     = { 1.17155, 1.00513 };
const double p0_err[2] = { 0.00477137, 0.028861 };
const double p1[2]     = { -0.000143935, -7.81881e-05 };
const double p1_err[2] = { 3.79477e-06, 1.20209e-05 };

double
AnalysisBase::get_ht_weight(DataStruct& data, const double& nSigmaHT)
{
  // Using method described by Julie Hogan:
  // https://indico.cern.ch/event/508384/contributions/2029874/attachments/1255336/1852975/JetRwtIssues_B2GWkshp_040816.pdf
  // Use linear functions calculated with scripts/CalcHTScaleFactors.C macro
  // linear function(s): p0 + p1 * HT

  // Calculate unscaled jet HT
  double ht = 0; for (const auto& pt : AK8_Pt) ht += pt;

  double w = 1.0;
  if (ht>=800&&ht<2000)
    w *= get_syst_weight(p0[0], p0_err[0]/p0[0], nSigmaHT) + get_syst_weight(p1[0], p1_err[0]/p1[0], nSigmaHT) * ht;
  else if (ht>=2000)
    w *= get_syst_weight(p0[1], p0_err[1]/p0[1], nSigmaHT) + get_syst_weight(p1[1], p1_err[1]/p1[1], nSigmaHT) * ht;

  return w;
}


//_______________________________________________________
//                  Get alpha_s weight
double
AnalysisBase::get_alphas_weight(const std::vector<float>& alphas_Weights, const double& nSigmaAlphaS, const int& LHA_PDF_ID)
{
  // A set of two weights corresponding to 
  // Powheg:  alpha_s = 0.118 -+ 0.002 
  // aMC@NLO: alpha_s = 0.118 -+ 0.001
  // Recommendation is to use +- 0.0015 --> rescale difference by 0.75 or 1.5
  // Treat weight as usual, gaussian, rescale to desired nSigma
  double w_alphas = 1;
  double w_alphas_up   = alphas_Weights[1];
  double w_alphas_down = alphas_Weights[0];
  double nSigma_0_0015 = nSigmaAlphaS;
  if (LHA_PDF_ID==260000||LHA_PDF_ID==260400) {
    // Powheg samples have -+ 0.001
    nSigma_0_0015 *= 1.5;
  } else {
    // aMC@NLO samples have -+ 0.002
    nSigma_0_0015 *= 0.75;
  }
  w_alphas = get_syst_weight(w_alphas, w_alphas_up, w_alphas_down, nSigma_0_0015);
  return w_alphas;
}


//_______________________________________________________
//                  Get scale weight
double
AnalysisBase::get_scale_weight(const std::vector<float>& scale_Weights, const std::vector<double>& scale_weight_norm, const double& nSigmaScale, const unsigned int& numScale)
{
  /*
    Typical LHE run info:
    <weightgroup combine="envelope" type="Central scale variation">
    <weight id="1"> mur=1 muf=1 </weight>
    <weight id="2"> mur=1 muf=2 </weight>     --> save [0]
    <weight id="3"> mur=1 muf=0.5 </weight>   --> save [1]
    <weight id="4"> mur=2 muf=1 </weight>     --> save [2]
    <weight id="5"> mur=2 muf=2 </weight>     --> save [3]
    <weight id="6"> mur=2 muf=0.5 </weight>
    <weight id="7"> mur=0.5 muf=1 </weight>   --> save [4]
    <weight id="8"> mur=0.5 muf=2 </weight>
    <weight id="9"> mur=0.5 muf=0.5 </weight> --> save [5]
    </weightgroup>

    SUSY GEN Lumi info:
    GEN:   LHE, id = 1, Central scale variation,  mur=1 muf=1                               
    GEN:   LHE, id = 2, Central scale variation,  mur=1 muf=2                               
    GEN:   LHE, id = 3, Central scale variation,  mur=1 muf=0.5                             
    GEN:   LHE, id = 4, Central scale variation,  mur=2 muf=1                               
    GEN:   LHE, id = 5, Central scale variation,  mur=2 muf=2                               
    GEN:   LHE, id = 6, Central scale variation,  mur=2 muf=0.5                             
    GEN:   LHE, id = 7, Central scale variation,  mur=0.5 muf=1                             
    GEN:   LHE, id = 8, Central scale variation,  mur=0.5 muf=2                             
    GEN:   LHE, id = 9, Central scale variation,  mur=0.5 muf=0.5

    https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L195-L202
    We save only ids: 2,3,4,5,7,9 (in this order)

    The idea here is to randomly choose to vary mu_f or mu_r or both simulataneously
    and rescale weight difference the usual way by desired nSigma
  */
  if (nSigmaScale==0) return 1; // No systematics
  if (scale_Weights.empty()) return 1; // ST samples are known to miss scale weights
  if (scale_weight_norm.empty()) {
    utils::error("AnalysisBase - Scale weight normalizations were not provided for this sample, rerun scripts/get_scaleweight_norm.py on unskimmed ntuple");
  }
  double w_scale = 1;
  double w_scale_up = 1;   // Corresponds to 0.5 (More signal events)
  double w_scale_down = 1; // Corresponds to 2.0 (Less signal events)
  if (numScale==1) {
    // Vary factorization scale
    // fix mu_r = 1.0, vary mu_f = 0,5, 2.0
    w_scale_up   = scale_Weights[1] / scale_weight_norm[1];
    w_scale_down = scale_Weights[0] / scale_weight_norm[0];
  } else if (numScale==2) {
    // Vary renormalization scale
    // fix mu_f = 1.0, vary mu_r = 0,5, 2.0
    w_scale_up   = scale_Weights[4] / scale_weight_norm[4];
    w_scale_down = scale_Weights[2] / scale_weight_norm[2];
  } else if (numScale==3) {
    // Vary both simulatneously
    // mu_r = mu_f = 0,5, 2.0
    w_scale_up   = scale_Weights[5] / scale_weight_norm[5];
    w_scale_down = scale_Weights[3] / scale_weight_norm[3];
  }
  w_scale = get_syst_weight(w_scale, w_scale_up, w_scale_down, nSigmaScale);
  return w_scale;
}

//_______________________________________________________
//                Benchmarking (batch) jobs

void
AnalysisBase::job_monitoring(const int& entry, const int& nevents, const std::string& curr_file, const float threshold=5)
{
  if (entry==0) {
    sw_1k_ ->Start(kFALSE);
    sw_10k_->Start(kFALSE);
    sw_job_->Start(kFALSE);
  } else {
    double time_1 = sw_1_->RealTime();
    sw_1_->Reset(); sw_1_->Start(kFALSE);
    if (time_1>threshold&&entry!=1) {
      ++bad_files[curr_file];
      //std::cout<<"Bad read - time threshold: "<<threshold<<"s, unresponsive time: "<<time_1<<" s, entry: "<<entry<<" occurence: "<<bad_files[curr_file]<<std::endl;
      //if(bad_files[curr_file]==5) {
      //  std::cout<<"Badly readable file found: "<<curr_file<<std::endl;
      //  if (crash_job) {
      //    std::cout<<"Reached "<<threshold<<" occurences, exiting the job and requesting new EOS copy"<<std::endl;
      //    exit(1);
      //  }
      //}
    }
    if (entry%1000==0) {
      double meas_1k = 1000/sw_1k_->RealTime();
      h_read_speed_1k->Fill(meas_1k);
      sw_1k_->Reset();
      sw_1k_->Start(kFALSE);
      //std::cout<<"Meas  1k: "<<meas_1k<<std::endl;
    }
    if (entry%10000==0) {
      double meas_10k = 10000/sw_10k_->RealTime();
      h_read_speed_10k->Fill(meas_10k);
      h_read_speed_vs_nevt_10k->Fill(entry, meas_10k);
      sw_10k_->Reset();
      sw_10k_->Start(kFALSE);
      //std::cout<<"Meas 10k: "<<meas_10k<<std::endl;
    }
    if (entry+1==nevents) {
      sw_job_->Stop();
      double meas_job = nevents/sw_job_->RealTime();
      h_read_speed_job->Fill(meas_job);
      h_read_speed_vs_nevt_job->Fill(nevents, meas_job);
      h_runtime_job->Fill(sw_job_->RealTime()/60.);
      h_runtime_vs_nevt_job->Fill(nevents, sw_job_->RealTime()/60.);
      std::cout<<"JobMonitoringReport RunTime(s): "<<sw_job_->RealTime()<<" Nevents: "<<nevents<<" Nevt/s: "<<meas_job<<std::endl;
      for (const auto& bad_file : bad_files)
	std::cout<<"Badly readable file found: "<<bad_file.first<<" N_occurence: "<<bad_file.second<<std::endl;
    }
  }
}

//_______________________________________________________
//                Calculate scale factors

TProfile* eff_btag_b_loose;
TProfile* eff_btag_c_loose;
TProfile* eff_btag_l_loose;
TProfile* eff_btag_b_medium;
TProfile* eff_btag_c_medium;
TProfile* eff_btag_l_medium;

TH2F* eff_full_ele_reco;
TH2F* eff_full_ele_vetoid;
TH2F* eff_full_ele_looseid;
TH2F* eff_full_ele_mediumid;
TH2F* eff_full_ele_mvalooseid_tightip2d;
TH2F* eff_full_ele_miniiso01;
TH2F* eff_full_ele_miniiso02;
TH2F* eff_full_ele_miniiso04;
TH2D* eff_fast_ele_vetoid;
TH2D* eff_fast_ele_looseid;
TH2D* eff_fast_ele_mediumid;
TH2D* eff_fast_ele_mvalooseid_tightip2d;
TH2D* eff_fast_ele_miniiso01;
TH2D* eff_fast_ele_miniiso02;
TH2D* eff_fast_ele_miniiso04;
TH2F* eff_full_ele_veto;
TGraphAsymmErrors* eff_full_muon_trk;
TGraphAsymmErrors* eff_full_muon_trk_veto;
TH2F* eff_full_muon_looseid;
TH2F* eff_full_muon_mediumid;
TH2F* eff_full_muon_miniiso04;
TH2F* eff_full_muon_miniiso02;
TH2F* eff_full_muon_looseip2d;
TH2F* eff_full_muon_tightip2d;
TH2D* eff_fast_muon_looseid;
TH2D* eff_fast_muon_mediumid;
TH2D* eff_fast_muon_miniiso04;
TH2D* eff_fast_muon_miniiso02;
TH2D* eff_fast_muon_looseip2d;
TH2D* eff_fast_muon_tightip2d;
TH2F* eff_full_muon_veto;

//TGraphAsymmErrors* eff_trigger;
TH2D* eff_trigger_veto;
TH2D* eff_trigger_veto_up;
TH2D* eff_trigger_veto_down;
TH2D* eff_trigger_ele;
TH2D* eff_trigger_ele_up;
TH2D* eff_trigger_ele_down;
TH2D* eff_trigger_pho;
TH2D* eff_trigger_pho_up;
TH2D* eff_trigger_pho_down;
TH2D* eff_trigger_mu;
TH2D* eff_trigger_mu_up;
TH2D* eff_trigger_mu_down;
TH2D* eff_trigger_F_met;
TH2D* eff_trigger_F_mu;
TH2D* eff_trigger_F_ele;
TH2D* eff_trigger_F_pho;

TGraphAsymmErrors* eff_full_fake_bW;
TGraphAsymmErrors* eff_full_fake_eW;
TGraphAsymmErrors* eff_full_fake_baW;
TGraphAsymmErrors* eff_full_fake_eaW;
TGraphAsymmErrors* eff_full_fake_bm0bW;
TGraphAsymmErrors* eff_full_fake_em0bW;
TGraphAsymmErrors* eff_full_fake_bTop;
TGraphAsymmErrors* eff_full_fake_eTop;
TGraphAsymmErrors* eff_full_fake_baTop;
TGraphAsymmErrors* eff_full_fake_eaTop;
TGraphAsymmErrors* eff_full_fake_bmTop;
TGraphAsymmErrors* eff_full_fake_emTop;
TGraphAsymmErrors* eff_full_fake_bm0bTop;
TGraphAsymmErrors* eff_full_fake_em0bTop;
TGraphAsymmErrors* eff_fast_bW;
TGraphAsymmErrors* eff_fast_eW;
TGraphAsymmErrors* eff_fast_bTop;
TGraphAsymmErrors* eff_fast_eTop;
TGraphAsymmErrors* eff_fast_fake_bW;
TGraphAsymmErrors* eff_fast_fake_eW;
TGraphAsymmErrors* eff_fast_fake_bTop;
TGraphAsymmErrors* eff_fast_fake_eTop;

void AnalysisBase::init_syst_input() {
  TString Sample(sample);
  
  // B-tagging
  // Efficiencies (Oct31 - test)
  TFile* f;
  if (Sample.Contains("FastSim"))
    f = TFile::Open("btag_eff/May19_withLepJets/FastSim_SMS-T5ttcc.root");
  else if (Sample.Contains("WJetsToLNu")) 
    f = TFile::Open("btag_eff/May19_withLepJets/WJetsToLNu.root");
  else if (Sample.Contains("TT")||Sample.Contains("ST")) 
    f = TFile::Open("btag_eff/May19_withLepJets/TT_powheg-pythia8.root");
  else 
    f = TFile::Open("btag_eff/May19_withLepJets/QCD.root");
  eff_btag_b_loose  = ((TH2D*)f->Get("btag_eff_b_loose"))->ProfileX();
  eff_btag_c_loose  = ((TH2D*)f->Get("btag_eff_c_loose"))->ProfileX();
  eff_btag_l_loose  = ((TH2D*)f->Get("btag_eff_l_loose"))->ProfileX();
  eff_btag_b_medium = ((TH2D*)f->Get("btag_eff_b_medium"))->ProfileX();
  eff_btag_c_medium = ((TH2D*)f->Get("btag_eff_c_medium"))->ProfileX();
  eff_btag_l_medium = ((TH2D*)f->Get("btag_eff_l_medium"))->ProfileX();
  eff_btag_b_loose  ->SetDirectory(0);
  eff_btag_c_loose  ->SetDirectory(0);
  eff_btag_l_loose  ->SetDirectory(0);
  eff_btag_b_medium ->SetDirectory(0);
  eff_btag_c_medium ->SetDirectory(0);
  eff_btag_l_medium ->SetDirectory(0);
  f->Close();
  // Moriond17 SFs
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco?rev=14#Supported_Algorithms_and_Operati
  // Summer16 FullSim
  btag_calib_full_ =  new BTagCalibration("csvv2", "scale_factors/btag/CSVv2_Moriond17_B_H.csv");
  // Loose WP
  btag_sf_full_loose_  = new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
  btag_sf_full_loose_->load(*btag_calib_full_, BTagEntry::FLAV_B,    "comb");
  btag_sf_full_loose_->load(*btag_calib_full_, BTagEntry::FLAV_C,    "comb");
  btag_sf_full_loose_->load(*btag_calib_full_, BTagEntry::FLAV_UDSG, "incl");
  // Medium WP
  btag_sf_full_medium_ = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
  btag_sf_full_medium_->load(*btag_calib_full_, BTagEntry::FLAV_B,    "comb");
  btag_sf_full_medium_->load(*btag_calib_full_, BTagEntry::FLAV_C,    "comb");
  btag_sf_full_medium_->load(*btag_calib_full_, BTagEntry::FLAV_UDSG, "incl");
  // Spring16 FastSim
  // This file needed minor formatting to be readable
  // sed 's;^";;;s; "\;;;;s;"";";g;' scale_factors/btag/fastsim_csvv2_ttbar_26_1_2017.csv
  btag_calib_fast_ =  new BTagCalibration("csvv2", "scale_factors/btag/fastsim_csvv2_ttbar_26_1_2017_fixed.csv");
  // Loose WP
  btag_sf_fast_loose_  = new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
  btag_sf_fast_loose_->load(*btag_calib_fast_, BTagEntry::FLAV_B,    "fastsim");
  btag_sf_fast_loose_->load(*btag_calib_fast_, BTagEntry::FLAV_C,    "fastsim");
  btag_sf_fast_loose_->load(*btag_calib_fast_, BTagEntry::FLAV_UDSG, "fastsim");
  // Medium WP
  btag_sf_fast_medium_ = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
  btag_sf_fast_medium_->load(*btag_calib_fast_, BTagEntry::FLAV_B,    "fastsim");
  btag_sf_fast_medium_->load(*btag_calib_fast_, BTagEntry::FLAV_C,    "fastsim");
  btag_sf_fast_medium_->load(*btag_calib_fast_, BTagEntry::FLAV_UDSG, "fastsim");

  // SoftDrop Mass correction for W tagging - Spring
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=43#Recipes_to_obtain_the_PUPPI_soft
  // Moriond17+ReReco
  TFile* file = TFile::Open("scale_factors/softdrop_mass_corr/puppiCorr.root");
  puppisd_corrGEN_      = (TF1*)((TF1*)file->Get("puppiJECcorr_gen"))->Clone();
  puppisd_corrRECO_cen_ = (TF1*)((TF1*)file->Get("puppiJECcorr_reco_0eta1v3"))->Clone();
  puppisd_corrRECO_for_ = (TF1*)((TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5"))->Clone();
  file->Close();

  // Lepton scale factors
  // Ele - Reconstruction  SF - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=38#Electron_efficiencies_and_scale
  eff_full_ele_reco                 = utils::getplot_TH2F("scale_factors/electron/reco/egammaEffi.txt_EGM2D.root","EGamma_SF2D", "ele1");
  // Ele - Data-FullSim    SF - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#Data_leading_order_FullSim_MC_co
  eff_full_ele_vetoid               = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","GsfElectronToCutBasedSpring15V", "ele2");
  eff_full_ele_looseid              = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","GsfElectronToCutBasedSpring15L", "ele3");
  eff_full_ele_mediumid             = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","GsfElectronToCutBasedSpring15M", "ele4");
  eff_full_ele_mvalooseid_tightip2d = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","GsfElectronToMVAVLooseTightIP2D","ele5");
  eff_full_ele_miniiso01            = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","MVAVLooseElectronToMini",        "ele6");
  eff_full_ele_miniiso02            = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","MVAVLooseElectronToMini2",       "ele7");
  eff_full_ele_miniiso04            = utils::getplot_TH2F("scale_factors/electron/fullsim/scaleFactors.root","MVAVLooseElectronToMini4",       "ele8");
  // Ele - FullSim-FastSim SF - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#FullSim_FastSim_TTBar_MC_compari
  eff_fast_ele_vetoid               = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_vetoCB.root",  "histo2D", "ele9");
  eff_fast_ele_looseid              = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_looseCB.root", "histo2D", "ele10");
  eff_fast_ele_mediumid             = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_mediumCB.root","histo2D", "ele11");
  eff_fast_ele_mvalooseid_tightip2d = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_vloose.root",  "histo2D", "ele12");
  eff_fast_ele_miniiso01            = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_mini01.root",  "histo2D", "ele13");  
  eff_fast_ele_miniiso02            = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_mini02.root",  "histo2D", "ele14");  
  eff_fast_ele_miniiso04            = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_mini04.root",  "histo2D", "ele15");  
  // Inclusive Razor Scale Factors
  eff_full_ele_veto                 = utils::getplot_TH2F("scale_factors/RazorRunAuxFiles_Expanded/"
							  "efficiency_results_VetoElectronSelectionEffDenominatorGen_2016_Rereco_Golden.root",
							  "ScaleFactor_VetoElectronSelectionEffDenominatorGen", "ele16");

  // Muon Tracking eff     SF - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#FullSim_FastSim_TTBar_MC_com_AN1
  eff_full_muon_trk   		    = utils::getplot_TGraphAsymmErrors("scale_factors/muon/tracking/Tracking_EfficienciesAndSF_BCDEFGH.root", "ratio_eff_eta3_dr030e030_corr", "mu1");
  eff_full_muon_trk_veto	    = utils::getplot_TGraphAsymmErrors("scale_factors/muon/tracking/Tracking_EfficienciesAndSF_BCDEFGH.root", "ratio_eff_eta3_tk0_dr030e030_corr", "mu2");
  // Muon Data-FullSim     SF - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#Data_leading_order_FullSim_M_AN1
  eff_full_muon_looseid		    = utils::getplot_TH2F("scale_factors/muon/fullsim/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root", "SF", "mu3");
  eff_full_muon_mediumid	    = utils::getplot_TH2F("scale_factors/muon/fullsim/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root","SF", "mu4");
  eff_full_muon_miniiso04	    = utils::getplot_TH2F("scale_factors/muon/fullsim/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root",  "SF", "mu5");
  eff_full_muon_miniiso02	    = utils::getplot_TH2F("scale_factors/muon/fullsim/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root", "SF", "mu6");
  eff_full_muon_looseip2d	    = utils::getplot_TH2F("scale_factors/muon/fullsim/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root",    "SF", "mu7");
  eff_full_muon_tightip2d	    = utils::getplot_TH2F("scale_factors/muon/fullsim/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root",    "SF", "mu8");
  // Muon FullSim-FastSim  SF - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#FullSim_FastSim_TTBar_MC_com_AN1
  eff_fast_muon_looseid		    = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_looseID.root",            "histo2D", "mu9");
  eff_fast_muon_mediumid	    = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_mediumID.root",           "histo2D", "mu10");
  eff_fast_muon_miniiso04	    = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_looseID_mini04.root",     "histo2D", "mu11");
  eff_fast_muon_miniiso02	    = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_mediumID_mini02.root",    "histo2D", "mu12");
  eff_fast_muon_looseip2d	    = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_mediumID_looseIP2D.root", "histo2D", "mu13");
  eff_fast_muon_tightip2d           = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_mediumID_tightIP2D.root", "histo2D", "mu14");
  // Inclusive Razor Scale Factors
  eff_full_muon_veto                = utils::getplot_TH2F("scale_factors/RazorRunAuxFiles_Expanded/"
							  "efficiency_results_VetoMuonSelectionEffDenominatorGen_2016_Rereco_Golden.root",
							  "ScaleFactor_VetoMuonSelectionEffDenominatorGen", "mu15");

  // 1D Trigger efficiency
  // TH1D* pass  = utils::getplot_TH1D("trigger_eff/Dec02_Golden_JSON/SingleLepton.root", "trigger_pass",  "trig1");
  // TH1D* total = utils::getplot_TH1D("trigger_eff/Dec02_Golden_JSON/SingleLepton.root", "trigger_total", "trig2");
  // eff_trigger = new TGraphAsymmErrors(pass, total);

  // 2D Trigger Efficiency (New) - Use combination of SingleElectron + MET datasets
#if USE_ISO_TRK_VETO > 0
  TH2D* veto_pass_2d  = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/MET.root",            "trigger2d_pass",   "trig1");
  TH2D* veto_total_2d = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/MET.root",            "trigger2d_total",  "trig2");
  TH2D* ele_pass_2d   = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/SingleElectron.root", "trigger2d_pass",   "trig3");
  TH2D* ele_total_2d  = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/SingleElectron.root", "trigger2d_total",  "trig4");
  TH2D* mu_pass_2d    = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/SingleMuon.root",     "trigger2d_pass",   "trig5");
  TH2D* mu_total_2d   = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/SingleMuon.root",     "trigger2d_total",  "trig6");
  TH2D* pho_pass_2d   = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/SingleElectron.root", "trigger2d_pass",   "trig7");
  TH2D* pho_total_2d  = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON_IsoTrkVeto/SingleElectron.root", "trigger2d_total",  "trig8");
#else
  TH2D* veto_pass_2d  = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/MET.root",            "trigger2d_pass",   "trig1");
  TH2D* veto_total_2d = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/MET.root",            "trigger2d_total",  "trig2");
  TH2D* ele_pass_2d   = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/SingleElectron.root", "trigger2d_pass",   "trig3");
  TH2D* ele_total_2d  = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/SingleElectron.root", "trigger2d_total",  "trig4");
  TH2D* mu_pass_2d    = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/SingleMuon.root",     "trigger2d_pass",   "trig5");
  TH2D* mu_total_2d   = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/SingleMuon.root",     "trigger2d_total",  "trig6");
  TH2D* pho_pass_2d   = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/SingleElectron.root", "trigger2d_pass",   "trig7");
  TH2D* pho_total_2d  = utils::getplot_TH2D("trigger_eff/Dec02_Golden_JSON/SingleElectron.root", "trigger2d_total",  "trig8");
#endif
  eff_trigger_veto      = (TH2D*)veto_total_2d->Clone("eff_trigger_veto");      eff_trigger_veto     ->Reset();
  eff_trigger_veto_up   = (TH2D*)veto_total_2d->Clone("eff_trigger_veto_up");   eff_trigger_veto_up  ->Reset();
  eff_trigger_veto_down = (TH2D*)veto_total_2d->Clone("eff_trigger_veto_down"); eff_trigger_veto_down->Reset();
  eff_trigger_ele       = (TH2D*)ele_total_2d ->Clone("eff_trigger_ele");       eff_trigger_ele      ->Reset();
  eff_trigger_ele_up    = (TH2D*)ele_total_2d ->Clone("eff_trigger_ele_up");    eff_trigger_ele_up   ->Reset();
  eff_trigger_ele_down  = (TH2D*)ele_total_2d ->Clone("eff_trigger_ele_down");  eff_trigger_ele_down ->Reset();
  eff_trigger_mu        = (TH2D*)mu_total_2d  ->Clone("eff_trigger_mu");        eff_trigger_mu       ->Reset();
  eff_trigger_mu_up     = (TH2D*)mu_total_2d  ->Clone("eff_trigger_mu_up");     eff_trigger_mu_up    ->Reset();
  eff_trigger_mu_down   = (TH2D*)mu_total_2d  ->Clone("eff_trigger_mu_down");   eff_trigger_mu_down  ->Reset();
  eff_trigger_pho       = (TH2D*)pho_total_2d ->Clone("eff_trigger_pho");       eff_trigger_pho      ->Reset();
  eff_trigger_pho_up    = (TH2D*)pho_total_2d ->Clone("eff_trigger_pho_up");    eff_trigger_pho_up   ->Reset();
  eff_trigger_pho_down  = (TH2D*)pho_total_2d ->Clone("eff_trigger_pho_down");  eff_trigger_pho_down ->Reset();
  for (int i=1; i<veto_total_2d->GetNbinsX()+1; i++) for (int j=1; j<veto_total_2d->GetNbinsY()+1; j++) {
    int veto_pass = veto_pass_2d->GetBinContent(i,j), veto_total = veto_total_2d->GetBinContent(i,j);
    if (veto_total>0) {
      TH1D p("p","",1,0,1); p.SetBinContent(1,veto_pass);  p.SetBinError(1,std::sqrt(veto_pass));
      TH1D t("t","",1,0,1); t.SetBinContent(1,veto_total); t.SetBinError(1,std::sqrt(veto_total));
      double eff = 0, err_up = 0, err_down = 0;
      utils::geteff_AE(TGraphAsymmErrors(&p,&t), 0, eff, err_up, err_down);
      //std::cout<<"Trigger efficiency: "<<i<<" "<<j<<" "<<eff-err_down<<" "<<eff<<" "<<eff+err_up<<std::endl;
      eff_trigger_veto     ->SetBinContent(i,j,eff);
      eff_trigger_veto_up  ->SetBinContent(i,j,eff+err_up);
      eff_trigger_veto_down->SetBinContent(i,j,eff-err_down);
      // SPECIAL: Set error to the total counts, so we know if a bin is not empty
      eff_trigger_veto     ->SetBinError(i,j,veto_total);
    }
    int ele_pass = ele_pass_2d->GetBinContent(i,j), ele_total = ele_total_2d->GetBinContent(i,j);
    if (ele_total>0) {
      TH1D p("p","",1,0,1); p.SetBinContent(1,ele_pass);  p.SetBinError(1,std::sqrt(ele_pass));
      TH1D t("t","",1,0,1); t.SetBinContent(1,ele_total); t.SetBinError(1,std::sqrt(ele_total));
      double eff = 0, err_up = 0, err_down = 0;
      utils::geteff_AE(TGraphAsymmErrors(&p,&t), 0, eff, err_up, err_down);
      //std::cout<<"Trigger efficiency: "<<i<<" "<<j<<" "<<eff-err_down<<" "<<eff<<" "<<eff+err_up<<std::endl;
      eff_trigger_ele     ->SetBinContent(i,j,eff);
      eff_trigger_ele_up  ->SetBinContent(i,j,eff+err_up);
      eff_trigger_ele_down->SetBinContent(i,j,eff-err_down);
      // SPECIAL: Set error to the total counts, so we know if a bin is not empty
      eff_trigger_ele     ->SetBinError(i,j,ele_total);
    }
    int pho_pass = pho_pass_2d->GetBinContent(i,j), pho_total = pho_total_2d->GetBinContent(i,j);
    if (pho_total>0) {
      TH1D p("p","",1,0,1); p.SetBinContent(1,pho_pass);  p.SetBinError(1,std::sqrt(pho_pass));
      TH1D t("t","",1,0,1); t.SetBinContent(1,pho_total); t.SetBinError(1,std::sqrt(pho_total));
      double eff = 0, err_up = 0, err_down = 0;
      utils::geteff_AE(TGraphAsymmErrors(&p,&t), 0, eff, err_up, err_down);
      //std::cout<<"Trigger efficiency: "<<i<<" "<<j<<" "<<eff-err_down<<" "<<eff<<" "<<eff+err_up<<std::endl;
      eff_trigger_pho     ->SetBinContent(i,j,eff);
      eff_trigger_pho_up  ->SetBinContent(i,j,eff+err_up);
      eff_trigger_pho_down->SetBinContent(i,j,eff-err_down);
      // SPECIAL: Set error to the total counts, so we know if a bin is not empty
      eff_trigger_pho     ->SetBinError(i,j,pho_total);
    }
    int mu_pass = mu_pass_2d->GetBinContent(i,j), mu_total = mu_total_2d->GetBinContent(i,j);
    if (mu_total>0) {
      TH1D p("p","",1,0,1); p.SetBinContent(1,mu_pass);  p.SetBinError(1,std::sqrt(mu_pass));
      TH1D t("t","",1,0,1); t.SetBinContent(1,mu_total); t.SetBinError(1,std::sqrt(mu_total));
      double eff = 0, err_up = 0, err_down = 0;
      utils::geteff_AE(TGraphAsymmErrors(&p,&t), 0, eff, err_up, err_down);
      //std::cout<<"Trigger efficiency: "<<i<<" "<<j<<" "<<eff-err_down<<" "<<eff<<" "<<eff+err_up<<std::endl;
      eff_trigger_mu     ->SetBinContent(i,j,eff);
      eff_trigger_mu_up  ->SetBinContent(i,j,eff+err_up);
      eff_trigger_mu_down->SetBinContent(i,j,eff-err_down);
      // SPECIAL: Set error to the total counts, so we know if a bin is not empty
      eff_trigger_mu     ->SetBinError(i,j,mu_total);
    }
  }
  // Same trigger efficiencies but in the F region (needed for fake rates)
  const char* fin = "trigger_eff/Dec02_Golden_JSON/F_Region.root";
  eff_trigger_F_met = utils::getplot_TH2D(fin, "met", "trig_f_met");
  eff_trigger_F_mu  = utils::getplot_TH2D(fin, "mu",  "trig_f_mu");
  eff_trigger_F_ele = utils::getplot_TH2D(fin, "ele", "trig_f_ele");
  eff_trigger_F_pho = utils::getplot_TH2D(fin, "pho", "trig_f_pho");

  // W/Top (anti-/mass-)tag (and fake rate) scale factors
  // From Changgi
  //    eff_full_fake_bW      = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "bW",      "full_fake_W_barrel");
  //    eff_full_fake_eW      = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "eW",      "full_fake_W_endcap");
  //    eff_full_fake_bm0bW   = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "bm0bW",   "full_fake_m0bW_barrel");
  //    eff_full_fake_em0bW   = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "em0bW",   "full_fake_m0bW_endcap");
  //    eff_full_fake_baW     = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "baW",     "full_fake_aW_barrel");
  //    eff_full_fake_eaW     = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "eaW",     "full_fake_aW_endcap");
  //    eff_full_fake_bTop    = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "bTop",    "full_fake_Top_barrel");
  //    eff_full_fake_eTop    = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "eTop",    "full_fake_Top_endcap");
  //    eff_full_fake_bmTop   = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "bmTop",   "full_fake_mTop_barrel");
  //    eff_full_fake_emTop   = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "emTop",   "full_fake_mTop_endcap");
  //    eff_full_fake_bm0bTop = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "bm0bTop", "full_fake_0bmTop_barrel");
  //    eff_full_fake_em0bTop = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "em0bTop", "full_fake_0bmTop_endcap");
  //    eff_full_fake_baTop   = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "baTop",   "full_fake_aTop_barrel");
  //    eff_full_fake_eaTop   = utils::getplot_TH1D("scale_factors/w_top_tag/WTopTagSF.root",                             "eaTop",   "full_fake_aTop_endcap");
  //    eff_fast_bW           = utils::getplot_TH1D("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap.root", "bWFF",    "fast_bW");
  //    eff_fast_eW           = utils::getplot_TH1D("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap.root", "eWFF",    "fast_eW");
  //    eff_fast_bTop         = utils::getplot_TH1D("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap.root", "bTopFF",  "fast_bTop");
  //    eff_fast_eTop         = utils::getplot_TH1D("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap.root", "eTopFF",  "fast_eTop");
  // From Janos
#if USE_ISO_TRK_VETO > 0
  eff_full_fake_bW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "bW",      "full_fake_W_barrel");
  eff_full_fake_eW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "eW",      "full_fake_W_endcap");
  eff_full_fake_bm0bW   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "bm0bW",   "full_fake_m0bW_barrel");
  eff_full_fake_em0bW   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "em0bW",   "full_fake_m0bW_endcap");
  eff_full_fake_baW     = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "baW",     "full_fake_aW_barrel");
  eff_full_fake_eaW     = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "eaW",     "full_fake_aW_endcap");
  eff_full_fake_bTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "bTop",    "full_fake_Top_barrel");
  eff_full_fake_eTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "eTop",    "full_fake_Top_endcap");
  eff_full_fake_bmTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "bmTop",   "full_fake_mTop_barrel");
  eff_full_fake_emTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "emTop",   "full_fake_mTop_endcap");
  eff_full_fake_bm0bTop = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "bm0bTop", "full_fake_0bmTop_barrel");
  eff_full_fake_em0bTop = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "em0bTop", "full_fake_0bmTop_endcap");
  eff_full_fake_baTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "baTop",   "full_fake_aTop_barrel");
  eff_full_fake_eaTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos_IsoTrkVeto.root",                             "eaTop",   "full_fake_aTop_endcap");
  eff_fast_bW           = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "bWFF",    "fast_bW");
  eff_fast_eW           = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "eWFF",    "fast_eW");
  eff_fast_bTop         = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "bTopFF",  "fast_bTop");
  eff_fast_eTop         = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "eTopFF",  "fast_eTop");
  eff_fast_fake_bW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "bMWFF",   "fast_fake_bW");
  eff_fast_fake_eW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "eMWFF",   "fast_fake_eW");
  eff_fast_fake_bTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "bMTopFF", "fast_fake_bTop");
  eff_fast_fake_eTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos_IsoTrkVeto.root", "eMTopFF", "fast_fake_eTop");
#else
  eff_full_fake_bW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "bW",      "full_fake_W_barrel");
  eff_full_fake_eW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "eW",      "full_fake_W_endcap");
  eff_full_fake_bm0bW   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "bm0bW",   "full_fake_m0bW_barrel");
  eff_full_fake_em0bW   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "em0bW",   "full_fake_m0bW_endcap");
  eff_full_fake_baW     = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "baW",     "full_fake_aW_barrel");
  eff_full_fake_eaW     = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "eaW",     "full_fake_aW_endcap");
  eff_full_fake_bTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "bTop",    "full_fake_Top_barrel");
  eff_full_fake_eTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "eTop",    "full_fake_Top_endcap");
  eff_full_fake_bmTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "bmTop",   "full_fake_mTop_barrel");
  eff_full_fake_emTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "emTop",   "full_fake_mTop_endcap");
  eff_full_fake_bm0bTop = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "bm0bTop", "full_fake_0bmTop_barrel");
  eff_full_fake_em0bTop = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "em0bTop", "full_fake_0bmTop_endcap");
  eff_full_fake_baTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "baTop",   "full_fake_aTop_barrel");
  eff_full_fake_eaTop   = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/WTopTagSF_Janos.root",                             "eaTop",   "full_fake_aTop_endcap");
  eff_fast_bW           = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "bWFF",    "fast_bW");
  eff_fast_eW           = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "eWFF",    "fast_eW");
  eff_fast_bTop         = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "bTopFF",  "fast_bTop");
  eff_fast_eTop         = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "eTopFF",  "fast_eTop");
  eff_fast_fake_bW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "bMWFF",   "fast_fake_bW");
  eff_fast_fake_eW      = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "eMWFF",   "fast_fake_eW");
  eff_fast_fake_bTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "bMTopFF", "fast_fake_bTop");
  eff_fast_fake_eTop    = utils::getplot_TGraphAsymmErrors("scale_factors/w_top_tag/fastsim/FullFastSimTagSF_BarrelEndcap_Janos.root", "eMTopFF", "fast_fake_eTop");
#endif
}


double AnalysisBase::calc_top_tagging_sf(DataStruct& data, const double& nSigmaTopTagSF, const double& nSigmaTopTagFastSimSF,
					 const double& nSigmaTopMisTagSF, const double& nSigmaTopMisTagFastSimSF, const bool& isFastSim) {
  double w = 1;
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (passHadTopTag[i]) {
      // Gen-matched tags
      if (hasGenTop[i]) {
	// Use POG scale factor for tag
	w *= get_syst_weight(TOP_TAG_SF, TOP_TAG_SF+TOP_TAG_SF_ERR_UP, TOP_TAG_SF-TOP_TAG_SF_ERR_DOWN, nSigmaTopTagSF);
	// Additionally use our scale factors for FastSim
	if (isFastSim) {
	  //double eff, err;
	  double eff, err_up, err_down;
	  if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	    //utils::geteff1D(eff_fast_bTop, data.jetsAK8.Pt[i], eff, err);
	    utils::geteff_AE(eff_fast_bTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
	  } else {
	    //utils::geteff1D(eff_fast_eTop, data.jetsAK8.Pt[i], eff, err);
	    utils::geteff_AE(eff_fast_eTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
	  }
	  //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaTopTagFastSimSF);
	  w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTopTagFastSimSF);
	}
      } else {
	// Mis-tags: Top tagging fake rate scale factor
	//double eff, err;
	double eff, err_up, err_down;
	if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	  //utils::geteff1D(eff_full_fake_bTop, data.jetsAK8.Pt[i], eff, err);
	  utils::geteff_AE(eff_full_fake_bTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
	} else {
	  //utils::geteff1D(eff_full_fake_eTop, data.jetsAK8.Pt[i], eff, err);	  
	  utils::geteff_AE(eff_full_fake_eTop, data.jetsAK8.Pt[i], eff, err_up, err_down);	  
	}
	//w *= get_syst_weight(eff, eff+err, eff-err, nSigmaTopMisTagSF);
	w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTopMisTagSF);
	if (isFastSim) {
	  if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	    utils::geteff_AE(eff_fast_fake_bTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
	  } else {
	    utils::geteff_AE(eff_fast_fake_eTop, data.jetsAK8.Pt[i], eff, err_up, err_down);	  
	  }
	  w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTopMisTagFastSimSF);
	}
      }
    }
  }

  return w;
}

double AnalysisBase::calc_fake_top_0b_mass_tagging_sf(DataStruct& data, const double& nSigmaTop0BMassTagSF) {
  double w = 1;
  if (nGenTop==0) while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (passHadTop0BMassTag[i]) {
      //double eff, err;
      double eff, err_up, err_down;
      if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	//utils::geteff1D(eff_full_fake_bm0bTop, data.jetsAK8.Pt[i], eff, err);
	utils::geteff_AE(eff_full_fake_bm0bTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
      } else {
	//utils::geteff1D(eff_full_fake_em0bTop, data.jetsAK8.Pt[i], eff, err);	
	utils::geteff_AE(eff_full_fake_em0bTop, data.jetsAK8.Pt[i], eff, err_up, err_down);	
      }
      //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaTop0BMassTagSF);
      w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTop0BMassTagSF);
    }
  }

  return w;
}

double AnalysisBase::calc_fake_top_mass_tagging_sf(DataStruct& data, const double& nSigmaTopMassTagSF) {
  double w = 1;
  if (nGenTop==0) while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (passHadTop0BMassTag[i]) {
      //double eff, err;
      double eff, err_up, err_down;
      if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	//utils::geteff1D(eff_full_fake_bmTop, data.jetsAK8.Pt[i], eff, err);
	utils::geteff_AE(eff_full_fake_bmTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
      } else {
	//utils::geteff1D(eff_full_fake_emTop, data.jetsAK8.Pt[i], eff, err);	
	utils::geteff_AE(eff_full_fake_emTop, data.jetsAK8.Pt[i], eff, err_up, err_down);	
      }
      //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaTopMassTagSF);
      w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTopMassTagSF);
    }
  }

  return w;
}

double AnalysisBase::calc_fake_top_anti_tagging_sf(DataStruct& data, const double& nSigmaTopAntiTagSF) {
  double w = 1;
  if (nGenTop==0) while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (passHadTop0BAntiTag[i]) {
      //double eff, err;
      double eff, err_up, err_down;
      if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	//utils::geteff1D(eff_full_fake_baTop, data.jetsAK8.Pt[i], eff, err);
	utils::geteff_AE(eff_full_fake_baTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
      } else {
	//utils::geteff1D(eff_full_fake_eaTop, data.jetsAK8.Pt[i], eff, err);
	utils::geteff_AE(eff_full_fake_eaTop, data.jetsAK8.Pt[i], eff, err_up, err_down);
      }
      //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaTopAntiTagSF);
      w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTopAntiTagSF);
    }
  }

  return w;
}

double AnalysisBase::calc_w_tagging_sf(DataStruct& data, const double& nSigmaWTagSF, const double& nSigmaWTagFastSimSF, 
				       const double& nSigmaWMisTagSF, const double& nSigmaWMisTagFastSimSF, const bool& isFastSim) {
  double w = 1.0;

  // New method
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (passTightWTag[i]) {
      if (hasGenW[i]) {
        // Use POG scale factor for efficiency scale factor
        w *= get_syst_weight(W_TAG_HP_SF, W_TAG_HP_SF_ERR, nSigmaWTagSF);
        // Additionally use our scale factors for FastSim
        if (isFastSim) {
	  //double eff, err;
	  double eff, err_up, err_down;
          if (std::abs(data.jetsAK8.Eta[i])<1.5) {
            //utils::geteff1D(eff_fast_bW, data.jetsAK8.Pt[i], eff, err);
            utils::geteff_AE(eff_fast_bW, data.jetsAK8.Pt[i], eff, err_up, err_down);
          } else {
            //utils::geteff1D(eff_fast_eW, data.jetsAK8.Pt[i], eff, err);
            utils::geteff_AE(eff_fast_eW, data.jetsAK8.Pt[i], eff, err_up, err_down);
          }
          //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaWTagFastSimSF);
          w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaWTagFastSimSF);
        }
      } else  {
	// W tagging fake rate scale factor
	//double eff, err;
	double eff, err_up, err_down;
        if (std::abs(data.jetsAK8.Eta[i])<1.5) {
          //utils::geteff1D(eff_full_fake_bW, data.jetsAK8.Pt[i], eff, err);
          utils::geteff_AE(eff_full_fake_bW, data.jetsAK8.Pt[i], eff, err_up, err_down);
        } else {
          //utils::geteff1D(eff_full_fake_eW, data.jetsAK8.Pt[i], eff, err);
          utils::geteff_AE(eff_full_fake_eW, data.jetsAK8.Pt[i], eff, err_up, err_down);
        }
        //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaWMisTagSF);
        w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaWMisTagSF);
        if (isFastSim) {
	  if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	    utils::geteff_AE(eff_fast_fake_bW, data.jetsAK8.Pt[i], eff, err_up, err_down);
	  } else {
	    utils::geteff_AE(eff_fast_fake_eW, data.jetsAK8.Pt[i], eff, err_up, err_down);
	  }
	  w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaWMisTagFastSimSF);
	}
      }
    } else if (passTightWAntiTag[i]) {
      // Also adding low purity SFs from the POG
      if (hasGenW[i]) {
        w *= get_syst_weight(W_TAG_LP_SF, W_TAG_LP_SF_ERR, nSigmaWTagSF);
      }
    }
  }
  
  return w;
}

double AnalysisBase::calc_fake_w_mass_tagging_sf(DataStruct& data, const double& nSigmaWMassTagSF) {
  double w = 1.0;

  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (nGenHadW==0) {
      if (passWMassTag[i]) {
	//double eff, err;
	double eff, err_up, err_down;
	if (std::abs(data.jetsAK8.Eta[i])<1.5) {
	  //utils::geteff1D(eff_full_fake_bm0bW, data.jetsAK8.Pt[i], eff, err);
          utils::geteff_AE(eff_full_fake_bm0bW, data.jetsAK8.Pt[i], eff, err_up, err_down);
	} else {
	  //utils::geteff1D(eff_full_fake_em0bW, data.jetsAK8.Pt[i], eff, err);
          utils::geteff_AE(eff_full_fake_em0bW, data.jetsAK8.Pt[i], eff, err_up, err_down);
	}
	//w *= get_syst_weight(eff, eff+err, eff-err, nSigmaWMassTagSF);
	w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaWMassTagSF);
      }
    }
  }
  
  return w;
}

double AnalysisBase::calc_fake_w_anti_tagging_sf(DataStruct& data, const double& nSigmaWAntiTagSF) {
  double w = 1.0;

  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    if (passTightWAntiTag[i]) {
      //double eff, err;
      double eff, err_up, err_down;
      if (std::abs(data.jetsAK8.Eta[i])<1.5) {      
	//utils::geteff1D(eff_full_fake_baW, data.jetsAK8.Pt[i], eff, err);
        utils::geteff_AE(eff_full_fake_baW, data.jetsAK8.Pt[i], eff, err_up, err_down);
      } else {
	//utils::geteff1D(eff_full_fake_eaW, data.jetsAK8.Pt[i], eff, err);
        utils::geteff_AE(eff_full_fake_eaW, data.jetsAK8.Pt[i], eff, err_up, err_down);
      }
      //w *= get_syst_weight(eff, eff+err, eff-err, nSigmaWAntiTagSF);
      w *= get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaWAntiTagSF);
    }
  }
  
  return w;
}


std::pair<double, double> AnalysisBase::calc_b_tagging_sf(DataStruct& data, const double& nSigmaBTagSF, const double& nSigmaBTagFastSimSF,  const bool& isFastSim) {

  double pMC_loose = 1, pData_loose = 1;
  double pMC_medium = 1, pData_medium = 1;
  while(data.jetsAK4.Loop()) {
    size_t i = data.jetsAK4.it;
    float pt = data.jetsAK4.Pt[i], eta = data.jetsAK4.Eta[i];
    // Jet ID
    if (passLooseJet[i]) {

      // Btag efficiencies (quark flavour dependent)
      BTagEntry::JetFlavor FLAV;
      double eff_medium = 1.0, eff_loose = 1.0;
      if (data.jetsAK4.HadronFlavour[i]==5) {
	FLAV = BTagEntry::FLAV_B;
	eff_loose  = utils::geteff1D(eff_btag_b_loose,  pt);
	eff_medium = utils::geteff1D(eff_btag_b_medium, pt);
      } else if (data.jetsAK4.HadronFlavour[i]==4) {
	FLAV = BTagEntry::FLAV_C;
	eff_loose  = utils::geteff1D(eff_btag_c_loose,  pt);
	eff_medium = utils::geteff1D(eff_btag_c_medium, pt);
      } else {
	FLAV = BTagEntry::FLAV_UDSG;
	eff_loose  = utils::geteff1D(eff_btag_l_loose,  pt);
	eff_medium = utils::geteff1D(eff_btag_l_medium, pt);
      }
      
      // Scale factors - FullSim
      double sf_loose_cen   = btag_sf_full_loose_ ->eval_auto_bounds("central", FLAV, eta, pt); 
      double sf_loose_up    = btag_sf_full_loose_ ->eval_auto_bounds("up",      FLAV, eta, pt);
      double sf_loose_down  = btag_sf_full_loose_ ->eval_auto_bounds("down",    FLAV, eta, pt); 
      double sf_medium_cen  = btag_sf_full_medium_->eval_auto_bounds("central", FLAV, eta, pt); 
      double sf_medium_up   = btag_sf_full_medium_->eval_auto_bounds("up",      FLAV, eta, pt);
      double sf_medium_down = btag_sf_full_medium_->eval_auto_bounds("down",    FLAV, eta, pt); 
      
      double sf_loose       = get_syst_weight(sf_loose_cen,  sf_loose_up,  sf_loose_down,  nSigmaBTagSF);
      double sf_medium      = get_syst_weight(sf_medium_cen, sf_medium_up, sf_medium_down, nSigmaBTagSF);
      
      // FastSim
      if (isFastSim) {
	sf_loose_cen   = btag_sf_fast_loose_ ->eval_auto_bounds("central", FLAV, eta, pt); 
	sf_loose_up    = btag_sf_fast_loose_ ->eval_auto_bounds("up",      FLAV, eta, pt);
	sf_loose_down  = btag_sf_fast_loose_ ->eval_auto_bounds("down",    FLAV, eta, pt); 
	sf_medium_cen  = btag_sf_fast_medium_->eval_auto_bounds("central", FLAV, eta, pt); 
	sf_medium_up   = btag_sf_fast_medium_->eval_auto_bounds("up",      FLAV, eta, pt);
	sf_medium_down = btag_sf_fast_medium_->eval_auto_bounds("down",    FLAV, eta, pt); 

	sf_loose      *= get_syst_weight(sf_loose_cen,  sf_loose_up,  sf_loose_down,  nSigmaBTagFastSimSF);
	sf_medium     *= get_syst_weight(sf_medium_cen, sf_medium_up, sf_medium_down, nSigmaBTagFastSimSF);
      }
      
      // Working points
      if (passLooseBTag[i]) {
	pMC_loose   *= eff_loose;
	pData_loose *= eff_loose * sf_loose;
      } else {
	pMC_loose   *= 1 - eff_loose;
	pData_loose *= 1 - eff_loose * sf_loose;
      }
      
      if (passMediumBTag[i]) {
	pMC_medium   *= eff_medium;
	pData_medium *= eff_medium * sf_medium;
      } else {
	pMC_medium   *= 1 - eff_medium;
	pData_medium *= 1 - eff_medium * sf_medium;
      }
    }
  }
  double weight_loose  = pData_loose/pMC_loose;
  double weight_medium = pData_medium/pMC_medium;
  return std::make_pair(weight_loose, weight_medium);
}

std::tuple<double, double, double> AnalysisBase::calc_ele_sf(DataStruct& data, const double& nSigmaEleRecoSF, const double& nSigmaEleIDSF, const double& nSigmaEleIsoSF, const double& nSigmaEleFastSimSF,const bool& isFastSim) {
  double reco_sf, reco_sf_err, sf, sf_err;
  double weight_veto  = 1.0, weight_loose = 1.0, weight_select = 1.0;
  while(data.ele.Loop()) {
    size_t i       = data.ele.it;
    double pt      = data.ele.Pt[i];
    double eta     = data.ele.Eta[i];
    float abseta   = std::abs(eta);
    float miniIso  = data.ele.MiniIso[i];
    float absd0    = std::abs(data.ele.Dxy[i]);
    float absdz    = std::abs(data.ele.Dz[i]);
#if USE_MVA_ID == 1
#if USE_POG_ID > 0
    // https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=210#Electrons
    //VLoose WP, pT > 10 GeV. A/B values are -0.48/-0.85, -0.67/-0.91, -0.49/-0.83 for |eta| [0, 0.8], [0.8, 1.479], [1.479, 2.5]
    //  In addition, for pT 5-10, the following values are used on the HZZ MVA: 0.46, -0.03, 0.06 for |eta| [0, 0.8], [0.8, 1.479], [1.479, 2.5] 
    int categ = data.ele.vidMvaHZZcateg[i];
    double mva = categ<3 ? data.ele.vidMvaHZZvalue[i] : data.ele.vidMvaGPvalue[i];
    std::vector<double> veto_cut = { 
      // pt<10
      0.46, -0.03, 0.06, 
      // pt>=10
      std::min(-0.48, std::max(-0.85 , -0.48 + (-0.85 - -0.48)/10.0*(pt-15))),
      std::min(-0.67, std::max(-0.91 , -0.67 + (-0.91 - -0.67)/10.0*(pt-15))),
      std::min(-0.49, std::max(-0.83 , -0.49 + (-0.83 - -0.49)/10.0*(pt-15)))
    };
    bool id_veto_noiso = (mva>veto_cut[categ]);
#endif
#else
    bool id_veto_noiso = (data.ele.vidVetonoiso[i] == 1.0);
#endif
    bool id_loose_noiso  = (data.ele.vidLoosenoiso[i] == 1.0);
    bool id_select_noiso = (data.ele.vidMediumnoiso[i] == 1.0);
    // Apply reconstruction scale factor - Warning! strange binning (pt vs eta)
    utils::geteff2D(eff_full_ele_reco, eta, pt, reco_sf, reco_sf_err);
    // If pt is below 20 or above 80 GeV increase error by 1%
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=38#Electron_efficiencies_and_scale
    if (pt<20||pt>=80) reco_sf_err = std::sqrt(reco_sf_err*reco_sf_err + 0.01*0.01);

    // For FullSim scale factors, we apply syst error from each bin separately
    // For FastSim scale factors, we apply a 2% error (per electron leg)

    // Veto Electrons
#if USE_MVA_ID == 1
#if USE_POG_ID > 0
    // Using previous POG ID scale factors
    // Apply ID + IP scale factor
    if ( id_veto_noiso &&
	 pt      >= ELE_VETO_PT_CUT &&
	 abseta  <  ELE_VETO_ETA_CUT //&& !(abseta>=1.442 && abseta< 1.556) &&
	 absd0   <  ELE_VETO_IP_D0_CUT &&
	 absdz   <  ELE_VETO_IP_DZ_CUT ) {
      utils::geteff2D(eff_full_ele_mvalooseid_tightip2d, pt, eta, sf, sf_err);
      weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      if (isFastSim) {
	utils::geteff2D(eff_fast_ele_mvalooseid_tightip2d, pt, eta, sf, sf_err);
	weight_veto *= sf;
      }
      if ( miniIso <  ELE_VETO_MINIISO_CUT ) {
	// Apply Iso scale factor
	if (ELE_VETO_MINIISO_CUT == 0.1)
	  utils::geteff2D(eff_full_ele_miniiso01, pt, eta, sf, sf_err);
	else if (ELE_VETO_MINIISO_CUT == 0.2)
	  utils::geteff2D(eff_full_ele_miniiso02, pt, eta, sf, sf_err);
	else if (ELE_VETO_MINIISO_CUT == 0.4)
	  utils::geteff2D(eff_full_ele_miniiso04, pt, eta, sf, sf_err);
	weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIsoSF);
	if (isFastSim) {
	  if (ELE_VETO_MINIISO_CUT == 0.1)
	    utils::geteff2D(eff_fast_ele_miniiso01, pt, eta, sf, sf_err);
	  else if (ELE_VETO_MINIISO_CUT == 0.2)
	    utils::geteff2D(eff_fast_ele_miniiso02, pt, eta, sf, sf_err);
	  else if (ELE_VETO_MINIISO_CUT == 0.4)
	    utils::geteff2D(eff_fast_ele_miniiso04, pt, eta, sf, sf_err);
	  weight_veto *= sf;
	  // Apply 2% error per electron leg
	  weight_veto *= get_syst_weight(1, 0.02, nSigmaEleFastSimSF);
	}
	// Apply the Reco SF
	weight_veto   *= get_syst_weight(reco_sf, reco_sf_err, nSigmaEleRecoSF);
      }
    }
#else
    // Using Razor Inclusive ID
    // Apply Ele Reco + RazorInclusive ID scale factors
    if ( passEleVeto[i] ) {
      utils::geteff2D(eff_full_ele_veto, pt, eta, sf, sf_err);
      weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      // Apply the Reco SF
      weight_veto *= get_syst_weight(reco_sf, reco_sf_err, nSigmaEleRecoSF);
    }
#endif
#else
    // Apply ID scale factor
    if ( id_veto_noiso &&
	 pt      >= ELE_VETO_PT_CUT &&
	 abseta  <  ELE_VETO_ETA_CUT //&& !(abseta>=1.442 && abseta< 1.556) ) {
      utils::geteff2D(eff_full_ele_vetoid, pt, eta, sf, sf_err);
      weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      if (isFastSim) {
	utils::geteff2D(eff_fast_ele_vetoid, pt, eta, sf, sf_err);
	weight_veto *= sf;
      }
      if ( miniIso <  ELE_VETO_MINIISO_CUT &&
           absd0   <  ELE_VETO_IP_D0_CUT &&
           absdz   <  ELE_VETO_IP_DZ_CUT ) {
	// Apply Iso scale factor
	if (ELE_VETO_MINIISO_CUT == 0.1)
	  utils::geteff2D(eff_full_ele_miniiso01, pt, eta, sf, sf_err);
	else if (ELE_VETO_MINIISO_CUT == 0.2)
	  utils::geteff2D(eff_full_ele_miniiso02, pt, eta, sf, sf_err);
	else if (ELE_VETO_MINIISO_CUT == 0.4)
	  utils::geteff2D(eff_full_ele_miniiso04, pt, eta, sf, sf_err);
	weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIsoSF);
	if (isFastSim) {
	  if (ELE_VETO_MINIISO_CUT == 0.1)
	    utils::geteff2D(eff_fast_ele_miniiso01, pt, eta, sf, sf_err);
	  else if (ELE_VETO_MINIISO_CUT == 0.2)
	    utils::geteff2D(eff_fast_ele_miniiso02, pt, eta, sf, sf_err);
	  else if (ELE_VETO_MINIISO_CUT == 0.4)
	    utils::geteff2D(eff_fast_ele_miniiso04, pt, eta, sf, sf_err);
	  weight_veto *= sf;
	  // Apply 2% error per electron leg
	  weight_veto *= get_syst_weight(1, 0.02, nSigmaEleFastSimSF);
	}
	// Apply the Reco SF
	weight_veto   *= get_syst_weight(reco_sf, reco_sf_err, nSigmaEleRecoSF);
      }
    }
#endif
    
    // Loose Electrons
    if ( id_loose_noiso &&
	 pt      >= ELE_LOOSE_PT_CUT &&
	 abseta  <  ELE_LOOSE_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_ele_looseid, pt, eta, sf, sf_err);
      weight_loose *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      if (isFastSim) {
	utils::geteff2D(eff_fast_ele_looseid, pt, eta, sf, sf_err);
	weight_loose *= sf;
      }
      if ( miniIso <  ELE_LOOSE_MINIISO_CUT &&
	   absd0   <  ELE_LOOSE_IP_D0_CUT &&
	   absdz   <  ELE_LOOSE_IP_DZ_CUT ) {
	// Apply Iso scale factor
	if (ELE_LOOSE_MINIISO_CUT == 0.1)
	  utils::geteff2D(eff_full_ele_miniiso01, pt, eta, sf, sf_err);
	else if (ELE_LOOSE_MINIISO_CUT == 0.2)
	  utils::geteff2D(eff_full_ele_miniiso02, pt, eta, sf, sf_err);
	else if (ELE_LOOSE_MINIISO_CUT == 0.4)
	  utils::geteff2D(eff_full_ele_miniiso04, pt, eta, sf, sf_err);
	weight_loose *= get_syst_weight(sf, sf_err, nSigmaEleIsoSF);
	if (isFastSim) {
	  if (ELE_LOOSE_MINIISO_CUT == 0.1)
	    utils::geteff2D(eff_fast_ele_miniiso01, pt, eta, sf, sf_err);
	  else if (ELE_LOOSE_MINIISO_CUT == 0.2)
	    utils::geteff2D(eff_fast_ele_miniiso02, pt, eta, sf, sf_err);
	  else if (ELE_LOOSE_MINIISO_CUT == 0.4)
	    utils::geteff2D(eff_fast_ele_miniiso04, pt, eta, sf, sf_err);
	  weight_loose *= sf;
	  // Apply 2% error per electron leg
	  weight_loose *= get_syst_weight(1, 0.02, nSigmaEleFastSimSF);
	}
	// Apply the Reco SF
	weight_loose   *= get_syst_weight(reco_sf, reco_sf_err, nSigmaEleRecoSF);
      }
    }

    // Selected Electrons
    if ( id_select_noiso &&
	 pt      >= ELE_SELECT_PT_CUT &&
	 abseta  <  ELE_SELECT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_ele_mediumid, pt, eta, sf, sf_err);
      weight_select *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      if (isFastSim) {
	utils::geteff2D(eff_fast_ele_mediumid, pt, eta, sf, sf_err);
	weight_select *= sf;
      }
      if ( miniIso <  ELE_SELECT_MINIISO_CUT &&
	   absd0   <  ELE_SELECT_IP_D0_CUT &&
	   absdz   <  ELE_SELECT_IP_DZ_CUT ) {
	// Apply Iso scale factor
	if (ELE_SELECT_MINIISO_CUT == 0.1)
	  utils::geteff2D(eff_full_ele_miniiso01, pt, eta, sf, sf_err);
	else if (ELE_SELECT_MINIISO_CUT == 0.2)
	  utils::geteff2D(eff_full_ele_miniiso02, pt, eta, sf, sf_err);
	else if (ELE_SELECT_MINIISO_CUT == 0.4)
	  utils::geteff2D(eff_full_ele_miniiso04, pt, eta, sf, sf_err);
	weight_select *= get_syst_weight(sf, sf_err, nSigmaEleIsoSF);
	if (isFastSim) {
	  if (ELE_SELECT_MINIISO_CUT == 0.1)
	    utils::geteff2D(eff_fast_ele_miniiso01, pt, eta, sf, sf_err);
	  else if (ELE_SELECT_MINIISO_CUT == 0.2)
	    utils::geteff2D(eff_fast_ele_miniiso02, pt, eta, sf, sf_err);
	  else if (ELE_SELECT_MINIISO_CUT == 0.4)
	    utils::geteff2D(eff_fast_ele_miniiso04, pt, eta, sf, sf_err);
	  weight_select *= sf;
	  // Apply 2% error per electron leg
	  weight_select *= get_syst_weight(1, 0.02, nSigmaEleFastSimSF);
	}
	// Apply the Reco SF
	weight_select   *= get_syst_weight(reco_sf, reco_sf_err, nSigmaEleRecoSF);
      }
    }

  }
  
  return std::make_tuple(weight_veto, weight_loose, weight_select);
}

std::tuple<double, double, double> AnalysisBase::calc_muon_sf(DataStruct& data, const double& nSigmaMuonTrkSF, const double& nSigmaMuonFullSimSF, const double& nSigmaMuonFastSimSF, const bool& isFastSim) {
  double trk_sf, trk_sf_err_up, trk_sf_err_down, trk_sf_veto, trk_sf_veto_err_up, trk_sf_veto_err_down, sf, sf_err;
  double weight_veto  = 1.0, weight_loose = 1.0, weight_select = 1.0;
  while(data.mu.Loop()) {
    size_t i       = data.mu.it;
    double pt      = data.mu.Pt[i];
    double eta     = data.mu.Eta[i];
    float abseta   = std::abs(eta);
    float miniIso  = data.mu.MiniIso[i];
#if USE_POG_ID > 0
    bool id_veto_noiso   = (data.mu.IsLooseMuon[i]  == 1.0);
#endif
    float absd0    = std::abs(data.mu.Dxy[i]);
    float absdz    = std::abs(data.mu.Dz[i]);
    bool id_loose_noiso  = (data.mu.IsLooseMuon[i]  == 1.0);
    bool id_select_noiso = (data.mu.IsMediumMuon[i] == 1.0);
    // Tacking efficiency scale factor
    utils::geteff_AE(eff_full_muon_trk,      eta, trk_sf,      trk_sf_err_up,      trk_sf_err_down);
    utils::geteff_AE(eff_full_muon_trk_veto, eta, trk_sf_veto, trk_sf_veto_err_up, trk_sf_veto_err_down);

    // Veto Muons
#if USE_POG_ID > 0
    // Using previous POG ID scale factors
    if ( id_veto_noiso &&
	 pt      >= MU_VETO_PT_CUT &&
	 abseta  <  MU_VETO_ETA_CUT &&
	 miniIso <  MU_VETO_MINIISO_CUT &&
	 absd0   <  MU_VETO_IP_D0_CUT &&
	 absdz   <  MU_VETO_IP_DZ_CUT ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_muon_looseid, pt, eta, sf, sf_err);
      weight_veto *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_looseid, pt, eta, sf, sf_err);
	weight_veto *= sf;
      }
      // Apply Isolation scale factor
      if (MU_VETO_MINIISO_CUT == 0.2)
	utils::geteff2D(eff_full_muon_miniiso02, pt, eta, sf, sf_err);
      else if (MU_VETO_MINIISO_CUT == 0.4)
	utils::geteff2D(eff_full_muon_miniiso04, pt, eta, sf, sf_err);
      weight_veto *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_miniiso02, pt, eta, sf, sf_err);
	weight_veto *= sf;
      }
      // Apply IP efficiency scale factor
      utils::geteff2D(eff_full_muon_looseip2d, pt, eta, sf, sf_err);
      weight_veto *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_looseip2d, pt, eta, sf, sf_err);
	weight_veto *= sf;
      }
      // Apply systematics
      weight_veto *= get_syst_weight(1, 0.03, nSigmaMuonFullSimSF);
      if (isFastSim) weight_veto *= get_syst_weight(1, 0.02, nSigmaMuonFastSimSF);
      // Apply Tracking scale factor here
      weight_veto *= get_syst_weight(trk_sf_veto, trk_sf_veto+trk_sf_veto_err_up, trk_sf_veto-trk_sf_veto_err_down, nSigmaMuonTrkSF);
    }
#else
    // Using Razor Inclusive ID
    // Apply Ele Reco + RazorInclusive ID scale factors
    if ( passMuVeto[i] ) {
      utils::geteff2D(eff_full_muon_veto, pt, eta, sf, sf_err);
      weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonFullSimSF);
      // Apply Tracking scale factor here
      if (MU_VETO_PT_CUT<10) {
	weight_veto *= get_syst_weight(trk_sf_veto, trk_sf_veto+trk_sf_veto_err_up, trk_sf_veto-trk_sf_veto_err_down, nSigmaMuonTrkSF);
      } else {
	weight_veto *= get_syst_weight(trk_sf, trk_sf+trk_sf_err_up, trk_sf-trk_sf_err_down, nSigmaMuonTrkSF);
      }
    }
#endif

    // Loose Muons
    if ( id_loose_noiso &&
	 pt      >= MU_LOOSE_PT_CUT &&
	 abseta  <  MU_LOOSE_ETA_CUT &&
	 miniIso <  MU_LOOSE_MINIISO_CUT &&
	 absd0   <  MU_LOOSE_IP_D0_CUT &&
	 absdz   <  MU_LOOSE_IP_DZ_CUT ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_muon_looseid, pt, eta, sf, sf_err);
      weight_loose *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_looseid, pt, eta, sf, sf_err);
	weight_loose *= sf;
      }
      // Apply Isolation scale factor
      if (MU_LOOSE_MINIISO_CUT == 0.2)
	utils::geteff2D(eff_full_muon_miniiso02, pt, eta, sf, sf_err);
      else if (MU_LOOSE_MINIISO_CUT == 0.4)
	utils::geteff2D(eff_full_muon_miniiso04, pt, eta, sf, sf_err);
      weight_loose *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_miniiso04, pt, eta, sf, sf_err);
	weight_loose *= sf;
      }
      // Apply IP efficiency scale factor
      utils::geteff2D(eff_full_muon_looseip2d, pt, eta, sf, sf_err);
      weight_loose *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_looseip2d, pt, eta, sf, sf_err);
	weight_loose *= sf;
      }
      // Apply systematics
      weight_loose *= get_syst_weight(1, 0.03, nSigmaMuonFullSimSF);
      if (isFastSim) weight_loose *= get_syst_weight(1, 0.02, nSigmaMuonFastSimSF);
      // Apply Tracking scale factor here
      if (MU_LOOSE_PT_CUT<10) {
	weight_loose *= get_syst_weight(trk_sf_veto, trk_sf_veto+trk_sf_veto_err_up, trk_sf_veto-trk_sf_veto_err_down, nSigmaMuonTrkSF);
      } else {
	weight_loose *= get_syst_weight(trk_sf, trk_sf+trk_sf_err_up, trk_sf-trk_sf_err_down, nSigmaMuonTrkSF);
      }
    }

    // Selected Muons
    if ( id_select_noiso &&
	 pt      >= MU_SELECT_PT_CUT &&
	 abseta  <  MU_SELECT_ETA_CUT &&
	 miniIso <  MU_SELECT_MINIISO_CUT &&
	 absd0   <  MU_SELECT_IP_D0_CUT &&
	 absdz   <  MU_SELECT_IP_DZ_CUT ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_muon_mediumid, pt, eta, sf, sf_err);
      weight_select *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_mediumid, pt, eta, sf, sf_err);
	weight_select *= sf;
      }
      // Apply Isolation scale factor
      if (MU_SELECT_MINIISO_CUT == 0.2)
	utils::geteff2D(eff_full_muon_miniiso02, pt, eta, sf, sf_err);
      else if (MU_SELECT_MINIISO_CUT == 0.4)
	utils::geteff2D(eff_full_muon_miniiso04, pt, eta, sf, sf_err);
      weight_select *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_miniiso02, pt, eta, sf, sf_err);
	weight_select *= sf;
      }
      // Apply IP efficiency scale factor
      utils::geteff2D(eff_full_muon_tightip2d, pt, eta, sf, sf_err);
      weight_select *= sf;
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_tightip2d, pt, eta, sf, sf_err);
	weight_select *= sf;
      }
      // Apply systematics
      weight_select *= get_syst_weight(1, 0.03, nSigmaMuonFullSimSF);
      if (isFastSim) weight_select *= get_syst_weight(1, 0.02, nSigmaMuonFastSimSF);
      // Apply Tracking scale factor here
      if (MU_SELECT_PT_CUT<10) {
	weight_select *= get_syst_weight(trk_sf_veto, trk_sf_veto+trk_sf_veto_err_up, trk_sf_veto-trk_sf_veto_err_down, nSigmaMuonTrkSF);
      } else {
	weight_select *= get_syst_weight(trk_sf, trk_sf+trk_sf_err_up, trk_sf-trk_sf_err_down, nSigmaMuonTrkSF);
      }
    }

  }

  return std::make_tuple(weight_veto, weight_loose, weight_select);
}


//  double AnalysisBase::calc_njet_weight(DataStruct& data, const double& nSigmaNJetWeight) {
//    // QCD Njet Reweighting
//    // First run this script to obtain the expo fit parameters
//    // root scripts/calc_njet_weight.C
//    // Up variation taken from Q region
//    // Down variation taken from G-1 region
//    double w_up = std::exp(-4.15259e-01 + 9.97141e-02 * nJet);
//    double w_dn = std::exp( 1.02951e+00 - 1.90498e-01 * nJet);
//    double weight = get_syst_weight(1.0, w_up, w_dn, nSigmaNJetWeight);
//    return weight;
//  }


double AnalysisBase::calc_lostlep_weight(DataStruct& data, const double& nSigmaLostLep) {
  // Lost Lepton event weight
  // First run this script to obtain the lost lepton uncertainties
  // python scripts/calc_lostlepsyst.py
  // Use the averages combination of all leptons for W/top
  double weight = 1;
  if (nEleVeto+nMuVeto+nTauVeto==0) {
    while(data.gen.Loop()) {
      size_t i = data.gen.it;
      // Check if the gen particle is a final state lepton from W mother
      if (passGenLepton[i]) {
	double unc = 0;
	double abseta = std::abs(data.gen.Eta[i]);
	if      (abseta<0.5) unc = 0.125;
	else if (abseta<1.0) unc = 0.126;
	else if (abseta<1.5) unc = 0.129;
	else if (abseta<2.0) unc = 0.143;
	else if (abseta<2.5) unc = 0.175;
	weight *= get_syst_weight(1.0, unc, nSigmaLostLep);
	//std::cout<<"Lost-lepton found: pt/eta/id = "<<data.gen.Pt[i]<<" "<<data.gen.Eta[i]<<" "<<data.gen.ID[i]<<" E/M/T="<<nEleVeto<<"/"<<nMuVeto<<"/"<<nTauVeto<<std::endl;
      }
    }
  }
  return weight;
}


double AnalysisBase::calc_trigger_efficiency(DataStruct& data, const double& nSigmaTrigger) {
  // 1D trigger efficiency method
  //  double eff, err_up, err_down;
  //  utils::geteff_AE(eff_trigger, AK4_Ht, eff, err_up, err_down);
  //  double w = get_syst_weight(eff, eff+err_up, eff-err_down, nSigmaTrigger);
  
  // Check the presence of a lepton/photon and apply different weights
  TH2D *h      = eff_trigger_veto;
  TH2D *h_up   = eff_trigger_veto_up;
  TH2D *h_down = eff_trigger_veto_down;
  TH2D *h_F    = eff_trigger_F_met;
  if (nEleVeto>=1) {
    h      = eff_trigger_ele;
    h_up   = eff_trigger_ele_up;
    h_down = eff_trigger_ele_down;
    h_F    = eff_trigger_F_ele;
  } else if (nMuVeto>=1) {
    h      = eff_trigger_mu;
    h_up   = eff_trigger_mu_up;
    h_down = eff_trigger_mu_down;
    h_F    = eff_trigger_F_mu;
  } else if (nPhotonPreSelect>=1) {
    h      = eff_trigger_pho;
    h_up   = eff_trigger_pho_up;
    h_down = eff_trigger_pho_down;
    h_F    = eff_trigger_F_pho;
  }

  // Trigger efficiency in the F region
  if (nJetAK8>0) {
    other_trigger_eff = utils::geteff2D(h_F, AK4_Ht, data.jetsAK8.Pt[iJetAK8[0]]);
  } else other_trigger_eff = 0.0;

  // 2D trigger efficiency (New)
  if (nJetAK8>0) {
    double eff = 0, total = 0;
    utils::geteff2D(h, AK4_Ht, data.jetsAK8.Pt[iJetAK8[0]], eff, total); // total was saved to histo error
    // For the time being only weight the measurable phase space
    // Rest is 0 --> Could weight with the TGraphAsymmErrors::Efficiency value (0.5+-0.5)
    if (total>0) {
      double eff_up   = utils::geteff2D(h_up,   AK4_Ht, data.jetsAK8.Pt[iJetAK8[0]]);
      double eff_down = utils::geteff2D(h_down, AK4_Ht, data.jetsAK8.Pt[iJetAK8[0]]);
      double w = get_syst_weight(eff, eff_up, eff_down, nSigmaTrigger);
      return w;
    } else return 0;
  } else return 0;
}


int AnalysisBase::calc_mrr2_bin(DataStruct& d, const char& region) {
  // Select the razor variable corresponding to the desired region
  double MR = d.evt.MR, R2 = d.evt.R2;
  if (region=='Z') {
    R2 = R2_ll;
  } else if (region=='L') {
    R2 = R2_1vl;
  } else if (region=='G') {
    R2 = R2_pho;
    MR = MR_pho;
  }
  int MRR2_bin = -1;
  if (MR>=800&&R2>=0.08) {
    MRR2_bin = 0;
    for (const auto& r2 : { 0.12, 0.16, 0.24, 0.4, 1.5 }) {
      if (R2<r2) break;
      ++MRR2_bin;
    }
    if (MRR2_bin<5) {
      for (const auto& mr : { 1000, 1200, 1600, 2000, 4000}) {
	if (MR<mr) break;
	MRR2_bin+=5;
      }
    } else MRR2_bin=9999;
  }
  return MRR2_bin;
}
