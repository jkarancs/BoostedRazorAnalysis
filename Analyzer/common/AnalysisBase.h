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
#include "TProfile.h"
#include "TStopwatch.h"
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
    sample(dirname)
  {
    sw_1_  = new TStopwatch;
    sw_1k_  = new TStopwatch;
    sw_10k_ = new TStopwatch;
    sw_job_ = new TStopwatch;
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

  void calculate_common_variables(DataStruct&, const unsigned int&);

  void init_common_histos();

  void fill_common_histos(DataStruct&, const unsigned int&, const double&);

  double get_xsec_from_ntuple(const std::vector<std::string>&, const std::string&);

  double get_xsec_from_txt_file(const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const std::string&);

  void calc_weightnorm_histo_from_ntuple(const std::vector<std::string>&, const double&, const std::vector<std::string>&,
					 const std::vector<std::string>&, bool);

  void init_pileup_reweightin(const std::string&, const std::string&, const std::vector<std::string>&);

  double get_pileup_weight(const int&, const double&);

  void rescale_smear_jet_met(DataStruct&, const bool&, const unsigned int&, const double&, const double&, const double&);

  double get_ht_weight(DataStruct&, const double&);

  double get_alphas_weight(const std::vector<float>&, const double&, const int&);

  double get_scale_weight(const std::vector<float>&, const double&, const unsigned int&);

  double get_syst_weight(const double&, const double&, const double&, const double&);

  double get_syst_weight(const double&, const double&, const double&);

  void job_monitoring(const int&, const int&, const std::string&, const float);

  void init_scale_factors();

  double calc_top_tagging_sf(DataStruct&, const double&);

  double calc_w_tagging_sf(DataStruct&, const double&);

  std::pair<double, double> calc_b_tagging_sf(DataStruct&, const double&, const bool&);

  std::pair<double, double> calc_ele_sf(DataStruct&, const double&, const double&, const double&, const double&, const bool&);

  std::pair<double, double> calc_muon_sf(DataStruct&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const bool&);
  
  double calc_trigger_efficiency(DataStruct&, const double&);

  std::map<char, std::vector<double> > scale_factors;
  std::map<char, double> sf_weight;

  const bool isData;
  const bool isSignal;
  const std::string sample;

private:
  TStopwatch *sw_1_, *sw_1k_, *sw_10k_, *sw_job_;
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

  std::map<char, std::vector<Cut> > analysis_cuts;

private:
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
  baseline_cuts.push_back({ .name="Clean_NGoodVtx",          .func = [&data] { return data.filter.goodVertices; } });
  
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
}

//_______________________________________________________
//                 Define common variables

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
  - pt > 200
  - |eta| < 2.4

*/
#define JET_AK4_PT_CUT  30
#define JET_AK4_ETA_CUT 2.4
#define JET_AK8_PT_CUT  200
#define JET_AK8_ETA_CUT 2.4

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
   Latest W-tagging Working points / scale factors (Spring16+PromptReco with ICHEP JEC)
   https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=43#Working_points_and_scale_factors
   
   - Using currently values for Spring16+ICHEP JEC
   TODO: Update to Moriond17+ReReco

   Choose:
   - Loose jet ID
   Tight Tag selection:
   - AK8 CHS jets
   - pt > 200
   - |eta| < 2.4
   - 65 <= Puppi SD Mass (corr) < 105
   - Puppi tau21 < 0.40

*/

#define W_PT_CUT            200
#define W_ETA_CUT           2.4
#define W_SD_MASS_CUT_LOW   65
#define W_SD_MASS_CUT_HIGH  105
#define W_TAU21_LOOSE_CUT   0.6 // DO NOT USE WITH PUPPI
#define W_TAU21_TIGHT_CUT   0.4

#define W_TAG_EFF_SF     1.03
#define W_TAG_EFF_SF_ERR 0.078
#define W_TAG_JMS_SF     1.00
#define W_TAG_JMS_SF_ERR 0.004
#define W_TAG_JMR_SF     1.08
#define W_TAG_JMR_SF_ERR 0.11

/*
  Top Tagging working points (No subjet B tag):
  https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging#13_TeV_working_points
  
  Latest WPs/SFs not yet on twiki (Angela Mc Lean, Mareike Meyer, Svenja Schumann):
  https://indico.cern.ch/event/523604/contributions/2147605/attachments/1263012/1868103/TopTaggingWp_v5.pdf
  https://indico.cern.ch/event/523604/contributions/2147605/attachments/1263012/1868207/TopTaggingSF_76X.pdf
  
  !! Warning, Scale factor to be updated for Puppi jets !!
  
  Choose:
  - Loose selection: e(S) = 51.2%, e(B) = 3% WP
  - AK8 Puppi jets
  - 105 < SD Mass < 200
  - tau32 < 0.67


  Top Tagging working point (With subjet B tag, 2015 Data)
  
  Latest WPs/SFs -  JME-16-003 PAS:
  http://cms.cern.ch/iCMS/analysisadmin/viewanalysis?id=1694&field=id&value=1694
  
  Choose:
  - Loose selection: e(S) ~= 55%, e(B) = 3.0% WP
  - AK8 Puppi jets
  - 105 <= SD Mass < 210
  - tau32 < 0.8
  - max subjet BTag CSV > 0.46
*/

#define USE_BTAG 1

#if USE_BTAG == 0

#define TOP_PT_CUT            400
#define TOP_SD_MASS_CUT_LOW   105 // prev 74X 110
#define TOP_SD_MASS_CUT_HIGH  200 // prev 74X 210
#define TOP_TAU32_CUT        0.67 // prev 74X 0.75

#define TOP_TAG_SF_LOW       0.97
#define TOP_TAG_SF_LOW_ERR   0.09
#define TOP_TAG_SF_HIGH      0.99
#define TOP_TAG_SF_HIGH_ERR  0.18

#else
#define TOP_PT_CUT            400
#define TOP_SD_MASS_CUT_LOW   105
#define TOP_SD_MASS_CUT_HIGH  210
#define TOP_TAU32_CUT        0.80
#define TOP_BTAG_CSV         0.46

#define TOP_TAG_SF_LOW       1.04
#define TOP_TAG_SF_LOW_ERR   0.11
#define TOP_TAG_SF_HIGH      1.05
#define TOP_TAG_SF_HIGH_ERR  0.26
#endif

/*
  Latest Electron IDs:
  [1] Cut Based  - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Working_points_for_2016_data_for
  [2] MVA        - https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2?rev=30#Recommended_MVA_recipes_for_2016
  [3] SUSY (Use) - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#Electrons

  Latest Isolation WPs:
  [4] SUSY MiniIso Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO_AN1

  Latest Impact Point Cut:
  [5] SUSY Loose/Tight IP2D (Use) - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO_AN1
  [6] POG  Tight - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Offline_selection_criteria

  For Veto Choose:
  - Spring15 Cut based Veto ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (Medium WP [4])
  - pt >= 5
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.2, |dz| < 0.5 (Loose IP cut Recommended by SUSY Group [5])

  For Selection Choose:
  - Spring15 Cut based Medium ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Loose IP cut Recommended by SUSY Group [5])

*/

#define ELE_VETO_PT_CUT        5
#define ELE_VETO_ETA_CUT       2.5
#define ELE_VETO_MINIISO_CUT   0.1
#define ELE_VETO_IP_D0_CUT     0.2
#define ELE_VETO_IP_DZ_CUT     0.5

#define ELE_SELECT_PT_CUT       10
#define ELE_SELECT_ETA_CUT      2.5
#define ELE_SELECT_MINIISO_CUT  0.1
#define ELE_SELECT_IP_D0_CUT   0.05
#define ELE_SELECT_IP_DZ_CUT    0.1

/*
  Latest Muon IDs (Loose/Medium):
  [1] POG Loose/Medium - https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=28#Short_Term_Instructions_for_Mori

  Latest Isolation WPs:
  [2] SUSY MiniISo Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO

  Latest Impact Point Cut (Loose/Tight):
  [3] SUSY Loose/Tight IP2D - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO
  
  For Veto Choose:
  - POG recommended Loose ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.4 (Loose WP [2])
  - pt >= 5
  - |eta| < 2.4
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [3])

  For Selection Choose:
  - POG recommended Medium ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.2 (tight WP [2])
  - pt >= 5
  - |eta| < 2.4
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [3])

*/

#define MU_VETO_PT_CUT         5
#define MU_VETO_ETA_CUT        2.4
#define MU_VETO_MINIISO_CUT    0.4
#define MU_VETO_IP_D0_CUT      0.2
#define MU_VETO_IP_DZ_CUT      0.5

#define MU_SELECT_PT_CUT        10
#define MU_SELECT_ETA_CUT       2.4
#define MU_SELECT_MINIISO_CUT   0.2
#define MU_SELECT_IP_D0_CUT     0.05
#define MU_SELECT_IP_DZ_CUT     0.1

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
std::vector<bool> passLooseBTag;
std::vector<bool> passMediumBTag;
std::vector<bool> passTightBTag;
unsigned int nJet;
unsigned int nLooseBTag;
unsigned int nMediumBTag;
unsigned int nTightBTag;
double AK4_Ht;
double minDeltaPhi; // Min(DeltaPhi(Jet_i, MET)), i=1,2,3

// AK8 jets
std::vector<size_t > iJetAK8;
std::vector<size_t > iWPreTag;
std::vector<size_t > iLooseWTag;
std::vector<size_t > iTightWTag;
std::vector<size_t > iTightWAntiTag;
std::vector<size_t > itJetAK8;
std::vector<size_t > itWPreTag;
std::vector<size_t > itLooseWTag;
std::vector<size_t > itTightWTag;
std::vector<size_t > itTightWAntiTag;
std::vector<double> tau21;
std::vector<double> tau31;
std::vector<double> tau32;
std::vector<float> softDropMassCorr; // Correction + uncertainties for W tagging
#if VER == 0
std::vector<double> maxSubjetCSV;
#endif
std::vector<bool> passSubjetBTag;
std::vector<bool> passLooseJetAK8;
std::vector<bool> passWPreTag;
std::vector<bool> passLooseWTag;
std::vector<bool> passTightWTag;
std::vector<bool> passTightWAntiTag;
std::vector<bool> passHadTopPreTag;
std::vector<bool> passHadTopTag;
unsigned int nJetAK8;
unsigned int nWPreTag;
unsigned int nLooseWTag;
unsigned int nTightWTag;
unsigned int nTightWAntiTag;
unsigned int nSubjetBTag;
unsigned int nHadTopTag;
unsigned int nHadTopPreTag;
double AK8_Ht;

// Event Letpons
std::vector<size_t > iEleSelect;
std::vector<size_t > iMuSelect;
std::vector<size_t > itEleSelect;
std::vector<size_t > itMuSelect;
std::vector<bool> passEleVeto;
std::vector<bool> passMuVeto;
std::vector<bool> passEleSelect;
std::vector<bool> passMuSelect;
unsigned int nEleVetoNoIso;
unsigned int nEleVeto;
unsigned int nEleSelect;
unsigned int nMuVetoNoIso;
unsigned int nMuVeto;
unsigned int nMuSelect;
unsigned int nLepVetoNoIso;
unsigned int nLepVeto;
unsigned int nLepSelect;
double MT;

void
AnalysisBase::calculate_common_variables(DataStruct& data, const unsigned int& syst_index)
{
  std::vector<TLorentzVector> veto_leptons_noiso, veto_leptons, selected_leptons;
  //std::vector<TLorentzVector> veto_muons_noiso, veto_muons, selected_muons;
  //std::vector<bool> veto_lep_in_jet;
  //std::vector<bool> veto_mu_in_jet, selected_mu_in_jet;

  // It only makes sense to calculate certain variables only once if they don't depend on jet energy
  if (syst_index == 0) {

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

    // Event Letpons
    iEleSelect   .clear();
    itEleSelect  .assign(data.ele.size, (size_t)-1);
    passEleVeto  .assign(data.ele.size, 0);
    passEleSelect.assign(data.ele.size, 0);
    nEleVetoNoIso = nEleVeto = nEleSelect = 0;
    while(data.ele.Loop()) {
      size_t i = data.ele.it;
      TLorentzVector ele_v4; ele_v4.SetPtEtaPhiE(data.ele.Pt[i], data.ele.Eta[i], data.ele.Phi[i], data.ele.E[i]);
      float pt = data.ele.Pt[i];
      float abseta = std::abs(data.ele.Eta[i]);
      float miniIso = data.ele.MiniIso[i]/data.ele.Pt[i];
      float absd0 = std::abs(data.ele.Dxy[i]);
      float absdz = std::abs(data.ele.Dz[i]);
      bool id_veto = (data.ele.vidVetonoiso[i] == 1.0);
      bool id_select = (data.ele.vidMediumnoiso[i] == 1.0);
      //bool id_veto = (data.ele.vidVeto[i] == 1.0);
      //bool id_select = (data.ele.vidTight[i] == 1.0);
      // Veto
      if (passEleVeto[i] = 
	  ( id_veto &&
	    pt      >= ELE_VETO_PT_CUT &&
	    abseta  <  ELE_VETO_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
	    absd0   <  ELE_VETO_IP_D0_CUT &&
	    absdz   <  ELE_VETO_IP_DZ_CUT) ) {
	veto_leptons_noiso.push_back(ele_v4);
	nEleVetoNoIso++;
	if (miniIso <  ELE_VETO_MINIISO_CUT) {
	  nEleVeto++;
	  veto_leptons.push_back(ele_v4);
	  //veto_lep_in_jet.push_back(data.ele.IsPartOfNearAK4Jet[i]);
	}
      }
      // Select
      if (passEleSelect[i] = 
	  ( id_select &&
	    pt        >= ELE_SELECT_PT_CUT &&
	    abseta    <  ELE_SELECT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
	    miniIso   <  ELE_SELECT_MINIISO_CUT &&
	    absd0     <  ELE_SELECT_IP_D0_CUT &&
	    absdz     <  ELE_SELECT_IP_DZ_CUT) ) {
	selected_leptons.push_back(ele_v4);
	iEleSelect.push_back(i);
	itEleSelect[i] = nEleSelect++;
      }
    }

    // Number of Veto/Select Muons
    iMuSelect    .clear();
    itMuSelect   .assign(data.mu.size,  (size_t)-1);
    passMuVeto   .assign(data.mu.size,  0);
    passMuSelect .assign(data.mu.size,  0);
    nMuVetoNoIso = nMuVeto = nMuSelect = 0;
    while(data.mu.Loop()) {
      size_t i = data.mu.it;
      TLorentzVector mu_v4; mu_v4.SetPtEtaPhiE(data.mu.Pt[i], data.mu.Eta[i], data.mu.Phi[i], data.mu.E[i]);
      float pt = data.mu.Pt[i];
      float abseta = std::abs(data.mu.Eta[i]);
      float miniIso = data.mu.MiniIso[i]/data.mu.Pt[i];
      float absd0 = std::abs(data.mu.Dxy[i]);
      float absdz = std::abs(data.mu.Dz[i]);
      bool id_veto = (data.mu.IsLooseMuon[i] == 1.0);
      bool id_select = (data.mu.IsMediumMuon[i] == 1.0);
      // Veto
      if (passMuVeto[i] =
	  (id_veto &&
	   pt      >= MU_VETO_PT_CUT &&
	   abseta  <  MU_VETO_ETA_CUT &&
	   absd0   <  MU_VETO_IP_D0_CUT &&
	   absdz   <  MU_VETO_IP_DZ_CUT) ) {
	veto_leptons_noiso.push_back(mu_v4);
	//veto_muons_noiso.push_back(mu_v4);
	nMuVetoNoIso++;
	if (miniIso <  MU_VETO_MINIISO_CUT) {
	  nMuVeto++;
	  veto_leptons.push_back(mu_v4);
	  //veto_muons.push_back(mu_v4);
	  //veto_lep_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
	  //veto_mu_in_jet.push_back(data.mu.IsPartOfNearAK4Jet[i]);
	}
      }
      // Select
      if (passMuSelect[i] =
	  ( id_select &&
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
    }

    nLepVetoNoIso = nEleVetoNoIso + nMuVetoNoIso;
    nLepVeto      = nEleVeto  + nMuVeto;
    nLepSelect    = nEleSelect + nMuSelect;

    // MT
    MT = 9999;
    if (nLepSelect==1) {
      if (nEleSelect==1) {
	MT = sqrt( 2*data.ele.Pt[iEleSelect[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.ele.Phi[iEleSelect[0]])) );
      } else if (nMuSelect==1) {
	MT = sqrt( 2*data.mu.Pt[iMuSelect[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.mu.Phi[iMuSelect[0]])) );
      }
    }
  }

  // Rest of the vairables need to be recalculated each time the jet energy is changed
  // eg. Jet selection, W/top tags, HT (obviously), etc. that depends on jet pt

  // AK4 jets
  iJet       .clear();
  iLooseBTag .clear();
  iMediumBTag.clear();
  iTightBTag .clear();
  itJet         .assign(data.jetsAK4.size, (size_t)-1);
  itLooseBTag   .assign(data.jetsAK4.size, (size_t)-1);
  itMediumBTag  .assign(data.jetsAK4.size, (size_t)-1);
  itTightBTag   .assign(data.jetsAK4.size, (size_t)-1);
  passLooseJet  .assign(data.jetsAK4.size, 0);
  passLooseBTag .assign(data.jetsAK4.size, 0);
  passMediumBTag.assign(data.jetsAK4.size, 0);
  passTightBTag .assign(data.jetsAK4.size, 0);
  nJet = 0;
  nLooseBTag  = 0;
  nMediumBTag = 0;
  nTightBTag  = 0;
  AK4_Ht = 0;
  minDeltaPhi = 9999;
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
      iJet.push_back(i);
      itJet[i] = nJet++;

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
      if (nJet<=3) {
	double dphi = std::abs(TVector2::Phi_mpi_pi(data.met.Phi[0] - data.jetsAK4.Phi[i]));
	if (dphi<minDeltaPhi) minDeltaPhi = dphi;
      }
      
    } // End Jet Selection
    
    // Online jet selection for HT (+ testing Additional Loose Jet ID)
    if ( //data.jetsAK4.looseJetID[i] == 1 &&
	 data.jetsAK4.Pt[i]         >  30 &&
	 std::abs(data.jetsAK4.Eta[i])  <  3.0 ) {
      AK4_Ht += data.jetsAK4.Pt[data.jetsAK4.it];

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
  
  // Add isolated leptons to HT computation
  //for (size_t ilep=0, nlep=veto_leptons.size(); ilep<nlep; ++ilep)
  //  if (add_lepton_to_ht[ilep]) AK4_Ht += veto_leptons[ilep].Pt();

  //for (size_t imu=0, nmu=selected_muons.size(); imu<nmu; ++imu)
  //  if (remove_muon_from_ht[imu]) AK4_Ht -= selected_muons[imu].Pt();
  
  // AK8 jets
  iJetAK8   .clear();
  iWPreTag  .clear();
  iLooseWTag.clear();
  iTightWTag.clear();
  iTightWAntiTag.clear();
  itJetAK8         .assign(data.jetsAK8.size, (size_t)-1);
  itWPreTag        .assign(data.jetsAK8.size, (size_t)-1);
  itLooseWTag      .assign(data.jetsAK8.size, (size_t)-1);
  itTightWTag      .assign(data.jetsAK8.size, (size_t)-1);
  itTightWAntiTag  .assign(data.jetsAK8.size, (size_t)-1);
  passLooseJetAK8  .assign(data.jetsAK8.size, 0);
  passWPreTag      .assign(data.jetsAK8.size, 0);
  passLooseWTag    .assign(data.jetsAK8.size, 0);
  passTightWTag    .assign(data.jetsAK8.size, 0);
  passTightWAntiTag.assign(data.jetsAK8.size, 0);
  passHadTopPreTag .assign(data.jetsAK8.size, 0);
  passHadTopTag    .assign(data.jetsAK8.size, 0);
  nJetAK8       = 0;
  nWPreTag      = 0;
  nLooseWTag    = 0;
  nTightWTag    = 0;
  nTightWAntiTag= 0;
  nSubjetBTag   = 0;
  nHadTopTag    = 0;
  nHadTopPreTag = 0;
  AK8_Ht   = 0;
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    // Jet ID
    if ( passLooseJetAK8[i] = 
	 ( data.jetsAK8.looseJetID[i] == 1 &&
	   data.jetsAK8.Pt[i]         >= JET_AK8_PT_CUT &&
	   std::abs(data.jetsAK8.Eta[i])  <  JET_AK8_ETA_CUT ) ) {
      iJetAK8.push_back(i);
      itJetAK8[i] = nJetAK8++;

      // Tagging Variables
      double pt      = data.jetsAK8.Pt[i];
      double abseta  = data.jetsAK8.Eta[i];
      double sd_mass = softDropMassCorr[i]; // Corrected, scaled, smeared
      double tau_21 = tau21[i];
      double tau_32 = tau32[i];

      // _______________________________________________________
      //                   Hadronic W Tag definition

      if (passWPreTag[i] = 
	  ( pt      >= W_PT_CUT &&
	    abseta  <  W_ETA_CUT &&
	    sd_mass >= W_SD_MASS_CUT_LOW && 
	    sd_mass <  W_SD_MASS_CUT_HIGH) ) {
	iWPreTag.push_back(i);
	itWPreTag[i] = nWPreTag++;
	// Loose/Tight W Tag Working points
	if (passLooseWTag[i] = (tau_21 < W_TAU21_LOOSE_CUT) ) {
	  iLooseWTag.push_back(i);
	  itLooseWTag[i] = nLooseWTag++;
	}
	if (passTightWTag[i] = (tau_21 < W_TAU21_TIGHT_CUT) ) {
	  iTightWTag.push_back(i);
	  itTightWTag[i] = nTightWTag++;
	} else {
	  iTightWAntiTag.push_back(i);
	  itTightWAntiTag[i] = nTightWAntiTag++;
	}
      }

      // _______________________________________________________
      //                  Hadronic Top Tag definition

      // New hadronic top tag
      if (passHadTopPreTag[i] = 
	  ( pt      >= TOP_PT_CUT && 
	    sd_mass >= TOP_SD_MASS_CUT_LOW &&
	    sd_mass <  TOP_SD_MASS_CUT_HIGH) ) {
#if USE_BTAG == 1
	if (passHadTopPreTag[i] = passSubjetBTag[i]) {
#endif
	  nHadTopPreTag++;
	  if (passHadTopTag[i] = (tau_32 < TOP_TAU32_CUT) ) nHadTopTag++;
#if USE_BTAG == 1
	}
#endif
      }

    } // End Jet Selection

    // Online jet selection for AK8 HT
    if ( data.jetsAK8.Pt[i]         > 150 &&
	 std::abs(data.jetsAK8.Eta[i])  <  2.5 ) {
      // Ht
      AK8_Ht += data.jetsAK8.Pt[i];
    }

  } // End AK8 Jet Loop
}


//_______________________________________________________
//                 List of Histograms

TH1D* h_totweight;
std::vector<TH2D*> vh_totweight_signal;
std::vector<TH2D*> vh_xsec_signal;
std::vector<TH2D*> vh_weightnorm_signal;
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
TH2D* h_btageff_b_loose;
TH2D* h_btageff_c_loose;
TH2D* h_btageff_l_loose;
TH2D* h_btageff_b_medium;
TH2D* h_btageff_c_medium;
TH2D* h_btageff_l_medium;

//_______________________________________________________
//              Define Histograms here
void
AnalysisBase::init_common_histos()
{
  // total weight
  h_totweight                  = new TH1D("totweight",          "MC;;Total (generator) event weight", 1,0,1);
  // signal weight
  vh_totweight_signal .push_back(new TH2D("totweight_T1tttt",   "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",        201,-12.5,5012.5, 201,-12.5,5012.5));
  vh_xsec_signal      .push_back(new TH2D("xsec_T1tttt",        "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Cross-section (pb)",  201,-12.5,5012.5, 201,-12.5,5012.5));
  vh_weightnorm_signal.push_back(new TH2D("weightnorm_T1tttt",  "T1tttt or T5ttcc or T5tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);weight norm. factor", 201,-12.5,5012.5, 201,-12.5,5012.5));
  vh_totweight_signal .push_back(new TH2D("totweight_T2tt",     "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",        401,-2.5,2002.5, 401,-2.5,2002.5));
  vh_xsec_signal      .push_back(new TH2D("xsec_T2tt",          "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Cross-section (pb)",  401,-2.5,2002.5, 401,-2.5,2002.5));
  vh_weightnorm_signal.push_back(new TH2D("weightnorm_T2tt",    "T2tt;M_{#tilde{s}} (GeV);M_{#tilde{#chi}^{0}} (GeV);weight norm. factor", 401,-2.5,2002.5, 401,-2.5,2002.5));
  // pileup
  h_pileup_data                = new TH1D("pileup_data",        "Pile-up distribution - Data (Nominal);Pile-up", 100,0,100);
  h_pileup_data_down           = new TH1D("pileup_data_down",   "Pile-up distribution - Data (down);Pile-up",    100,0,100);
  h_pileup_data_up             = new TH1D("pileup_data_up",     "Pile-up distribution - Data (up);Pile-up",      100,0,100);
  h_pileup_mc                  = new TH1D("pileup_mc",          "Pile-up distribution - MC;Pile-up",             100,0,100);
  h_pileup_weight              = new TH1D("pileup_weight",      "Pile-up weights - Nominal MB X-sec (69 mb);Pile-up;Weight",    100,0,100);
  h_pileup_weight_down         = new TH1D("pileup_weight_down", "Pile-up weights - MB X-sec up 5% (72.45 mb);Pile-up;Weight",   100,0,100);
  h_pileup_weight_up           = new TH1D("pileup_weight_up",   "Pile-up weights - MB X-sec down 5% (65.55 mb);Pile-up;Weight", 100,0,100);
  h_nvtx                       = new TH1D("nvtx",               "Number of vertices - Nominal;N_{Vertices}",                      100,0,100);
  h_nvtx_rw                    = new TH1D("nvtx_rw",            "Number of vertices - Pile-up reweighted (MC only);N_{Vertices}", 100,0,100);
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
  double btagbins[3] = { -0.5,0.5,1.5 };
  h_btageff_b_loose            = new TH2D("btageff_b_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,btagbins);
  h_btageff_c_loose            = new TH2D("btageff_c_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,btagbins);
  h_btageff_l_loose            = new TH2D("btageff_l_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,btagbins);
  h_btageff_b_medium           = new TH2D("btageff_b_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,btagbins);
  h_btageff_c_medium           = new TH2D("btageff_c_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,btagbins);
  h_btageff_l_medium           = new TH2D("btageff_l_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,btagbins);
}

//_______________________________________________________
//               Fill Histograms here
void
AnalysisBase::fill_common_histos(DataStruct& d, const unsigned int& syst_index, const double& weight)
{
  if (syst_index == 0) {
    while(d.jetsAK4.Loop()) {
      size_t i = d.jetsAK4.it;
      if (passLooseJet[i]) {
	if (d.jetsAK4.HadronFlavour[i]==5) {
	  h_btageff_b_loose ->Fill(d.jetsAK4.Pt[i], passLooseBTag[i]);
	  h_btageff_b_medium->Fill(d.jetsAK4.Pt[i], passMediumBTag[i]);
	} else if (d.jetsAK4.HadronFlavour[i]==4) {
	  h_btageff_c_loose ->Fill(d.jetsAK4.Pt[i], passLooseBTag[i]);
	  h_btageff_c_medium->Fill(d.jetsAK4.Pt[i], passMediumBTag[i]);
	} else {
	  h_btageff_l_loose ->Fill(d.jetsAK4.Pt[i], passLooseBTag[i]);
	  h_btageff_l_medium->Fill(d.jetsAK4.Pt[i], passMediumBTag[i]);
	}
      }
    }
  }
}

//_______________________________________________________
//           Read cross-section from ntuple
double
AnalysisBase::get_xsec_from_ntuple(const std::vector<std::string>& filenames, const std::string& treename)
{
  float evt_XSec=0, prev_XSec=0;
  for (const auto& filename : filenames) {
    TFile *f = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)f->Get(treename.c_str());
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
double
AnalysisBase::get_xsec_from_txt_file(const std::string& txt_file)
{
  double evt_XSec = 0;
  std::ifstream xsecFile(txt_file.c_str());
  if ( !xsecFile.good() ) {
    return -9999.0;
    std::cout<<"Unable to open cross-section file: "<<txt_file<<std::endl;
    utils::error("Please provide the correct txt file for Cross-sections in settings.h!");
  } else {

    // Read all nSigmas, nums
    std::string line;
    std::string shortname;
    std::string primary_dataset;
    double xsec;
    while ( std::getline(xsecFile, line) ) {
      std::stringstream nth_line;
      nth_line<<line;
      nth_line>>shortname;
      nth_line>>primary_dataset;
      nth_line>>xsec;
      if (sample==shortname) evt_XSec = xsec;
    }
  }
  if (evt_XSec == 0) {
    std::cout<<"No crossection found for "<<sample<<" in cross section file: "<<txt_file<<std::endl;
    utils::error("Please fix the cross-section file in settings.h!");
  }

  return evt_XSec;
}

//_______________________________________________________
//          Read total weight from ntuple histos
double
AnalysisBase::get_totweight_from_ntuple(const std::vector<std::string>& filenames, const std::string& histoname)
{
  // Merging totweight histos
  for (const auto& filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    h_totweight->Add((TH1D*)f->Get(histoname.c_str()));
    f->Close();
  }
  return h_totweight->GetBinContent(1);
}

//_______________________________________________________
//       Calculate weight normalization for signal
void
AnalysisBase::calc_weightnorm_histo_from_ntuple(const std::vector<std::string>& filenames, const double& intLumi, const std::vector<std::string>& vname_signal,
						const std::vector<std::string>& vname_totweight, bool verbose=1)
{
  // Find the index of the current signal
  int signal_index = -1;
  std::string signal_name = "";
  if (filenames.size()>0) for (size_t i=0, n=vname_signal.size(); i<n; ++i) 
    if (filenames[0].find(vname_signal[i])!=std::string::npos&&signal_index==-1) {
      signal_index = i;
      signal_name = vname_signal[i];
    }
  signal_index = (signal_index>=3); // 0: Mlsp vs Mgluino - T1tttt, T5ttcc, T5tttt; 1: Mlsp vs Mstop - T2tt

  // Merge totweight histos
  std::map<int, double> xsec_mother;
  for (const auto& filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    // Get total weight
    TH2D* totweight = (TH2D*)f->Get(vname_totweight[signal_index].c_str());
    vh_totweight_signal[signal_index]->Add(totweight);
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
  if (verbose) {
    std::cout<<"- Signal: "<<signal_name<<std::endl;
    for (int binx=1, nbinx=vh_xsec_signal[signal_index]->GetNbinsX(); binx<=nbinx; ++binx) 
      for (int biny=1, nbiny=vh_xsec_signal[signal_index]->GetNbinsY(); biny<=nbiny; ++biny) {
        double mMother = vh_xsec_signal[signal_index]->GetXaxis()->GetBinCenter(binx);
        double mLSP = vh_xsec_signal[signal_index]->GetYaxis()->GetBinCenter(biny);
        double xsec  = vh_xsec_signal[signal_index]      ->GetBinContent(binx, biny);
        double totw  = vh_totweight_signal[signal_index] ->GetBinContent(binx, biny);
        double wnorm = vh_weightnorm_signal[signal_index]->GetBinContent(binx, biny);
        if (totw>0) std::cout<<(signal_index?"  Bin: M(s~)=":"  Bin: M(g~)=")<<mMother<<" M(LSP)="<<mLSP<<":   xsec="<<xsec<<" totweight="<<totw<<" weightnorm="<<wnorm<<std::endl;
      }
    std::cout<<std::endl;
  }
}


//_______________________________________________________
//             Load pile-up reweighting infos
void
AnalysisBase::init_pileup_reweightin(const std::string& pileupDir, const std::string& mcPileupHistoName, const std::vector<std::string>& filenames)
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
  //   h_pileup_mc->Add((TH1D*)f_pileup_mc->Get(mcPileupHistoName.c_str()));
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
//                  Get pile-up weight
double
AnalysisBase::get_pileup_weight(const int& NtrueInt, const double& nSigmaPU)
{
  int pu_bin = NtrueInt+1; // eg. pileup 0, is filled in bin 1
  double w_pileup = h_pileup_weight->GetBinContent(pu_bin);
  double w_pileup_up = h_pileup_weight_up->GetBinContent(pu_bin);
  double w_pileup_down = h_pileup_weight_down->GetBinContent(pu_bin);
  double w = get_syst_weight(w_pileup, w_pileup_up, w_pileup_down, nSigmaPU);
  return w;
}


//_______________________________________________________
//              Rescale jet 4-momenta

std::vector<float> AK4_E, AK4_Pt;
std::vector<float> AK8_E, AK8_Pt, AK8_softDropMass;//, AK8_trimmedMass, AK8_prunedMass, AK8_filteredMass;
std::vector<float> AK8_softDropMassCorr; // Correction for W tagging

std::vector<float> AK4_JERSmearFactor,     AK8_JERSmearFactor,     AK8_JMRSmearFactor;  
std::vector<float> AK4_JERSmearFactorUp,   AK8_JERSmearFactorUp,   AK8_JMRSmearFactorUp;
std::vector<float> AK4_JERSmearFactorDown, AK8_JERSmearFactorDown, AK8_JMRSmearFactorDown;

TVector3 met;
TVector3 dmet_JESUp,  dmet_JESDown;
TVector3 dmet_JERUp,  dmet_JERDown;
TVector3 dmet_RestUp, dmet_RestDown;

void
AnalysisBase::rescale_smear_jet_met(DataStruct& data, const bool& applySmearing, const unsigned int& syst_index,
				    const double& nSigmaJES, const double& nSigmaJER, const double& nSigmaRestMET)
{
  // Apply Jet Energy Scale (JES) and Jet Energy Resolution (JER) corrections
  // For AK8 jets which are used for W tagging (only):
  // - Additionally apply jet mass scale (JMS) and jet mass resolutin (JMR) corrections
  //   (We use the combination measured by the JMAR group)

  // Initialization (needed for later variations
  if (syst_index==0) {
    if (applySmearing) {
      // Calculate the smear factors
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
      AK8_JMRSmearFactor    .clear();
      AK8_JMRSmearFactorUp  .clear();
      AK8_JMRSmearFactorDown.clear();
      for (size_t i=0; i<data.jetsAK8.size; ++i) {
        double JERSmear     = data.jetsAK8.SmearedPt[i]/data.jetsAK8.Pt[i];
        double JERSmearUp   = 1 + (JERSmear-1) * (data.jetsAK8.JERSFUp[i]  -1) / (data.jetsAK8.JERSF[i]-1);
        double JERSmearDown = 1 + (JERSmear-1) * (data.jetsAK8.JERSFDown[i]-1) / (data.jetsAK8.JERSF[i]-1);
        AK8_JERSmearFactor    .push_back(JERSmear);
        AK8_JERSmearFactorUp  .push_back(JERSmearUp);
        AK8_JERSmearFactorDown.push_back(JERSmearDown);
        // For the JMR Smearing apply the same procedure, but swap to the JMR scale factor
        double JMRSmear     = 1 + (JERSmear-1) * (W_TAG_JMR_SF                 -1) / (data.jetsAK8.JERSF[i]-1);
        double JMRSmearUp   = 1 + (JERSmear-1) * (W_TAG_JMR_SF+W_TAG_JMR_SF_ERR-1) / (data.jetsAK8.JERSF[i]-1);
        double JMRSmearDown = 1 + (JERSmear-1) * (W_TAG_JMR_SF-W_TAG_JMR_SF_ERR-1) / (data.jetsAK8.JERSF[i]-1);
        AK8_JMRSmearFactor    .push_back(JMRSmear);
        AK8_JMRSmearFactorUp  .push_back(JMRSmearUp);
        AK8_JMRSmearFactorDown.push_back(JMRSmearDown);      
      }
    }
    // Save the original values for later (before applying any systematics)
    AK4_E            = data.jetsAK4.E;
    AK4_Pt           = data.jetsAK4.Pt;
    AK8_E            = data.jetsAK8.E;
    AK8_Pt           = data.jetsAK8.Pt;
#if VER == 0
    AK8_softDropMass = data.jetsAK8.softDropMass;
#else
    AK8_softDropMass = data.jetsAK8.softDropMassPuppi;
#endif
    //AK8_trimmedMass  = data.jetsAK8.trimmedMass;
    //AK8_prunedMass   = data.jetsAK8.prunedMass;
    //AK8_filteredMass = data.jetsAK8.filteredMass;

    // Correction for Puppi SoftDrop Mass
    // (Needed for W tagging)
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=43#Working_points_and_scale_factors
    // TODO: Use uncorrected SoftDrop Mass (also for the next ntuple production)
    // - Fix Crash with Eval() functions - TF1s cannot correctly be loaded currently
    if (puppisd_corrGEN_==0) {
      AK8_softDropMassCorr = AK8_softDropMass;
    } else {
      AK8_softDropMassCorr.clear();
      for (size_t i=0; i<data.jetsAK8.size; ++i) {
#if VER == 0
	double puppi_pt  = data.jetsAK8.Pt[i];
	double puppi_eta = data.jetsAK8.Eta[i];
	double puppi_sd_mass = data.jetsAK8.softDropMass[i]/data.jetsAK8.jecFactor0[i];
#else
	double puppi_pt  = data.jetsAK8.PtPuppi[i];
	double puppi_eta = data.jetsAK8.EtaPuppi[i];
	double puppi_sd_mass = data.jetsAK8.softDropMassPuppi[i];
#endif
	double corr = puppisd_corrGEN_->Eval(puppi_pt);
	if(std::abs(puppi_eta)<=1.3) corr *= puppisd_corrRECO_cen_->Eval(puppi_pt);
	else corr *= puppisd_corrRECO_for_->Eval(puppi_pt);
	
	AK8_softDropMassCorr.push_back(puppi_sd_mass * corr);
      }
    }

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
    // Use the phi direction of the largest remaining systematc
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
	  phi_up    = dmet.Phi();
	  ptsum_up  = std::sqrt(ptsum_up*ptsum_up + dmet.Perp2());
	}
      } else {
	// Rest Down
	if (dmet.Pt()>maxdpt_down) {
	  maxdpt_down = dmet.Pt();
	  phi_down    = dmet.Phi();
	  ptsum_down  = std::sqrt(ptsum_down*ptsum_down + dmet.Perp2());
	}
      }
    }
    dmet_RestUp.  SetPtEtaPhi(ptsum_up,   0, phi_up);
    dmet_RestDown.SetPtEtaPhi(ptsum_down, 0, phi_down);
#endif
  }


  // Apply systematic variations
  // Even if Sigmas=0, we still smear jets!
  // AK4 jets
  //AK4_Ht = 0;
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
    //AK4_Ht += data.jetsAK4.Pt[i];
  }
  // AK8 jets
  //AK8_Ht = 0;
  softDropMassCorr.clear();
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
    //AK8_Ht += data.jetsAK8.Pt[i];

    // For Top jet mass, similarly apply only JES + JER for now
    // (Since there's no other recommendation)
#if VER == 0
    data.jetsAK8.softDropMass[i]    = AK8_softDropMass[i] * scaleJES;
    if (applySmearing) data.jetsAK8.softDropMass[i] *= scaleJER;
#else
    data.jetsAK8.softDropMassPuppi[i]    = AK8_softDropMass[i] * scaleJES;
    if (applySmearing) data.jetsAK8.softDropMassPuppi[i] *= scaleJER;
#endif
    //data.jetsAK8.trimmedMass[i]  = AK8_trimmedMass[i]  * scaleJES * scaleJER;
    //data.jetsAK8.prunedMass[i]   = AK8_prunedMass[i]   * scaleJES * scaleJER;
    //data.jetsAK8.filteredMass[i] = AK8_filteredMass[i] * scaleJES * scaleJER;
    
    // For W jet mass apply combination of both JES+JMS and JER+JMR
    // Measurements of the combined uncertainties are provided by the JMAR group
    // We use those instead of only applying the JES and JER
    double scaleJMS = get_syst_weight(W_TAG_JMS_SF, W_TAG_JMS_SF_ERR, nSigmaJES);
    double scaled_corrected_mass = AK8_softDropMassCorr[i] * scaleJMS;
    if (applySmearing) {
      double scaleJMR = get_syst_weight(AK8_JMRSmearFactor[i], AK8_JMRSmearFactorUp[i], AK8_JMRSmearFactorDown[i], nSigmaJER);
      scaled_corrected_mass *= scaleJMR;
    }
    softDropMassCorr.push_back(scaled_corrected_mass);
  }

  TVector3 dmet(0,0,0);
#if VER > 0
  // MET Uncertainties
  if      (nSigmaJES    >0) dmet += std::abs(nSigmaJES) * dmet_JESUp;
  else if (nSigmaJES    <0) dmet += std::abs(nSigmaJES) * dmet_JESDown;
  if (applySmearing) {
    if      (nSigmaJER    >0) dmet += std::abs(nSigmaJES) * dmet_JERUp;
    else if (nSigmaJER    <0) dmet += std::abs(nSigmaJES) * dmet_JERDown;
  }
  if      (nSigmaRestMET>0) dmet += std::abs(nSigmaJES) * dmet_RestUp;
  else if (nSigmaRestMET<0) dmet += std::abs(nSigmaJES) * dmet_RestDown;
#endif
  TVector3 shifted_met = met + dmet;
  data.met.Pt[0]  = shifted_met.Pt();
  data.met.Phi[0] = shifted_met.Phi();

  if (applySmearing||syst_index!=0) {
    // Recalculation of Razor variables
    // Has to be done after jet uncertainties (including smearing, even if no systematic)
    // Get selected AK4 jets (input for megajets)
    std::vector<TLorentzVector> selected_jets_AK4;
    while(data.jetsAK4.Loop()) {
      size_t i = data.jetsAK4.it;
      TLorentzVector jet_v4; jet_v4.SetPtEtaPhiE(data.jetsAK4.Pt[i], data.jetsAK4.Eta[i], data.jetsAK4.Phi[i], data.jetsAK4.E[i]);
      // Pass jet selection criteria
      if ( data.jetsAK4.looseJetID[i] == 1 &&
           data.jetsAK4.Pt[i]             >= JET_AK4_PT_CUT &&
           std::abs(data.jetsAK4.Eta[i])  <  JET_AK4_ETA_CUT ) {
        selected_jets_AK4.push_back(jet_v4);
      }
    }
    // Razor variables
    if (selected_jets_AK4.size() < 2) {
      data.evt.MR  = -9999;
      data.evt.MTR = -9999;
      data.evt.R   = -9999;
      data.evt.R2  = -9999;
    } else {
      std::vector<TLorentzVector> hemis_AK4 = Razor::CombineJets(selected_jets_AK4);
      data.evt.MR  = Razor::CalcMR(hemis_AK4[0], hemis_AK4[1]);
      data.evt.MTR = Razor::CalcMTR(hemis_AK4[0], hemis_AK4[1], shifted_met);
      data.evt.R   = data.evt.MTR/data.evt.MR;
      data.evt.R2  = data.evt.R*data.evt.R;
    }
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
AnalysisBase::get_scale_weight(const std::vector<float>& scale_Weights, const double& nSigmaScale, const unsigned int& numScale)
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
  double w_scale = 1;
  double w_scale_up = 1;   // Corresponds to 0.5 (More signal events)
  double w_scale_down = 1; // Corresponds to 2.0 (Less signal events)
  if (numScale==1) {
    // fix mu_r = 1.0, vary mu_f = 0,5, 2.0
    w_scale_up   = scale_Weights[1];
    w_scale_down = scale_Weights[0];
  } else if (numScale==2) {
    // fix mu_f = 1.0, vary mu_r = 0,5, 2.0
    w_scale_up   = scale_Weights[4];
    w_scale_down = scale_Weights[2];
  } else if (numScale==3) {
    // vary simulataneously mu_r = mu_f = 0,5, 2.0
    w_scale_up   = scale_Weights[5];
    w_scale_down = scale_Weights[3];
  }
  w_scale = get_syst_weight(w_scale, w_scale_up, w_scale_down, nSigmaScale);
  return w_scale;
}

//_______________________________________________________
//    Apply analysis cuts in the specified search region

bool
Analysis::apply_all_cuts(char region) {
  return apply_ncut(region, analysis_cuts[region].size());
}

bool
Analysis::apply_ncut(char region, unsigned int ncut) {
  if (ncut>analysis_cuts[region].size()) return 0;
  for (unsigned int i=0; i<ncut; ++i) if ( ! analysis_cuts[region][i].func() ) return 0;
  return 1;
}

// Cuts to apply/exclude by cut name
bool
Analysis::apply_cut(char region, std::string cut_name) {
  for (const auto& cut : analysis_cuts[region]) if (cut_name == cut.name) return cut.func();
  return 0;
}

bool
Analysis::apply_cuts(char region, std::vector<std::string> cuts) {
  for (const auto& cut_in_region : analysis_cuts[region]) for (const auto& cut : cuts) 
    if (cut == cut_in_region.name) if (!cut_in_region.func()) return 0;
  return 1;
}

bool
Analysis::apply_all_cuts_except(char region, std::string cut_to_skip) {
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
    std::cout<<"No cut to be skipped exsists in seaerch region \""<<region<<"\" with name: \""<<cut_to_skip<<"\""<<std::endl;
    utils::error("Analysis - the second argument for apply_all_cuts_except() is a non-sensical cut");
  }
  return result;
}

bool
Analysis::apply_all_cuts_except(char region, std::vector<std::string> cuts_to_skip) {
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
    utils::error("Analysis - the second argument for apply_all_cuts_except() contains at least one non-sensical cut");
  }
  return result;
}


// Same functions but with cut index which is faster (can use an enum, to make it nicer)
bool
Analysis::apply_cut(char region, unsigned int cut_index) { return analysis_cuts[region][cut_index].func(); }

bool
Analysis::apply_cuts(char region, std::vector<unsigned int> cuts) {
  for (const unsigned int& cut : cuts) if ( ! analysis_cuts[region][cut].func() ) return 0;
  return 1;
}

bool
Analysis::apply_all_cuts_except(char region, unsigned int cut_to_skip) {
  if (cut_to_skip>=analysis_cuts[region].size()) {
    std::cout<<"Index ("<<cut_to_skip<<") is too high for the cut to be skipped in search region '"<<region<<"'"<<std::endl;
    utils::error("Analysis::apply_all_cuts_except(char region, unsigned int cut_to_skip)");
  }
  for (unsigned int i=0, n=analysis_cuts[region].size(); i<n; ++i) {
    if (i==cut_to_skip) continue;
    if ( ! analysis_cuts[region][i].func() ) return 0;
  }
  return 1;
}

bool
Analysis::apply_all_cuts_except(char region, std::vector<unsigned int> cuts_to_skip) {
  for (unsigned int i=0, n=analysis_cuts[region].size(); i<n; ++i) {
    for (const unsigned int& cut_to_skip : cuts_to_skip) if (i!=cut_to_skip) 
      if ( ! analysis_cuts[region][i].func() ) return 0;
  }
  return 1;
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
TH2F* eff_full_ele_mediumid;
TH2F* eff_full_ele_miniiso01;
TH2D* eff_fast_ele_vetoid_miniiso01;
TH2D* eff_fast_ele_mediumid_miniiso01;
TGraphAsymmErrors* eff_full_muon_trk;
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

TGraphAsymmErrors* eff_trigger;

void AnalysisBase::init_scale_factors() {
  TString Sample(sample);
  
  // B-tagging
  // Efficiencies (Oct31 - test)
  TFile* f;
  if (Sample.Contains("FastSim"))
    f = TFile::Open("btageff/Oct31/FastSim_SMS-T5ttcc.root");
  else if (Sample.Contains("TT")||Sample.Contains("ST")) 
    f = TFile::Open("btageff/Oct31/TT_powheg-pythia8_ext4.root");
  else 
    f = TFile::Open("btageff/Oct31/QCD.root");
  eff_btag_b_loose  = ((TH2D*)f->Get("btageff_b_loose"))->ProfileX();
  eff_btag_c_loose  = ((TH2D*)f->Get("btageff_c_loose"))->ProfileX();
  eff_btag_l_loose  = ((TH2D*)f->Get("btageff_l_loose"))->ProfileX();
  eff_btag_b_medium = ((TH2D*)f->Get("btageff_b_medium"))->ProfileX();
  eff_btag_c_medium = ((TH2D*)f->Get("btageff_c_medium"))->ProfileX();
  eff_btag_l_medium = ((TH2D*)f->Get("btageff_l_medium"))->ProfileX();
  eff_btag_b_loose  ->SetDirectory(0);
  eff_btag_c_loose  ->SetDirectory(0);
  eff_btag_l_loose  ->SetDirectory(0);
  eff_btag_b_medium ->SetDirectory(0);
  eff_btag_c_medium ->SetDirectory(0);
  eff_btag_l_medium ->SetDirectory(0);
  f->Close();
  // Moriond17 SFs
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco?rev=5#Supported_Algorithms_and_Operati
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
  // Moriond17 - Reconstruction efficiency. Scale factors for 80X
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=38#Electron_efficiencies_and_scale
  // Muon tracking efficiency SFs
  // https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2?rev=21#Tracking_efficiency_provided_by
  // Moriond17 FullSim - SUSY
  // https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#Data_leading_order_FullSim_MC_co
  // https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#Data_leading_order_FullSim_M_AN1
  // Spring16 FastSim to FullSim
  // https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=179#FullSim_FastSim_TTBar_MC_compari
  // https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=179#FullSim_FastSim_TTBar_MC_com_AN1
  eff_full_ele_reco               = utils::getplot_TH2F("scale_factors/electron/egammaEffi.txt_EGM2D.root","EGamma_SF2D", "ele1");
  eff_full_ele_vetoid             = utils::getplot_TH2F("scale_factors/electron/scaleFactors.root","GsfElectronToCutBasedSpring15V", "ele2");
  eff_full_ele_mediumid           = utils::getplot_TH2F("scale_factors/electron/scaleFactors.root","GsfElectronToCutBasedSpring15M", "ele3");
  eff_full_ele_miniiso01          = utils::getplot_TH2F("scale_factors/electron/scaleFactors.root","MVAVLooseElectronToMini"       , "ele4");
  eff_fast_ele_vetoid_miniiso01   = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_vetoCB_mini01.root",  "histo2D"      , "ele5");
  eff_fast_ele_mediumid_miniiso01 = utils::getplot_TH2D("scale_factors/electron/fastsim/sf_el_mediumCB_mini01.root","histo2D"      , "ele6");
  eff_full_muon_trk   		  = utils::getplot_TGraphAsymmErrors("scale_factors/muon/ratios.root", "ratio_eta", "mu1");
  eff_full_muon_looseid		  = utils::getplot_TH2F("scale_factors/muon/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root", "SF", "mu2");
  eff_full_muon_mediumid	  = utils::getplot_TH2F("scale_factors/muon/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root","SF", "mu3");
  eff_full_muon_miniiso04	  = utils::getplot_TH2F("scale_factors/muon/TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta.root",  "SF", "mu4");
  eff_full_muon_miniiso02	  = utils::getplot_TH2F("scale_factors/muon/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root", "SF", "mu5");
  eff_full_muon_looseip2d	  = utils::getplot_TH2F("scale_factors/muon/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root",    "SF", "mu6");
  eff_full_muon_tightip2d	  = utils::getplot_TH2F("scale_factors/muon/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root",    "SF", "mu7");
  eff_fast_muon_looseid		  = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_loose.root",          "histo2D", "mu8");
  eff_fast_muon_mediumid	  = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_medium.root",         "histo2D", "mu9");
  eff_fast_muon_miniiso04	  = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_looseID_mini04.root", "histo2D", "mu10");
  eff_fast_muon_miniiso02	  = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_mediumID_mini02.root","histo2D", "mu11");
  eff_fast_muon_looseip2d	  = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_looseIP2D.root",      "histo2D", "mu12");
  eff_fast_muon_tightip2d         = utils::getplot_TH2D("scale_factors/muon/fastsim/sf_mu_tightIP2D.root",      "histo2D", "mu13");

  // Trigger efficiency
  // TODO: Update for latest ntuples (Jan12 or later)
  TH1D* num = utils::getplot_TH1D("trigger_eff/Oct21_Golden_JSON/SingleLepton.root", "h_HT_pre_pass", "trig1");
  TH1D* den = utils::getplot_TH1D("trigger_eff/Oct21_Golden_JSON/SingleLepton.root", "h_HT_pre",      "trig2");
  eff_trigger = new TGraphAsymmErrors(num, den);
}


double AnalysisBase::calc_top_tagging_sf(DataStruct& data, const double& nSigmaTopTagSF) {
  double w = 1.0;

  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    double pt = data.jetsAK8.Pt[i];
    if (passHadTopTag[i]) {
      // Top-tagged AK8 jets
      if (pt >= 400 && pt < 550)
	w *= get_syst_weight(TOP_TAG_SF_LOW, TOP_TAG_SF_LOW_ERR, nSigmaTopTagSF);
      else if (pt >= 550)
	w *= get_syst_weight(TOP_TAG_SF_HIGH, TOP_TAG_SF_HIGH_ERR, nSigmaTopTagSF);
    }
  }

  return w;
}


double AnalysisBase::calc_w_tagging_sf(DataStruct& data, const double& nSigmaWTagSF) {
  double w = 1.0;

  while(data.jetsAK8.Loop()) if (passTightWTag[data.jetsAK8.it])
    w *= get_syst_weight(W_TAG_EFF_SF, W_TAG_EFF_SF_ERR, nSigmaWTagSF);

  return w;
}

std::pair<double, double> AnalysisBase::calc_b_tagging_sf(DataStruct& data, const double& nSigmaBTagSF, const bool& isFastSim) {

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

	sf_loose      *= get_syst_weight(sf_loose_cen,  sf_loose_up,  sf_loose_down,  nSigmaBTagSF);
	sf_medium     *= get_syst_weight(sf_medium_cen, sf_medium_up, sf_medium_down, nSigmaBTagSF);
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

std::pair<double, double> AnalysisBase::calc_ele_sf(DataStruct& data, const double& nSigmaEleRecoSF, const double& nSigmaEleIDSF, const double& nSigmaEleIsoSF, const double& nSigmaEleFastSimSF, const bool& isFastSim) {
  double sf, sf_err;
  double weight_veto  = 1.0, weight_select = 1.0;
  while(data.ele.Loop()) {
    size_t i       = data.ele.it;
    double pt      = data.ele.Pt[i];
    double eta     = data.ele.Eta[i];
    float abseta   = std::abs(eta);
    float miniIso  = data.ele.MiniIso[i]/data.ele.Pt[i];
    float absd0    = std::abs(data.ele.Dxy[i]);
    float absdz    = std::abs(data.ele.Dz[i]);
    bool id_veto   = (data.ele.vidVetonoiso[i] == 1.0);
    bool id_select = (data.ele.vidMediumnoiso[i] == 1.0);
    // Apply reconstruction scale factor - Warning! strange binning (pt vs eta)
    utils::geteff2D(eff_full_ele_reco, eta, pt, sf, sf_err);
    // If pt is below 20 or above 80 GeV increase error by 1%
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2?rev=38#Electron_efficiencies_and_scale
    if (pt<20||pt>=80) sf_err = std::sqrt(sf_err*sf_err + 0.01+0.01);
    weight_veto   *= get_syst_weight(sf, sf_err, nSigmaEleRecoSF);
    weight_select *= get_syst_weight(sf, sf_err, nSigmaEleRecoSF);

    // Veto Electrons
    if ( id_veto &&
	 pt      >= ELE_VETO_PT_CUT &&
	 abseta  <  ELE_VETO_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_ele_vetoid, pt, eta, sf, sf_err);
      weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      if ( miniIso <  ELE_VETO_MINIISO_CUT ) {
	// Apply Isolation scale factor (No separate SF for IP cut which is included)
	if ( absd0   <  ELE_VETO_IP_D0_CUT &&
	     absdz   <  ELE_VETO_IP_DZ_CUT ) {
	  utils::geteff2D(eff_full_ele_miniiso01, pt, eta, sf, sf_err);
	  weight_veto *= get_syst_weight(sf, sf_err, nSigmaEleIsoSF);
	}
	// Apply FastSim to FullSim scale factor
	// Error in histo is statistical, apply 2% instead
	if (isFastSim) {
	  utils::geteff2D(eff_fast_ele_vetoid_miniiso01, pt, eta, sf, sf_err);
	  weight_veto *= get_syst_weight(sf, 0.02, nSigmaEleFastSimSF);
	}
      }
    }


    // Selected Electrons
    if ( id_select &&
	 pt        >= ELE_SELECT_PT_CUT &&
	 abseta    <  ELE_SELECT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_ele_mediumid, pt, eta, sf, sf_err);
      weight_select *= get_syst_weight(sf, sf_err, nSigmaEleIDSF);
      if ( miniIso   <  ELE_SELECT_MINIISO_CUT ) {
	// Apply FastSim to FullSim scale factor on top
	// Error in histo is statistical, apply 2% instead
	if (isFastSim) {
	  utils::geteff2D(eff_fast_ele_mediumid_miniiso01, pt, eta, sf, sf_err);
	  weight_select *= get_syst_weight(sf, 0.02, nSigmaEleFastSimSF);
	}
	// Apply Isolation scale factor (No separate SF for IP cut)
	if ( absd0     <  ELE_SELECT_IP_D0_CUT &&
	     absdz     <  ELE_SELECT_IP_DZ_CUT ) {
	  utils::geteff2D(eff_full_ele_miniiso01, pt, eta, sf, sf_err);
	  weight_select *= get_syst_weight(sf, sf_err, nSigmaEleIsoSF);
	}
      }
    }

  }
  
  return std::make_pair(weight_veto, weight_select);
}

std::pair<double, double> AnalysisBase::calc_muon_sf(DataStruct& data, const double& nSigmaMuonTrkSF, const double& nSigmaMuonIDSF, const double& nSigmaMuonIsoSF, const double& nSigmaMuonIPSF, const double& nSigmaMuonFastSimIDSF, const double& nSigmaMuonFastSimIsoSF, const double& nSigmaMuonFastSimIPSF, const bool& isFastSim) {
  double sf, sf_err, sf_err_down, sf_err_up;
  double weight_veto  = 1.0, weight_select = 1.0;
  while(data.mu.Loop()) {
    size_t i       = data.mu.it;
    double pt      = data.mu.Pt[i];
    double eta     = data.mu.Eta[i];
    float abseta   = std::abs(eta);
    float miniIso  = data.mu.MiniIso[i]/data.mu.Pt[i];
    float absd0    = std::abs(data.mu.Dxy[i]);
    float absdz    = std::abs(data.mu.Dz[i]);
    bool id_veto   = (data.mu.IsLooseMuon[i] == 1.0);
    bool id_select = (data.mu.IsMediumMuon[i] == 1.0);
    // Apply tracking efficiency scale factor
    if (pt>=10) {
      utils::geteff_AE(eff_full_muon_trk, eta, sf, sf_err_down, sf_err_up);
      weight_veto   *= get_syst_weight(sf, sf-sf_err_down, sf+sf_err_up, nSigmaMuonTrkSF);
      weight_select *= get_syst_weight(sf, sf-sf_err_down, sf+sf_err_up, nSigmaMuonTrkSF);
    }
    
    // Veto Muons
    if ( id_veto &&
	 pt      >= MU_VETO_PT_CUT &&
	 abseta  <  MU_VETO_ETA_CUT) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_muon_looseid, pt, eta, sf, sf_err);
      weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonIDSF);
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_looseid, pt, eta, sf, sf_err);
	weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonFastSimIDSF);
      }
      // Apply Isolation scale factor
      if ( miniIso <  MU_VETO_MINIISO_CUT ) {
	utils::geteff2D(eff_full_muon_miniiso04, pt, eta, sf, sf_err);
	weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonIsoSF);
	if (isFastSim) {
	  utils::geteff2D(eff_fast_muon_miniiso04, pt, eta, sf, sf_err);
	  weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonFastSimIsoSF);
	}
      }
      // Apply IP efficiency scale factor
      if ( absd0   <  MU_VETO_IP_D0_CUT &&
	   absdz   <  MU_VETO_IP_DZ_CUT ) {
	utils::geteff2D(eff_full_muon_looseip2d, pt, eta, sf, sf_err);
	weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonIPSF);
	if (isFastSim) {
	  utils::geteff2D(eff_fast_muon_looseip2d, pt, eta, sf, sf_err);
	  weight_veto *= get_syst_weight(sf, sf_err, nSigmaMuonFastSimIPSF);
	}
      }
    }


    // Selected muons
    if ( id_select &&
	 pt      >= MU_SELECT_PT_CUT &&
	 abseta  <  MU_SELECT_ETA_CUT ) {
      // Apply ID scale factor
      utils::geteff2D(eff_full_muon_mediumid, pt, eta, sf, sf_err);
      weight_select *= get_syst_weight(sf, sf_err, nSigmaMuonIDSF);
      if (isFastSim) {
	utils::geteff2D(eff_fast_muon_mediumid, pt, eta, sf, sf_err);
	weight_select *= get_syst_weight(sf, sf_err, nSigmaMuonFastSimIDSF);
      }
      // Apply Isolation scale factor
      if ( miniIso <  MU_SELECT_MINIISO_CUT ) {
	utils::geteff2D(eff_full_muon_miniiso02, pt, eta, sf, sf_err);
	weight_select *= get_syst_weight(sf, sf_err, nSigmaMuonIsoSF);
	if (isFastSim) {
	  utils::geteff2D(eff_fast_muon_miniiso02, pt, eta, sf, sf_err);
	  weight_select *= get_syst_weight(sf, sf_err, nSigmaMuonFastSimIsoSF);
	}
      }
      // Apply IP efficiency scale factor
      if ( absd0   <  MU_SELECT_IP_D0_CUT &&
	   absdz   <  MU_SELECT_IP_DZ_CUT ) {
	utils::geteff2D(eff_full_muon_tightip2d, pt, eta, sf, sf_err);
	weight_select *= get_syst_weight(sf, sf_err, nSigmaMuonIPSF);
	if (isFastSim) {
	  utils::geteff2D(eff_fast_muon_tightip2d, pt, eta, sf, sf_err);
	  weight_select *= get_syst_weight(sf, sf_err, nSigmaMuonFastSimIPSF);
	}
      }
    }

  }

  return std::make_pair(weight_veto, weight_select);
}


double AnalysisBase::calc_trigger_efficiency(DataStruct& data, const double& nSigmaTrigger) {
  double eff, err_down, err_up;
  utils::geteff_AE(eff_trigger, AK4_Ht, eff, err_down, err_up);
  double w = get_syst_weight(eff, eff-err_down, eff+err_up, nSigmaTrigger);
  return w;
}
