#include <iostream>
#include <functional>
#include <map>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TStopwatch.h"
#include <thread>
#include <chrono>

#include "utils.h"
#include "DataStruct.h"
#include "GluinoXSec.h"
#include "StopXSec.h"

// _____________________________________________________________
//        AnalysisBase: Methods common in all analysis

class AnalysisBase
{
public:
  AnalysisBase();
  ~AnalysisBase();

  typedef struct Cut { std::string name; std::function<bool()> func; } Cut;
  std::vector<Cut> baseline_cuts;

  // Functions used by the Analyzer
  void define_preselections(const DataStruct&, const bool&, const bool&);

  void calculate_common_variables(DataStruct&, const unsigned int&);

  void init_common_histos();

  double get_xsec_from_ntuple(const std::vector<std::string>&, const std::string&);

  double get_xsec_from_txt_file(const std::string&, const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const std::string&);

  void calc_weightnorm_histo_from_ntuple(const std::vector<std::string>&, const double&, const std::vector<std::string>&,
					 const std::vector<std::string>&, bool);

  void init_pileup_reweightin(const std::string&, const std::string&, const std::vector<std::string>&);

  double get_pileup_weight(const int&, const double&);

  void rescale_jets(DataStruct&, const unsigned int&, const double&);

  double get_w_tagging_sf(DataStruct&, const double&);

  double get_b_tagging_sf(DataStruct&, const double&);

  double get_top_tagging_sf(DataStruct&, const double&);

  double get_ht_weight(DataStruct&, const double&);

  double get_alphas_weight(const std::vector<float>&, const double&, const int&);

  double get_scale_weight(const std::vector<float>&, const double&, const unsigned int&);

  double get_syst_weight(const double&, const double&, const double&, const double&);

  double get_syst_weight(const double&, const double&, const double&);

  void benchmarking(const int&, const int&, const bool);

private:

  TStopwatch *sw_1k_, *sw_10k_, *sw_job_;

  void moderate_job_(TH1D* h, double, int, int);

};

// _____________________________________________________________
//         Analysis: Analysis specific methods/histos

class Analysis : public AnalysisBase
{
public:
  Analysis();
  ~Analysis();

  void calculate_variables(DataStruct&, const unsigned int&);

  double get_analysis_weight(DataStruct&);

  bool pass_skimming(DataStruct&);

  void define_selections(const DataStruct&);

  virtual bool signal_selection(const DataStruct&);

  void define_histo_options(const double&, const DataStruct&, const unsigned int&, const unsigned int&, std::string, bool);

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
//                       Constructor
AnalysisBase::AnalysisBase() {
  sw_1k_  = new TStopwatch;
  sw_10k_ = new TStopwatch;
  sw_job_ = new TStopwatch;
}


//_______________________________________________________
//                       Destructor
AnalysisBase::~AnalysisBase() {
  delete sw_1k_;
  delete sw_10k_;
  delete sw_job_;
}


//_______________________________________________________
//                 Define baseline cuts
void
AnalysisBase::define_preselections(const DataStruct& data, const bool& isData, const bool& isSignal)
{ 
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
  baseline_cuts.push_back({ .name="Clean_CSC_Halo_Tight",    .func = [&data,isSignal] { return isSignal ? 1 : data.filter.globalTightHalo2016Filter; } });
  baseline_cuts.push_back({ .name="Clean_HBHE_Noise",        .func = [&data] { return data.filter.HBHENoiseFilter; } });
  baseline_cuts.push_back({ .name="Clean_HBHE_IsoNoise",     .func = [&data] { return data.filter.HBHENoiseIsoFilter; } });
  baseline_cuts.push_back({ .name="Clean_Ecal_Dead_Cell_TP", .func = [&data] { return data.filter.EcalDeadCellTriggerPrimitiveFilter; } });
  baseline_cuts.push_back({ .name="Clean_EE_Bad_Sc",         .func = [&data,isData] { return isData ? data.filter.eeBadScFilter : 1; } });
  // Not in MiniAODv2 (producer added)
  baseline_cuts.push_back({ .name="Clean_Bad_Muon",          .func = [&data] { return data.filter.BadPFMuonFilter; } });
  baseline_cuts.push_back({ .name="Clean_Bad_Charged",       .func = [&data] { return data.filter.BadChargedCandidateFilter; } });
}

//_______________________________________________________
//                 Define common variables

/*
  Jet ID:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

  For AK4 Jet Selection Choose:
  - Loose jet ID
  - pt > 30
  - |eta| < 2.4

  For AK4 Jet Selection Choose:
  - Tight jet ID
  - pt > 200
  - |eta| < 2.4

*/
#define JET_AK4_PT_CUT  30
#define JET_AK4_ETA_CUT 2.4
#define JET_AK8_PT_CUT  200
#define JET_AK8_ETA_CUT 2.4

/*
  B tagging working points:
  https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80X

  Latest WPs/SFs:
  https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80X?rev=18#Supported_Algorithms_and_Operati

  Choose:
  - CombinedSecondaryVertex v2
  - CSVv2 >= 0.46 (Loose - for Veto)
  - CSVv2 >= 0.8  (Medium - for Tag)

*/
#define B_SUBJET_CSV_LOOSE_CUT 0.460
#define B_CSV_LOOSE_CUT        0.460
#define B_CSV_MEDIUM_CUT       0.800
#define B_CSV_TIGHT_CUT        0.935

/* 
   W tagging working points:
   https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=36#Working_points_and_scale_factors

   Latest WPs/SFs:
   https://indico.cern.ch/event/559594/contributions/2258668/attachments/1317079/1973381/JMAR_Meeting_Gelli.pdf

   Choose:
   - Tight jet ID
   Tight Tag selection e(S) = 66.5%:
   - AK8 Puppi jets
   - pt > 200
   - eta < 2.4
   - 65 <= SD Mass < 105
   - tau21 < 0.45

*/

#define W_PT_CUT            200
#define W_ETA_CUT           2.4
#define W_SD_MASS_CUT_LOW   65
#define W_SD_MASS_CUT_HIGH  105
#define W_TAU21_LOOSE_CUT   0.6
#define W_TAU21_TIGHT_CUT   0.45

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
*/
/*
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
  Latest Electron IDs (Veto/Tight):
  [1] POG Veto/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=38#Working_points_for_2016_data_for

  Latest Isolation WPs:
  [2] SUSY MiniIso Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=166#ID_IP_ISO_AN1

  Latest Impact Point Cut (Loose/Tight):
  [3] SUSY Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=166#ID_IP_ISO
  [4] POG  Tight - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=38#Offline_selection_criteria

  For Veto Choose:
  - Cut based Veto ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (tight WP [2])
  - pt >= 10    (Added in B2G ntuple level - I can loosen it if needed in the next production)
  - |eta| < 2.5
  - |d0| < 0.2, |dz| < 0.5 (Loose IP cut Recommended by SUSY Group [3])

  For Selection Choose:
  - Cut based Tight ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (tight WP [2])
  - pt >= 10
  - |eta| < 2.5
  - |d0| < 0.05 (0.1), |dz| < 0.1 (0.2) (Tight IP cut Recommended by EGamma POG [4])

*/
/*
  Latest Muon IDs (Loose/Tight):
  [1] POG Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=26#Muon_Identification

  Latest Isolation WPs:
  [2] SUSY MiniISo Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=166#ID_IP_ISO

  Latest Impact Point Cut (Loose/Tight):
  [3] SUSY Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=166#ID_IP_ISO
  [4] POG Tight  - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=38#Offline_selection_criteria
  
  For Veto Choose:
  - POG recommended Loose ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.2 (tight WP [2])
  - pt >= 5
  - |eta| < 2.4
  - |d0| < 0.2, |dz| < 0.5 (Loose IP cut Recommended by SUSY Group [3])

  For Selection Choose:
  - POG recommended Tight ID (No Iso, Loose IP)
  - Mini-Isolation (EA)/pt < 0.2 (tight WP [2])
  - pt >= 5
  - |eta| < 2.4
  - |d0| < 0.05, |dz| < 0.1 (Tight IP cut Recommended by EGamma POG [4])

  Not (yet) used - variable needs to be added for next production:
  - Impact point: |d0| < 0.2, |dz| < 0.5

*/

#define ELE_VETO_PT_CUT        10  // Same cut as in B2G ntuples
#define ELE_VETO_ETA_CUT       2.5
#define ELE_VETO_MINIISO_CUT   0.1
#define ELE_VETO_IP_D0_CUT     0.2
#define ELE_VETO_IP_DZ_CUT     0.5

#define ELE_TIGHT_PT_CUT       10  // Same cut as in B2G ntuples
#define ELE_TIGHT_ETA_CUT      2.5
#define ELE_TIGHT_MINIISO_CUT  0.1
#define ELE_TIGHT_IP_EB_D0_CUT 0.05
#define ELE_TIGHT_IP_EB_DZ_CUT 0.1
#define ELE_TIGHT_IP_EE_D0_CUT 0.1
#define ELE_TIGHT_IP_EE_DZ_CUT 0.2

#define MU_VETO_PT_CUT         5
#define MU_VETO_ETA_CUT        2.4
#define MU_VETO_MINIISO_CUT    0.2
#define MU_VETO_IP_D0_CUT      0.2
#define MU_VETO_IP_DZ_CUT      0.5

#define MU_TIGHT_PT_CUT        10
#define MU_TIGHT_ETA_CUT       2.4
#define MU_TIGHT_MINIISO_CUT   0.2
#define MU_TIGHT_IP_D0_CUT     0.2
#define MU_TIGHT_IP_DZ_CUT     0.5

// AK4 jets
/*
  convention:

  iObject  -  gives the index of the nth selected object in the original collection

  example:
  for (size_t i=0; i<nJet; ++i) h_pt->Fill( data.jetsAK4Puppi.Pt[iJet[i]] );  
  or
  if (nJet>0) vh_pt[0] -> Fill( data.jetsAK4Puppi.Pt[iJet[0]] );
  if (nJet>1) vh_pt[1] -> Fill( data.jetsAK4Puppi.Pt[iJet[1]] );


  itObject  -  gives the index in the selected collection

  example:
  for (size_t it=0; it<data.jetsAK4Puppi.size; ++it)
    if (passLooseJet[it]) vh_pt[itJet[it]]->Fill( data.jetsAK4Puppi.Pt[it] );

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
double AK4Puppi_Ht;
double minDeltaPhi; // Min(DeltaPhi(Jet_i, MET)), i=1,2,3

// AK8 jets
std::vector<size_t > iJetAK8;
std::vector<size_t > iWPreTag;
std::vector<size_t > iLooseWTag;
std::vector<size_t > iTightWTag;
std::vector<size_t > itJetAK8;
std::vector<size_t > itWPreTag;
std::vector<size_t > itLooseWTag;
std::vector<size_t > itTightWTag;
std::vector<double> tau21;
std::vector<double> tau31;
std::vector<double> tau32;
std::vector<double> maxSubjetCSV;
std::vector<bool> passSubjetBTag;
std::vector<bool> passTightJetAK8;
std::vector<bool> passWPreTag;
std::vector<bool> passLooseWTag;
std::vector<bool> passTightWTag;
std::vector<bool> passHadTopPreTag;
std::vector<bool> passHadTopTag;
unsigned int nJetAK8;
unsigned int nWPreTag;
unsigned int nLooseWTag;
unsigned int nTightWTag;
unsigned int nSubjetBTag;
unsigned int nHadTopTag;
unsigned int nHadTopPreTag;
double AK8Puppi_Ht;

// Event Letpons
std::vector<size_t > iEleTight;
std::vector<size_t > iMuTight;
std::vector<size_t > itEleTight;
std::vector<size_t > itMuTight;
std::vector<bool> passEleVeto;
std::vector<bool> passMuVeto;
std::vector<bool> passEleTight;
std::vector<bool> passMuTight;
unsigned int nEleVeto;
unsigned int nEleTight;
unsigned int nMuVeto;
unsigned int nMuTight;
unsigned int nLepVeto;
unsigned int nLepTight;
double MT;

void
AnalysisBase::calculate_common_variables(DataStruct& data, const unsigned int& syst_index)
{
  // It only makes sense to calculate certain variables only once if they don't depend on jet energy
  if (syst_index == 0) {

    // Loop on AK8 Puppi jets
    tau21         .assign(data.jetsAK8Puppi.size, 9999);
    tau31         .assign(data.jetsAK8Puppi.size, 9999);
    tau32         .assign(data.jetsAK8Puppi.size, 9999);
    maxSubjetCSV  .assign(data.jetsAK8Puppi.size, 0);
    passSubjetBTag.assign(data.jetsAK8Puppi.size, 0);
    nSubjetBTag = 0;
    while(data.jetsAK8Puppi.Loop()) {
      size_t i = data.jetsAK8Puppi.it;
      // N-subjettiness
      if (data.jetsAK8Puppi.tau1[i]>0) tau21[i] = data.jetsAK8Puppi.tau2[i]/data.jetsAK8Puppi.tau1[i];
      if (data.jetsAK8Puppi.tau1[i]>0) tau31[i] = data.jetsAK8Puppi.tau3[i]/data.jetsAK8Puppi.tau1[i];
      if (data.jetsAK8Puppi.tau2[i]>0) tau32[i] = data.jetsAK8Puppi.tau3[i]/data.jetsAK8Puppi.tau2[i];
      // Maximum Subjet btag discriminator
      maxSubjetCSV[i] = -9999;
      int i_sj0 = data.jetsAK8Puppi.vSubjetIndex0[i], i_sj1 = data.jetsAK8Puppi.vSubjetIndex1[i];
      if (i_sj0 != -1) if (data.subjetsAK8Puppi.CSVv2[i_sj0] > maxSubjetCSV[i]) maxSubjetCSV[i] = data.subjetsAK8Puppi.CSVv2[i_sj0];
      if (i_sj1 != -1) if (data.subjetsAK8Puppi.CSVv2[i_sj1] > maxSubjetCSV[i]) maxSubjetCSV[i] = data.subjetsAK8Puppi.CSVv2[i_sj1];
#if USE_BTAG == 1
      if (passSubjetBTag[i] = (maxSubjetCSV[i] >= TOP_BTAG_CSV) ) nSubjetBTag++;
#else
      if (passSubjetBTag[i] = (maxSubjetCSV[i] >= B_SUBJET_CSV_LOOSE_CUT) ) nSubjetBTag++;
#endif
    }

    // Event Letpons
    iEleTight   .clear();
    itEleTight  .assign(data.ele.size, (size_t)-1);
    passEleVeto .assign(data.ele.size, 0);
    passEleTight.assign(data.ele.size, 0);
    nEleVeto = nEleTight = 0;
    while(data.ele.Loop()) {
      size_t i = data.ele.it;
      float pt = data.ele.Pt[i];
      float abseta = fabs(data.ele.Eta[i]);
      float miniIso = data.ele.MiniIso[i]/data.ele.Pt[i];
      float absd0 = fabs(data.ele.Dxy[i]);
      float absdz = fabs(data.ele.Dz[i]);
      bool id_veto = (data.ele.vidVetonoiso[i] == 1.0);
      bool id_tight = (data.ele.vidTightnoiso[i] == 1.0);
      // Veto
      if (passEleVeto[i] = 
	  ( id_veto &&
	    pt      >= ELE_VETO_PT_CUT &&
	    abseta  <  ELE_VETO_ETA_CUT &&
	    miniIso <  ELE_VETO_MINIISO_CUT &&
	    absd0   <  ELE_VETO_IP_D0_CUT &&
	    absdz   <  ELE_VETO_IP_DZ_CUT) )
	nEleVeto++;
      // Tight
      if (passEleTight[i] = 
	  ( id_tight &&
	    pt        >= ELE_TIGHT_PT_CUT &&
	    abseta    <  ELE_TIGHT_ETA_CUT &&
	    miniIso   <  ELE_TIGHT_MINIISO_CUT &&
	    (abseta<1.479 ?
	     (absd0   <  ELE_TIGHT_IP_EB_D0_CUT &&
	      absdz   <  ELE_TIGHT_IP_EB_DZ_CUT) :
	     (absd0   <  ELE_TIGHT_IP_EE_D0_CUT &&
	      absdz   <  ELE_TIGHT_IP_EE_DZ_CUT) ) ) ) {
	iEleTight.push_back(i);
	itEleTight[i] = nEleTight++;
      }
    }

    // Number of Veto/Tight Muons
    iMuTight    .clear();
    itMuTight   .assign(data.mu.size,  (size_t)-1);
    passMuVeto  .assign(data.mu.size,  0);
    passMuTight .assign(data.mu.size,  0);
    nMuVeto = nMuTight = 0;
    while(data.mu.Loop()) {
      size_t i = data.mu.it;
      float pt = data.mu.Pt[i];
      float abseta = fabs(data.mu.Eta[i]);
      float miniIso = data.mu.MiniIso[i]/data.mu.Pt[i];
      float absd0 = fabs(data.mu.Dxy[i]);
      float absdz = fabs(data.mu.Dz[i]);
      bool id_veto = (data.mu.IsLooseMuon[i] == 1.0);
      bool id_tight = (data.mu.IsTightMuon[i] == 1.0);
      // Veto
      if (passMuVeto[i] =
	  (id_veto &&
	   pt      >= MU_VETO_PT_CUT &&
	   abseta  <  MU_VETO_ETA_CUT &&
	   miniIso <  MU_VETO_MINIISO_CUT &&
	   absd0   <  MU_VETO_IP_D0_CUT &&
	   absdz   <  MU_VETO_IP_DZ_CUT) )
	nMuVeto++;
      // Tight
      if (passMuTight[i] =
	  ( id_tight &&
	    pt      >= MU_TIGHT_PT_CUT &&
	    abseta  <  MU_TIGHT_ETA_CUT &&
	    miniIso <  MU_TIGHT_MINIISO_CUT &&
	    absd0   <  MU_TIGHT_IP_D0_CUT &&
	    absdz   <  MU_TIGHT_IP_DZ_CUT) ) {
	iMuTight.push_back(i);
	itMuTight[i] = nMuTight++;
      }
    }

    nLepVeto  = nEleVeto  + nMuVeto;
    nLepTight = nEleTight + nMuTight;

    // MT
    MT = 9999;
    if (nLepTight==1) {
      if (nEleTight==1) {
	MT = sqrt( 2*data.ele.Pt[iEleTight[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.ele.Phi[iEleTight[0]])) );
      } else if (nMuTight==1) {
	MT = sqrt( 2*data.mu.Pt[iMuTight[0]]*data.met.Pt[0] * (1 - std::cos(data.met.Phi[0]-data.mu.Phi[iMuTight[0]])) );
      }
    }
  }

  // Rest of the vairables need to be recalculated each time the jet energy is changed
  // eg. Jet selection, W/top tags, HT (obviously), etc. that depends on jet pt

  // AK4 Puppi jets
  iJet       .clear();
  iLooseBTag .clear();
  iMediumBTag.clear();
  iTightBTag .clear();
  itJet         .assign(data.jetsAK4Puppi.size, (size_t)-1);
  itLooseBTag   .assign(data.jetsAK4Puppi.size, (size_t)-1);
  itMediumBTag  .assign(data.jetsAK4Puppi.size, (size_t)-1);
  itTightBTag   .assign(data.jetsAK4Puppi.size, (size_t)-1);
  passLooseJet  .assign(data.jetsAK4Puppi.size, 0);
  passLooseBTag .assign(data.jetsAK4Puppi.size, 0);
  passMediumBTag.assign(data.jetsAK4Puppi.size, 0);
  passTightBTag .assign(data.jetsAK4Puppi.size, 0);
  nJet = 0;
  nLooseBTag  = 0;
  nMediumBTag = 0;
  nTightBTag  = 0;
  AK4Puppi_Ht = 0;
  minDeltaPhi = 9999;
  while(data.jetsAK4Puppi.Loop()) {
    size_t i = data.jetsAK4Puppi.it;
    // Jet ID
    if ( passLooseJet[i] = 
	 ( data.jetsAK4Puppi.looseJetID[i] == 1 &&
	   data.jetsAK4Puppi.Pt[i]         >= JET_AK4_PT_CUT &&
	   fabs(data.jetsAK4Puppi.Eta[i])  <  JET_AK4_ETA_CUT ) ) {
      iJet.push_back(i);
      itJet[i] = nJet++;

      // B tagging
      if (passLooseBTag[i]  = (data.jetsAK4Puppi.CSVv2[i] >= B_CSV_LOOSE_CUT ) ) {
	iLooseBTag.push_back(i);
	itLooseBTag[i] = nLooseBTag++;
      }
      if (passMediumBTag[i] = (data.jetsAK4Puppi.CSVv2[i] >= B_CSV_MEDIUM_CUT) ) {
	iMediumBTag.push_back(i);
	itMediumBTag[i] = nMediumBTag++;
      }
      if (passTightBTag[i]  = (data.jetsAK4Puppi.CSVv2[i] >= B_CSV_TIGHT_CUT ) ) {
	iTightBTag.push_back(i);
	itTightBTag[i] = nTightBTag++;
      }

      // Ht
      AK4Puppi_Ht += data.jetsAK4Puppi.Pt[data.jetsAK4Puppi.it];    

      // minDeltaPhi
      if (nJet<=3) {
	double dphi = fabs(TVector2::Phi_mpi_pi(data.met.Phi[0] - data.jetsAK4Puppi.Phi[i]));
	if (dphi<minDeltaPhi) minDeltaPhi = dphi;
      }



    } // End Jet Selection
  } // End AK4 Jet Loop
  
  // AK8 Puppi jets
  iJetAK8   .clear();
  iWPreTag  .clear();
  iLooseWTag.clear();
  iTightWTag.clear();
  itJetAK8   .assign(data.jetsAK4Puppi.size, (size_t)-1);
  itWPreTag  .assign(data.jetsAK4Puppi.size, (size_t)-1);
  itLooseWTag.assign(data.jetsAK4Puppi.size, (size_t)-1);
  itTightWTag.assign(data.jetsAK4Puppi.size, (size_t)-1);
  passTightJetAK8 .assign(data.jetsAK8Puppi.size, 0);
  passWPreTag     .assign(data.jetsAK8Puppi.size, 0);
  passLooseWTag   .assign(data.jetsAK8Puppi.size, 0);
  passTightWTag   .assign(data.jetsAK8Puppi.size, 0);
  passHadTopPreTag.assign(data.jetsAK8Puppi.size, 0);
  passHadTopTag   .assign(data.jetsAK8Puppi.size, 0);
  nJetAK8       = 0;
  nWPreTag      = 0;
  nLooseWTag    = 0;
  nTightWTag    = 0;
  nSubjetBTag   = 0;
  nHadTopTag    = 0;
  nHadTopPreTag = 0;
  AK8Puppi_Ht   = 0;
  while(data.jetsAK8Puppi.Loop()) {
    size_t i = data.jetsAK8Puppi.it;
    // Jet ID
    if ( passTightJetAK8[i] = 
	 ( data.jetsAK8Puppi.tightJetID[i] == 1 &&
	   data.jetsAK8Puppi.Pt[i]         >= JET_AK8_PT_CUT &&
	   fabs(data.jetsAK8Puppi.Eta[i])  <  JET_AK8_ETA_CUT ) ) {
      iJetAK8.push_back(i);
      itJetAK8[i] = nJetAK8++;

      // Tagging Variables
      double pt      = data.jetsAK8Puppi.Pt[i];
      double abseta  = data.jetsAK8Puppi.Eta[i];
      double sd_mass = data.jetsAK8Puppi.softDropMass[i];
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

      // Ht
      AK8Puppi_Ht += data.jetsAK8Puppi.Pt[i];

    } // End Jet Selection
  } // End AK4 Jet Loop
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
  // benchmarking histos
  h_read_speed_1k              = new TH1D("read_speed_1k",          ";Read speed (Events/s);Measurement/1k Event",  1000,0,10000);
  h_read_speed_10k             = new TH1D("read_speed_10k",         ";Read speed (Events/s);Measurement/10k Event", 1000,0,10000);
  h_read_speed_job             = new TH1D("read_speed_job",         ";Read speed (Events/s);Measurement/Job",       1000,0,10000);
  h_read_speed_vs_nevt_10k     = new TH2D("read_speed_vs_nevt_10k", ";Entry;Read speed (Events/s)/10k Event",       100,0,10000000, 200,0,10000);
  h_read_speed_vs_nevt_job     = new TH2D("read_speed_vs_nevt_job", ";Total Entries;Read speed (Events/s)/Job",     100,0,10000000, 200,0,10000);
  h_runtime_job                = new TH1D("runtime_job",            ";Total job run-time (min)",                    600,0,600);
  h_runtime_vs_nevt_job        = new TH2D("runtime_vs_nevt_job",    ";Total Entries;Total job run-time (min)",      100,0,10000000, 600,0,600);
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
AnalysisBase::get_xsec_from_txt_file(const std::string& txt_file, const std::string& dirname)
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
      if (dirname==shortname) evt_XSec = xsec;
    }
  }
  if (evt_XSec == 0) {
    std::cout<<"No crossection found for "<<dirname<<" in cross section file: "<<txt_file<<std::endl;
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
  for (auto filename : filenames) {
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
  for (auto filename : filenames) {
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
  // for (auto filename : filenames) {
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

std::vector<float> AK4Puppi_E, AK4Puppi_Pt;
std::vector<float> AK8Puppi_E, AK8Puppi_Pt, AK8Puppi_softDropMass;//, AK8Puppi_trimmedMass, AK8Puppi_prunedMass, AK8Puppi_filteredMass;

void
AnalysisBase::rescale_jets(DataStruct& data, const unsigned int& syst_index, const double& nSigmaJEC)
{
  if (syst_index==0) {
    AK4Puppi_E            = data.jetsAK4Puppi.E;
    AK4Puppi_Pt           = data.jetsAK4Puppi.Pt;
    AK8Puppi_E            = data.jetsAK8Puppi.E;
    AK8Puppi_Pt           = data.jetsAK8Puppi.Pt;
    AK8Puppi_softDropMass = data.jetsAK8Puppi.softDropMass;
    //AK8Puppi_trimmedMass  = data.jetsAK8Puppi.trimmedMass;
    //AK8Puppi_prunedMass   = data.jetsAK8Puppi.prunedMass;
    //AK8Puppi_filteredMass = data.jetsAK8Puppi.filteredMass;
  }
  if (nSigmaJEC != 0) {
    // AK4 Puppi jets
    AK4Puppi_Ht = 0;
    while(data.jetsAK4Puppi.Loop()) {
      size_t i = data.jetsAK4Puppi.it;
      double scale = get_syst_weight(1.0, data.jetsAK4Puppi.jecUncertainty[i], nSigmaJEC);
      data.jetsAK4Puppi.Pt[i] = AK4Puppi_Pt[i] * scale;
      data.jetsAK4Puppi.E[i]  = AK4Puppi_E[i]  * scale;
      AK4Puppi_Ht += data.jetsAK4Puppi.Pt[i];    
    }
    // AK8 Puppi jets
    AK8Puppi_Ht = 0;
    while(data.jetsAK8Puppi.Loop()) {
      size_t i = data.jetsAK4Puppi.it;
      double scale = get_syst_weight(1.0, data.jetsAK8Puppi.jecUncertainty[i], nSigmaJEC);
      data.jetsAK8Puppi.Pt[i]           = AK8Puppi_Pt[i]           * scale;
      data.jetsAK8Puppi.E[i]            = AK8Puppi_E[i]            * scale;
      data.jetsAK8Puppi.softDropMass[i] = AK8Puppi_softDropMass[i] * scale;
      //data.jetsAK8Puppi.trimmedMass[i]  = AK8Puppi_trimmedMass[i]  * scale;
      //data.jetsAK8Puppi.prunedMass[i]   = AK8Puppi_prunedMass[i]   * scale;
      //data.jetsAK8Puppi.filteredMass[i] = AK8Puppi_filteredMass[i] * scale;
      AK8Puppi_Ht += data.jetsAK8Puppi.Pt[i];
    }
  }
}


//____________________________________________________
//                  Scale factors

double
AnalysisBase::get_w_tagging_sf(DataStruct& data, const double& nSigmaWTagSF)
{
  double w = 1.0;
  return w;
}

double
AnalysisBase::get_b_tagging_sf(DataStruct& data, const double& nSigmaBTagSF)
{
  double w = 1.0;
  return w;
}

double
AnalysisBase::get_top_tagging_sf(DataStruct& data, const double& nSigmaHadTopTagSF)
{
  double w = 1.0;

  while(data.jetsAK8Puppi.Loop()) {
    size_t i = data.jetsAK8Puppi.it;
    double pt = data.jetsAK8Puppi.Pt[i];
    double sd_mass = data.jetsAK8Puppi.softDropMass[i];
    double tau32 = 9999;
    if (data.jetsAK8Puppi.tau2[i]>0) tau32 = data.jetsAK8Puppi.tau3[i]/data.jetsAK8Puppi.tau2[i];
    if (pt >= TOP_PT_CUT && sd_mass>=TOP_SD_MASS_CUT_LOW && sd_mass<TOP_SD_MASS_CUT_HIGH && tau32 < TOP_TAU32_CUT) {
#if USE_BTAG == 1
      if (passSubjetBTag[i]) {
#endif
	// Top-tagged AK8 jets
	if (data.jetsAK8Puppi.Pt[i] >= 400 && data.jetsAK8Puppi.Pt[i] < 550)
	  w *= get_syst_weight(TOP_TAG_SF_LOW, TOP_TAG_SF_LOW_ERR, nSigmaHadTopTagSF);
	else if (data.jetsAK8Puppi.Pt[i] >= 550)
	  w *= get_syst_weight(TOP_TAG_SF_HIGH, TOP_TAG_SF_HIGH_ERR, nSigmaHadTopTagSF);
#if USE_BTAG == 1
      }
#endif
    }
  }

  return w;
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
  double ht = 0; for (auto pt : AK8Puppi_Pt) ht += pt;

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
    for (auto cut : cuts_to_skip) std::cout<<cut<<", "; std::cout<<std::endl;
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
AnalysisBase::benchmarking(const int& entry, const int& nevents, const bool moderate=0)
{
  if (entry==0) {
    sw_1k_ ->Start(kFALSE);
    sw_10k_->Start(kFALSE);
    sw_job_->Start(kFALSE);
  } else {
    if (entry%1000==0) {
      double meas_1k = 1000/sw_1k_->RealTime();
      h_read_speed_1k->Fill(meas_1k);
      if (moderate) moderate_job_(h_read_speed_1k, meas_1k, 2, 5);
      sw_1k_->Reset();
      sw_1k_->Start(kFALSE);
      //std::cout<<"Meas  1k: "<<meas_1k<<std::endl;
    }
    if (entry%10000==0) {
      double meas_10k = 10000/sw_10k_->RealTime();
      h_read_speed_10k->Fill(meas_10k);
      h_read_speed_vs_nevt_10k->Fill(entry, meas_10k);
      if (moderate) moderate_job_(h_read_speed_10k, meas_10k, 2, 10);
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
    }
  }
}

void
AnalysisBase::moderate_job_(TH1D* h, double speed, int nrms_threshold, int sleep_s) {
  // Median - 2 sigma
  int n = h->GetXaxis()->GetNbins(); 
  std::vector<double> x(n);
  h->GetXaxis()->GetCenter(&x[0]);
  const double * y = h->GetArray();
  double median = TMath::Median(n, &x[0], &y[1]);
  //double mean = h->GetMean();
  //double threshold = mean -nrms_threshold*h->GetRMS();
  double threshold = median -nrms_threshold*h->GetRMS();
  if (speed<threshold) {
    sw_1k_->Stop(); sw_10k_->Stop();
    std::this_thread::sleep_for (std::chrono::seconds(sleep_s));
    sw_1k_->Continue(); sw_10k_->Continue();
  }
  //std::cout<<"- Moderating for "<<sleep_s<<" sec, median: "<<median<<" rms: "<<h->GetRMS()<<" threshold"<<threshold<<" speed: "<<speed<<std::endl;
}
