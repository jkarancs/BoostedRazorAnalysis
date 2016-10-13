#include <iostream>
#include <functional>
#include <map>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

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
  void define_preselections(const DataStruct&);

  void calculate_common_variables(DataStruct&, const unsigned int&);

  void init_common_histos();

  double get_xsec_from_ntuple(const std::vector<std::string>&, const std::string&, const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const std::string&);

  void calc_weightnorm_histo_from_ntuple(const std::vector<std::string>&, const double&, const std::vector<std::string>&,
					 const std::vector<std::string>&, bool);

  void init_pileup_reweightin(const std::string&, const std::string&, const std::vector<std::string>&);

  double get_pileup_weight(const int&, const double&);

  void rescale_jets(DataStruct&, const unsigned int&, const double&);

  double get_top_tagging_sf(DataStruct&, const double&);

  double get_ht_weight(DataStruct&, const double&);

  double get_alphas_weight(const std::vector<float>&, const double&, const int&);

  double get_scale_weight(const std::vector<float>&, const double&, const unsigned int&);

  double get_syst_weight(const double&, const double&, const double&, const double&);

  double get_syst_weight(const double&, const double&, const double&);

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

  std::vector<Cut> analysis_cuts;

private:
  bool _apply_cut(std::string);
  bool _apply_ncut(size_t);

  typedef struct Sample { std::string postfix; std::string legend; std::string color; std::vector<std::string> dirs; } Sample;
  typedef struct PostfixOptions { size_t index; std::string postfixes; std::string legends; std::string colors; } PostfixOptions;
  PostfixOptions get_pf_opts_(std::vector<std::vector<Sample> > lists, std::string);
};

//_______________________________________________________
//                       Constructor
AnalysisBase::AnalysisBase() { }


//_______________________________________________________
//                       Destructor
AnalysisBase::~AnalysisBase() { }


//_______________________________________________________
//                 Define baseline cuts
void
AnalysisBase::define_preselections(const DataStruct& data)
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
  
  // Recommended event filters by MET group - Updated for 76X!
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_76X_v2_produced_with_the
  // 
  // Select at least one good vertex (z<24, rho<2, ndof>=4)
  // NGoodVtx defined in:
  // https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L272-L275
  // baseline_cuts.push_back({ .name="met_filter_NGoodVtx",          .func = [&data](){ return data.evt.NGoodVtx>0; } });
  baseline_cuts.push_back({ .name="met_filter_NGoodVtx",          .func = [&data](){ return data.filter.goodVertices; } }); // Now works in 76X MiniAOD
  
  // Other filters (From MiniAOD)
  baseline_cuts.push_back({ .name="met_filter_EE_Bad_Sc",         .func = [&data](){ return data.filter.eeBadScFilter; } });
  baseline_cuts.push_back({ .name="met_filter_Ecal_Dead_Cell_TP", .func = [&data](){ return data.filter.EcalDeadCellTriggerPrimitiveFilter; } });
  baseline_cuts.push_back({ .name="met_filter_HBHE_Noise",        .func = [&data](){ return data.filter.HBHENoiseFilter; } });
  baseline_cuts.push_back({ .name="met_filter_HBHE_IsoNoise",     .func = [&data](){ return data.filter.HBHENoiseIsoFilter; } });
  baseline_cuts.push_back({ .name="met_filter_CSC_Halo_Tight",    .func = [&data](){ return data.filter.CSCTightHalo2015Filter; } });
  //baseline_cuts.push_back({ .name="met_filter_Muon_Bad_Track",    .func = [&data](){ return data.filter.muonBadTrackFilter; } });
  //baseline_cuts.push_back({ .name="met_filter_CH_Track_Resol",    .func = [&data](){ return data.filter.chargedHadronTrackResolutionFilter; } });
}

//_______________________________________________________
//                 Define common variables

/*
  Jet ID:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

  Choose:
  - Loose jet ID
*/

/*
   Top Tagging working points:
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

/*
  Latest TOP Tagging working point (2015 Data)

   Latest WPs/SFs -  JME-16-003 PAS:
   http://cms.cern.ch/iCMS/analysisadmin/viewanalysis?id=1694&field=id&value=1694
   
   Choose:
   - Loose selection: e(S) ~= 55%, e(B) = 3.0% WP
   - AK8 Puppi jets
   - 105 <= SD Mass < 210
   - tau32 < 0.8
   - max subjet BTag CSV > 0.46
*/

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
   W tagging working points:
   https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging

   Latest WPs/SFs not yet on twiki (Thea Aarrestad, slide 77):
   https://indico.cern.ch/event/530683/contributions/2164853/attachments/1271780/1884879/WtagSF_JMAR_TAarrestad.pdf

   Choose:
   - Very loose selection e(S) = 93.3%
   - AK8 Puppi jets
   - 65 < SD Mass < 105
   - tau21 <= 0.56
*/

#define W_PT_CUT            200
#define W_SD_MASS_CUT_LOW    65
#define W_SD_MASS_CUT_HIGH  105
#define W_TAU21_CUT        0.56

// Further analysis cuts
#define R_CUT             0.4
#define R_CUT_LOW         0.2
#define DPHI_CUT          2.7

// AK8 jets
std::vector<int> passSubjetBTag;
std::vector<int> passHadTopTag;
std::vector<int> passHadTopPreTag;
std::vector<int> passHadWTag;
std::vector<double> subjetBTagDiscr;

// Event - Jets
unsigned int nLooseJet;
unsigned int nSubjetBTag;
unsigned int nHadTopTag;
unsigned int nHadTopPreTag;
unsigned int nHadWTag;
double AK4Puppi_Ht, AK8Puppi_Ht;

/*
  Electron Veto ID:
  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2

  ID used in AugXX/SepXX ntuple production:
  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns

  Latest ID (upcoming production):
  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for

  Choose:
  - Cut based Veto ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1
  - pt >= 5
  - |eta| < 2.5

-------------------------------
  
  Muon Loose ID:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  
  Latest ID used:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Loose_Muon
  
  Choose:
  - POG recommended Loose ID
  - Mini-Isolation (EA)/pt < 0.2
  - pt >= 5
  - |eta| < 2.4

  Not (yet) used - variable needs to be added for next production:
  - Impact point: |d0| < 0.2, |dz| < 0.5

*/

#define ELE_PT_CUT 5
#define ELE_ABSETA_CUT 2.5
#define ELE_MINIISO_CUT 0.1
#define MU_PT_CUT 5
#define MU_ABSETA_CUT 2.4
#define MU_MINIISO_CUT 0.2
//#define MU_ABSD0_CUT 0.2
//#define MU_ABSDZ_CUT 0.5

// Event Letpons
unsigned int nEleVeto;
unsigned int nMuLoose;

// Min(DeltaPhi(Jet_i, MET)), i=1,2,3
double minDeltaPhi;

void
AnalysisBase::calculate_common_variables(DataStruct& data, const unsigned int& syst_index)
{
  // It only makes sense to calculate certain variables only once if they don't depend on jet energy
  if (syst_index == 0) {
    nLooseJet = 0;
    nSubjetBTag = 0;
    passSubjetBTag.clear();
    subjetBTagDiscr.clear();

    // Loop on AK8 Puppi jets
    while(data.jetsAK8Puppi.Loop()) {
      if (data.jetsAK8Puppi.looseJetID[data.jetsAK8Puppi.it]) nLooseJet++;

      // Maximum Subjet btag discriminator
      // Puppi jets doesn't store them correctly, get them from CHS jets for now
      double max_subjet_btag_discr = -9999;
      int subjet0_index = data.jetsAK8Puppi.vSubjetIndex0[data.jetsAK8Puppi.it], subjet1_index = data.jetsAK8Puppi.vSubjetIndex1[data.jetsAK8Puppi.it];
      if (subjet0_index != -1) if (data.subjetsAK8Puppi.CSVv2[subjet0_index] > max_subjet_btag_discr) max_subjet_btag_discr = data.subjetsAK8Puppi.CSVv2[subjet0_index];
      if (subjet1_index != -1) if (data.subjetsAK8Puppi.CSVv2[subjet1_index] > max_subjet_btag_discr) max_subjet_btag_discr = data.subjetsAK8Puppi.CSVv2[subjet1_index];
#if USE_BTAG == 1
      passSubjetBTag.push_back((max_subjet_btag_discr > TOP_BTAG_CSV));
      if (max_subjet_btag_discr > TOP_BTAG_CSV) ++nSubjetBTag;
#else
      passSubjetBTag.push_back((max_subjet_btag_discr> 0.46));
      if (max_subjet_btag_discr > 0.46) ++nSubjetBTag;
#endif
      subjetBTagDiscr.push_back(max_subjet_btag_discr);
    }

    // AK4 Puppi jets
    minDeltaPhi = 9999;
    while(data.jetsAK4Puppi.Loop()) {
      // minDeltaPhi
      if (data.jetsAK4Puppi.it<3) {
	double dphi = std::abs(TVector2::Phi_mpi_pi(data.met.Phi[0] - data.jetsAK4Puppi.Phi[data.jetsAK4Puppi.it]));
	if (dphi<minDeltaPhi) minDeltaPhi = dphi;
      }
    }

    // Number of Veto Electrons
    nEleVeto = 0;
    while(data.ele.Loop()) {
      if (data.ele.Pt[data.ele.it] < ELE_PT_CUT) continue;
      if (fabs(data.ele.Eta[data.ele.it]) >= ELE_ABSETA_CUT) continue;
      if (!data.ele.IDVeto_NoIso[data.ele.it]) continue;
      if (data.ele.MiniIso[data.ele.it]/data.ele.Pt[data.ele.it] >= ELE_MINIISO_CUT) continue;
      ++nEleVeto;
    }

    // Number of Loose Muons
    nMuLoose = 0;
    while(data.mu.Loop()) {
      if (data.mu.Pt[data.mu.it] < MU_PT_CUT) continue;
      if (fabs(data.mu.Eta[data.mu.it]) >= MU_ABSETA_CUT) continue;
      if (!data.mu.IsLooseMuon[data.mu.it]) continue;
      if (data.mu.MiniIso[data.mu.it]/data.mu.Pt[data.mu.it] >= MU_MINIISO_CUT) continue;
      //if (fabs(data.mu.Dxy[data.mu.it]) >= MU_ABSD0_CUT) continue;
      //if (fabs(data.mu.Dz[data.mu.it]) >= MU_ABSDZ_CUT) continue;
      ++nMuLoose;
    }
  }

  // Rest of the vairables need to be recalculated each time the jet energy is changed
  // eg. top/W tags and HT (obviously) depends on jet pt
  AK4Puppi_Ht = 0;
  while(data.jetsAK4Puppi.Loop()) AK4Puppi_Ht += data.jetsAK4Puppi.Pt[data.jetsAK4Puppi.it];    

  AK8Puppi_Ht = 0;
  nHadTopTag = nHadTopPreTag = nHadWTag = 0;
  passHadTopTag.clear();
  passHadTopPreTag.clear();
  passHadWTag.clear();

  // Loop on AK8 Puppi jets
  while(data.jetsAK8Puppi.Loop()) {
    AK8Puppi_Ht += data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it];

    // _______________________________________________________
    //                  Boosted Objects

    double pt = data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it];
    double sd_mass = data.jetsAK8Puppi.softDropMass[data.jetsAK8Puppi.it];
    double tau21 = data.jetsAK8Puppi.tau2[data.jetsAK8Puppi.it];
    if (data.jetsAK8Puppi.tau1[data.jetsAK8Puppi.it]!=0) tau21 /= data.jetsAK8Puppi.tau1[data.jetsAK8Puppi.it];
    else tau21 = 9999;
    double tau32 = data.jetsAK8Puppi.tau3[data.jetsAK8Puppi.it];
    if (data.jetsAK8Puppi.tau2[data.jetsAK8Puppi.it]!=0) tau32 /= data.jetsAK8Puppi.tau2[data.jetsAK8Puppi.it];
    else tau32 = 9999;
    
    // _______________________________________________________
    //                  Hadronic Top Tag definition
    
    passHadTopTag.push_back(0);
    passHadTopPreTag.push_back(0);
    // New hadronic top tag
    if (pt >= TOP_PT_CUT && sd_mass>=TOP_SD_MASS_CUT_LOW && sd_mass<TOP_SD_MASS_CUT_HIGH) {
      passHadTopPreTag[data.jetsAK8Puppi.it] = 1;
      ++nHadTopPreTag;
      if (tau32 < TOP_TAU32_CUT) {
#if USE_BTAG == 1
	if (passSubjetBTag[data.jetsAK8Puppi.it]) {
#endif
	  passHadTopTag[data.jetsAK8Puppi.it] = 1;
	  ++nHadTopTag;
#if USE_BTAG == 1
	}
#endif
      }
    }

    // _______________________________________________________
    //                  Hadronic W Tag definition
    passHadWTag.push_back(0);
    // New hadronic top tag
    if (pt >= W_PT_CUT && sd_mass>=W_SD_MASS_CUT_LOW && sd_mass<W_SD_MASS_CUT_HIGH) {
      if (tau21 < W_TAU21_CUT) {
	passHadWTag[data.jetsAK8Puppi.it] = 1;
	++nHadWTag;
      }
    }
  } // end of AK8 (Puppi) jet loop
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
}


//_______________________________________________________
//           Read cross-section from ntuple
double
AnalysisBase::get_xsec_from_ntuple(const std::vector<std::string>& filenames, const std::string& treename, const std::string& dirname)
{
  /*
    Aug17 ntuple production did not correctly fill cross-sections

  float evt_XSec=0, prev_XSec=0;
  for (auto filename : filenames) {
    TFile *f = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)f->Get(treename.c_str());
    tree->GetBranch("evt_XSec")->SetAddress(&evt_XSec);
    tree->GetEntry(0);
    f->Close();
    if (prev_XSec!=0&&prev_XSec!=evt_XSec) {
      cout << "!! Error !! Analysis - Files added with different cross-sections. Please, add them separately!" << endl;
      return 0;
    }
    prev_XSec = evt_XSec;
  }
  */

  // Temporarily read cross-sections from txt file
  double evt_XSec = 0;
  std::ifstream xsecFile("common/BackGroundXSec.txt");
  if ( !xsecFile.good() ) {
    return -9999.0;
    std::cout<<"Unable to open cross-section file: common/BackGroundXSec.txt"<<std::endl;
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
      if (dirname==shortname) {
	evt_XSec = xsec;
	std::cout<<"( TEMPORARY FIX: Cross section successfully read from common/BackGroundXSec.txt )"<<std::endl;
      }
    }
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
      double scale = get_syst_weight(1.0, data.jetsAK4Puppi.jecUncertainty[data.jetsAK4Puppi.it], nSigmaJEC);
      data.jetsAK4Puppi.Pt[data.jetsAK4Puppi.it] = AK4Puppi_Pt[data.jetsAK4Puppi.it] * scale;
      data.jetsAK4Puppi.E[data.jetsAK4Puppi.it]  = AK4Puppi_E[data.jetsAK4Puppi.it]  * scale;
      AK4Puppi_Ht += data.jetsAK4Puppi.Pt[data.jetsAK4Puppi.it];    
    }
    // AK8 Puppi jets
    AK8Puppi_Ht = 0;
    while(data.jetsAK8Puppi.Loop()) {
      double scale = get_syst_weight(1.0, data.jetsAK8Puppi.jecUncertainty[data.jetsAK8Puppi.it], nSigmaJEC);
      data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]           = AK8Puppi_Pt[data.jetsAK8Puppi.it]           * scale;
      data.jetsAK8Puppi.E[data.jetsAK8Puppi.it]            = AK8Puppi_E[data.jetsAK8Puppi.it]            * scale;
      data.jetsAK8Puppi.softDropMass[data.jetsAK8Puppi.it] = AK8Puppi_softDropMass[data.jetsAK8Puppi.it] * scale;
      //data.jetsAK8Puppi.trimmedMass[data.jetsAK8Puppi.it]  = AK8Puppi_trimmedMass[data.jetsAK8Puppi.it]  * scale;
      //data.jetsAK8Puppi.prunedMass[data.jetsAK8Puppi.it]   = AK8Puppi_prunedMass[data.jetsAK8Puppi.it]   * scale;
      //data.jetsAK8Puppi.filteredMass[data.jetsAK8Puppi.it] = AK8Puppi_filteredMass[data.jetsAK8Puppi.it] * scale;
      AK8Puppi_Ht += data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it];
    }
  }
}


//____________________________________________________
//               Top-Tagging Scale factor
double
AnalysisBase::get_top_tagging_sf(DataStruct& data, const double& nSigmaHadTopTagSF)
{
  double w = 1.0;

  while(data.jetsAK8Puppi.Loop()) {
    double pt = data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it];
    double sd_mass = data.jetsAK8Puppi.softDropMass[data.jetsAK8Puppi.it];
    double tau32 = data.jetsAK8Puppi.tau3[data.jetsAK8Puppi.it];
    if (data.jetsAK8Puppi.tau2[data.jetsAK8Puppi.it]!=0) tau32 /= data.jetsAK8Puppi.tau2[data.jetsAK8Puppi.it];
    else tau32 = 9999;
    if (pt >= TOP_PT_CUT && sd_mass>=TOP_SD_MASS_CUT_LOW && sd_mass<TOP_SD_MASS_CUT_HIGH && tau32 < TOP_TAU32_CUT) {
#if USE_BTAG == 1
      if (passSubjetBTag[data.jetsAK8Puppi.it]) {
#endif
	// Top-tagged AK8 jets
	if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it] >= 400 && data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it] < 550)
	  w *= get_syst_weight(TOP_TAG_SF_LOW, TOP_TAG_SF_LOW_ERR, nSigmaHadTopTagSF);
	else if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it] >= 550)
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
