#include <iostream>
#include <functional>
#include <map>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

#include "DataStruct.h"
#include "GluinoXSec.h"

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

  double get_xsec_from_ntuple(const std::vector<std::string>&, const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const std::string&);

  void calc_weightnorm_histo_from_ntuple(const std::vector<std::string>&, const double&, const std::vector<std::string>&,
					 const std::vector<std::string>&, const std::vector<std::string>&, bool);

  void init_pileup_reweightin(const std::string&, const std::string&, const std::vector<std::string>&);

  double get_pileup_weight(const int&, const double&);

  void rescale_jets(DataStruct&, const unsigned int&, const double&);

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

   Latest WPs not yet on twiki (Gregor Kasieczka, Mareike Meyer):
   https://indico.cern.ch/event/518509/contributions/2032850/attachments/1256008/1854115/TopTagging11_04.pdf

   Scale factors (?):
   

   Choose:
   - Loose selection: e(B) = 10% WP
   - AK8 Puppi jets
   - 60 < SD Mass < 190
   - tau32 < 0.76
*/

#define TOP_PT_CUT            400
#define TOP_SD_MASS_CUT_LOW    60 // prev 110
#define TOP_SD_MASS_CUT_HIGH  190 // prev 210
#define TOP_TAU32_CUT        0.76 // prev 0.75


/* 
   W tagging working points:
   https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging

   Latest WPs not yet on twiki (Thea Aarrestad, slide 77):
   https://indico.cern.ch/event/530683/contributions/2164853/attachments/1271780/1884879/WtagSF_JMAR_TAarrestad.pdf

   Scale factors (Thea Aarrestad, slide 77):
   https://indico.cern.ch/event/530683/contributions/2164853/attachments/1271780/1884879/WtagSF_JMAR_TAarrestad.pdf

   Choose:
   - Loose selection e(S) = 93.3%
   - AK8 Puppi jets
   - 65 < SD Mass < 105
   - tau21 <= 0.56
*/

#define W_PT_CUT            200
#define W_SD_MASS_CUT_LOW    65
#define W_SD_MASS_CUT_HIGH  105
#define W_TAU21_CUT        0.56

// AK8 jets
std::vector<int> passLooseJetID;
std::vector<int> passHadTopTag;
std::vector<int> passHadTopPreTag;
std::vector<int> passHadWTag;

// Event
unsigned int nLooseJet;
unsigned int nHadTopTag;
unsigned int nHadTopPreTag;
unsigned int nHadWTag;
double AK4_Ht, AK4Puppi_Ht, AK8Puppi_Ht;

void
AnalysisBase::calculate_common_variables(DataStruct& data, const unsigned int& syst_index)
{
  // It only makes sense to calculate certain variables only once if they don't depend on jet energy
  if (syst_index == 0) {
    nLooseJet = 0;

    // Jet ID - Loose (Fractions should not depend on JEC)
    // Loop on AK8 Puppi jets
    while(data.jetsAK8Puppi.Loop()) {
      double eta = data.jetsAK8Puppi.Eta[data.jetsAK8Puppi.it];
      double nhf = data.jetsAK8Puppi.neutralHadronEnergy[data.jetsAK8Puppi.it] / data.jetsAK8Puppi.E[data.jetsAK8Puppi.it];
      double nemf = data.jetsAK8Puppi.neutralEmEnergy[data.jetsAK8Puppi.it] / data.jetsAK8Puppi.E[data.jetsAK8Puppi.it];
      double chf = data.jetsAK8Puppi.chargedHadronEnergy[data.jetsAK8Puppi.it]/data.jetsAK8Puppi.E[data.jetsAK8Puppi.it];
      double cemf = data.jetsAK8Puppi.chargedEmEnergy[data.jetsAK8Puppi.it]/data.jetsAK8Puppi.E[data.jetsAK8Puppi.it];
      int NumConst = data.jetsAK8Puppi.chargedMultiplicity[data.jetsAK8Puppi.it] + data.jetsAK8Puppi.neutralMultiplicity[data.jetsAK8Puppi.it];
      int chm = data.jetsAK8Puppi.chargedMultiplicity[data.jetsAK8Puppi.it];
      int NumNeutralParticle = data.jetsAK8Puppi.neutralMultiplicity[data.jetsAK8Puppi.it];
      bool pass_Loose_ID = ( (nhf<0.99 && nemf<0.99 && NumConst>1) && ((fabs(eta)<=2.4 && chf>0 && chm>0 && cemf<0.99) || fabs(eta)>2.4) && fabs(eta)<=3.0 )
        || ( nemf<0.9 && NumNeutralParticle>10 && fabs(eta)>3.0 );
      passLooseJetID.push_back(pass_Loose_ID);
      if (pass_Loose_ID) nLooseJet++;
    }
  }

  // Rest of the vairables need to be recalculated each time the jet energy is changed
  // eg. top/W tags and HT (obviously) depends on jet pt
  AK4_Ht = 0;
  while(data.jetsAK4.Loop()) AK4_Ht += data.jetsAK4.Pt[data.jetsAK4.it];    

  AK4Puppi_Ht = 0;
  while(data.jetsAK4Puppi.Loop()) AK4Puppi_Ht += data.jetsAK4Puppi.Pt[data.jetsAK4Puppi.it];    

  AK8Puppi_Ht = 0;
  nHadTopTag = nHadTopPreTag = nHadWTag = 0;

  // Loop on AK8 Puppi jets
  while(data.jetsAK8Puppi.Loop()) {
    AK8Puppi_Ht += data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it];

    // _______________________________________________________
    //                  Boosted Objects

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
    if (sd_mass>=TOP_SD_MASS_CUT_LOW && sd_mass<TOP_SD_MASS_CUT_HIGH) {
      passHadTopPreTag[data.jetsAK8Puppi.it] = 1;
      ++nHadTopPreTag;
      if (tau32 < TOP_TAU32_CUT) {
	passHadTopTag[data.jetsAK8Puppi.it] = 1;
	++nHadTopTag;
      }
    }
    
    // _______________________________________________________
    //                  Hadronic W Tag definition
    passHadWTag.push_back(0);
    // New hadronic top tag
    if (sd_mass>=W_SD_MASS_CUT_LOW && sd_mass<W_SD_MASS_CUT_HIGH) {
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
  // signal weight and xsec
  vh_totweight_signal .push_back(new TH2D("totweight_T1tttt",   "T1tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Total Weight",        201,-12.5,5012.5, 201,-12.5,5012.5));
  vh_xsec_signal      .push_back(new TH2D("xsec_T1tttt",        "T1tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);Cross-section (pb)",  201,-12.5,5012.5, 201,-12.5,5012.5));
  vh_weightnorm_signal.push_back(new TH2D("weightnorm_T1tttt",  "T1tttt;M_{#tilde{g}} (GeV);M_{#tilde{#chi}^{0}} (GeV);weight norm. factor", 201,-12.5,5012.5, 201,-12.5,5012.5));
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
AnalysisBase::get_xsec_from_ntuple(const std::vector<std::string>& filenames, const std::string& treename)
{
  float evt_XSec=0, prev_XSec=0;
  for (auto filename : filenames) {
    TTree* tree = (TTree*)TFile::Open(filename.c_str())->Get(treename.c_str());
    tree->GetBranch("evt_XSec")->SetAddress(&evt_XSec);
    tree->GetEntry(0);
    if (prev_XSec!=0&&prev_XSec!=evt_XSec) {
      cout << "!! Error !! Analysis - Files added with different cross-sections. Please, add them separately!" << endl;
      return 0;
    }
    prev_XSec = evt_XSec;
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
						const std::vector<std::string>& vname_xsec, const std::vector<std::string>& vname_totweight, bool verbose=1)
{
  // Get gluino cross sections and merge totweight histos
  std::map<int, double> xsec_glu;
  for (auto filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    for (size_t i=0, n=vname_xsec.size(); i<n; ++i) {
      // Get cross section (set only once for each new gluino weight)
      /*          
		  The filling of xsec histo has some bug in B2GTTree level
		  
      TH2D* xsec = (TH2D*)f->Get(vname_xsec[i].c_str());
      for (int binx=1, nbinx=xsec->GetNbinsX(); binx<=nbinx; ++binx) 
	for (int biny=1, nbiny=xsec->GetNbinsY(); biny<=nbiny; ++biny) {
	  double xs = xsec->GetBinContent(binx, biny);
	  if (xs!=0&&!xsec_glu.count(binx)) xsec_glu[binx] = xs;
        }
      */
      // Get total weight
      TH2D* totweight = (TH2D*)f->Get(vname_totweight[i].c_str());
      vh_totweight_signal[i]->Add(totweight);
    }
    f->Close();
  }
  // Set xsec for each gluino mass bin
  // This way we fix temporary bug with xsec for some bins
  // Read gluino xsec from same file used in TTree step
  for (size_t i=0, n=vname_xsec.size(); i<n; ++i) 
    for (int binx=1, nbinx=vh_xsec_signal[i]->GetNbinsX(); binx<=nbinx; ++binx) {
      double mGlu = vh_xsec_signal[i]->GetBinCenter(binx);
      xsec_glu[binx] = GetGluinoXSec(mGlu).first; // first: mean xsec, second: error
      for (int biny=1, nbiny=vh_xsec_signal[i]->GetNbinsY(); biny<=nbiny; ++biny)
	vh_xsec_signal[i]->SetBinContent(binx, biny, xsec_glu[binx]);
    }
  // Calculate weight normalization
  for (size_t i=0, n=vname_xsec.size(); i<n; ++i) {
    // weightnorm = (settings.intLumi*xsec)/totweight;
    // Divide(h1,h2,c1,c2) --> c1*h1/(c2*h2)
    vh_weightnorm_signal[i]->Divide(vh_xsec_signal[i], vh_totweight_signal[i], intLumi);
    if (verbose) {
      std::cout<<"- Signal: "<<vname_signal[i]<<std::endl;
      for (int binx=1, nbinx=vh_xsec_signal[i]->GetNbinsX(); binx<=nbinx; ++binx) 
	for (int biny=1, nbiny=vh_xsec_signal[i]->GetNbinsY(); biny<=nbiny; ++biny) {
	  double mGlu = vh_xsec_signal[i]->GetXaxis()->GetBinCenter(binx);
	  double mLSP = vh_xsec_signal[i]->GetYaxis()->GetBinCenter(biny);
	  double xsec  = vh_xsec_signal[i]      ->GetBinContent(binx, biny);
	  double totw  = vh_totweight_signal[i] ->GetBinContent(binx, biny);
	  double wnorm = vh_weightnorm_signal[i]->GetBinContent(binx, biny);
	  if (totw>0) std::cout<<"  Bin: M(g~)="<<mGlu<<" M(LSP)="<<mLSP<<":   xsec="<<xsec<<" totweight="<<totw<<" weightnorm="<<wnorm<<std::endl;
	}
      std::cout<<std::endl;
    }
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
  if (nSigma == 0) {
    return w;
  } else {
    // Compute the weight according to the systematic variation considered
    // Use difference between nominal and up/down as 1 sigma variation 
    double dw_up = weight_nominal * uncertainty;
    double dw_down = -1.0 * dw_up;
    if (nSigma >= 0.) {
      w += nSigma*dw_up; 
    } else {
      w += nSigma*dw_down;
    }
    return w; 
  }
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
std::vector<float> AK4Puppi_E, AK4Puppi_Pt;
std::vector<float> AK8_E, AK8_Pt, AK8_softDropMass, AK8_trimmedMass, AK8_prunedMass, AK8_filteredMass;
std::vector<float> AK8Puppi_E, AK8Puppi_Pt, AK8Puppi_softDropMass, AK8Puppi_trimmedMass, AK8Puppi_prunedMass, AK8Puppi_filteredMass;

void
AnalysisBase::rescale_jets(DataStruct& data, const unsigned int& syst_index, const double& nSigmaJEC)
{
  if (syst_index==0) {
    AK4_E            = data.jetsAK4.E;
    AK4_Pt           = data.jetsAK4.Pt;
    AK8_E            = data.jetsAK8.E;
    AK8_Pt           = data.jetsAK8.Pt;
    AK8_softDropMass = data.jetsAK8.softDropMass;
    AK8_trimmedMass  = data.jetsAK8.trimmedMass;
    AK8_prunedMass   = data.jetsAK8.prunedMass;
    AK8_filteredMass = data.jetsAK8.filteredMass;
    AK4Puppi_E            = data.jetsAK4Puppi.E;
    AK4Puppi_Pt           = data.jetsAK4Puppi.Pt;
    AK8Puppi_E            = data.jetsAK8Puppi.E;
    AK8Puppi_Pt           = data.jetsAK8Puppi.Pt;
    AK8Puppi_softDropMass = data.jetsAK8Puppi.softDropMass;
    AK8Puppi_trimmedMass  = data.jetsAK8Puppi.trimmedMass;
    AK8Puppi_prunedMass   = data.jetsAK8Puppi.prunedMass;
    AK8Puppi_filteredMass = data.jetsAK8Puppi.filteredMass;
  }
  if (nSigmaJEC != 0) {
    // AK4 CHS jets
    while(data.jetsAK4.Loop()) {
      double scale = get_syst_weight(1.0, data.jetsAK4.jecUncertainty[data.jetsAK4.it], nSigmaJEC);
      data.jetsAK4.Pt[data.jetsAK4.it] = AK4_Pt[data.jetsAK4.it] * scale;
      data.jetsAK4.E[data.jetsAK4.it]  = AK4_E[data.jetsAK4.it]  * scale;
    }
    // AK4 Puppi jets
    while(data.jetsAK4Puppi.Loop()) {
      double scale = get_syst_weight(1.0, data.jetsAK4Puppi.jecUncertainty[data.jetsAK4Puppi.it], nSigmaJEC);
      data.jetsAK4Puppi.Pt[data.jetsAK4Puppi.it] = AK4Puppi_Pt[data.jetsAK4Puppi.it] * scale;
      data.jetsAK4Puppi.E[data.jetsAK4Puppi.it]  = AK4Puppi_E[data.jetsAK4Puppi.it]  * scale;
    }
    // AK8 CHS jets
    while(data.jetsAK8.Loop()) {
      double scale = get_syst_weight(1.0, data.jetsAK8.jecUncertainty[data.jetsAK8.it], nSigmaJEC);
      data.jetsAK8.Pt[data.jetsAK8.it]           = AK8_Pt[data.jetsAK8.it]           * scale;
      data.jetsAK8.E[data.jetsAK8.it]            = AK8_E[data.jetsAK8.it]            * scale;
      data.jetsAK8.softDropMass[data.jetsAK8.it] = AK8_softDropMass[data.jetsAK8.it] * scale;
      data.jetsAK8.trimmedMass[data.jetsAK8.it]  = AK8_trimmedMass[data.jetsAK8.it]  * scale;
      data.jetsAK8.prunedMass[data.jetsAK8.it]   = AK8_prunedMass[data.jetsAK8.it]   * scale;
      data.jetsAK8.filteredMass[data.jetsAK8.it] = AK8_filteredMass[data.jetsAK8.it] * scale;
    }
    // AK8 Puppi jets
    while(data.jetsAK8Puppi.Loop()) {
      double scale = get_syst_weight(1.0, data.jetsAK8Puppi.jecUncertainty[data.jetsAK8Puppi.it], nSigmaJEC);
      data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]           = AK8Puppi_Pt[data.jetsAK8Puppi.it]           * scale;
      data.jetsAK8Puppi.E[data.jetsAK8Puppi.it]            = AK8Puppi_E[data.jetsAK8Puppi.it]            * scale;
      data.jetsAK8Puppi.softDropMass[data.jetsAK8Puppi.it] = AK8Puppi_softDropMass[data.jetsAK8Puppi.it] * scale;
      data.jetsAK8Puppi.trimmedMass[data.jetsAK8Puppi.it]  = AK8Puppi_trimmedMass[data.jetsAK8Puppi.it]  * scale;
      data.jetsAK8Puppi.prunedMass[data.jetsAK8Puppi.it]   = AK8Puppi_prunedMass[data.jetsAK8Puppi.it]   * scale;
      data.jetsAK8Puppi.filteredMass[data.jetsAK8Puppi.it] = AK8Puppi_filteredMass[data.jetsAK8Puppi.it] * scale;
    }
  }
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

    https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L195-L202
    We save only ids: 2,3,4,5,7,9 (in this order)

    The idea here is to randomly choose to vary mu_f or mu_r or both simulataneously
    and rescale weight difference the usual way by desired nSigma
  */
  double w_scale = 1;
  double w_scale_up = 1;
  double w_scale_down = 1;
  if (numScale==1) {
    // fix mu_r = 1.0, vary mu_f = 2.0, 0.5
    w_scale_up   = scale_Weights[0];
    w_scale_down = scale_Weights[1];
  } else if (numScale==2) {
    // fix mu_f = 1.0, vary mu_r = 2.0, 0.5
    w_scale_up   = scale_Weights[2];
    w_scale_down = scale_Weights[4];
  } else if (numScale==3) {
    // vary simulataneously mu_r = mu_f = 2.0, 0.5
    w_scale_up   = scale_Weights[3];
    w_scale_down = scale_Weights[5];
  }
  w_scale = get_syst_weight(w_scale, w_scale_up, w_scale_down, nSigmaScale);
  return w_scale;
}
