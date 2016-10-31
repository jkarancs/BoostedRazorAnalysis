#include "TLorentzVector.h"
#include "common/AnalysisBase.h"

//_______________________________________________________
//                       Constructor
Analysis::Analysis() : AnalysisBase() { }


//_______________________________________________________
//                       Destructor
Analysis::~Analysis() { }

//_______________________________________________________
//                  Calculate variables

void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{
}

//_______________________________________________________
//          Define Analysis specific weights
double
Analysis::get_analysis_weight(DataStruct& data)
{
  double w = 1;
    
  return w;
}

//_______________________________________________________
//          Define Analysis event selection cuts

bool
Analysis::pass_skimming(DataStruct& data)
{
  if (!(nJet>=2)) return 0;
  if (!(nJetAK8>=1)) return 0;
  return 1;

  // Signal skim
  //return _apply_ncut(analysis_cuts.size());
}

void
Analysis::define_selections(const DataStruct& d)
{
  analysis_cuts.clear();

  // MET Filters, etc. are already applied in AnalysisBase.h, See baseline_cuts
  analysis_cuts.push_back({ .name="2Jet",     .func = []    { return nJet>=2;                       }});
  analysis_cuts.push_back({ .name="3Jet",     .func = []    { return nJet>=3;                       }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts.push_back({ .name="1JetAK8",  .func = []    { return nJetAK8>=1;                    }}); // Similar to pt>200, one AK8 jet has pt>170
  analysis_cuts.push_back({ .name="MR_R2",    .func = [&d]  { return d.evt.MR>800 && d.evt.R2>0.08; }});
  analysis_cuts.push_back({ .name="0Ele",     .func = []    { return nEleVeto==0;                   }});
  analysis_cuts.push_back({ .name="0Mu",      .func = []    { return nMuVeto==0;                    }});
  //analysis_cuts.push_back({ .name="0TauTrk",  .func = []   { return;  }});
  analysis_cuts.push_back({ .name="1b",       .func = []    { return nMediumBTag>=1;                }});
  analysis_cuts.push_back({ .name="1W",       .func = []    { return nTightWTag>=1;                 }});
  //analysis_cuts.push_back({ .name="mDPhiHat", .func = []   { return;  }});
  analysis_cuts.push_back({ .name="mDPhi",    .func = []    { return minDeltaPhi>=0.4;              }}); // Decreased it to the AK4 cone size (from 0.5)
}

bool
Analysis::_apply_cut(std::string cut_name) {
  for (auto cut : analysis_cuts) if (cut_name == cut.name) return cut.func();
  return 0;
}

bool
Analysis::_apply_ncut(size_t ncut) {
  if (ncut>analysis_cuts.size()) return 0;
  for (size_t i=0; i<ncut; ++i) if ( ! analysis_cuts[i].func() ) return 0;
  return 1;
}

//_______________________________________________________
//                 Signal Region
// Must define at some point, because we need to blind this region in data

bool
Analysis::signal_selection(const DataStruct& data) {
  // This will blind the data in the Signal selection
  return _apply_ncut(analysis_cuts.size());
}

//_______________________________________________________
//                 List of Histograms
TH1D* h_njet;
TH1D* h_nb;
TH1D* h_nw;
TH1D* h_ht_gen;
TH1D* h_ht_AK4Puppi;
TH1D* h_ht_AK8Puppi;
TH1D* h_jet1_pt;
TH1D* h_jet2_pt;
TH1D* h_jet3_pt;
std::vector<TH1D*> vh_jet1_pt;

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  h_njet         = new TH1D("njet",         ";N_{jet}",                20, 0,  20);
  h_nw           = new TH1D("nw",           ";N_{W tag}",              20, 0,  20);
  h_nb           = new TH1D("nb",           ";N_{b tag}",              20, 0,  20);
  h_ht_gen       = new TH1D("ht_gen",       ";H_{T}^{gen}",            200, 0,2000);
  h_ht_AK4Puppi  = new TH1D("ht_AK4Puppi",  ";H_{T}",                  200, 0,2000);
  h_ht_AK8Puppi  = new TH1D("ht_AK8Puppi",  ";H_{T}^{AK8}",            200, 0,2000);
  h_jet1_pt      = new TH1D("jet1_pt",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt      = new TH1D("jet2_pt",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt      = new TH1D("jet3_pt",      ";p_{T, jet3}",            200, 0,2000);
  
  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    histoname<<"jet1_pt_syst"<<i;
    title<<"Systematic variation #="<<i;
    vh_jet1_pt.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    vh_jet1_pt[i]->Sumw2();
  }
}


//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_analysis_histos(DataStruct& data, const unsigned int& syst_index, const double& weight)
{
  if (syst_index == 0) {
    // syst_index should only be non-0 if settings.varySystematics is true
    // in case of studying systematics, one should fill a different histogram for each syst_index
    // this variable can be used to chose the correct vector element in case there is a vector of histograms
    // It makes sense, to cut on syst_index == 0, for all ordinary plots
    // syst_index == 0 always guarantees, there are no variations in any systematics
    
    // Check what common variables are available in AnalysisBase.h
    // There a good chance a lot of stuff is already calculated!
    // Especially common object selections or variables to cut on in Analysis
    h_njet   ->Fill(nJet,        weight);
    h_nb     ->Fill(nMediumBTag, weight);
    h_nw     ->Fill(nTightWTag,  weight);
    
    h_ht_gen->Fill(data.evt.Gen_Ht,  weight);  // in ntuple
    h_ht_AK4Puppi->Fill(AK4Puppi_Ht, weight); // Calculated in AnalysisBase.h
    h_ht_AK8Puppi->Fill(AK8Puppi_Ht, weight); // Calculated in AnalysisBase.h
    
    if (_apply_ncut(3)) {
      h_jet1_pt->Fill(data.jetsAK4Puppi.Pt[iJet[0]], weight);
      h_jet2_pt->Fill(data.jetsAK4Puppi.Pt[iJet[1]], weight);
      h_jet3_pt->Fill(data.jetsAK4Puppi.Pt[iJet[2]], weight);
    }
  }
  
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  if (_apply_ncut(3)) vh_jet1_pt[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[0]], weight);
}

// Methods used by SmartHistos (Plotter)
// Can leave them empty
void
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst, 
			       const unsigned int& syst_index, std::string dirname, bool runOnSkim=false)
{
}

void
Analysis::load_analysis_histos(std::string inputfile)
{
}

void
Analysis::save_analysis_histos(bool draw=0)
{
}
