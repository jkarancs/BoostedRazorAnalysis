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

unsigned int nLooseIDHadTopTagJets;
unsigned int nLooseIDHadWTagJets;

void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{
  // Jet variables (initialize)
  nLooseIDHadTopTagJets = nLooseIDHadWTagJets = 0;

  // Loop on jets
  while(data.jetsAK8Puppi.Loop()) {
    if (passLooseJetID[data.jetsAK8Puppi.it]) {
      if (passHadTopTag[data.jetsAK8Puppi.it]) ++nLooseIDHadTopTagJets;
      if (passHadWTag[data.jetsAK8Puppi.it]) ++nLooseIDHadWTagJets;
    }
  } // end of AK8 jet loop
}

//_______________________________________________________
//          Define Analysis specific weights
double
Analysis::get_analysis_weight(DataStruct& data)
{
  double w = 1;
  
  /*
  //____________________________________________________
  //               Jet pT reweighting
  
  // Use linear function calculated by scripts/JetPtScaleFactors.C script
  // reweight event corresponding to the product of SFs for each jet
  // linear function: p0 + p1 * jet pt
  
  const double p0 = 0.970718;
  const double p1 = -0.000331894;
  while(data.jetsAK8.Loop()) {
    w *= p0 + p1 * data.jetsAK8.Pt[data.jetsAK8.it];
  }
  */
  
  return w;
}

//_______________________________________________________
//          Define Analysis event selection cuts

bool
Analysis::pass_skimming(DataStruct& data)
{
  float pt_threshold = 300;
  int N_CHS = 0, N_Puppi = 0;
  while(data.jetsAK8.Loop()) if (data.jetsAK8.Pt[data.jetsAK8.it]>=pt_threshold) ++N_CHS;
  while(data.jetsAK8Puppi.Loop()) if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=pt_threshold) ++N_Puppi;
  return (N_CHS >= 1 || N_Puppi >= 1);
}

void
Analysis::define_selections(const DataStruct& data)
{
  // cut1: njet >= 1
  analysis_cuts.push_back({ .name="1jet",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0;
			      return 1;
			    } });

  // cut2: jet 1 pass loose jet id
  analysis_cuts.push_back({ .name="jet1_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      return passLooseJetID[0];
			    } });

  // cut3: jet 1 eta < 2.4
  analysis_cuts.push_back({ .name="jet1_eta",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (fabs(data.jetsAK8Puppi.Eta[0])>=2.4) return 0;
			      return 1;
			    } });

  // cut4: jet 1 pt >= 400
  analysis_cuts.push_back({ .name="jet1_pt",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (data.jetsAK8Puppi.Pt[0]<400) return 0;
			      return 1;
			    } });

  // cut5: Full-hadronic trigger
  analysis_cuts.push_back({ .name="hlt_ak8ht700_mass50", .func = [&data](){
			      // Define cut function here:
			      return data.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; 
			    } });

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
  return 0;
}

//_______________________________________________________
//                 List of Histograms
TH1D* h_njet;
TH1D* h_nhadtop;
TH1D* h_nhadw;
TH1D* h_ht_gen;
TH1D* h_ht_AK4;
TH1D* h_ht_AK4Puppi;
TH1D* h_ht_AK8;
TH1D* h_ht_AK8Puppi;
TH1D* h_jet1_pt;
std::vector<TH1D*> vh_jet1_pt;

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  h_njet         = new TH1D("njet",         ";N_{AK8 (Puppi), loose ID}",  20, 0,  20);
  h_nhadtop      = new TH1D("nhadtop",      ";N_{top tag}",                20, 0,  20);
  h_nhadw        = new TH1D("nhadw",        ";N_{W tag}",                  20, 0,  20);
  h_ht_gen       = new TH1D("ht_gen",       ";H_{T}^{gen}",            200, 0,2000);
  h_ht_AK4       = new TH1D("ht_AK4",       ";H_{T}^{AK4 (CHS)}",      200, 0,2000);
  h_ht_AK4Puppi  = new TH1D("ht_AK4Puppi",  ";H_{T}^{AK4 (Puppi)}",    200, 0,2000);
  h_ht_AK8       = new TH1D("ht_AK8",       ";H_{T}^{AK8 (CHS)}",      200, 0,2000);
  h_ht_AK8Puppi  = new TH1D("ht_AK8Puppi",  ";H_{T}^{AK8 (Puppi)}",    200, 0,2000);
  h_jet1_pt      = new TH1D("jet1_pt",      ";p_{T, jet1}",            200, 0,2000);
  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    histoname<<"jet1_pt_syst"<<i;
    title<<"Systematic variation #="<<i<<";p_{T, jet1}";
    vh_jet1_pt.push_back(new TH1D(histoname.str().c_str(), title.str().c_str(), 200, 0,2000));
    vh_jet1_pt[i]->Sumw2();
  }
  
  h_njet->Sumw2();
  h_nhadtop->Sumw2();
  h_nhadw->Sumw2();
  h_ht_gen->Sumw2();
  h_ht_AK4->Sumw2();
  h_ht_AK4Puppi->Sumw2();
  h_ht_AK8->Sumw2();
  h_ht_AK8Puppi->Sumw2();
  h_jet1_pt->Sumw2();
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
    
    h_njet   ->Fill(nLooseJet, weight);
    h_nhadtop->Fill(nLooseIDHadTopTagJets, weight);
    h_nhadw  ->Fill(nLooseIDHadWTagJets, weight);
    
    h_ht_gen->Fill(data.evt.Gen_Ht, weight);  // in ntuple
    h_ht_AK4->Fill(AK4_Ht, weight);           // Calculated in AnalysisBase.h
    h_ht_AK4Puppi->Fill(AK4Puppi_Ht, weight); // Calculated in AnalysisBase.h
    h_ht_AK8->Fill(data.evt.Ht, weight);      // in ntuple, AK8 CHS is default, will switch to Puppi
    h_ht_AK8Puppi->Fill(AK8Puppi_Ht, weight); // Calculated in AnalysisBase.h
    
    if (_apply_ncut(2)) h_jet1_pt->Fill(data.jetsAK8Puppi.Pt[0], weight);
  }
  
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  if (_apply_ncut(2)) vh_jet1_pt[syst_index]->Fill(data.jetsAK8Puppi.Pt[0], weight);
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
