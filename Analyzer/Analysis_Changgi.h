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
    if (data.jetsAK8Puppi.looseJetID[data.jetsAK8Puppi.it]) {
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
  int N_Puppi = 0;
  while(data.jetsAK8Puppi.Loop()) if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=pt_threshold) ++N_Puppi;
  return (N_Puppi >= 1);
}

void
Analysis::define_selections(const DataStruct& data)
{
  // cut1: njet >= 2
  analysis_cuts.push_back({ .name="2jet",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0;
			      return 1;
			    } });

  // cut2: jet 1 pass loose jet id
  analysis_cuts.push_back({ .name="jet1_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      return data.jetsAK8Puppi.looseJetID[0]; 
			    } });

  // cut3: jet 2 pass loose jet id
  analysis_cuts.push_back({ .name="jet2_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      return data.jetsAK8Puppi.looseJetID[1];
			    } });

  // cut4: jet 1 eta < 2.4
  analysis_cuts.push_back({ .name="jet1_eta",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (fabs(data.jetsAK8Puppi.Eta[0])>=2.4) return 0;
			      return 1;
			    } });

  // cut5: jet 2 eta < 2.4
  analysis_cuts.push_back({ .name="jet2_eta",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (fabs(data.jetsAK8Puppi.Eta[1])>=2.4) return 0;
			      return 1;
			    } });

  // cut6: jet 1 pt >= 400
  analysis_cuts.push_back({ .name="jet1_pt",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (data.jetsAK8Puppi.Pt[0]<TOP_PT_CUT) return 0;
			      return 1;
			    } });

  // cut7: jet 2 pt >= 400
  analysis_cuts.push_back({ .name="jet2_pt",   .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (data.jetsAK8Puppi.Pt[1]<TOP_PT_CUT) return 0;
			      return 1;
			    } });

  // cut8: 110 <= jet 1 mass (softdrop) < 210
  analysis_cuts.push_back({ .name="jet1_mass", .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (data.jetsAK8Puppi.softDropMass[0]<TOP_SD_MASS_CUT_LOW) return 0;
			      if (data.jetsAK8Puppi.softDropMass[0]>=TOP_SD_MASS_CUT_HIGH) return 0;
			      return 1;
			    } });

  // cut9: 110 <= jet 2 mass (softdrop) < 210
  analysis_cuts.push_back({ .name="jet2_mass", .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (data.jetsAK8Puppi.softDropMass[1]<TOP_SD_MASS_CUT_LOW) return 0;
			      if (data.jetsAK8Puppi.softDropMass[1]>=TOP_SD_MASS_CUT_HIGH) return 0;
			      return 1;
			    } });

/*
  // cut10: Full-hadronic trigger
  analysis_cuts.push_back({ .name="hlt_ak8ht700_mass50", .func = [&data](){
			      // Define cut function here:
			      return data.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; 
			    } });
*/

  // cut11: | DeltaPhi | < DPHI_CUT
  analysis_cuts.push_back({ .name="delta_phi", .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      TLorentzVector jet1, jet2;
			      jet1.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[0], data.jetsAK8Puppi.Eta[0], data.jetsAK8Puppi.Phi[0], data.jetsAK8Puppi.E[0]);
			      jet2.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[1], data.jetsAK8Puppi.Eta[1], data.jetsAK8Puppi.Phi[1], data.jetsAK8Puppi.E[1]);
			      double dPhi = fabs(jet1.DeltaPhi(jet2));
			      if (dPhi>=DPHI_CUT) return 0;
			      return 1;
			    } });

  // cut12: jet 1 tau32 < TOP_TAU32_CUT
  analysis_cuts.push_back({ .name="jet1_tau32",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      double tau32 = data.jetsAK8Puppi.tau3[0];
			      if (data.jetsAK8Puppi.tau2[0]!=0) tau32 /= data.jetsAK8Puppi.tau2[0];
			      else tau32 = 9999;
			      if (tau32>=TOP_TAU32_CUT) return 0;
			      return 1;
			    } });

  // cut13: jet 2 tau32 < TOP_TAU32_CUT
  analysis_cuts.push_back({ .name="jet2_tau32",   .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      double tau32 = data.jetsAK8Puppi.tau3[1];
			      if (data.jetsAK8Puppi.tau2[1]!=0) tau32 /= data.jetsAK8Puppi.tau2[1];
			      else tau32 = 9999;
			      if (tau32>=TOP_TAU32_CUT) return 0;
			      return 1;
			    } });
/*
  // cut14: R >= R_CUT
  analysis_cuts.push_back({ .name="R", .func = [&data](){
			      // Define cut function here:
			      if (data.evt.R>=R_CUT) return 0;
			      return 1;
			    } });
*/
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
TH1D* h_ht_AK4Puppi;
TH1D* h_ht_AK8Puppi;
TH1D* h_AK4Puppi_jet1_pt;
TH1D* h_AK8Puppi_jet1_pt;
TH1D* h_AK4Puppi_jet2_pt;
TH1D* h_AK8Puppi_jet2_pt;
std::vector<TH1D*> vh_AK4Puppi_jet1_pt;
std::vector<TH1D*> vh_AK8Puppi_jet1_pt;
TH1D* h_AK8Puppi_MR;
TH1D* h_AK8Puppi_MTR;
TH1D* h_AK8Puppi_R;
TH1D* h_AK8Puppi_R2;
TH1D* h_AK8Puppi_tau32;
TH1D* h_AK8Puppi_tau31;
TH1D* h_AK8Puppi_tau21;
TH1D* h_MET;
TH1D* h_softDropMass;
TH2D* h_GluinoLSPMass;
std::vector<TH2D*> vh_abcd;
//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  h_njet         = new TH1D("njet",         ";N_{AK8 (Puppi), loose ID}",  20, 0,  20);
  h_nhadtop      = new TH1D("nhadtop",      ";N_{top tag}",                20, 0,  20);
  h_nhadw        = new TH1D("nhadw",        ";N_{W tag}",                  20, 0,  20);
  h_ht_gen       = new TH1D("ht_gen",       ";H_{T}^{gen}",            800, 0,4000);
  h_ht_AK4Puppi  = new TH1D("ht_AK4Puppi",  ";H_{T}^{AK4 (Puppi)}",    800, 0,4000);
  h_ht_AK8Puppi  = new TH1D("ht_AK8Puppi",  ";H_{T}^{AK8 (Puppi)}",    800, 0,4000);
  h_AK4Puppi_jet1_pt = new TH1D("jet1_AK4Puppi_pt",";p_{T, jet1}",           200, 0,2000);
  h_AK8Puppi_jet1_pt = new TH1D("jet1_AK8Puppi_pt",";p_{T, jet1}",           200, 0,2000);
  h_AK4Puppi_jet2_pt = new TH1D("jet2_AK4Puppi_pt",";p_{T, jet2}",           200, 0,2000);
  h_AK8Puppi_jet2_pt = new TH1D("jet2_AK8Puppi_pt",";p_{T, jet2}",           200, 0,2000);

  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream h1, h2, h3, h4, title;
    h1<<"jet1_pt_syst"<<i;
    h2<<"jet1_AK8Puppi_pt_syst"<<i;
    h3<<"jet1_AK4_pt_syst"<<i;
    h4<<"jet1_AK4Puppi_pt_syst"<<i;
    title<<"Systematic variation #="<<i<<";p_{T, jet1}";
    vh_AK8Puppi_jet1_pt.push_back(new TH1D(h2.str().c_str(), title.str().c_str(), 200, 0,2000));
    vh_AK4Puppi_jet1_pt.push_back(new TH1D(h4.str().c_str(), title.str().c_str(), 200, 0,2000));
  }

  h_AK8Puppi_MR      = new TH1D("AK8Puppi_MR",   ";MR_{AK8Puppi}",         200, 0,4000);
  h_AK8Puppi_MTR     = new TH1D("AK8Puppi_MTR",  ";MTR_{AK8Puppi}",        200, 0,2000);
  h_AK8Puppi_R       = new TH1D("AK8Puppi_R",    ";R_{AK8Puppi}",          500, 0,1);
  h_AK8Puppi_R2      = new TH1D("AK8Puppi_R2",   ";R2_{AK8Puppi}",         500, 0,1);

  h_AK8Puppi_tau32 = new TH1D("tau32", "", 200,0,1);
  h_AK8Puppi_tau31 = new TH1D("tau31", "", 200,0,1);
  h_AK8Puppi_tau21 = new TH1D("tau21", "", 200,0,1);
  h_MET = new TH1D("met", "", 400,0,2000);
  h_softDropMass = new TH1D("softDropMass", "", 100,0,500);
  h_GluinoLSPMass = new TH2D("GluinoLSPMass","",34,600,2300,64,0,1600);


  // ABCD background estimation
  Double_t R_bins[3] = { 0.2, 0.4, 100 };
  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream h_abcd_name, h_abcd_title;
    h_abcd_name<<"abcd_syst"<<i;
    h_abcd_title<<"Systematic variation #="<<i;
    vh_abcd.push_back(new TH2D(h_abcd_name.str().c_str(),  (h_abcd_title.str()+";R (AK8 Puppi);Both jets pass tau32 cuts").c_str(), 2,R_bins, 2,0,2));
  }

h_njet->Sumw2();
h_nhadtop->Sumw2();
h_nhadw->Sumw2();
h_ht_gen->Sumw2();
h_ht_AK4Puppi->Sumw2();
h_ht_AK8Puppi->Sumw2();
h_AK4Puppi_jet1_pt->Sumw2();
h_AK8Puppi_jet1_pt->Sumw2();
h_AK4Puppi_jet2_pt->Sumw2();
h_AK8Puppi_jet2_pt->Sumw2();

h_AK8Puppi_MR->Sumw2();
h_AK8Puppi_MTR->Sumw2();
h_AK8Puppi_R->Sumw2();
h_AK8Puppi_R2->Sumw2();
h_AK8Puppi_tau32->Sumw2();
h_AK8Puppi_tau31->Sumw2();
h_AK8Puppi_tau21->Sumw2();
h_MET->Sumw2();
h_softDropMass->Sumw2();
h_GluinoLSPMass->Sumw2();

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
    
    if (_apply_ncut(analysis_cuts.size())){
    h_GluinoLSPMass->Fill(data.evt.SUSY_Gluino_Mass,data.evt.SUSY_LSP_Mass,weight);
    h_ht_gen->Fill(data.evt.Gen_Ht, weight);  // in ntuple
    h_ht_AK4Puppi->Fill(AK4Puppi_Ht, weight); // Calculated in AnalysisBase.h
    h_ht_AK8Puppi->Fill(AK8Puppi_Ht, weight); // Calculated in AnalysisBase.h
    
    h_AK4Puppi_jet1_pt->Fill(data.jetsAK4Puppi.Pt[0], weight);
    h_AK8Puppi_jet1_pt->Fill(data.jetsAK8Puppi.Pt[0], weight);
    h_AK4Puppi_jet2_pt->Fill(data.jetsAK4Puppi.Pt[1], weight);
    h_AK8Puppi_jet2_pt->Fill(data.jetsAK8Puppi.Pt[1], weight);

    h_AK8Puppi_MR->Fill(data.evt.MR, weight);
    h_AK8Puppi_MTR->Fill(data.evt.MTR, weight);
    h_AK8Puppi_R->Fill(data.evt.R, weight);
    h_AK8Puppi_R2->Fill(data.evt.R2, weight);

    h_AK8Puppi_tau32->Fill(data.jetsAK8Puppi.tau3.at(0)/data.jetsAK8Puppi.tau2.at(0),weight);
    h_AK8Puppi_tau31->Fill(data.jetsAK8Puppi.tau3.at(0)/data.jetsAK8Puppi.tau1.at(0),weight);
    h_AK8Puppi_tau21->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),weight);
    h_MET->Fill(data.met.Pt.at(0),weight);
    h_softDropMass->Fill(data.jetsAK8Puppi.softDropMass.at(0),weight);
   }
  }

  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  if (_apply_ncut(analysis_cuts.size())){
     vh_AK8Puppi_jet1_pt[syst_index]->Fill(data.jetsAK8Puppi.Pt[0], weight);
     vh_AK4Puppi_jet1_pt[syst_index]->Fill(data.jetsAK4Puppi.Pt[0], weight);
  }

  // ABCD background estimation
  if (_apply_ncut(analysis_cuts.size()-3)) vh_abcd[syst_index]->Fill(data.evt.R, passHadTopTag[0]&&passHadTopTag[1], weight);
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
