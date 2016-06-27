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
  // cut1: njet >= 3 // cut from thesis
  analysis_cuts.push_back({ .name="1jet",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<3) return 0;
			      return 1;
			    } });

  // cut2: jet 1 pass loose jet id
  analysis_cuts.push_back({ .name="jet1_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      return passLooseJetID[0];
			    } });

  // cut3: jet 2 pass loose jet id
  analysis_cuts.push_back({ .name="jet2_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      return passLooseJetID[1];
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
  // cut6: one of jet pt >= 200  // cut from thesis
  analysis_cuts.push_back({ .name="jet1_pt",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (data.jetsAK8Puppi.Pt[0]<200) return 0;
			      return 1;
			    } });

  // cut7: 110 <= jet 1 mass (softdrop) < 210
  // cut7: 70  <= jet 1 mass (softdrop) < 100 // cut from thesis
  analysis_cuts.push_back({ .name="jet1_mass", .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<1) return 0; // for safety
            if (data.jetsAK8Puppi.softDropMass.at(0)< 100.) return 0;
            if (data.jetsAK8Puppi.softDropMass.at(0)>= 70.) return 0;
            return 1;
          } });

  // cut8: 110 <= jet 2 mass (softdrop) < 210
  // cut8: 70  <= jet 2 mass (softdrop) < 100 // cut from thesis
  analysis_cuts.push_back({ .name="jet2_mass", .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            if (data.jetsAK8Puppi.softDropMass.at(1)< 100.) return 0;
            if (data.jetsAK8Puppi.softDropMass.at(1)>= 70.) return 0;
            return 1;
          } });

  // cut9: Full-hadronic trigger
  analysis_cuts.push_back({ .name="hlt_ak8ht700_mass50", .func = [&data](){
  //analysis_cuts.push_back({ .name="hlt_pfht800", .func = [&data](){
  //analysis_cuts.push_back({ .name="hlt_pfjet450", .func = [&data](){
			      // Define cut function here:
			      return data.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; 
			      //return data.hlt.PFHT800==1; 
			      //return data.hlt.PFJet450==1; 
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
TH1D* h_AK4_jet1_pt;
std::vector<TH1D*> vh_AK4_jet1_pt;
TH1D* h_AK4Puppi_jet1_pt;
std::vector<TH1D*> vh_AK4Puppi_jet1_pt;
TH1D* h_AK8Puppi_jet1_pt;
std::vector<TH1D*> vh_AK8Puppi_jet1_pt;

TH1D* h_MR;
TH1D* h_MTR;
TH1D* h_R;
TH1D* h_R2;
TH1D* h_AK8Puppi_MR;
TH1D* h_AK8Puppi_MTR;
TH1D* h_AK8Puppi_R;
TH1D* h_AK8Puppi_R2;
TH1D* h_AK4_MR;
TH1D* h_AK4_MTR;
TH1D* h_AK4_R;
TH1D* h_AK4_R2;
TH1D* h_AK8_tau1;
TH1D* h_AK8_tau2;
TH1D* h_AK8_tau3;
/*
TH2D* h_jet1pt_jet1pt_AK4;
TH2D* h_jet1pt_jet1pt_AK4Puppi;
TH2D* h_jet1pt_jet1pt_AK8Puppi;
TH2D* h_R_R2;
TH2D* h_MR_MR_AK4;
TH2D* h_MTR_MTR_AK4;
TH2D* h_R_R_AK4;
TH2D* h_MR_MR_AK8Puppi;
TH2D* h_MTR_MTR_AK8Puppi;
TH2D* h_R_R_AK8Puppi;
*/
TH1D* h_AK8_tau32;
TH1D* h_AK8_tau31;
TH1D* h_AK8_tau21;
TH1D* h_MET;
TH1D* h_softDropMass;
TH1D* h_trimmedMass;
TH1D* h_prunedMass;
TH1D* h_filteredMass;

TH2D* h_AK8R_DPhi;
TH2D* h_AK8PuppiR_DPhi;
TH2D* h_AK4R_DPhi;
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
  h_AK4_jet1_pt = new TH1D("jet1_AK4_pt",";p_{T, jet1}",           200, 0,2000);
  h_AK4Puppi_jet1_pt = new TH1D("jet1_AK4Puppi_pt",";p_{T, jet1}",           200, 0,2000);
  h_AK8Puppi_jet1_pt = new TH1D("jet1_AK8Puppi_pt",";p_{T, jet1}",           200, 0,2000);

  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    histoname<<"jet1_pt_syst"<<i;
    title<<"Systematic variation #="<<i<<";p_{T, jet1}";
    vh_jet1_pt.push_back(new TH1D(histoname.str().c_str(), title.str().c_str(), 200, 0,2000));
    vh_AK8Puppi_jet1_pt.push_back(new TH1D(histoname.str().c_str(), title.str().c_str(), 200, 0,2000));
    vh_AK4_jet1_pt.push_back(new TH1D(histoname.str().c_str(), title.str().c_str(), 200, 0,2000));
    vh_AK4Puppi_jet1_pt.push_back(new TH1D(histoname.str().c_str(), title.str().c_str(), 200, 0,2000));
    vh_jet1_pt[i]->Sumw2();
    vh_AK8Puppi_jet1_pt[i]->Sumw2();
    vh_AK4_jet1_pt[i]->Sumw2();
    vh_AK4Puppi_jet1_pt[i]->Sumw2();
  }

  h_MR      = new TH1D("MR",     ";MR_{AK8?}",             200, 0,2000);
  h_MTR     = new TH1D("MTR",    ";MTR_{AK8?}",            200, 0,2000);
  h_R       = new TH1D("R",      ";R_{AK8?}",              200, -2,2);
  h_R2      = new TH1D("R2",     ";R2_{AK8?}",             200, -2,2);

  h_AK8Puppi_MR      = new TH1D("AK8Puppi_MR",   ";MR_{AK8Puppi}",         200, 0,2000);
  h_AK8Puppi_MTR     = new TH1D("AK8Puppi_MTR",  ";MTR_{AK8Puppi}",        200, 0,2000);
  h_AK8Puppi_R       = new TH1D("AK8Puppi_R",    ";R_{AK8Puppi}",          200, -2,2);
  h_AK8Puppi_R2      = new TH1D("AK8Puppi_R2",   ";R2_{AK8Puppi}",         200, -2,2);

  h_AK4_MR      = new TH1D("AK4_MR",     ";MR_{AK4}",              200, 0,2000);
  h_AK4_MTR     = new TH1D("AK4_MTR",    ";MTR_{AK4}",             200, 0,2000);
  h_AK4_R       = new TH1D("AK4_R",      ";R_{AK4}",               200, -2,2);
  h_AK4_R2      = new TH1D("AK4_R2",     ";R2_{AK4}",              200, -2,2);

  h_AK8_tau1    = new TH1D("tau1",";tau1_{AK8}",           200,-1,1);
  h_AK8_tau2    = new TH1D("tau2",";tau1_{AK8}",           200,-1,1);
  h_AK8_tau3    = new TH1D("tau3",";tau1_{AK8}",           200,-1,1);
/*
  h_jet1pt_jet1pt_AK4 = new TH2D("jet1pt_jet1pt_AK4", "", 200,0,2000,200,0,2000);
  h_jet1pt_jet1pt_AK8Puppi = new TH2D("jet1pt_jet1pt_AK8Puppi", "", 200,0,2000,200,0,2000);
  h_jet1pt_jet1pt_AK4Puppi = new TH2D("jet1pt_jet1pt_AK4Puppi", "", 200,0,2000,200,0,2000);
  h_R_R2 = new TH2D("R_R2", "", 200,-2,2,200,-2,2);
  h_MR_MR_AK4 = new TH2D("MR_MR_AK4", "", 200,0,2000,200,0,2000);
  h_MTR_MTR_AK4 = new TH2D("MTR_MTR_AK4", "", 200,0,2000,200,0,2000);
  h_R_R_AK4 = new TH2D("R_R_AK4", "", 200,-2,2,200,-2,2);
  h_MR_MR_AK8Puppi = new TH2D("MR_MR_AK8Puppi", "", 200,0,2000,200,0,2000);
  h_MTR_MTR_AK8Puppi = new TH2D("MTR_MTR_AK8Puppi", "", 200,0,2000,200,0,2000);
  h_R_R_AK8Puppi = new TH2D("R_R_AK8Puppi", "", 200,-2,2,200,-2,2);
*/
  h_AK8_tau32 = new TH1D("tau32", "", 200,-1,1);
  h_AK8_tau31 = new TH1D("tau31", "", 200,-1,1);
  h_AK8_tau21 = new TH1D("tau21", "", 200,-1,1);
  h_MET = new TH1D("met", "", 200,0,2000);
  h_softDropMass = new TH1D("softDropMass", "", 200,0,2000);
  h_trimmedMass = new TH1D("trimmedMass", "", 200,0,2000);
  h_prunedMass = new TH1D("prunedMass", "", 200,0,2000);
  h_filteredMass = new TH1D("filteredMass", "", 200,0,2000);

  h_AK8R_DPhi = new TH2D("h_AK8R_DPhi", "", 350, 0, 3.5, 150, 0, 1.5);
  h_AK8PuppiR_DPhi = new TH2D("h_AK8PuppiR_DPhi", "", 350, 0, 3.5, 150, 0, 1.5);
  h_AK4R_DPhi = new TH2D("h_AK4R_DPhi", "", 350, 0, 3.5, 150, 0, 1.5);

  h_njet->Sumw2();
  h_nhadtop->Sumw2();
  h_nhadw->Sumw2();
  h_ht_gen->Sumw2();
  h_ht_AK4->Sumw2();
  h_ht_AK4Puppi->Sumw2();
  h_ht_AK8->Sumw2();
  h_ht_AK8Puppi->Sumw2();
  h_jet1_pt->Sumw2();
  h_AK4_jet1_pt->Sumw2();
  h_MR->Sumw2();
  h_MTR->Sumw2();
  h_R->Sumw2();
  h_R2->Sumw2();
  h_AK8Puppi_MR->Sumw2();
  h_AK8Puppi_MTR->Sumw2();
  h_AK8Puppi_R->Sumw2();
  h_AK8Puppi_R2->Sumw2();
  h_AK4_MR->Sumw2();
  h_AK4_MTR->Sumw2();
  h_AK4_R->Sumw2();
  h_AK4_R2->Sumw2();
  h_AK8_tau1->Sumw2();
  h_AK8_tau2->Sumw2();
  h_AK8_tau3->Sumw2();
/*
  h_jet1pt_jet1pt_AK4->Sumw2();
  h_jet1pt_jet1pt_AK4Puppi->Sumw2();
  h_jet1pt_jet1pt_AK8Puppi->Sumw2();
  h_R_R2->Sumw2();
  h_MR_MR_AK4->Sumw2();
  h_MTR_MTR_AK4->Sumw2();
  h_R_R_AK4->Sumw2();
  h_MR_MR_AK8Puppi->Sumw2();
  h_MTR_MTR_AK8Puppi->Sumw2();
  h_R_R_AK8Puppi->Sumw2();
*/
  h_AK8_tau32->Sumw2();
  h_AK8_tau31->Sumw2();
  h_AK8_tau21->Sumw2();
  h_MET->Sumw2();
  h_softDropMass->Sumw2();
  h_trimmedMass->Sumw2();
  h_prunedMass->Sumw2();
  h_filteredMass->Sumw2();

  h_AK8R_DPhi->Sumw2();
  h_AK8PuppiR_DPhi->Sumw2();
  h_AK4R_DPhi->Sumw2();
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
    
    if (_apply_ncut(9)){
    h_ht_gen->Fill(data.evt.Gen_Ht, weight);  // in ntuple
    h_ht_AK4->Fill(AK4_Ht, weight);           // Calculated in AnalysisBase.h
    h_ht_AK4Puppi->Fill(AK4Puppi_Ht, weight); // Calculated in AnalysisBase.h
    h_ht_AK8->Fill(data.evt.Ht, weight);      // in ntuple, AK8 CHS is default, will switch to Puppi
    h_ht_AK8Puppi->Fill(AK8Puppi_Ht, weight); // Calculated in AnalysisBase.h
    
      h_jet1_pt->Fill(data.jetsAK8Puppi.Pt[0], weight);
    h_jet1_pt->Fill(data.jetsAK8.Pt[0], weight);
    h_AK4_jet1_pt->Fill(data.jetsAK4.Pt[0], weight);
    h_AK4Puppi_jet1_pt->Fill(data.jetsAK4Puppi.Pt[0], weight);
    h_AK8Puppi_jet1_pt->Fill(data.jetsAK8Puppi.Pt[0], weight);

    h_MR->Fill(data.evt.MR, weight);
    h_MTR->Fill(data.evt.MTR, weight);
    h_R->Fill(data.evt.R, weight);
    h_R2->Fill(data.evt.R2, weight);

    h_AK8Puppi_MR->Fill(data.evt.AK8Puppi_MR, weight);
    h_AK8Puppi_MTR->Fill(data.evt.AK8Puppi_MTR, weight);
    h_AK8Puppi_R->Fill(data.evt.AK8Puppi_R, weight);
    h_AK8Puppi_R2->Fill(data.evt.AK8Puppi_R2, weight);

    h_AK4_MR->Fill(data.evt.AK4_MR, weight);
    h_AK4_MTR->Fill(data.evt.AK4_MTR, weight);
    h_AK4_R->Fill(data.evt.AK4_R, weight);
    h_AK4_R2->Fill(data.evt.AK4_R2, weight);

    h_AK8_tau1->Fill(data.jetsAK8.tau1.at(0),weight);
    h_AK8_tau2->Fill(data.jetsAK8.tau2.at(0),weight);
    h_AK8_tau3->Fill(data.jetsAK8.tau3.at(0),weight);
/*
  h_jet1pt_jet1pt_AK4->Fill(data.jetsAK8.Pt[0],data.jetsAK4.Pt[0],weight);
  h_jet1pt_jet1pt_AK4Puppi->Fill(data.jetsAK8.Pt[0],data.jetsAK4Puppi.Pt[0],weight);
  h_jet1pt_jet1pt_AK8Puppi->Fill(data.jetsAK8.Pt[0],data.jetsAK8Puppi.Pt[0],weight);
  h_R_R2->Fill(data.evt.R,data.evt.R2,weight);
  h_MR_MR_AK4->Fill(data.evt.MR,data.evt.AK4_MR,weight);
  h_MTR_MTR_AK4->Fill(data.evt.MTR,data.evt.AK4_MTR,weight);
  h_R_R_AK4->Fill(data.evt.R,data.evt.AK4_R,weight);
  h_MR_MR_AK8Puppi->Fill(data.evt.MR,data.evt.AK8Puppi_MR,weight);
  h_MTR_MTR_AK8Puppi->Fill(data.evt.MTR,data.evt.AK8Puppi_MTR,weight);
  h_R_R_AK8Puppi->Fill(data.evt.R,data.evt.AK8Puppi_R,weight);
*/
    h_AK8_tau32->Fill(data.jetsAK8.tau3.at(0)/data.jetsAK8.tau2.at(0),weight);
    h_AK8_tau31->Fill(data.jetsAK8.tau3.at(0)/data.jetsAK8.tau1.at(0),weight);
    h_AK8_tau21->Fill(data.jetsAK8.tau2.at(0)/data.jetsAK8.tau1.at(0),weight);
    h_MET->Fill(data.met.Pt.at(0),weight);
    h_softDropMass->Fill(data.jetsAK8.softDropMass.at(0),weight);
    h_trimmedMass->Fill(data.jetsAK8.trimmedMass.at(0),weight);
    h_prunedMass->Fill(data.jetsAK8.prunedMass.at(0),weight);
    h_filteredMass->Fill(data.jetsAK8.filteredMass.at(0),weight);

     float DPhi = -9999.;
     TLorentzVector jet1, jet2;

     jet1.SetPtEtaPhiE(data.jetsAK8.Pt[0], data.jetsAK8.Eta[0], data.jetsAK8.Phi[0], data.jetsAK8.E[0]);
     jet2.SetPtEtaPhiE(data.jetsAK8.Pt[1], data.jetsAK8.Eta[1], data.jetsAK8.Phi[1], data.jetsAK8.E[1]);
     DPhi = fabs(jet1.DeltaPhi(jet2));
     h_AK8R_DPhi->Fill(DPhi, data.evt.R, weight);
     jet1.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[0], data.jetsAK8Puppi.Eta[0], data.jetsAK8Puppi.Phi[0], data.jetsAK8Puppi.E[0]);
     jet2.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[1], data.jetsAK8Puppi.Eta[1], data.jetsAK8Puppi.Phi[1], data.jetsAK8Puppi.E[1]);
     DPhi = fabs(jet1.DeltaPhi(jet2));
     h_AK8PuppiR_DPhi->Fill(DPhi, data.evt.AK8Puppi_R, weight);
     jet1.SetPtEtaPhiE(data.jetsAK4.Pt[0], data.jetsAK4.Eta[0], data.jetsAK4.Phi[0], data.jetsAK4.E[0]);
     jet2.SetPtEtaPhiE(data.jetsAK4.Pt[1], data.jetsAK4.Eta[1], data.jetsAK4.Phi[1], data.jetsAK4.E[1]);
     DPhi = fabs(jet1.DeltaPhi(jet2));
     h_AK4R_DPhi->Fill(DPhi, data.evt.AK4_R, weight);
   }
  }
  
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  if (_apply_ncut(9)){
     vh_jet1_pt[syst_index]->Fill(data.jetsAK8.Pt[0], weight);
     vh_AK8Puppi_jet1_pt[syst_index]->Fill(data.jetsAK8Puppi.Pt[0], weight);
     vh_AK4_jet1_pt[syst_index]->Fill(data.jetsAK4.Pt[0], weight);
     vh_AK4Puppi_jet1_pt[syst_index]->Fill(data.jetsAK4Puppi.Pt[0], weight);
  }
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
