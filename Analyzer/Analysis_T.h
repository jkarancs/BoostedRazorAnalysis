#include "TLorentzVector.h"
#include "common/AnalysisBase.h"
#include "common/SmartHistos.h"

SmartHistos sh;


//_______________________________________________________
//                       Constructor
Analysis::Analysis() : AnalysisBase() { }


//_______________________________________________________
//                       Destructor
Analysis::~Analysis() { }

//_______________________________________________________
//                  Calculate variables


unsigned int nLooseIDHadTopTagJets;  


void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{

 nLooseIDHadTopTagJets = 0;

 // Loop on AK8Puppi jets

  while(data.jetsAK8Puppi.Loop()) {
    if (data.jetsAK8Puppi.looseJetID[data.jetsAK8Puppi.it]) {
      if (passHadTopTag[data.jetsAK8Puppi.it]) ++nLooseIDHadTopTagJets; 
      
    }
  } // end of AK8 jet loop
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
//                Define Skimming cuts
//   (Not needed, unless you want to skim the ntuple)

bool
Analysis::pass_skimming(DataStruct& data)
{
  if (!(nJetAK8>=1)) return 0;
  if (!(nJet>=2)) return 0;
  return 1;

  float pt_threshold = 300;
  int N_Puppi = 0;
  while(data.jetsAK8Puppi.Loop()) if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=pt_threshold) ++N_Puppi;
  return (N_Puppi >= 1);

  // Signal skim
  //return apply_all_cuts("S");
}

//_______________________________________________________
//          Define Analysis event selection cuts
//     Can define all sorts of Signal/Control regions

void
Analysis::define_selections(const DataStruct& data)
{
 analysis_cuts.clear();

  // Define here cuts that are common in all Signal/Control regions
  // MET Filters, etc. are already applied in AnalysisBase.h, See baseline_cuts
  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                    }}); // Similar to pt>200, one AK8 jet has pt>170
  baseline_cuts.push_back({ .name="Skim_2Jet",       .func = []    { return nJet>=2;                       }});
  baseline_cuts.push_back({ .name="Baseline_3Jet",   .func = []    { return nJet>=3;                       }}); // Separate cut, so one can exclude (N-1)
 // baseline_cuts.push_back({ .name="Baseline_MR_R2",  .func = [&d]  { return d.evt.MR>800 && d.evt.R2>0.08; }});
  
  // S: Signal region
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                   }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                    }});
  //analysis_cuts["S"].push_back({ .name="0TauTrk",    .func = []    { return;  }});
  analysis_cuts['S'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                }});
  analysis_cuts['S'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                 }});
  //analysis_cuts['S'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['S'].push_back({ .name="mDPhi>=0p4", .func = []    { return minDeltaPhi>=0.4;              }}); // Decreased it to the AK4 cone size (from 0.5)
  
  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepTight==1;                  }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==1;                 }});
  analysis_cuts['W'].push_back({ .name="1Wpre",      .func = []    { return nWPreTag>=1;                   }});
  //analysis_cuts['W'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['W'].push_back({ .name="mDPhi>=0p4",  .func = []   { return minDeltaPhi>=0.4;              }}); // Decreased it to the AK4 cone size (from 0.5)
  analysis_cuts['W'].push_back({ .name="30<=MT<100",  .func = []   { return MT>=30 && MT<100;              }});
  
  // T: Top enriched control sample
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepTight==1;                  }});
  analysis_cuts['T'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                }});
  analysis_cuts['T'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                 }});
  //analysis_cuts['T'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['T'].push_back({ .name="mDPhi>=0p4", .func = []    { return minDeltaPhi>=0.4;              }}); // Decreased it to the AK4 cone size (from 0.5)
  analysis_cuts['T'].push_back({ .name="MT<100",     .func = []    { return MT<100;                        }});
  
  // Q: QCD enriched control sample
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                   }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                    }});
  //analysis_cuts['Q'].push_back({ .name="0TauTrk",    .func = []    { return;  }});
  analysis_cuts['Q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                 }});
  analysis_cuts['Q'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                 }});
  //analysis_cuts['Q'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['Q'].push_back({ .name="mDPhi<0.25", .func = []    { return minDeltaPhi<0.25;              }}); // Decreased it to 0.25 (from 0.3)
 


// cut1: njet >= 2

  analysis_cuts['J'].push_back({ .name="2jet",   .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0;
            return 1;
          } });

  // cut2: jet 1 pass loose jet id
 analysis_cuts['J'].push_back({ .name="jet1_id",   .func = [&data](){
            // Define 'cut function here:
            if (data.jetsAK8Puppi.size<1) return 0; // for safety
            return data.jetsAK8Puppi.looseJetID[0]; 
          } });

  // cut3: jet 2 pass loose jet id
  analysis_cuts['J'].push_back({ .name="jet2_id",   .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            return data.jetsAK8Puppi.looseJetID[1];
          } });

  // cut4: jet 1 eta < 2.4
  analysis_cuts['J'].push_back({ .name="jet1_eta",   .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<1) return 0; // for safety
            if (fabs(data.jetsAK8Puppi.Eta[0])>=2.4) return 0;
            return 1;
          } });

  // cut5: jet 2 eta < 2.4
  analysis_cuts['J'].push_back({ .name="jet2_eta",   .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            if (fabs(data.jetsAK8Puppi.Eta[1])>=2.4) return 0;
            return 1;
          } });

  // cut6: jet 1 pt >= 400
  analysis_cuts['J'].push_back({ .name="jet1_pt",   .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<1) return 0; // for safety
            if (data.jetsAK8Puppi.Pt[0]<TOP_PT_CUT) return 0;
            return 1;
          } });

  // cut7: jet 2 pt >= 400
  analysis_cuts['J'].push_back({ .name="jet2_pt",   .func = [&data](){ 
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            if (data.jetsAK8Puppi.Pt[1]<TOP_PT_CUT) return 0;
            return 1;
          } });

  // cut8: 105 <= jet 1 mass (softdrop) < 210
  analysis_cuts['J'].push_back({ .name="jet1_mass", .func = [&data](){ 
            // Define cut function here:
            if (data.jetsAK8Puppi.size<1) return 0; // for safety
            if (data.jetsAK8Puppi.softDropMass[0]<TOP_SD_MASS_CUT_LOW) return 0;
            if (data.jetsAK8Puppi.softDropMass[0]>=TOP_SD_MASS_CUT_HIGH) return 0;
            return 1;
          } });

  // cut9: 105 <= jet 2 mass (softdrop) < 210
 analysis_cuts['J'].push_back({ .name="jet2_mass", .func = [&data](){ 
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            if (data.jetsAK8Puppi.softDropMass[1]<TOP_SD_MASS_CUT_LOW) return 0;
            if (data.jetsAK8Puppi.softDropMass[1]>=TOP_SD_MASS_CUT_HIGH) return 0;
            return 1;
          } });
 
  // cut10: | DeltaPhi | < DPHI_CUT
 analysis_cuts['J'].push_back({ .name="delta_phi", .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            TLorentzVector jet1, jet2;
            jet1.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[0], data.jetsAK8Puppi.Eta[0], data.jetsAK8Puppi.Phi[0], data.jetsAK8Puppi.E[0]);
            jet2.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[1], data.jetsAK8Puppi.Eta[1], data.jetsAK8Puppi.Phi[1], data.jetsAK8Puppi.E[1]);
            double dPhi = fabs(jet1.DeltaPhi(jet2)); 
            if (dPhi>=DPHI_CUT) return 0;            
            return 1;
          } });    

  // cut11: jet 1 tau32 < TOP_TAU32_CUT
 analysis_cuts['J'].push_back({ .name="jet1_tau32",   .func = [&data](){
            // Define cut function here:
            if (data.jetsAK8Puppi.size<1) return 0; // for safety
            double tau32 = data.jetsAK8Puppi.tau3[0];
            if (data.jetsAK8Puppi.tau2[0]!=0) tau32 /= data.jetsAK8Puppi.tau2[0];
            else tau32 = 9999;
            if (tau32>=TOP_TAU32_CUT) return 0;
            return 1;
          } });

  // cut12: jet 2 tau32 < TOP_TAU32_CUT
  analysis_cuts['J'].push_back({ .name="jet2_tau32",   .func = [&data](){ 
            // Define cut function here:
            if (data.jetsAK8Puppi.size<2) return 0; // for safety
            double tau32 = data.jetsAK8Puppi.tau3[1];
            if (data.jetsAK8Puppi.tau2[1]!=0) tau32 /= data.jetsAK8Puppi.tau2[1];
            else tau32 = 9999;
            if (tau32>=TOP_TAU32_CUT) return 0;
            return 1;
          } }); 

}//End of define_selections


//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& data) {
  return apply_all_cuts('S');
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

TH1D* h_R2;
TH1D* h_MR;
TH2D* h_R2_MR;

TH1D* h_R2_Scaled;
TH2D* h_R2_MR_Scaled;
TH1D* h_MR_Scaled;

TH1D* h_gen_toppt;

TH1D* h_NtopMult;

TH1D* h_nhadtop; 

TH1D* h_AK8Puppi_tau32;
TH1D* h_AK8Puppi_tau31;
TH1D* h_AK8Puppi_tau21;



std::vector<TH1D*> vh_jet1_pt;

Analysis::PostfixOptions
Analysis::get_pf_opts_(std::vector<std::vector<Sample> > lists, std::string dirname) {
  std::vector<Sample> samples;
  for (auto list : lists) samples.insert(samples.end(), list.begin(), list.end());
  PostfixOptions opt{ (size_t)-1, "", "", "" };
  for (size_t i=0; i<samples.size(); ++i) {
    // Find index of matching directory
    for (size_t j=0; j<samples[i].dirs.size(); ++j)
      if (samples[i].dirs[j] == dirname) opt.index = i;
    opt.postfixes += samples[i].postfix;
    opt.legends += samples[i].legend;
    opt.colors += samples[i].color;
    if (i+1!=samples.size()) {
      opt.postfixes +=  ";";
      opt.legends += ";";
      opt.colors += ",";
    }
  }
  return opt;
}

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  h_njet           = new TH1D("njet",         ";N_{jet}",                20, 0,  20);
  h_nw             = new TH1D("nw",           ";N_{W tag}",              20, 0,  20);
  h_nb             = new TH1D("nb",           ";N_{b tag}",              20, 0,  20);
  h_ht_gen         = new TH1D("ht_gen",       ";H_{T}^{gen}",            200, 0,2000);
  h_ht_AK4Puppi    = new TH1D("ht_AK4Puppi",  ";H_{T}",                  200, 0,2000);
  h_ht_AK8Puppi    = new TH1D("ht_AK8Puppi",  ";H_{T}^{AK8}",            200, 0,2000);
  h_jet1_pt        = new TH1D("jet1_pt",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt        = new TH1D("jet2_pt",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt        = new TH1D("jet3_pt",      ";p_{T, jet3}",            200, 0,2000);
  
  h_R2             = new TH1D("R2 ",";R2",100,0,1);
  h_MR             = new TH1D("MR ",";MR",100,0,5000);
  h_R2_MR          = new TH2D("R2_MR ",";MR" ,100,0,5000, 100,0,1);
 
  h_R2_Scaled      = new TH1D("R2_Scaled ",";R2",100,0,1);
  h_R2_MR_Scaled   = new TH2D("R2_MR_Scaled",";R2", 100,0,1 ,100,0,5000);
  h_MR_Scaled      = new TH1D("MR_Scaled ",";MR",100,0,5000);

  h_gen_toppt      = new TH1D("h_gen_toppt", "h_gen_toppt", 50, 0, 1000);
  h_NtopMult       = new TH1D("h_NtopMult", "h_NtopMult", 10, 0, 10);

  h_nhadtop        = new TH1D("nhadtop", ";N_{top tag}", 5, 0, 5);
  h_AK8Puppi_tau32 = new TH1D("tau32", "", 200,0,1);
  h_AK8Puppi_tau31 = new TH1D("tau31", "", 200,0,1);
  h_AK8Puppi_tau21 = new TH1D("tau21", "", 200,0,1);

  

  
 for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    histoname<<"jet1_pt_syst"<<i;
    title<<"Systematic variation #="<<i;
    vh_jet1_pt.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    vh_jet1_pt[i]->Sumw2();
  }

h_nhadtop->Sumw2(); 
h_AK8Puppi_tau32->Sumw2();
h_AK8Puppi_tau31->Sumw2();
h_AK8Puppi_tau21->Sumw2();

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


    h_R2->Fill(data.evt.R2, weight);
    h_MR->Fill(data.evt.MR, weight);
    h_R2_MR->Fill(data.evt.MR, data.evt.R2, weight);

    h_R2_Scaled->Fill(data.evt.R2);
    h_R2_MR_Scaled->Fill(data.evt.R2, data.evt.MR);
    h_MR_Scaled->Fill(data.evt.MR);

if (apply_all_cuts('J')) h_AK8Puppi_tau32->Fill(data.jetsAK8Puppi.tau3.at(0)/data.jetsAK8Puppi.tau2.at(0),weight); 
if (apply_all_cuts('J')) h_AK8Puppi_tau31->Fill(data.jetsAK8Puppi.tau3.at(0)/data.jetsAK8Puppi.tau1.at(0),weight); 
if (apply_all_cuts('J')) h_AK8Puppi_tau21->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),weight); 
if (apply_all_cuts('J')) h_nhadtop->Fill(nLooseIDHadTopTagJets, weight);



  if (h_R2_Scaled->Integral()!=0)
  h_R2_Scaled->Scale(1/h_R2_Scaled->Integral());   

  if (h_R2_MR_Scaled->Integral()!=0)
   h_R2_MR_Scaled->Scale(1/h_R2_MR_Scaled->Integral());

  if (h_MR_Scaled->Integral()!=0)
   h_MR_Scaled->Scale(1/ h_MR_Scaled->Integral());

   

    // For example this applies the first three cuts in signal region
    // ele/mu veto


//GenSize 

int NtopMult=0;
 
for (unsigned int i=0; i<data.gen.size; i++) 
{
  

   // if (data.gen.Status[i] != 3) continue;
    if (fabs(data.gen.ID[i]) == 6) 
    {

  
        if ((fabs(data.gen.Dau1ID[i]) == 5 && fabs(data.gen.Dau0ID[i]) == 24) ||
           (fabs(data.gen.Dau1ID[i]) == 24 && fabs(data.gen.Dau0ID[i]) == 5))
       { 
  
          h_gen_toppt->Fill(data.gen.Pt[i]);

          NtopMult++;
  
       }
    }
    

}

        h_NtopMult->Fill(NtopMult);

     
    if (apply_ncut('S', 2)) {
      h_jet1_pt->Fill(data.jetsAK4Puppi.Pt[iJet[0]], weight);
      h_jet2_pt->Fill(data.jetsAK4Puppi.Pt[iJet[1]], weight);
      h_jet3_pt->Fill(data.jetsAK4Puppi.Pt[iJet[2]], weight);

    }
    /* 
      Other examples to use analysis_cuts object

      if (apply_cut("S","1W"))                          --> 1 Cut from S region
      if (apply_cut("W","1Wpre"))                       --> 1 Cut from W region
      if (apply_all_cuts("T"))                          --> All cuts in T region
      if (apply_all_cuts_except("Q", "mDPhi<0.25"))     --> N-1 cut
      if (apply_all_cuts_except("S", {"0Ele", "0Mu" })) --> S without Lep veto

      But be aware: Whatever is defined in the baseline_cuts will apply to all histograms
      Also if you use skimmed ntuples (very likely) then those cuts are already applied
      This is because unskimmed ntuple is 4.3 TB in size, and we cannot have them on EOS
    */
  }
  
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  if (apply_all_cuts('S')) vh_jet1_pt[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[0]], weight);
}

////>>>>>>>>>>>>>>>>>> Methods used by SmartHistos (Plotter)>>>>>>>>>>>>>>>>>>

void
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst, 
             const unsigned int& syst_index, std::string dirname, bool runOnSkim=false)
{

  std::vector<Sample> signal_all, signal_selected, signal_fastsim, signal_gluino, signal_stop;
  signal_all.push_back({ .postfix="T5ttcc",       .legend="T5ttcc",      .color="12", /*DGrey*/ .dirs={ "FastSim_SMS-T5ttcc" } });
  signal_all.push_back({ .postfix="T5tttt",       .legend="T5tttt",      .color="862",/*Azure*/ .dirs={ "FastSim_SMS-T5tttt" } });
  signal_all.push_back({ .postfix="T1tttt",       .legend="T1tttt",      .color="841",/*Teal*/  .dirs={ "FastSim_SMS-T1tttt" } });
  signal_all.push_back({ .postfix="T2tt",         .legend="T2tt",        .color="403",/*DYell*/ .dirs={ 
         "FastSim_SMS-T2tt_mStop-150to250", "FastSim_SMS-T2tt_mStop-250to350",
         "FastSim_SMS-T2tt_mStop-350to400", "FastSim_SMS-T2tt_mStop-400to1200" 
       } });
  signal_all.push_back({ .postfix="T2tt_FullSim", .legend="T2tt (FullSim)", .color="804",/*DOran*/ .dirs={
         "FullSim_SMS-T2tt_mStop-425_mLSP-325", "FullSim_SMS-T2tt_mStop-500_mLSP-325",
         "FullSim_SMS-T2tt_mStop-850_mLSP-100" 
       } });
  signal_selected.push_back(signal_all[0]);
  for (int i=0; i<4; ++i) signal_fastsim.push_back(signal_all[i]);
  for (int i=0; i<3; ++i) signal_gluino .push_back(signal_all[i]);
  for (int i=3; i<5; ++i) signal_stop .push_back(signal_all[i]);



static const PostfixOptions signals_background_opt = get_pf_opts_({signal_all}, dirname);
  sh.AddNewPostfix("Signals,Background",  [&d] { 
         // Select gluino/stop mass to give ~1k events with 40 fb^-1
         if (signals_background_opt.index==1) {
           if (d.evt.SUSY_Gluino_Mass == 1200 && d.evt.SUSY_LSP_Mass == 200) return (size_t)-1; // T5tttt
         } else if (signals_background_opt.index==3) {
           if (d.evt.SUSY_Stop_Mass  == 800 && d.evt.SUSY_LSP_Mass == 100) return (size_t)-1; // T2tt - Same as FullSim point
         }
         return signals_background_opt.index; 
}, signals_background_opt.postfixes, signals_background_opt.legends, signals_background_opt.colors);

}

void
Analysis::load_analysis_histos(std::string inputfile)
{
}

void
Analysis::save_analysis_histos(bool draw=0)
{

}

