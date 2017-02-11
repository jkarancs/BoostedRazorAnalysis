#include "TLorentzVector.h"
#include "TMath.h"
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

  int mGluino[581] ={800,800,800,800,800,800,800,850,900,900,900,900,900,900,900,900,950,950,950,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1050,1050,1050,1050,1050,1100,1100,1100,1100,1100,1100,1100,1100,1100,1100,1100,1100,1150,1150,1150,1150,1150,1150,1150,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1200,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1250,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1300,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1350,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1400,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1450,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1550,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1650,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1700,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1750,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1850,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1900,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,1950,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2050,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2100,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2150,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2200,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2250,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300,2300};           

  int mLSP[581] = {0,100,200,300,400,500,525,
                   575,
                   0,100,200,300,400,500,600,625,
                   600,650,675,
                   0,100,200,300,400,500,600,650,700,725,
                   600,650,700,750,775,
                   0,100,200,300,400,500,600,650,700,750,800,825,
                   600,650,700,750,800,850,875,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,925,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,975,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1025,    
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1075,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1100,1125,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1175,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1225,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1275,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1325,
                   0,50,100,150,200,300,400,500,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1375,
                   0,50,100,150,200,300,400,500,600,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1425,
                   0,50,100,150,200,300,400,500,600,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1475,
                   0,50,100,150,200,300,400,500,600,700,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1525,
                   0,50,100,150,200,300,400,500,600,700,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1575,
                   0,50,100,150,200,300,400,500,600,700,800,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1100,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1100,1200,1250,1300,1350,1400,1450,1500,1550,1600,
                   0,50,100,150,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1350,1400,1450,1500,1550,1600};


//_______________________________________________________
//          Define Analysis specific weights
double
Analysis::get_analysis_weight(DataStruct& data)
{
  double w = 1;
  //w *=(TMath::Erf((AK4Puppi_Ht-784.036)/83.7739)/2.+0.5); //This isn't worked. It has a problem.
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

  // Signal skim
  //return apply_all_cuts('S');
}

Double_t* getVariableBinEdges(int num_entries, Double_t* tmp_array)
{ 
  Double_t* my_array = new Double_t[num_entries];
  for (int i = 0; i != num_entries; ++i) {
    my_array[i] = tmp_array[i];
    
    cout << "bin edge " << i << " : " << my_array[i] << endl;
  }
   return my_array;
}


//_______________________________________________________
//          Define Analysis event selection cuts
//     Can define all sorts of Signal/Control regions

void
Analysis::define_selections(const DataStruct& d)
{
  analysis_cuts.clear();

  // Define here cuts that are common in all Signal/Control regions
  // MET Filters, etc. are already applied in AnalysisBase.h, See baseline_cuts

/*
    // cut0: signal mass region
  baseline_cuts.push_back({ .name="signal_mass_selection",   .func = [&d]{ 
            int num = 1;
            return d.evt.SUSY_Gluino_Mass==mGluino[num] && d.evt.SUSY_LSP_Mass==mLSP[num];
          } });
*/

  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                    }}); // Similar to pt>200, one AK8 jet has pt>170
  baseline_cuts.push_back({ .name="Skim_2Jet",       .func = []    { return nJet>=2;                       }});
  baseline_cuts.push_back({ .name="Baseline_3Jet",   .func = []    { return nJet>=3;                       }}); // Separate cut, so one can exclude (N-1)
  baseline_cuts.push_back({ .name="Baseline_MR_R2",  .func = [&d]  { return d.evt.MR>800 && d.evt.R2>0.08; }});
  //baseline_cuts.push_back({ .name="Trigger",   .func = []    { return passTrigger==1;}});
//temporaray
  //baseline_cuts.push_back({ .name="One Electron",       .func = []    { return nEleTight>=1;                  }});
  //baseline_cuts.push_back({ .name="Electron Trigger",   .func = [&d]    { return d.hlt.Ele27_WPTight_Gsf==1;}});
  //baseline_cuts.push_back({ .name="One Muon",       .func = []    { return nMuTight>=1;                  }});
  //baseline_cuts.push_back({ .name="Muon Trigger",   .func = [&d]    { return d.hlt.IsoMu24==1;        }});

  // S: Signal region
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                   }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                    }});
  //analysis_cuts['S'].push_back({ .name="0TauTrk",    .func = []    { return;  }});
  analysis_cuts['S'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                }});
  analysis_cuts['S'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                 }});
  //analysis_cuts['S'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['S'].push_back({ .name="mDPhi>=0p4", .func = []    { return minDeltaPhi>=0.4;              }}); // Decreased it to the AK4 cone size (from 0.5)
  
  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepTight==1;                  }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                 }});
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
  analysis_cuts['Q'].push_back({ .name="mDPhi<0p25", .func = []    { return minDeltaPhi<0.25;              }}); // Decreased it to 0.25 (from 0.3)
  
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& data) {
  return apply_all_cuts('S');
  //return 1;
}

//_______________________________________________________
//                 List of Histograms
TH1D* h_njet;
TH1D* h_nb;
TH1D* h_nw;
TH1D* h_ht_gen;
TH1D* h_ht_AK4Puppi;
TH1D* h_ht_AK8Puppi;

TH1D* h_ht_AK4Puppi_S;
TH1D* h_ht_AK8Puppi_S;
TH1D* h_jet1_pt_S;
TH1D* h_jet2_pt_S;
TH1D* h_jet3_pt_S;
TH1D* h_MR_S;
TH1D* h_MTR_S;
TH1D* h_R_S;
TH1D* h_R2_S;
TH1D* h_tau21_S;
TH1D* h_MET_S;

TH1D* h_ht_AK4Puppi_W;
TH1D* h_ht_AK8Puppi_W;
TH1D* h_jet1_pt_W;
TH1D* h_jet2_pt_W;
TH1D* h_jet3_pt_W;
TH1D* h_MR_W;
TH1D* h_MTR_W;
TH1D* h_R_W;
TH1D* h_R2_W;
TH1D* h_tau21_W;
TH1D* h_MET_W;

TH1D* h_ht_AK4Puppi_T;
TH1D* h_ht_AK8Puppi_T;
TH1D* h_jet1_pt_T;
TH1D* h_jet2_pt_T;
TH1D* h_jet3_pt_T;
TH1D* h_MR_T;
TH1D* h_MTR_T;
TH1D* h_R_T;
TH1D* h_R2_T;
TH1D* h_tau21_T;
TH1D* h_MET_T;

TH1D* h_ht_AK4Puppi_Q;
TH1D* h_ht_AK8Puppi_Q;
TH1D* h_jet1_pt_Q;
TH1D* h_jet2_pt_Q;
TH1D* h_jet3_pt_Q;
TH1D* h_MR_Q;
TH1D* h_MTR_Q;
TH1D* h_R_Q;
TH1D* h_R2_Q;
TH1D* h_tau21_Q;
TH1D* h_MET_Q;


TH1D* h_softDropMass;
TH1D* h_StopMass;
TH1D* h_GluinoMass;
TH1D* h_LSPMass;
TH1D* h_minDeltaPhi;
TH2D* h_R2_MR;
TH2D* h_GluinoLSPMass;

//std::vector<TH1D*> vh_jet1_pt;

TH1D* h_njet_pre;
TH1D* h_nb_pre;
TH1D* h_nw_pre;
TH1D* h_j2_pt_pre;
TH1D* h_j3_pt_pre;
TH1D* h_MR_pre;
TH1D* h_R2_pre;
TH1D* h_tau21_pre;
TH1D* h_MET_pre;
TH1D* h_softDropMass_pre;

TH1D* h_njet_pre_pass;
TH1D* h_nb_pre_pass;
TH1D* h_nw_pre_pass;
TH1D* h_j2_pt_pre_pass;
TH1D* h_j3_pt_pre_pass;
TH1D* h_MR_pre_pass;
TH1D* h_R2_pre_pass;
TH1D* h_tau21_pre_pass;
TH1D* h_MET_pre_pass;
TH1D* h_softDropMass_pre_pass;

TH1D* h_njet_pre_passj1pt;
TH1D* h_nb_pre_passj1pt;
TH1D* h_nw_pre_passj1pt;
TH1D* h_j2_pt_pre_passj1pt;
TH1D* h_j3_pt_pre_passj1pt;
TH1D* h_MR_pre_passj1pt;
TH1D* h_R2_pre_passj1pt;
TH1D* h_tau21_pre_passj1pt;
TH1D* h_MET_pre_passj1pt;
TH1D* h_softDropMass_pre_passj1pt;

TH1D* h_njet_pre_passHT;
TH1D* h_nb_pre_passHT;
TH1D* h_nw_pre_passHT;
TH1D* h_j2_pt_pre_passHT;
TH1D* h_j3_pt_pre_passHT;
TH1D* h_MR_pre_passHT;
TH1D* h_R2_pre_passHT;
TH1D* h_tau21_pre_passHT;
TH1D* h_MET_pre_passHT;
TH1D* h_softDropMass_pre_passHT;

TH1D *h_HT_pre;
TH1D *h_HT_pre_pass;
TH1D *h_HT_pre_passj1pt;
TH1D *h_HT_pre_passHT;
TH1D *h_HT_Q;
TH1D *h_HT_Q_pass;
TH1D *h_HT_Q_passj1pt;
TH1D *h_HT_Q_passHT;
TH1D *h_HT_S;
TH1D *h_HT_S_pass;
TH1D *h_HT_S_passj1pt;
TH1D *h_HT_S_passHT;

TH1D *h_j1_pt_pre;
TH1D *h_j1_pt_pre_pass;
TH1D *h_j1_pt_pre_passHT;
TH1D *h_j1_pt_pre_passj1pt;
TH1D *h_j1_pt_Q;
TH1D *h_j1_pt_Q_pass;
TH1D *h_j1_pt_Q_passHT;
TH1D *h_j1_pt_Q_passj1pt;
TH1D *h_j1_pt_S;
TH1D *h_j1_pt_S_pass;
TH1D *h_j1_pt_S_passHT;
TH1D *h_j1_pt_S_passj1pt;

TH2D *h_HT_j1pt_pre;
TH2D *h_HT_j1pt_pre_pass;
TH2D *h_HT_j1pt_pre_passHT;
TH2D *h_HT_j1pt_pre_passj1pt;
TH2D *h_HT_j1pt_Q;
TH2D *h_HT_j1pt_Q_pass;
TH2D *h_HT_j1pt_Q_passHT;
TH2D *h_HT_j1pt_Q_passj1pt;
TH2D *h_HT_j1pt_S;
TH2D *h_HT_j1pt_S_pass;
TH2D *h_HT_j1pt_S_passHT ;
TH2D *h_HT_j1pt_S_passj1pt;
TH2D *h_HT_j1pt_T;
TH2D *h_HT_j1pt_T_pass;
TH2D *h_HT_j1pt_T_passHT;
TH2D *h_HT_j1pt_T_passj1pt;
TH2D *h_HT_j1pt_W;
TH2D *h_HT_j1pt_W_pass;
TH2D *h_HT_j1pt_W_passHT;
TH2D *h_HT_j1pt_W_passj1pt;

TH2D* h_HT_weight;
//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  h_HT_weight = new TH2D("HT_weight","",40,0,2000,40,0,2);
  h_njet         = new TH1D("njet",         ";N_{jet}",                20, 0,  20);
  h_nw           = new TH1D("nw",           ";N_{W tag}",              20, 0,  20);
  h_nb           = new TH1D("nb",           ";N_{b tag}",              20, 0,  20);
  h_ht_gen       = new TH1D("ht_gen",       ";H_{T}^{gen}",            300, 0,3000);

  h_ht_AK4Puppi  = new TH1D("ht_AK4Puppi",  ";H_{T}",                  300, 0,3000);
  h_ht_AK8Puppi  = new TH1D("ht_AK8Puppi",  ";H_{T}^{AK8}",            300, 0,3000);

  h_ht_AK4Puppi_S = new TH1D("ht_AK4Puppi_S",  ";H_{T}",                  300, 0,3000);
  h_ht_AK8Puppi_S = new TH1D("ht_AK8Puppi_S",  ";H_{T}^{AK8}",            300, 0,3000);
  h_jet1_pt_S = new TH1D("jet1_pt_S",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_S = new TH1D("jet2_pt_S",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_S = new TH1D("jet3_pt_S",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_S = new TH1D("MR_S",   ";MR_{AK4}",         200, 0,4000);
  h_MTR_S = new TH1D("MTR_S",  ";MTR_{AK4}",        200, 0,2000);
  h_R_S = new TH1D("R_S",    ";R_{AK4}",          500, 0,1);
  h_R2_S = new TH1D("R2_S",   ";R2_{AK4}",         500, 0,1);
  h_tau21_S = new TH1D("tau21_S", ";tau21", 200,0,1);
  h_MET_S = new TH1D("MET_S", ";MET", 400,0,2000);

  h_ht_AK4Puppi_W = new TH1D("ht_AK4Puppi_W",  ";H_{T}",                  300, 0,3000);
  h_ht_AK8Puppi_W = new TH1D("ht_AK8Puppi_W",  ";H_{T}^{AK8}",            300, 0,3000);
  h_jet1_pt_W = new TH1D("jet1_pt_W",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_W = new TH1D("jet2_pt_W",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_W = new TH1D("jet3_pt_W",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_W = new TH1D("MR_W",   ";MR_{AK4}",         200, 0,4000);
  h_MTR_W = new TH1D("MTR_W",  ";MTR_{AK4}",        200, 0,2000);
  h_R_W = new TH1D("R_W",    ";R_{AK4}",          500, 0,1);
  h_R2_W = new TH1D("R2_W",   ";R2_{AK4}",         500, 0,1);
  h_tau21_W = new TH1D("tau21_W", ";tau21", 200,0,1);
  h_MET_W = new TH1D("MET_W", ";MET", 400,0,2000);

  h_ht_AK4Puppi_T = new TH1D("ht_AK4Puppi_T",  ";H_{T}",                  300, 0,3000);
  h_ht_AK8Puppi_T = new TH1D("ht_AK8Puppi_T",  ";H_{T}^{AK8}",            300, 0,3000);
  h_jet1_pt_T = new TH1D("jet1_pt_T",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_T = new TH1D("jet2_pt_T",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_T = new TH1D("jet3_pt_T",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_T = new TH1D("MR_T",   ";MR_{AK4}",         200, 0,4000);
  h_MTR_T = new TH1D("MTR_T",  ";MTR_{AK4}",        200, 0,2000);
  h_R_T = new TH1D("R_T",    ";R_{AK4}",          500, 0,1);
  h_R2_T = new TH1D("R2_T",   ";R2_{AK4}",         500, 0,1);
  h_tau21_T = new TH1D("tau21_T", ";tau21", 200,0,1);
  h_MET_T = new TH1D("MET_T", ";MET", 400,0,2000);

  h_ht_AK4Puppi_Q = new TH1D("ht_AK4Puppi_Q",  ";H_{T}",                  300, 0,3000);
  h_ht_AK8Puppi_Q = new TH1D("ht_AK8Puppi_Q",  ";H_{T}^{AK8}",            300, 0,3000);
  h_jet1_pt_Q = new TH1D("jet1_pt_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Q = new TH1D("jet2_pt_Q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Q = new TH1D("jet3_pt_Q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Q = new TH1D("MR_Q",   ";MR_{AK4}",         200, 0,4000);
  h_MTR_Q = new TH1D("MTR_Q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Q = new TH1D("R_Q",    ";R_{AK4}",          500, 0,1);
  h_R2_Q = new TH1D("R2_Q",   ";R2_{AK4}",         500, 0,1);
  h_tau21_Q = new TH1D("tau21_Q", ";tau21", 200,0,1);
  h_MET_Q = new TH1D("MET_Q", ";MET", 400,0,2000);


  h_minDeltaPhi = new TH1D("minDeltaPhi", ";min_#Delta#phi(AK4_{1~3}, MET)", 314,0,3.14);
  h_softDropMass = new TH1D("softDropMass", "", 100,0,500);
  h_GluinoLSPMass = new TH2D("GluinoLSPMass","",34,600,2300,64,0,1600);
  h_StopMass = new TH1D("StopMass", "", 92,0,2300);
  h_GluinoMass = new TH1D("GluinoMass", "", 34,600,2300);
  h_LSPMass = new TH1D("LSPMass", "", 64,0,1600);
  h_R2_MR = new TH2D("R2_MR", ";MR_{AK4};R2_{AK4}", 20,0,4000,10,0,1);
  
/*
  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    histoname<<"jet1_pt_syst"<<i;
    title<<"Systematic variation #="<<i;
    vh_jet1_pt.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
  }
*/

  // Histograms:

  double htbn = 20;
  double htmn = 0;
  double htmx = 1500;

  double j1ptbn = 20;
  double j1ptmn = 0;
  double j1ptmx = 1000;

  //double MRbn = 7;
  //double MRmn = 0;
  //double MRmx = 4000;
  //double Rbn = 7;
  //double Rmn = 0;
  //double Rmx = 1;

  // Variable binning

  int nbn_HT = 19;
  int nbn_j1pt = 12;
  int nbn_MR = 6;
  int nbn_R = 6;
  Double_t bn_HT_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,900.,1000.,2500.};
  Double_t* bn_HT = 0;
  bn_HT = getVariableBinEdges(nbn_HT+1,bn_HT_tmp);
  Double_t bn_j1pt_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,700.,1000.};
  Double_t* bn_j1pt = 0;
  bn_j1pt = getVariableBinEdges(nbn_j1pt+1,bn_j1pt_tmp);
  Double_t bn_MR_tmp[] = {0.,800.,1000.,1200.,1600.,2000.,4000.};
  Double_t* bn_MR = 0;
  bn_MR = getVariableBinEdges(nbn_MR+1,bn_MR_tmp);
  Double_t bn_R_tmp[] = {0.,0.08,.12,.16,0.24,0.5,1.};
  Double_t* bn_R = 0;
  bn_R = getVariableBinEdges(nbn_R+1,bn_R_tmp);
  

  h_njet_pre = new TH1D("h_njet_pre",         ";N_{jet}",                20, 0,  20);
  h_nw_pre = new TH1D("h_nw_pre",           ";N_{W tag}",              20, 0,  20);
  h_nb_pre = new TH1D("h_nb_pre",           ";N_{b tag}",              20, 0,  20);
  h_j2_pt_pre = new TH1D("h_j2_pt_pre",      ";p_{T, jet2}",            j1ptbn,j1ptmn,j1ptmx);
  h_j3_pt_pre = new TH1D("h_j3_pt_pre",      ";p_{T, jet3}",            j1ptbn,j1ptmn,j1ptmx);
  h_MR_pre = new TH1D("h_MR_pre",   ";MR_{AK4}",        nbn_MR,bn_MR);
  h_R2_pre = new TH1D("h_R2_pre",    ";R2_{AK4}",          nbn_R,bn_R);
  h_tau21_pre = new TH1D("h_tau21_pre", ";tau21", 20,0,1);
  h_MET_pre = new TH1D("h_MET_pre", ";MET", j1ptbn,j1ptmn,j1ptmx);
  h_softDropMass_pre = new TH1D("h_softDropMass_pre", "", 20,0,500);

  h_njet_pre_pass = new TH1D("h_njet_pre_pass",         ";N_{jet}",                20, 0,  20);
  h_nw_pre_pass = new TH1D("h_nw_pre_pass",           ";N_{W tag}",              20, 0,  20);
  h_nb_pre_pass = new TH1D("h_nb_pre_pass",           ";N_{b tag}",              20, 0,  20);
  h_j2_pt_pre_pass = new TH1D("h_j2_pt_pre_pass",      ";p_{T, jet2}",            j1ptbn,j1ptmn,j1ptmx);
  h_j3_pt_pre_pass = new TH1D("h_j3_pt_pre_pass",      ";p_{T, jet3}",            j1ptbn,j1ptmn,j1ptmx);
  h_MR_pre_pass = new TH1D("h_MR_pre_pass",   ";MR_{AK4}",        nbn_MR,bn_MR);
  h_R2_pre_pass = new TH1D("h_R2_pre_pass",    ";R2_{AK4}",          nbn_R,bn_R);
  h_tau21_pre_pass = new TH1D("h_tau21_pre_pass", ";tau21", 20,0,1);
  h_MET_pre_pass = new TH1D("h_MET_pre_pass", ";MET", j1ptbn,j1ptmn,j1ptmx);
  h_softDropMass_pre_pass = new TH1D("h_softDropMass_pre_pass", "", 20,0,500);

  h_njet_pre_passj1pt = new TH1D("h_njet_pre_passj1pt",         ";N_{jet}",                20, 0,  20);
  h_nw_pre_passj1pt = new TH1D("h_nw_pre_passj1pt",           ";N_{W tag}",              20, 0,  20);
  h_nb_pre_passj1pt = new TH1D("h_nb_pre_passj1pt",           ";N_{b tag}",              20, 0,  20);
  h_j2_pt_pre_passj1pt = new TH1D("h_j2_pt_pre_passj1pt",      ";p_{T, jet2}",            j1ptbn,j1ptmn,j1ptmx);
  h_j3_pt_pre_passj1pt = new TH1D("h_j3_pt_pre_passj1pt",      ";p_{T, jet3}",            j1ptbn,j1ptmn,j1ptmx);
  h_MR_pre_passj1pt = new TH1D("h_MR_pre_passj1pt",   ";MR_{AK4}",        nbn_MR,bn_MR);
  h_R2_pre_passj1pt = new TH1D("h_R2_pre_passj1pt",    ";R2_{AK4}",          nbn_R,bn_R);
  h_tau21_pre_passj1pt = new TH1D("h_tau21_pre_passj1pt", ";tau21", 20,0,1);
  h_MET_pre_passj1pt = new TH1D("h_MET_pre_passj1pt", ";MET", j1ptbn,j1ptmn,j1ptmx);
  h_softDropMass_pre_passj1pt = new TH1D("h_softDropMass_pre_passj1pt", "", 20,0,500);

  h_njet_pre_passHT = new TH1D("h_njet_pre_passHT",         ";N_{jet}",                20, 0,  20);
  h_nw_pre_passHT = new TH1D("h_nw_pre_passHT",           ";N_{W tag}",              20, 0,  20);
  h_nb_pre_passHT = new TH1D("h_nb_pre_passHT",           ";N_{b tag}",              20, 0,  20);
  h_j2_pt_pre_passHT = new TH1D("h_j2_pt_pre_passHT",      ";p_{T, jet2}",            j1ptbn,j1ptmn,j1ptmx);
  h_j3_pt_pre_passHT = new TH1D("h_j3_pt_pre_passHT",      ";p_{T, jet3}",            j1ptbn,j1ptmn,j1ptmx);
  h_MR_pre_passHT = new TH1D("h_MR_pre_passHT",   ";MR_{AK4}",        nbn_MR,bn_MR);
  h_R2_pre_passHT = new TH1D("h_R2_pre_passHT",    ";R2_{AK4}",          nbn_R,bn_R);
  h_tau21_pre_passHT = new TH1D("h_tau21_pre_passHT", ";tau21", 20,0,1);
  h_MET_pre_passHT = new TH1D("h_MET_pre_passHT", ";MET", j1ptbn,j1ptmn,j1ptmx);
  h_softDropMass_pre_passHT = new TH1D("h_softDropMass_pre_passHT", "", 20,0,500);

  // HT
  h_HT_pre = new TH1D("h_HT_pre", "h_HT_pre", htbn, htmn, htmx);
  h_HT_pre_pass = new TH1D("h_HT_pre_pass", "h_HT_pre_pass", htbn, htmn, htmx);
  h_HT_pre_passj1pt = new TH1D("h_HT_pre_passj1pt", "h_HT_pre_passj1pt", htbn, htmn, htmx);
  h_HT_pre_passHT = new TH1D("h_HT_pre_passHT", "h_HT_pre_passHT", htbn, htmn, htmx);

  h_HT_Q = new TH1D("h_HT_Q", "h_HT_Q", htbn, htmn, htmx);
  h_HT_Q_pass = new TH1D("h_HT_Q_pass", "h_HT_Q_pass", htbn, htmn, htmx);
  h_HT_Q_passj1pt = new TH1D("h_HT_Q_passj1pt", "h_HT_Q_passj1pt", htbn, htmn, htmx);
  h_HT_Q_passHT = new TH1D("h_HT_Q_passHT", "h_HT_Q_passHT", htbn, htmn, htmx);

  h_HT_S = new TH1D("h_HT_S", "h_HT_S", htbn, htmn, htmx);
  h_HT_S_pass = new TH1D("h_HT_S_pass", "h_HT_S_pass", htbn, htmn, htmx);
  h_HT_S_passj1pt = new TH1D("h_HT_S_passj1pt", "h_HT_S_passj1pt", htbn, htmn, htmx);
  h_HT_S_passHT = new TH1D("h_HT_S_passHT", "h_HT_S_passHT", htbn, htmn, htmx);


  h_j1_pt_pre = new TH1D("h_j1_pt_pre", ";p_{T, jet1}_AK4_pre", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_pre_pass = new TH1D("h_j1_pt_pre_pass", ";p_{T, jet1}_AK4_pre_pass", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_pre_passHT = new TH1D("h_j1_pt_pre_passHT", ";p_{T, jet1}_AK4_pre_passHT", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_pre_passj1pt = new TH1D("h_j1_pt_pre_passj1pt", ";p_{T, jet1}_AK4_pre_passj1pt", j1ptbn, j1ptmn, j1ptmx);

  h_j1_pt_Q = new TH1D("h_j1_pt_Q", ";p_{T, jet1}_AK4_Q", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_Q_pass = new TH1D("h_j1_pt_Q_pass", ";p_{T, jet1}_AK4_Q_pass", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_Q_passHT = new TH1D("h_j1_pt_Q_passHT", ";p_{T, jet1}_AK4_Q_passHT", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_Q_passj1pt = new TH1D("h_j1_pt_Q_passj1pt", ";p_{T, jet1}_AK4_Q_passj1pt", j1ptbn, j1ptmn, j1ptmx);

  h_j1_pt_S = new TH1D("h_j1_pt_S", ";p_{T, jet1}_AK4_S", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_S_pass = new TH1D("h_j1_pt_S_pass", ";p_{T, jet1}_AK4_S_pass", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_S_passHT = new TH1D("h_j1_pt_S_passHT", ";p_{T, jet1}_AK4_S_passHT", j1ptbn, j1ptmn, j1ptmx);
  h_j1_pt_S_passj1pt = new TH1D("h_j1_pt_S_passj1pt", ";p_{T, jet1}_AK4_S_passj1pt", j1ptbn, j1ptmn, j1ptmx);


  h_HT_j1pt_pre = new TH2D("h_HT_j1pt_pre", "h_HT_j1pt_pre", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_pre_pass = new TH2D("h_HT_j1pt_pre_pass", "h_HT_j1pt_pre_pass", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_pre_passHT = new TH2D("h_HT_j1pt_pre_passHT", "h_HT_j1pt_pre_passHT", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_pre_passj1pt = new TH2D("h_HT_j1pt_pre_passj1pt", "h_HT_j1pt_pre_passj1pt", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

  h_HT_j1pt_Q = new TH2D("h_HT_j1pt_Q", "h_HT_j1pt_Q", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_Q_pass = new TH2D("h_HT_j1pt_Q_pass", "h_HT_j1pt_Q_pass", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_Q_passHT = new TH2D("h_HT_j1pt_Q_passHT", "h_HT_j1pt_Q_passHT", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_Q_passj1pt = new TH2D("h_HT_j1pt_Q_passj1pt", "h_HT_j1pt_Q_passj1pt", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

  h_HT_j1pt_S = new TH2D("h_HT_j1pt_S", "h_HT_j1pt_S", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_S_pass = new TH2D("h_HT_j1pt_S_pass", "h_HT_j1pt_S_pass", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_S_passHT = new TH2D("h_HT_j1pt_S_passHT", "h_HT_j1pt_S_passHT", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_S_passj1pt = new TH2D("h_HT_j1pt_S_passj1pt", "h_HT_j1pt_S_passj1pt", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

  h_HT_j1pt_T = new TH2D("h_HT_j1pt_T", "h_HT_j1pt_T", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_T_pass = new TH2D("h_HT_j1pt_T_pass", "h_HT_j1pt_T_pass", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_T_passHT = new TH2D("h_HT_j1pt_T_passHT", "h_HT_j1pt_T_passHT", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_T_passj1pt = new TH2D("h_HT_j1pt_T_passj1pt", "h_HT_j1pt_T_passj1pt", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

  h_HT_j1pt_W = new TH2D("h_HT_j1pt_W", "h_HT_j1pt_W", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_W_pass = new TH2D("h_HT_j1pt_W_pass", "h_HT_j1pt_W_pass", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_W_passHT = new TH2D("h_HT_j1pt_W_passHT", "h_HT_j1pt_W_passHT", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);
  h_HT_j1pt_W_passj1pt = new TH2D("h_HT_j1pt_W_passj1pt", "h_HT_j1pt_W_passj1pt", nbn_HT, bn_HT, nbn_j1pt, bn_j1pt);

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

   bool pass = data.hlt.AK8PFJet360_TrimMass30+data.hlt.PFHT800 == 0 ? false : true;

Analysis ana;
    h_HT_weight->Fill(AK4Puppi_Ht,ana.get_analysis_weight(data),weight);
    float reweight=1;
  //if(AK4Puppi_Ht > 1080) cout << "HT : " << AK4Puppi_Ht << " weight : " << w << endl;
   // if(AK4Puppi_Ht > 1080) cout << "HT : " << AK4Puppi_Ht << " weight : " << ana.get_analysis_weight(data) << endl;
    //reweight = weight*(TMath::Erf((AK4Puppi_Ht-784.036)/83.7739)/2.+0.5);
    //reweight = weight;

    if(AK4Puppi_Ht<450) reweight=weight*0;
    if(AK4Puppi_Ht>=450 && AK4Puppi_Ht<525) reweight=weight*0.0132325;
    if(AK4Puppi_Ht>=525 && AK4Puppi_Ht<600) reweight=weight*0.0304311;
    if(AK4Puppi_Ht>=600 && AK4Puppi_Ht<675) reweight=weight*0.0544904;
    if(AK4Puppi_Ht>=675 && AK4Puppi_Ht<750) reweight=weight*0.1318;
    if(AK4Puppi_Ht>=750 && AK4Puppi_Ht<825) reweight=weight*0.477929;
    if(AK4Puppi_Ht>=825 && AK4Puppi_Ht<900) reweight=weight*0.920323;
    if(AK4Puppi_Ht>=900 && AK4Puppi_Ht<975) reweight=weight*0.99669;
    if(AK4Puppi_Ht>=975 && AK4Puppi_Ht<1050) reweight=weight*0.998733;
    if(AK4Puppi_Ht>=1050&& AK4Puppi_Ht<1125) reweight=weight*0.999553;
    if(AK4Puppi_Ht>=1125&& AK4Puppi_Ht<1200) reweight=weight*1;
    if(AK4Puppi_Ht>=1200&& AK4Puppi_Ht<1275) reweight=weight*0.997535;
    if(AK4Puppi_Ht>=1275&& AK4Puppi_Ht<1350) reweight=weight*0.998922;
    if(AK4Puppi_Ht>=1350) reweight=weight*1.;
    
  
    h_njet   ->Fill(nJet,        reweight);
    h_nb     ->Fill(nMediumBTag, reweight);
    h_nw     ->Fill(nTightWTag,  reweight);
    h_ht_gen->Fill(data.evt.Gen_Ht,  reweight);  // in ntuple
    h_ht_AK4Puppi->Fill(AK4Puppi_Ht, reweight); // Calculated in AnalysisBase.h
    h_ht_AK8Puppi->Fill(AK8Puppi_Ht, reweight); // Calculated in AnalysisBase.h
    
    if (apply_all_cuts('W')) {
      h_ht_AK4Puppi_W->Fill(AK4Puppi_Ht, reweight);
      h_ht_AK8Puppi_W->Fill(AK8Puppi_Ht, reweight);
      h_jet1_pt_W->Fill(data.jetsAK4Puppi.Pt[iJet[0]], reweight);
      h_jet2_pt_W->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_jet3_pt_W->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_W->Fill(data.evt.MR, reweight);
      h_MTR_W->Fill(data.evt.MTR, reweight);
      h_R_W->Fill(data.evt.R, reweight);
      h_R2_W->Fill(data.evt.R2, reweight);
      h_tau21_W->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_W->Fill(data.met.Pt.at(0),reweight);

      h_HT_j1pt_W->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      if(data.hlt.PFHT800) h_HT_j1pt_W_passHT->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      if(data.hlt.AK8PFJet360_TrimMass30) h_HT_j1pt_W_passj1pt->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      if(pass) h_HT_j1pt_W_pass->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }

    if (apply_all_cuts('T')) {
      h_ht_AK4Puppi_T->Fill(AK4Puppi_Ht, reweight);
      h_ht_AK8Puppi_T->Fill(AK8Puppi_Ht, reweight);
      h_jet1_pt_T->Fill(data.jetsAK4Puppi.Pt[iJet[0]], reweight);
      h_jet2_pt_T->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_jet3_pt_T->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_T->Fill(data.evt.MR, reweight);
      h_MTR_T->Fill(data.evt.MTR, reweight);
      h_R_T->Fill(data.evt.R, reweight);
      h_R2_T->Fill(data.evt.R2, reweight);
      h_tau21_T->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_T->Fill(data.met.Pt.at(0),reweight);

      h_HT_j1pt_T->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      if(data.hlt.PFHT800) h_HT_j1pt_T_passHT->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      if(data.hlt.AK8PFJet360_TrimMass30) h_HT_j1pt_T_passj1pt->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      if(pass) h_HT_j1pt_T_pass->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }

    if (apply_all_cuts('Q')) {
      h_ht_AK4Puppi_Q->Fill(AK4Puppi_Ht, reweight);
      h_ht_AK8Puppi_Q->Fill(AK8Puppi_Ht, reweight);
      h_jet1_pt_Q->Fill(data.jetsAK4Puppi.Pt[iJet[0]], reweight);
      h_jet2_pt_Q->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_jet3_pt_Q->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_Q->Fill(data.evt.MR, reweight);
      h_MTR_Q->Fill(data.evt.MTR, reweight);
      h_R_Q->Fill(data.evt.R, reweight);
      h_R2_Q->Fill(data.evt.R2, reweight);
      h_tau21_Q->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_Q->Fill(data.met.Pt.at(0),reweight);

      h_HT_Q->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_Q->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_Q->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);

    if(data.hlt.PFHT800){
      h_HT_Q_passHT->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_Q_passHT->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_Q_passHT->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    if(data.hlt.AK8PFJet360_TrimMass30){
      h_HT_Q_passj1pt->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_Q_passj1pt->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_Q_passj1pt->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    if(pass){
      h_HT_Q_pass->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_Q_pass->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_Q_pass->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    }

    //if (apply_all_cuts_except('S', "mDPhi>=0p4")) {
    if (apply_all_cuts('S')) {
      h_ht_AK4Puppi_S->Fill(AK4Puppi_Ht, reweight);
      h_ht_AK8Puppi_S->Fill(AK8Puppi_Ht, reweight);
      h_jet1_pt_S->Fill(data.jetsAK4Puppi.Pt[iJet[0]], reweight);
      h_jet2_pt_S->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_jet3_pt_S->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_S->Fill(data.evt.MR, reweight);
      h_R_S->Fill(data.evt.R, reweight);
      h_MTR_S->Fill(data.evt.MTR, reweight);
      h_R2_S->Fill(data.evt.R2, reweight);
      h_tau21_S->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_S->Fill(data.met.Pt.at(0),reweight);

      h_GluinoLSPMass->Fill(data.evt.SUSY_Gluino_Mass,data.evt.SUSY_LSP_Mass,reweight);
      h_R2_MR->Fill(data.evt.MR, data.evt.R2, reweight);
      h_minDeltaPhi->Fill(minDeltaPhi, reweight);
      h_softDropMass->Fill(data.jetsAK8Puppi.softDropMass.at(0),reweight);
      h_StopMass->Fill(data.evt.SUSY_Stop_Mass,reweight);
      h_GluinoMass->Fill(data.evt.SUSY_Gluino_Mass,reweight);
      h_LSPMass->Fill(data.evt.SUSY_LSP_Mass,reweight);

      h_HT_S->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_S->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_S->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
 
      if(data.hlt.PFHT800){
        h_HT_S_passHT->Fill(AK4Puppi_Ht,reweight);
        h_j1_pt_S_passHT->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
        h_HT_j1pt_S_passHT->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      }
    if(data.hlt.AK8PFJet360_TrimMass30){
      h_HT_S_passj1pt->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_S_passj1pt->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_S_passj1pt->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    if(pass){
      h_HT_S_pass->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_S_pass->Fill(data.jetsAK8Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_S_pass->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    }

    h_njet_pre->Fill(nJet,        reweight);
    h_nb_pre->Fill(nMediumBTag, reweight);
    h_nw_pre->Fill(nTightWTag,  reweight);
    h_j2_pt_pre->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
    h_j3_pt_pre->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
    h_MR_pre->Fill(data.evt.MR, reweight);
    h_R2_pre->Fill(data.evt.R2, reweight);
    h_tau21_pre->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
    h_MET_pre->Fill(data.met.Pt.at(0),reweight);
    h_softDropMass_pre->Fill(data.jetsAK8Puppi.softDropMass.at(0),reweight);
    h_HT_pre->Fill(AK4Puppi_Ht,reweight);
    h_j1_pt_pre->Fill(data.jetsAK4Puppi.Pt[iJet[0]],reweight);
    h_HT_j1pt_pre->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);

    if(data.hlt.PFHT800){
      h_njet_pre_passHT->Fill(nJet,        reweight);
      h_nb_pre_passHT->Fill(nMediumBTag, reweight);
      h_nw_pre_passHT->Fill(nTightWTag,  reweight);
      h_j2_pt_pre_passHT->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_j3_pt_pre_passHT->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_pre_passHT->Fill(data.evt.MR, reweight);
      h_R2_pre_passHT->Fill(data.evt.R2, reweight);
      h_tau21_pre_passHT->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_pre_passHT->Fill(data.met.Pt.at(0),reweight);
      h_softDropMass_pre_passHT->Fill(data.jetsAK8Puppi.softDropMass.at(0),reweight);
      h_HT_pre_passHT->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_pre_passHT->Fill(data.jetsAK4Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_pre_passHT->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    if(data.hlt.AK8PFJet360_TrimMass30){
      h_njet_pre_passj1pt->Fill(nJet,        reweight);
      h_nb_pre_passj1pt->Fill(nMediumBTag, reweight);
      h_nw_pre_passj1pt->Fill(nTightWTag,  reweight);
      h_j2_pt_pre_passj1pt->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_j3_pt_pre_passj1pt->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_pre_passj1pt->Fill(data.evt.MR, reweight);
      h_R2_pre_passj1pt->Fill(data.evt.R2, reweight);
      h_tau21_pre_passj1pt->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_pre_passj1pt->Fill(data.met.Pt.at(0),reweight);
      h_softDropMass_pre_passj1pt->Fill(data.jetsAK8Puppi.softDropMass.at(0),reweight);
      h_HT_pre_passj1pt->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_pre_passj1pt->Fill(data.jetsAK4Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_pre_passj1pt->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }
    if(pass){
      h_njet_pre_pass->Fill(nJet,        reweight);
      h_nb_pre_pass->Fill(nMediumBTag, reweight);
      h_nw_pre_pass->Fill(nTightWTag,  reweight);
      h_j2_pt_pre_pass->Fill(data.jetsAK4Puppi.Pt[iJet[1]], reweight);
      h_j3_pt_pre_pass->Fill(data.jetsAK4Puppi.Pt[iJet[2]], reweight);
      h_MR_pre_pass->Fill(data.evt.MR, reweight);
      h_R2_pre_pass->Fill(data.evt.R2, reweight);
      h_tau21_pre_pass->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),reweight);
      h_MET_pre_pass->Fill(data.met.Pt.at(0),reweight);
      h_softDropMass_pre_pass->Fill(data.jetsAK8Puppi.softDropMass.at(0),reweight);
      h_HT_pre_pass->Fill(AK4Puppi_Ht,reweight);
      h_j1_pt_pre_pass->Fill(data.jetsAK4Puppi.Pt[iJet[0]],reweight);
      h_HT_j1pt_pre_pass->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],reweight);
    }

    /*
      Other examples to use analysis_cuts object

      if (apply_cut('S',"1W"))                          --> 1 Cut from S region
      if (apply_cut('W',"1Wpre"))                       --> 1 Cut from W region
      if (apply_all_cuts('T'))                          --> All cuts in T region
      if (apply_all_cuts_except('Q', "mDPhi<0.25"))     --> N-1 cut
      if (apply_all_cuts_except('S', {"0Ele", "0Mu" })) --> S without Lep veto

      But be aware: Whatever is defined in the baseline_cuts will apply to all histograms
      Also if you use skimmed ntuples (very likely) then those cuts are already applied
      This is because unskimmed ntuple is 4.3 TB in size, and we cannot have them on EOS
    */
  }
  
/*
    h_njet_pre[syst_index]->Fill(nJet,        weight);
    h_nb_pre[syst_index]->Fill(nMediumBTag, weight);
    h_nw_pre[syst_index]->Fill(nTightWTag,  weight);
    h_j1_pt_pre[syst_index]->Fill(data.jetsAK8Puppi.Pt[iJet[0]],weight);
    h_j2_pt_pre[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[1]], weight);
    h_j3_pt_pre[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[2]], weight);
    h_MR_pre[syst_index]->Fill(data.evt.MR, weight);
    h_R2_pre[syst_index]->Fill(data.evt.R2, weight);
    h_tau21_pre[syst_index]->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),weight);
    h_MET_pre[syst_index]->Fill(data.met.Pt.at(0),weight);
    h_softDropMass_pre[syst_index]->Fill(data.jetsAK8Puppi.softDropMass.at(0),weight);
    h_HT_pre[syst_index]->Fill(AK4Puppi_Ht,weight);
    h_HT_j1pt_pre[syst_index]->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],weight);

    if(pass){
      h_njet_pre_pass[syst_index]->Fill(nJet,        weight);
      h_nb_pre_pass[syst_index]->Fill(nMediumBTag, weight);
      h_nw_pre_pass[syst_index]->Fill(nTightWTag,  weight);
      h_j1 pt_pre_pass[syst_index]->Fill(data.jetsAK8Puppi.Pt[iJet[0]],weight);
      h_j2_pt_pre_pass[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[1]], weight);
      h_j3_pt_pre_pass[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[2]], weight);
      h_MR_pre_pass[syst_index]->Fill(data.evt.MR, weight);
      h_R2_pre_pass[syst_index]->Fill(data.evt.R2, weight);
      h_tau21_pre_pass[syst_index]->Fill(data.jetsAK8Puppi.tau2.at(0)/data.jetsAK8Puppi.tau1.at(0),weight);
      h_MET_pre_pass[syst_index]->Fill(data.met.Pt.at(0),weight);
      h_softDropMass_pre_pass[syst_index]->Fill(data.jetsAK8Puppi.softDropMass.at(0),weight);
      h_HT_pre_pass[syst_index]->Fill(AK4Puppi_Ht,weight);
      h_HT_j1pt_pre_pass[syst_index]->Fill(AK4Puppi_Ht,data.jetsAK8Puppi.Pt[iJet[0]],weight);
    }
*/
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  //if (apply_all_cuts('S')) vh_jet1_pt[syst_index]->Fill(data.jetsAK4Puppi.Pt[iJet[0]], weight);
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
