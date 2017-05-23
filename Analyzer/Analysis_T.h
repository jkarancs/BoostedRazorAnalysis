#ifndef VER
#define VER 0
#endif

#include "TLorentzVector.h"
#include "common/AnalysisBase.h"



//_______________________________________________________
//                  Calculate variables


//unsigned int nLooseIDHadTopTagJets;  


void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{

  //nLooseIDHadTopTagJets = 0;
  //
  //// Loop on AK8 jets
  //
  //while(data.jetsAK8.Loop()) {
  //  if (data.jetsAK8.looseJetID[data.jetsAK8.it]) {
  //    if (passHadTopTag[data.jetsAK8.it]) ++nLooseIDHadTopTagJets;       
  //  }
  //} // end of AK8 jet loop
}

//_______________________________________________________
//                Define Skimming cuts
//   (Not needed, unless you want to skim the ntuple)

bool
Analysis::pass_skimming(DataStruct& data)
{
  //int NJetAK8 = 0;
  //while(data.jetsAK8.Loop()) {
  //  size_t i = data.jetsAK8.it;
  //  // pt cut intentionally removed to accept all jets for systematics
  //  if ( data.jetsAK8.looseJetID[i] == 1 &&
  //       std::abs(data.jetsAK8.Eta[i])  <  JET_AK8_ETA_CUT ) {
  //    NJetAK8++;
  //  }
  //}
  //if (!(NJetAK8>=1)) return 0;
  //if (!(data.evt.R2>=0.04)) return 0;

  float pt_threshold = 300;
  int N = 0;
  while(data.jetsAK8.Loop()) if (data.jetsAK8.Pt[data.jetsAK8.it]>=pt_threshold) ++N;
  return (N >= 1);

  // Signal skim
  //return apply_all_cuts("S");
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

  baseline_cuts.push_back({ .name="signal_mass_selection",  .func = [this,&d]  { return isSignal ? d.evt.SUSY_Gluino_Mass==1200 && d.evt.SUSY_LSP_Mass == 300 : 1;     }});
  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  baseline_cuts.push_back({ .name="Skim_R2",         .func = [&d]  { return d.evt.R2>=0.04;                  }}); // New skim cut introduced in 2017 february
  
  // Q: QCD enriched control sample
  analysis_cuts['Q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
  analysis_cuts['Q'].push_back({ .name="MR_R2", .func = [&d] { return d.evt.MR>=800 && d.evt.R2>=0.04; }});
#if VER != 0
  analysis_cuts['Q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif

  analysis_cuts['Q'].push_back({ .name="1topAntiBtag",         .func = []    { return nHadTop0BAntiTag>=1;                   }});
  analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                 }});


  // Q': Dphi Control region of QCD enriched sample
  analysis_cuts['q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
  analysis_cuts['q'].push_back({ .name="MR_R2", .func = [&d] { return d.evt.MR>=800 && d.evt.R2>=0.04; }});
#if VER != 0
  analysis_cuts['q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['q'].push_back({ .name="1aTop",        .func = []    { return nHadTop0BAntiTag>=1;               }});
  analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});


  // T: Region Top enriched control sample 
  analysis_cuts['T'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['T'].push_back({ .name="MR_R2", .func = [&d] { return d.evt.MR>=800 && d.evt.R2>=0.04; }});
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepSelect==1;                   }});
  analysis_cuts['T'].push_back({ .name="1top",         .func = []    { return nHadTopTag>=1;                  }});
  analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});
  analysis_cuts['T'].push_back({ .name="MT<100",     .func = []    { return MT<100;                          }});


  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['W'].push_back({ .name="MR_R2", .func = [&d] { return d.evt.MR>=800 && d.evt.R2>=0.04; }});
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepSelect==1;                   }});
  analysis_cuts['W'].push_back({ .name="1topBpretag",         .func = []    { return nHadTop0BPreTag>=1;                   }});
  analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});
  analysis_cuts['W'].push_back({ .name="30<=MT<100", .func = []    { return MT>=30 && MT<100;                }});


  // t: >=1 top Signal region
  analysis_cuts['S'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
  analysis_cuts['S'].push_back({ .name="MR_R2", .func = [&d] { return d.evt.MR>=800 && d.evt.R2>=0.04; }});
#if VER != 0
  analysis_cuts['S'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['S'].push_back({ .name="1top",       .func = []    { return nHadTopTag>=1;                   }});
  analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});

 // S': DPhi Control region of Signal region
  analysis_cuts['s'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet360_TrimMass30==1 || d.hlt.PFHT800==1 : 1; }});
  analysis_cuts['s'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['s'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
  analysis_cuts['s'].push_back({ .name="MR_R2", .func = [&d] { return d.evt.MR>=800 && d.evt.R2>=0.04; }});
#if VER != 0
  analysis_cuts['s'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
 analysis_cuts['s'].push_back({ .name="1top",       .func = []    { return nHadTopTag>=1;                   }});
 analysis_cuts['s'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi<0.5;                }});
  

}//End of define_selections


void
Analysis::apply_scale_factors(DataStruct& d, const unsigned int& s, const std::vector<std::vector<double> >& nSigmaSFs)
{
  bool isFastSim = TString(sample).Contains("FastSim");
  size_t i = 0;

  // Don't forget to specify the total number of sigmas you use in settings_*.h !
  // Here it is 12 currently (4 mu, 7 ele, 1 top)

  // Electron SFs (4 sigmas - reco, iso, id, fastsim)
  std::pair<double, double> sf_ele = calc_ele_sf(d, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s], isFastSim);
  double sf_ele_veto = sf_ele.first, sf_ele_medium = sf_ele.second;
  i+=4;

  // Muon SFs (7 sigmas - tracking, fullsim id/iso/ip, fastsim id/iso/ip)
  std::pair<double, double> sf_muon = calc_muon_sf(d, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s],
						   nSigmaSFs[i+4][s], nSigmaSFs[i+5][s], nSigmaSFs[i+6][s], isFastSim);
  double sf_muon_veto = sf_muon.first, sf_muon_medium = sf_muon.second;
  i+=7;

  //W+b  // W tagging SF  (1 sigma - efficiency)
  //W+b  double sf_w = calc_w_tagging_sf(d, nSigmaSFs[i][s]);
  //W+b  i+=1;
  //W+b  
  //W+b  // b tagging SFs (1 sigma)
  //W+b  std::pair<double, double> sf_btag = calc_b_tagging_sf(d, nSigmaSFs[i][s], isFastSim);
  //W+b  double sf_btag_loose = sf_btag.first, sf_btag_medium = sf_btag.second;
  //W+b  i+=1;

  // top tagging SF (1 sigma)
  double sf_top = calc_top_tagging_sf(d, nSigmaSFs[i++][s]);
  i+=1;

  // Select scale factors to use
  for (auto& sf : scale_factors) sf.second.clear();

  /*
    W+b analysis scale factors

  scale_factors['S'].push_back(sf_ele_veto);
  scale_factors['S'].push_back(sf_muon_veto);
  scale_factors['S'].push_back(sf_btag_medium);
  scale_factors['S'].push_back(sf_w);
*/
  scale_factors['s'] = scale_factors['S'];

  scale_factors['Q'].push_back(sf_ele_veto);
  scale_factors['Q'].push_back(sf_muon_veto);
  scale_factors['Q'].push_back(sf_top);
  //scale_factors['Q'].push_back(sf_w);

  scale_factors['q'] = scale_factors['Q'];

  scale_factors['T'].push_back(sf_ele_medium);
  scale_factors['T'].push_back(sf_muon_medium);
  scale_factors['T'].push_back(sf_top);
 // scale_factors['T'].push_back(sf_w);

  scale_factors['W'].push_back(sf_ele_medium);
  scale_factors['W'].push_back(sf_muon_medium);
  scale_factors['W'].push_back(sf_top);

  

  // >=1 top analysis
  scale_factors['S'].push_back(sf_ele_veto);
  scale_factors['S'].push_back(sf_muon_veto);
  scale_factors['S'].push_back(sf_top);
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& d) {
  return apply_all_cuts('S');
}

//_______________________________________________________
//                 List of Histograms


TH1D* h_S_MR;
TH1D* h_s_MR;
TH1D* h_S_R2;
TH1D* h_s_R2;

TH1D* h_S_MR_nj35;
TH1D* h_S_R2_nj35;
TH1D* h_s_MR_nj35;
TH1D* h_s_R2_nj35;

TH1D* h_S_MR_nj5;
TH1D* h_S_R2_nj5;
TH1D* h_s_MR_nj5;
TH1D* h_s_R2_nj5;

TH1D* h_Q_MR;
TH1D* h_Q_R2;
TH1D* h_Q_MR_nj35;
TH1D* h_Q_R2_nj35;
TH1D* h_Q_MR_nj5;
TH1D* h_Q_R2_nj5;
TH1D* h_q_MR;
TH1D* h_q_R2;
TH1D* h_q_MR_nj35;
TH1D* h_q_R2_nj35;
TH1D* h_q_MR_nj5;
TH1D* h_q_R2_nj5;

TH1D* h_W_MR;
TH1D* h_W_R2;
TH1D* h_W_MR_nj35;
TH1D* h_W_R2_nj35;
TH1D* h_W_MR_nj5;
TH1D* h_W_R2_nj5;


TH1D* h_T_MR;
TH1D* h_T_R2;
TH1D* h_T_MR_nj35;
TH1D* h_T_R2_nj35;
TH1D* h_T_MR_nj5;
TH1D* h_T_R2_nj5;


std::vector<TH1D*> vh_jet1_pt;


//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
 // h_gen_toppt      = new TH1D("h_gen_toppt", "h_gen_toppt", 50, 0, 1000);
  
 // int nbn_HT = 19;
  int nbn_MR = 6;
  int nbn_R2 = 7;
 
  Double_t bn_MR_tmp[] = {0.,800.,1000.,1200.,1600.,2000.,4000.};
  Double_t* bn_MR = 0;
  bn_MR = utils::getVariableBinEdges(nbn_MR+1,bn_MR_tmp);
  Double_t bn_R2_tmp[] = {0.,0.04,0.08,0.12,0.16,0.24,0.5,1.};
  Double_t* bn_R2 = 0;
  bn_R2 = utils:: getVariableBinEdges(nbn_R2+1,bn_R2_tmp);

//S Region; S, S34, Sge5--> nJet, HT, MET, MR, R2, MRvsR2...15 plots
  
 
  h_S_MR           = new TH1D("S_MR ",";S_MR", nbn_MR, bn_MR);
  h_S_R2           = new TH1D("S_R2 ",";S_R2", nbn_R2, bn_R2 );
  h_s_MR           = new TH1D("s_MR ",";s_MR", nbn_MR, bn_MR);
  h_s_R2           = new TH1D("s_R2 ",";s_R2", nbn_R2, bn_R2 );
 
  //nJet>=3 && nJet<5 
 
 	 
  h_S_MR_nj35      = new TH1D("S_MR_nj35",";S_MR", nbn_MR, bn_MR);
  h_S_R2_nj35      = new TH1D("S_R2_nj35",";S_R2", nbn_R2, bn_R2 );
  h_s_MR_nj35      = new TH1D("s_MR_nj35",";s_MR", nbn_MR, bn_MR);
  h_s_R2_nj35      = new TH1D("s_R2_nj35",";s_R2", nbn_R2, bn_R2 );
  
  //nJet>=5

 
  h_S_MR_nj5       = new TH1D("S_MR_nj5",";S_MR", nbn_MR, bn_MR);
  h_S_R2_nj5       = new TH1D("S_R2_nj5",";S_R2", nbn_R2, bn_R2 );
  h_s_MR_nj5       = new TH1D("s_MR_nj5",";s_MR", nbn_MR, bn_MR);
  h_s_R2_nj5       = new TH1D("s_R2_nj5",";s_R2", nbn_R2, bn_R2 ); 
  

//Q Region; Q, Q34, Qge5--> MR, R2, MRvsR2 ...9 Plots
  h_Q_MR           = new TH1D("Q_MR ",";Q_MR", nbn_MR, bn_MR);
  h_Q_R2           = new TH1D("Q_R2 ",";Q_R2", nbn_R2, bn_R2 );
  
  h_Q_MR_nj35      = new TH1D("Q_MR_nj35",";Q_MR", nbn_MR, bn_MR);         
  h_Q_R2_nj35      = new TH1D("Q_R2_nj35",";Q_R2", nbn_R2, bn_R2 );

  h_q_MR           = new TH1D("q_MR ",";q_MR", nbn_MR, bn_MR);
  h_q_R2           = new TH1D("q_R2 ",";q_R2", nbn_R2, bn_R2 );
  
  h_q_MR_nj35      = new TH1D("q_MR_nj35",";q_MR", nbn_MR, bn_MR);         
  h_q_R2_nj35      = new TH1D("q_R2_nj35",";q_R2", nbn_R2, bn_R2 );
  
  //nJet>=5
  h_Q_MR_nj5       = new TH1D("Q_MR_nj5",";Q_MR", nbn_MR, bn_MR);
  h_Q_R2_nj5       = new TH1D("Q_R2_nj5",";Q_R2", nbn_R2, bn_R2 );

  h_q_MR_nj5       = new TH1D("q_MR_nj5",";q_MR", nbn_MR, bn_MR);
  h_q_R2_nj5       = new TH1D("q_R2_nj5",";q_R2", nbn_R2, bn_R2 );
 

//W Region; W, W34, Wge5--> MR, R2, MRvsR2 ...9 Plots
  h_W_MR           = new TH1D("W_MR ",";W_MR", nbn_MR, bn_MR);   
  h_W_R2           = new TH1D("W_R2 ",";W_R2", nbn_R2, bn_R2 );
  
  //nJet>=3 && nJet<5
  h_W_MR_nj35      = new TH1D("W_MR_nj35",";W_MR", nbn_MR, bn_MR);
  h_W_R2_nj35      = new TH1D("W_R2_nj35",";W_R2", nbn_R2, bn_R2 );
  
  //nJet>=5
  h_W_MR_nj5       = new TH1D("W_MR_nj5",";W_MR", nbn_MR, bn_MR);
  h_W_R2_nj5       = new TH1D("W_R2_nj5",";W_R2", nbn_R2, bn_R2 ); 
 
//T Region; T, T34, Tge5--> MR, R2, MRvsR2 ...9 Plots 
  h_T_MR           = new TH1D("T_MR ",";T_MR", nbn_MR, bn_MR);
  h_T_R2           = new TH1D("T_R2 ",";T_R2", nbn_R2, bn_R2 );
  
  //nJet>=3 && nJet<5
  h_T_MR_nj35      = new TH1D("T_MR_nj35",";T_MR", nbn_MR, bn_MR);
  h_T_R2_nj35      = new TH1D("T_R2_nj35",";T_R2", nbn_R2, bn_R2);
 
  //nJet>=5
  h_T_MR_nj5       = new TH1D("T_MR_nj5",";T_MR", nbn_MR, bn_MR);
  h_T_R2_nj5       = new TH1D("T_R2_nj5",";T_R2", nbn_R2, bn_R2 );
 	

  
  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    histoname<<"jet1_pt_syst"<<i;
    title<<"Systematic variation #="<<i;
    vh_jet1_pt.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
}
}


//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_analysis_histos(DataStruct& d, const unsigned int& syst_index, const double& weight)
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

    /*
      Weight:
      They now include trigger efficiencies for MC by default
      w is the event weight without any scale factor applied
      Because scale factors are region dependend, then
      in order to apply them, one has to use the sf_weight[region] variable instead
      eg. sf_weight['t']
     */

    // Baseline cuts 
    // Additionally, let's apply the trigger selection
    // since the weight contains the trigger scaling in MC
    // no specific region, so don't apply scale factors
    // Especially for comparison plots with Changgi
    // Alternatively, could apply SF of the Signal regio

    //double w = weight; // No scale factor applied
    //double w = sf_weight['t']; // Scale factors applied for the Signal region
    
    

 
    //  // For example this applies the first three cuts in signal region
    //  // HLT, ele/mu veto
    //  if (apply_ncut('t', 3)) {
    //    h_jet1_pt->Fill(d.jetsAK4.Pt[iJet[0]], w);
    //    h_jet2_pt->Fill(d.jetsAK4.Pt[iJet[1]], w);
    //    h_jet3_pt->Fill(d.jetsAK4.Pt[iJet[2]], w);
    //  }

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
   
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
   double w = sf_weight['S'];
   if (apply_all_cuts('S')) {
	vh_jet1_pt[syst_index]->Fill(d.jetsAK4.Pt[iJet[0]], w);
	h_S_MR -> Fill(d.evt.MR, w);	
 	h_S_R2 -> Fill(d.evt.R2, w);
	h_s_MR -> Fill(d.evt.MR, w);	
 	h_s_R2 -> Fill(d.evt.R2, w);
		
	
	if(nJet>=3 && nJet<5){
		h_S_MR_nj35 -> Fill(d.evt.MR, w);								
		h_S_R2_nj35 -> Fill(d.evt.R2, w);
		h_s_MR_nj35 -> Fill(d.evt.MR, w);								
		h_s_R2_nj35 -> Fill(d.evt.R2, w);
		
	}
	
	if(nJet>5){
		h_S_MR_nj5 -> Fill(d.evt.MR, w);
		h_S_R2_nj5 -> Fill(d.evt.R2, w);
		h_s_MR_nj5 -> Fill(d.evt.MR, w);
		h_s_R2_nj5 -> Fill(d.evt.R2, w);
		
	}
   }

    w = sf_weight['W'];
    if (apply_all_cuts('W')){
	h_W_MR -> Fill(d.evt.MR, w);	  
	h_W_R2 -> Fill(d.evt.R2, w);
	
	
	if(nJet>=3 && nJet<5){
		h_W_MR_nj35 -> Fill(d.evt.MR, w);	
		h_W_R2_nj35 -> Fill(d.evt.R2, w);
		
	}
	
	if(nJet>5){
		h_W_MR_nj5 -> Fill(d.evt.MR, w);
		h_W_R2_nj5 -> Fill(d.evt.R2, w);
		
	}	
   }

    w = sf_weight['Q'];
    if (apply_all_cuts('Q')){
        h_Q_MR -> Fill(d.evt.MR, w);
        h_Q_R2 -> Fill(d.evt.R2, w);
	h_q_MR -> Fill(d.evt.MR, w);
        h_q_R2 -> Fill(d.evt.R2, w);
       

        if(nJet>=3 && nJet<5){
                h_Q_MR_nj35 -> Fill(d.evt.MR, w);               
                h_Q_R2_nj35 -> Fill(d.evt.R2, w);
		h_q_MR_nj35 -> Fill(d.evt.MR, w);               
                h_q_R2_nj35 -> Fill(d.evt.R2, w);
                
        }

        if(nJet>5){
                h_Q_MR_nj5 -> Fill(d.evt.MR, w);
                h_Q_R2_nj5 -> Fill(d.evt.R2, w);
                h_q_MR_nj5 -> Fill(d.evt.MR, w);
                h_q_R2_nj5 -> Fill(d.evt.R2, w);
                
        }
   }

    w = sf_weight['T'];
    if (apply_all_cuts('T')){
        h_T_MR -> Fill(d.evt.MR, w);
        h_T_R2 -> Fill(d.evt.R2, w);
       

        if(nJet>=3 && nJet<5){
                h_T_MR_nj35 -> Fill(d.evt.MR, w);               
                h_T_R2_nj35 -> Fill(d.evt.R2, w);
                
        }

        if(nJet>5){
                h_T_MR_nj5 -> Fill(d.evt.MR, w);
                h_T_R2_nj5 -> Fill(d.evt.R2, w);
                
        }
   }

}
}

////>>>>>>>>>>>>>>>>>> Methods used by SmartHistos (Plotter)>>>>>>>>>>>>>>>>>>

void
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst, const unsigned int& syst_index, bool runOnSkim=false)
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

