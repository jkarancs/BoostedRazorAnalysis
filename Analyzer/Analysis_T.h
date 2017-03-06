#ifndef VER
#define VER 0
#endif

#include "TLorentzVector.h"
#include "common/AnalysisBase.h"
#include "common/SmartHistos.h"

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


  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  baseline_cuts.push_back({ .name="Skim_R2",         .func = [&d]  { return d.evt.R2>=0.04;                  }}); // New skim cut introduced in 2017 february
  baseline_cuts.push_back({ .name="Baseline_3Jet",   .func = []    { return nJet>=3;                       }}); // Separate cut, so one can exclude (N-1)
  baseline_cuts.push_back({ .name="Baseline_MR_R2",  .func = [&d]  { return d.evt.MR>800 && d.evt.R2>0.08; }});

  /*
    W+b analysis selection

  // S: Signal region
  analysis_cuts['S'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
#if VER != 0
  analysis_cuts['S'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['S'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                  }});
  analysis_cuts['S'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                   }});
  analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});

  // S': DPhi Control region of Signal region
  analysis_cuts['s'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['s'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['s'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
#if VER != 0
  analysis_cuts['s'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['s'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                  }});
  analysis_cuts['s'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                   }});
  analysis_cuts['s'].push_back({ .name="InvmDPhi",   .func = []    { return minDeltaPhi<0.5;                 }});

  // Q: QCD enriched control sample
  analysis_cuts['Q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
#if VER != 0
  analysis_cuts['Q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['Q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                   }});
  analysis_cuts['Q'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;               }});
  analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                 }});

  // Q': Dphi Control region of QCD enriched sample
  analysis_cuts['q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
#if VER != 0
  analysis_cuts['q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                   }});
  analysis_cuts['q'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;               }});
  analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});

  // T: Top enriched control sample
  analysis_cuts['T'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepSelect==1;                   }});
  analysis_cuts['T'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                  }});
  analysis_cuts['T'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                   }});
  analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});
  analysis_cuts['T'].push_back({ .name="MT<100",     .func = []    { return MT<100;                          }});

  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepSelect==1;                   }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                   }});
  analysis_cuts['W'].push_back({ .name="1mW",        .func = []    { return nWPreTag>=1;                     }});
  analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});
  analysis_cuts['W'].push_back({ .name="30<=MT<100", .func = []    { return MT>=30 && MT<100;                }});

  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  baseline_cuts.push_back({ .name="Skim_2Jet",       .func = []    { return nJet>=2;                         }});
  baseline_cuts.push_back({ .name="Baseline_3Jet",   .func = []    { return nJet>=3;                       }}); // Separate cut, so one can exclude (N-1)
  baseline_cuts.push_back({ .name="Baseline_MR_R2",  .func = [&d]  { return d.evt.MR>800 && d.evt.R2>0.08; }});
  
  */


  // t: >=1 top Signal region
  analysis_cuts['t'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['t'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['t'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
#if VER != 0
  analysis_cuts['t'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                }});
#endif
  analysis_cuts['t'].push_back({ .name="1top",       .func = []    { return nHadTopTag>=1;                   }});
  analysis_cuts['t'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                }});


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
  double sf_ele_veto = sf_ele.first/*, sf_ele_medium = sf_ele.second*/;
  i+=4;

  // Muon SFs (7 sigmas - tracking, fullsim id/iso/ip, fastsim id/iso/ip)
  std::pair<double, double> sf_muon = calc_muon_sf(d, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s],
						   nSigmaSFs[i+4][s], nSigmaSFs[i+5][s], nSigmaSFs[i+6][s], isFastSim);
  double sf_muon_veto = sf_muon.first/*, sf_muon_medium = sf_muon.second*/;
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

  scale_factors['s'] = scale_factors['S'];

  scale_factors['Q'].push_back(sf_ele_veto);
  scale_factors['Q'].push_back(sf_muon_veto);
  scale_factors['Q'].push_back(sf_btag_loose);
  scale_factors['Q'].push_back(sf_w);

  scale_factors['q'] = scale_factors['Q'];

  scale_factors['T'].push_back(sf_ele_medium);
  scale_factors['T'].push_back(sf_muon_medium);
  scale_factors['T'].push_back(sf_btag_medium);
  scale_factors['T'].push_back(sf_w);

  scale_factors['W'].push_back(sf_ele_medium);
  scale_factors['W'].push_back(sf_muon_medium);
  scale_factors['W'].push_back(sf_btag_loose);

  */

  // >=1 top analysis
  scale_factors['t'].push_back(sf_ele_veto);
  scale_factors['t'].push_back(sf_muon_veto);
  scale_factors['t'].push_back(sf_top);
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& d) {
  return apply_all_cuts('t');
}

//_______________________________________________________
//                 List of Histograms
TH1D* h_njet;
TH1D* h_nb;
TH1D* h_nw;
TH1D* h_ht_gen;
TH1D* h_ht_AK4;
TH1D* h_ht_AK8;
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

TH1D* h_AK8_tau32;
TH1D* h_AK8_tau31;
TH1D* h_AK8_tau21;



std::vector<TH1D*> vh_jet1_pt;

//Analysis::PostfixOptions
//Analysis::get_pf_opts_(std::vector<std::vector<Sample> > lists, std::string dirname) {
//  std::vector<Sample> samples;
//  for (auto list : lists) samples.insert(samples.end(), list.begin(), list.end());
//  PostfixOptions opt{ (size_t)-1, "", "", "" };
//  for (size_t i=0; i<samples.size(); ++i) {
//    // Find index of matching directory
//    for (size_t j=0; j<samples[i].dirs.size(); ++j)
//      if (samples[i].dirs[j] == dirname) opt.index = i;
//    opt.postfixes += samples[i].postfix;
//    opt.legends += samples[i].legend;
//    opt.colors += samples[i].color;
//    if (i+1!=samples.size()) {
//      opt.postfixes +=  ";";
//      opt.legends += ";";
//      opt.colors += ",";
//    }
//  }
//  return opt;
//}

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  h_njet           = new TH1D("njet",         ";N_{jet}",                20, 0,  20);
  h_nw             = new TH1D("nw",           ";N_{W tag}",              20, 0,  20);
  h_nb             = new TH1D("nb",           ";N_{b tag}",              20, 0,  20);
  h_ht_gen         = new TH1D("ht_gen",       ";H_{T}^{gen}",            200, 0,2000);
  h_ht_AK4         = new TH1D("ht_AK4",  ";H_{T}",                  200, 0,2000);
  h_ht_AK8         = new TH1D("ht_AK8",  ";H_{T}^{AK8}",            200, 0,2000);
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
  h_AK8_tau32      = new TH1D("tau32", "", 200,0,1);
  h_AK8_tau31      = new TH1D("tau31", "", 200,0,1);
  h_AK8_tau21      = new TH1D("tau21", "", 200,0,1);


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

    double w = weight; // No scale factor applied
    //double w = sf_weight['t']; // Scale factors applied for the Signal region

    if (apply_cut('t',"HLT")) {

      h_njet   ->Fill(nJet,        w);
      h_nb     ->Fill(nMediumBTag, w);
      h_nw     ->Fill(nTightWTag,  w);
    
      h_ht_gen->Fill(d.evt.Gen_Ht,  w);  // in ntuple
      h_ht_AK4->Fill(AK4_Ht, w); // Calculated in AnalysisBase.h
      h_ht_AK8->Fill(AK8_Ht, w); // Calculated in AnalysisBase.h


      h_R2->Fill(d.evt.R2, w);
      h_MR->Fill(d.evt.MR, w);
      h_R2_MR->Fill(d.evt.MR, d.evt.R2, w);

      // Janos Did you want to apply scale factors?
      // (This is the new way to do it)
      w = sf_weight['t'];
      h_R2_Scaled->Fill(d.evt.R2, w);
      h_R2_MR_Scaled->Fill(d.evt.R2, d.evt.MR, w);
      h_MR_Scaled->Fill(d.evt.MR, w);
    }




    // >=1 top Signal region
    // Apply scale factors corresponding to this region
    w = sf_weight['t'];
    if (apply_all_cuts('t')) h_AK8_tau32->Fill(tau21.at(0),w); 
    if (apply_all_cuts('t')) h_AK8_tau31->Fill(tau21.at(0),w); 
    if (apply_all_cuts('t')) h_AK8_tau21->Fill(tau21.at(0),w); 
    if (apply_all_cuts('t')) h_nhadtop->Fill(nHadTopTag,w);

    // Janos: what is this?
    // I think this should be moved to the
    // save_analysis_histos method (a final operation before writing them to files
    // But I usually do such operations after producing the plots in a separate process
    // If you instead wanted to apply scale factor I already did it above, and you can delete this)
    /*
    if (h_R2_Scaled->Integral()!=0)
      h_R2_Scaled->Scale(1/h_R2_Scaled->Integral());   

    if (h_R2_MR_Scaled->Integral()!=0)
      h_R2_MR_Scaled->Scale(1/h_R2_MR_Scaled->Integral());

    if (h_MR_Scaled->Integral()!=0)
      h_MR_Scaled->Scale(1/ h_MR_Scaled->Integral());
    */
   


    //GenSize 

    int NtopMult=0;
 
    for (unsigned int i=0; i<d.gen.size; i++) 
      {
	// if (d.gen.Status[i] != 3) continue;
	if (fabs(d.gen.ID[i]) == 6) 
	  {
	    if ((fabs(d.gen.Dau1ID[i]) == 5 && fabs(d.gen.Dau0ID[i]) == 24) ||
		(fabs(d.gen.Dau1ID[i]) == 24 && fabs(d.gen.Dau0ID[i]) == 5))
	      { 
		h_gen_toppt->Fill(d.gen.Pt[i]);

		NtopMult++;
	      }
	  }
      }
    h_NtopMult->Fill(NtopMult);

     
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
  }
  
  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  double w = sf_weight['t'];
  if (apply_all_cuts('t')) vh_jet1_pt[syst_index]->Fill(d.jetsAK4.Pt[iJet[0]], w);
}

////>>>>>>>>>>>>>>>>>> Methods used by SmartHistos (Plotter)>>>>>>>>>>>>>>>>>>

void
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst, 
			       const unsigned int& syst_index, bool runOnSkim=false)
{

  //   std::vector<Sample> signal_all, signal_selected, signal_fastsim, signal_gluino, signal_stop;
  //   signal_all.push_back({ .postfix="T5ttcc",       .legend="T5ttcc",      .color="12", /*DGrey*/ .dirs={ "FastSim_SMS-T5ttcc" } });
  //   signal_all.push_back({ .postfix="T5tttt",       .legend="T5tttt",      .color="862",/*Azure*/ .dirs={ "FastSim_SMS-T5tttt" } });
  //   signal_all.push_back({ .postfix="T1tttt",       .legend="T1tttt",      .color="841",/*Teal*/  .dirs={ "FastSim_SMS-T1tttt" } });
  //   signal_all.push_back({ .postfix="T2tt",         .legend="T2tt",        .color="403",/*DYell*/ .dirs={ 
  //          "FastSim_SMS-T2tt_mStop-150to250", "FastSim_SMS-T2tt_mStop-250to350",
  //          "FastSim_SMS-T2tt_mStop-350to400", "FastSim_SMS-T2tt_mStop-400to1200" 
  //        } });
  //   signal_all.push_back({ .postfix="T2tt_FullSim", .legend="T2tt (FullSim)", .color="804",/*DOran*/ .dirs={
  //          "FullSim_SMS-T2tt_mStop-425_mLSP-325", "FullSim_SMS-T2tt_mStop-500_mLSP-325",
  //          "FullSim_SMS-T2tt_mStop-850_mLSP-100" 
  //        } });
  //   signal_selected.push_back(signal_all[0]);
  //   for (int i=0; i<4; ++i) signal_fastsim.push_back(signal_all[i]);
  //   for (int i=0; i<3; ++i) signal_gluino .push_back(signal_all[i]);
  //   for (int i=3; i<5; ++i) signal_stop .push_back(signal_all[i]);
  // 
  // 
  // 
  //   static const PostfixOptions signals_background_opt = get_pf_opts_({signal_all}, dirname);
  //   sh.AddNewPostfix("Signals,Background",  [&d] { 
  // 		     // Select gluino/stop mass to give ~1k events with 40 fb^-1
  // 		     if (signals_background_opt.index==1) {
  // 		       if (d.evt.SUSY_Gluino_Mass == 1200 && d.evt.SUSY_LSP_Mass == 200) return (size_t)-1; // T5tttt
  // 		     } else if (signals_background_opt.index==3) {
  // 		       if (d.evt.SUSY_Stop_Mass  == 800 && d.evt.SUSY_LSP_Mass == 100) return (size_t)-1; // T2tt - Same as FullSim point
  // 		     }
  // 		     return signals_background_opt.index; 
  // 		   }, signals_background_opt.postfixes, signals_background_opt.legends, signals_background_opt.colors);
}

void
Analysis::load_analysis_histos(std::string inputfile)
{
}

void
Analysis::save_analysis_histos(bool draw=0)
{

}

