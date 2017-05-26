#ifndef VER
#define VER 0
#endif

#include "TLorentzVector.h"
#include "common/AnalysisBase.h"
#include "common/SmartHistos.h"

SmartHistos sh;

//_______________________________________________________
//                  Calculate variables

// Cut variables
static size_t cut_index;
std::map<char, unsigned int> cutbits;
std::map<char, bool> pass_all_cuts;

// N-1 weights
std::map<char, std::vector<double> > w_nm1;

void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{
  cut_index = -1;

  // It only makes sense to calculate certain variables only once if they don't depend on jet energy
  //if (syst_index == 0) {
  //  
  //}

  // Calculate decision of each individual cut
  for (const auto& region : analysis_cuts) {
    cutbits[region.first] = 0;
    for (size_t i=0, n=analysis_cuts[region.first].size(); i<n; ++i)
      if (analysis_cuts[region.first][i].func()) cutbits[region.first] += 1<<i;
    pass_all_cuts[region.first] = (cutbits[region.first]==(unsigned int)((1<<analysis_cuts[region.first].size())-1));
  }
}

//_______________________________________________________
//                Define Skimming cuts
//   (Not needed, unless you want to skim the ntuple)

bool
Analysis::pass_skimming(DataStruct& data)
{
  int NJetAK8 = 0;
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    // pt cut intentionally removed to allow studying systematics
    if ( data.jetsAK8.looseJetID[i] == 1 &&
         std::abs(data.jetsAK8.Eta[i])  <  JET_AK8_ETA_CUT ) {
      NJetAK8++;
    }
  }
  if (!(NJetAK8>=1)) return 0;
  if (!(data.evt.R2>=0.04)) return 0;
  return 1;
}

//_______________________________________________________
//          Define Analysis event selection cuts
//     Can define all sorts of Signal/Control regions

enum SCuts { S_3Jet, S_MR_R2, S_HLT, S_0Ele, S_0Mu, S_0IsoTrk, S_1b, S_1W,  S_mDPhi };
enum QCuts { Q_3Jet, Q_MR_R2, Q_HLT, Q_0Ele, Q_0Mu, Q_0IsoTrk, Q_0b, Q_1aW, Q_InvmDPhi0p3 };
enum TCuts { T_3Jet, T_MR_R2, T_HLT, T_1Lep,                   T_1b, T_1W,  T_mDPhi,    T_MT};
enum WCuts { W_3Jet, W_MR_R2, W_HLT, W_1Lep,                   W_0b, W_1mW, W_mDPhi,    W_MT};

void
Analysis::define_selections(const DataStruct& d)
{
  analysis_cuts.clear();

  // Define here cuts that are common in all Signal/Control regions
  // MET Filters, etc. are already applied in AnalysisBase.h, See baseline_cuts
  //baseline_cuts.push_back({ .name="Skim_R2",         .func = [&d]  { return d.evt.R2>=0.04;                   }}); // New skim cut introduced in 2017 february
  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                       }}); // Similar to pt>200, one AK8 jet has pt>200

  // Remove baseline cuts for btag efficiencies or skimming
  baseline_cuts.clear();

  // S: Signal region
  analysis_cuts['S'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['S'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['S'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
#if VER != 0
  analysis_cuts['S'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
#endif
  analysis_cuts['S'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['S'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});

  // S': DPhi Control region of Signal region
  analysis_cuts['s'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['s'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['s'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['s'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['s'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
#if VER != 0
  analysis_cuts['s'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
#endif
  analysis_cuts['s'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['s'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  analysis_cuts['s'].push_back({ .name="InvmDPhi",   .func = []    { return minDeltaPhi<0.5;                  }});

  // Q: QCD enriched control sample
  analysis_cuts['Q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['Q'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['Q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
#if VER != 0
  analysis_cuts['Q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
#endif
  analysis_cuts['Q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['Q'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;                }});
  analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                  }});

  // Q': Dphi Control region of QCD enriched sample
  analysis_cuts['q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['q'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
#if VER != 0
  analysis_cuts['q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
#endif
  analysis_cuts['q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['q'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;                }});
  analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});

  // T: Top enriched control sample
  analysis_cuts['T'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['T'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['T'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['T'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['T'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['T'].push_back({ .name="MT",         .func = []    { return MT_vetolep<100;                   }});

  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['W'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['W'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['W'].push_back({ .name="1mW",        .func = []    { return nWPreTag>=1;                      }});
  analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['W'].push_back({ .name="MT",         .func = []    { return MT_vetolep>=30 && MT_vetolep<100; }});

  // Z: Z->ll enriched control sample
  analysis_cuts['Z'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['Z'].push_back({ .name="MR_R2ll",    .func = [&d]  { return d.evt.MR>=800 && R2_ll>=0.08;     }});
  analysis_cuts['Z'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Z'].push_back({ .name="2Lep",       .func = []    { return (nEleSelect==2&&nMuVeto==0)||(nMuSelect==2&&nEleVeto==0); }});
  analysis_cuts['Z'].push_back({ .name="OppCharge",  .func = [&d]  { 
				   if (nEleSelect==2) return (d.ele.Charge[iEleSelect[0]] + d.ele.Charge[iEleSelect[1]])==0;
				   else if (nMuSelect==2) return (d.mu.Charge[iMuSelect[0]] + d.mu.Charge[iMuSelect[1]])==0;
				   return false;
				 }});
  analysis_cuts['Z'].push_back({ .name="1mW",        .func = []    { return nWPreTag>=1;                      }});
  analysis_cuts['Z'].push_back({ .name="mDPhill",    .func = []    { return minDeltaPhi_ll>=0.5;              }});
  analysis_cuts['Z'].push_back({ .name="Mll",        .func = []    { return std::abs(M_ll-91.2)<10;           }});

  // t: Boosted Top Signal region
  analysis_cuts['t'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['t'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['t'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['t'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['t'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
#if VER != 0
  analysis_cuts['t'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
#endif
  analysis_cuts['t'].push_back({ .name="1Top",       .func = []    { return nHadTopTag>=1;                    }});
  analysis_cuts['t'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});

}

//____________________________________________________
//          Analysis Specific Scale factors
//    (Defined for each search region separately)

void
Analysis::apply_scale_factors(DataStruct& data, const unsigned int& s, const std::vector<std::vector<double> >& nSigmaSFs)
{
  bool isFastSim = TString(sample).Contains("FastSim");
  size_t i = 0;

  // Don't forget to specify the total number of sigmas you use in settings_*.h !

  // Electron SFs (5 sigmas - reco, fullsim id/iso, fastsim)
  double sf_ele_veto, sf_ele_loose, sf_ele_medium;
  std::tie(sf_ele_veto, sf_ele_loose, sf_ele_medium) = calc_ele_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s], isFastSim);
  i+=4;

  // Muon SFs (3 sigmas - tracking, fullsim, fastsim)
  double sf_muon_veto, sf_muon_loose, sf_muon_medium;
  std::tie(sf_muon_veto, sf_muon_loose, sf_muon_medium) =  calc_muon_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], isFastSim);
  i+=3;

  // W tagging SF  (1 sigma - efficiency)
  double sf_w = calc_w_tagging_sf(data, nSigmaSFs[i][s]);
  i+=1;

  // b tagging SFs (1 sigma)
  std::pair<double, double> sf_btag = calc_b_tagging_sf(data, nSigmaSFs[i][s], isFastSim);
  double sf_btag_loose = sf_btag.first, sf_btag_medium = sf_btag.second;
  i+=1;

  // top tagging SF (1 sigma)
  double sf_top = calc_top_tagging_sf(data, nSigmaSFs[i][s]);
  i+=1;

  // Select scale factors to use
  for (auto& sf : scale_factors) sf.second.clear();
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

  scale_factors['T'].push_back(sf_ele_veto);
  scale_factors['T'].push_back(sf_muon_veto);
  scale_factors['T'].push_back(sf_btag_medium);
  scale_factors['T'].push_back(sf_w);

  scale_factors['W'].push_back(sf_ele_veto);
  scale_factors['W'].push_back(sf_muon_veto);
  scale_factors['W'].push_back(sf_btag_loose);
  scale_factors['W'].push_back(sf_w);

  scale_factors['Z'].push_back(sf_ele_medium);
  scale_factors['Z'].push_back(sf_muon_medium);
  scale_factors['Z'].push_back(sf_w);

  // Top analysis
  scale_factors['t'].push_back(sf_ele_veto);
  scale_factors['t'].push_back(sf_muon_veto);
  scale_factors['t'].push_back(sf_btag_medium);
  scale_factors['t'].push_back(sf_top);

  // N-1 weights
  // Calculate weight for all search regions, but without a specific weight
  if (!isData) for (const auto& region : analysis_cuts) {
    size_t n=all_weights.size()+scale_factors[region.first].size();
    for (size_t i=0; i<n; ++i) {
      w_nm1[region.first][i] = 1;
      for (size_t j=0; j<n; ++j) if (j!=i) {
	if (j<all_weights.size()) w_nm1[region.first][i] *= all_weights[j];
	else w_nm1[region.first][i] *= scale_factors[region.first][j-all_weights.size()];
      }
    }
  }
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& data) {
  return 0;
  //return apply_all_cuts('S');
}

//_______________________________________________________
//      Define Histo options: Filling/Postfixes
std::vector<std::string> all_cuts;

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

void
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst,
			       const unsigned int& syst_index, bool runOnSkim=false)
{
  const int debug = 0;

  if (debug) std::cout<<"Analysis::define_histo_options: start"<<std::endl;

  sh.SetHistoWeights({ [&weight] { return weight; } });

  // Keep this to be able to use analysis cuts
  define_preselections(d);
  define_selections(d);
  if (debug) std::cout<<"Analysis::define_histo_options: weight, selections ok"<<std::endl;

  // Initialize containers for N-1 weights
  for (const auto& region : analysis_cuts) w_nm1[region.first].resize(20,1);

  if (debug) std::cout<<"Analysis::define_histo_options: set containers ok"<<std::endl;

  // --------------------------------------------------------------------
  //                            Colors
  // --------------------------------------------------------------------

  // Common Histo colorings
  // 400 kYellow  800 kOrange
  // 416 kGreen	  820 kSpring
  // 432 kCyan	  840 kTeal
  // 600 kBlue	  860 kAzure
  // 616 kMagenta 880 kViolet
  // 632 kRed     900 kPink

  std::string col3_red_to_blue = "633,618,601,"; // red, purple, blue
  std::string col4_red_to_cyan = "633,618,601,434,"; // Red, purple, blue, cyan
  std::string col4_cyan_to_red = "434,601,618,633,"; // Cyan, blue, purple, red
  std::string col5_green_to_red = "418,434,601,618,633,"; // green, cyan, blue, purple, red
  std::string col5_red_to_green = "633,618,601,434,418,"; // red, , purple, blue, cyan, green
  std::string col6_rainbow_dark = "601,434,418,402,633,618,"; // blue, cyan, green, yellow, red, purple
  std::string col8 = "1,601,434,418,402,807,633,618,"; // above plus black and orange
  std::string col10 = "4,6,2,800,402,417,433,9,618,633,";
  std::string col12 = "1,4,6,2,800,402,417,433,9,618,633,924,"; // black, blue, magenta, red, orange, darker yellow, darker green, darker cyan, blue-purple, dark purple, dark red
  std::string col12_rainbow = "402,416,433,600,617,632,802,813,833,863,883,892,"; // Go around first bright and then dark colors

  // --------------------------------------------------------------------
  //                            Cuts
  // --------------------------------------------------------------------

  // Pass each cut
  //sh.AddNewCut("PassAnaSelection",    [this] { return analysis_cuts['S'][cut_index].func(); });

  //__________________________________
  //            Postfixes
  //__________________________________

  // Postfixes are vector definitions for histograms
  // They attach _<string> after histogram names
  // where <string> is chosen from a vector of strings
  // you need to give a natural number as a vector index
  // for the histogram to choose which histo to fill

  // Notation:
  // AddNewPostfix("Name of Postfix", lambda function returning non-negative integer, "postfix 1;2;3; ...", "legend text 1;2;3; ...", "ROOT color 1,2,3, ...")



  // Sample postfixes
  // Determine them from the directory names in which the input files are
  // Map directory names to postfix name, legend entry and color
  std::vector<Sample> bkg_ttbars;
  bkg_ttbars.push_back({ .postfix="TTJets_madgraph_HT",       .legend="t#bar{t} (madgraphMLM, HT600toInf)", .color="634",/*DRed*/   .dirs={ 
			   "TTJets_HT-600to800", "TTJets_HT-800to1200", "TTJets_HT-1200to2500", "TTJets_HT-2500toInf" 
			 } });
  
  //bkg_ttbars.push_back({ .postfix="TTJets_madgraph_FullSim",  .legend="t#bar{t} (madgraphMLM, FullSim)", .color="901",/*Pink*/   .dirs={ "TTJets_madgraphMLM-pythia8" } });
  //bkg_ttbars.push_back({ .postfix="TTJets_madgraph_FastSim",  .legend="t#bar{t} (madgraphMLM, FastSim)", .color="903",/*DPink*/  .dirs={ "TTJets_madgraphMLM_FastSim" } });
  bkg_ttbars.push_back({ .postfix="TTJets_amcatnlo",          .legend="t#bar{t} (aMC@NLO FxFx)",         .color="617",/*Magent*/ .dirs={ "TTJets_amcatnloFXFX-pythia8" } });

  //bkg_ttbars.push_back({ .postfix="TT_amcatnlo",              .legend="t#bar{t} (aMC@NLO)",              .color="619",/*DMagen*/ .dirs={ "TT_amcatnlo-pythia8_ext1" } });
  bkg_ttbars.push_back({ .postfix="TT_powheg_pythia8",        .legend="t#bar{t}",                        .color="633",/*Red*/    .dirs={ "TT_powheg-pythia8", "TT_powheg-pythia8_backup" } });
  bkg_ttbars.push_back({ .postfix="TT_powheg_herwigpp",       .legend="t#bar{t} (powheg, herwigpp)",     .color="803",/*DOran*/  .dirs={ "TT_powheg-herwigpp", "TT_powheg-herwigpp_ext2", "TT_powheg-herwigpp_ext3" } });

  if (debug) std::cout<<"Analysis::define_histo_options: ok1"<<std::endl;
  std::vector<Sample> bkg_nonttbars;
  /*
    Signal colors:
    "12"  //DarkGrey
    "862" //Azure
    "841" //Teal
    "804" //DarkOrange
    "403" //DarkYellow

    Free colors:
    "435" //DarkCyan
  */

  if (debug) std::cout<<"Analysis::define_histo_options: ok2"<<std::endl;
  bkg_nonttbars.push_back({ .postfix="Multijet",   .legend="Multijet",                                .color="619",/*DMagen*/ .dirs={ 
			      "QCD_HT100to200",  "QCD_HT200to300",   "QCD_HT300to500",   "QCD_HT500to700",
			      "QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf",
			      "QCD_HT300to500_ext1", "QCD_HT500to700_ext1", "QCD_HT700to1000_ext1",
			      "QCD_HT1000to1500_ext1", "QCD_HT1500to2000_ext1", "QCD_HT2000toInf_ext1",
			      "ZJetsToQQ_HT600toInf", "WJetsToQQ_HT180", "DYJetsToQQ_HT180", // V Multi
			      "WWTo4Q", "ZZTo4Q" // VV Multi
			    } });
  bkg_nonttbars.push_back({ .postfix="WToLNu",     .legend="W(#rightarrowl#nu)",                      .color="418",/*Green*/  .dirs={ 
			      "WJetsToLNu_HT-70To100", "WJetsToLNu_HT-100To200", "WJetsToLNu_HT-200To400", "WJetsToLNu_HT-400To600",
			      "WJetsToLNu_HT-600To800", "WJetsToLNu_HT-800To1200", "WJetsToLNu_HT-1200To2500", "WJetsToLNu_HT-2500ToInf"
			    } });
  bkg_nonttbars.push_back({ .postfix="ZToNuNu",    .legend="Z(#rightarrow#nu#nu)",                    .color="401",/*Yellow*/ .dirs={ 
			      "ZJetsToNuNu_HT-100To200", "ZJetsToNuNu_HT-200To400", "ZJetsToNuNu_HT-400To600", "ZJetsToNuNu_HT-600To800", 
			      "ZJetsToNuNu_HT-800To1200", "ZJetsToNuNu_HT-1200To2500", "ZJetsToNuNu_HT-2500ToInf"
			    } });
  bkg_nonttbars.push_back({ .postfix="Multiboson", .legend="VV(V)+t#bar{t}X",                         .color="601",/*Blue*/   .dirs={
        		      "WWTo2L2Nu", "WWToLNuQQ", /*"WWToLNuQQ_ext1",*/
			      "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2L2Q", "WZTo2Q2Nu", "WZTo3LNu",
        		      "ZZTo2L2Nu", "ZZTo2L2Q", "ZZTo2Q2Nu", "ZZTo4L", 
        		      "WWW", "WWZ", "WZZ", "ZZZ"
  //      		    } });
  //bkg_nonttbars.push_back({ .postfix="TTX",        .legend="t#bar{t} + W/Z/#gamma/t#bar{t}",          .color="843",/*DTeal*/  .dirs={ 
			      "TTWJetsToLNu", "TTWJetsToQQ",
			      "TTZToLLNuNu", "TTZToQQ",
			      "TTGJets",
			      "TTTT"
			    } });
  bkg_nonttbars.push_back({ .postfix="Top",        .legend="Top",                                     .color="433",/*Cyan*/   .dirs={ 
			      "ST_s-channel_4f_leptonDecays", "ST_s-channel_4f_InclusiveDecays",
			      "ST_t-channel_top_4f_inclusiveDecays", "ST_t-channel_antitop_4f_inclusiveDecays",
			      "ST_tW_top_5f_inclusiveDecays",        "ST_tW_antitop_5f_inclusiveDecays"
			    } });
  bkg_nonttbars.push_back({ .postfix="DYToLL",     .legend="Drell-Yan",                               .color="803",/*Brown*/  .dirs={ 
			      "DYJetsToLL_M-50_HT-100to200",    "DYJetsToLL_M-50_HT-200to400",    "DYJetsToLL_M-50_HT-2500toInf",   "DYJetsToLL_M-50_HT-400to600",
			      "DYJetsToLL_M-50_HT-600to800",    "DYJetsToLL_M-50_HT-70to100",     "DYJetsToLL_M-50_HT-800to1200",   "DYJetsToLL_M-50_HT-1200to2500",
			      "DYJetsToLL_M-5to50_HT-100to200", "DYJetsToLL_M-5to50_HT-200to400", "DYJetsToLL_M-5to50_HT-400to600", "DYJetsToLL_M-5to50_HT-600toInf"
			    } });

  std::vector<Sample> bkg_all, bkg_selected, ttbar_selected;
  bkg_all.insert(bkg_all.end(), bkg_ttbars.begin(), bkg_ttbars.end());
  bkg_all.insert(bkg_all.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());
  ttbar_selected.push_back(bkg_ttbars[2]); // powheg
  bkg_selected  .push_back(bkg_ttbars[2]); // powheg
  bkg_selected.insert(bkg_selected.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());

  if (debug) std::cout<<"Analysis::define_histo_options: ok3"<<std::endl;
  std::vector<Sample> data_all, data_selected, single_ele, single_mu, met;
  data_all.push_back({ .postfix="Data",      .legend="Data",             .color="1", .dirs={
			 "JetHT_Run2016B_03Feb2017_v2", "JetHT_Run2016C_03Feb2017", "JetHT_Run2016D_03Feb2017", "JetHT_Run2016E_03Feb2017",
			 "JetHT_Run2016F_03Feb2017", "JetHT_Run2016G_03Feb2017", "JetHT_Run2016H_03Feb2017_v2", "JetHT_Run2016H_03Feb2017_v3",
			 "JetHT_Run2016B_RRv3", "JetHT_Run2016C_RRv1", "JetHT_Run2016D_RRv1", "JetHT_Run2016E_RRv1",
			 "JetHT_Run2016F_RRv1", "JetHT_Run2016G_RRv1", "JetHT_Run2016H_PRv2", "JetHT_Run2016H_PRv3",
			 "JetHT_Run2016C_RRv1_recovery", "JetHT_Run2016D_RRv1_recovery", "JetHT_Run2016F_RRv1_recovery",
			 "JetHT_Run2016G_RRv1_recovery", "JetHT_Run2016H_PRv2_recovery"
		       } });
  data_all.push_back({ .postfix="SingleEle", .legend="Data (SingleEle)", .color="1", .dirs={
			 "SingleElectron_Run2016B_03Feb2017_v2", "SingleElectron_Run2016C_03Feb2017", "SingleElectron_Run2016D_03Feb2017", "SingleElectron_Run2016E_03Feb2017",
			 "SingleElectron_Run2016F_03Feb2017", "SingleElectron_Run2016G_03Feb2017", "SingleElectron_Run2016H_03Feb2017_v2", "SingleElectron_Run2016H_03Feb2017_v3",
			 "SingleElectron_Run2016B_RRv3", "SingleElectron_Run2016C_RRv1", "SingleElectron_Run2016D_RRv1", "SingleElectron_Run2016E_RRv1",
			 "SingleElectron_Run2016F_RRv1", "SingleElectron_Run2016G_RRv1", "SingleElectron_Run2016H_PRv2", "SingleElectron_Run2016H_PRv3",
			 "SingleElectron_Run2016C_RRv1_recovery", "SingleElectron_Run2016D_RRv1_recovery", "SingleElectron_Run2016E_RRv1_recovery",
			 "SingleElectron_Run2016F_RRv1_recovery", "SingleElectron_Run2016G_RRv1_recovery", "SingleElectron_Run2016H_PRv2_recovery",
		       } });
  data_all.push_back({ .postfix="SingleMu",  .legend="Data (SingleMu)",  .color="1", .dirs={
			 "SingleMuon_Run2016B_03Feb2017_v2", "SingleMuon_Run2016C_03Feb2017", "SingleMuon_Run2016D_03Feb2017", "SingleMuon_Run2016E_03Feb2017",
			 "SingleMuon_Run2016F_03Feb2017", "SingleMuon_Run2016G_03Feb2017", "SingleMuon_Run2016H_03Feb2017_v2", "SingleMuon_Run2016H_03Feb2017_v3",
			 "SingleMuon_Run2016B_RRv3", "SingleMuon_Run2016C_RRv1", "SingleMuon_Run2016D_RRv1", "SingleMuon_Run2016E_RRv1",
			 "SingleMuon_Run2016F_RRv1", "SingleMuon_Run2016G_RRv1", "SingleMuon_Run2016H_PRv2", "SingleMuon_Run2016H_PRv3",
			 "SingleMuon_Run2016B_RRv3_recovery", "SingleMuon_Run2016C_RRv1_recovery", "SingleMuon_Run2016E_RRv1_recovery",
			 "SingleMuon_Run2016F_RRv1_recovery", "SingleMuon_Run2016G_RRv1_recovery", "SingleMuon_Run2016H_PRv2_recovery",
		       } });
  data_all.push_back({ .postfix="MET",       .legend="Data (MET)",       .color="1", .dirs={
			 "MET_Run2016B_03Feb2017_v2", "MET_Run2016C_03Feb2017", "MET_Run2016D_03Feb2017", "MET_Run2016E_03Feb2017",
			 "MET_Run2016F_03Feb2017", "MET_Run2016G_03Feb2017", "MET_Run2016H_03Feb2017_v2", "MET_Run2016H_03Feb2017_v3",
			 "MET_Run2016B_RRv3", "MET_Run2016C_RRv1", "MET_Run2016D_RRv1", "MET_Run2016E_RRv1",
			 "MET_Run2016F_RRv1", "MET_Run2016G_RRv1", "MET_Run2016H_PRv2", "MET_Run2016H_PRv3",
			 "MET_Run2016B_RRv3_recovery", "MET_Run2016D_RRv1_recovery", "MET_Run2016E_RRv1_recovery",
			 "MET_Run2016H_PRv2_recovery"
		       } });
  data_selected.push_back(data_all[0]);
  single_ele.push_back(data_all[1]);
  single_mu.push_back(data_all[2]);
  met.push_back(data_all[3]);

  if (debug) std::cout<<"Analysis::define_histo_options: ok4"<<std::endl;
  std::vector<Sample> signal_all, signal_selected, signal_fastsim, signal_gluino, signal_stop;
  std::vector<Sample> T5ttcc, T5tttt, T1tttt, T2tt;
  signal_all.push_back({ .postfix="T5ttcc",       .legend="T5ttcc",      .color="12", /*DGrey*/ .dirs={ "FastSim_SMS-T5ttcc", "FastSim_SMS-T5ttcc_mGluino1750to2300" } });
  signal_all.push_back({ .postfix="T5tttt",       .legend="T5tttt",      .color="862",/*Azure*/ .dirs={ "FastSim_SMS-T5tttt" } });
  signal_all.push_back({ .postfix="T1tttt",       .legend="T1tttt",      .color="841",/*Teal*/  .dirs={ "FastSim_SMS-T1tttt" } });
  signal_all.push_back({ .postfix="T2tt",         .legend="T2tt",        .color="403",/*DYell*/ .dirs={ 
			   "FastSim_SMS-T2tt_mStop-150to250", "FastSim_SMS-T2tt_mStop-250to350",
			   "FastSim_SMS-T2tt_mStop-350to400", "FastSim_SMS-T2tt_mStop-400to1200" 
			 } });
  //signal_all.push_back({ .postfix="T2tt_FullSim", .legend="T2tt (FullSim)", .color="804",/*DOran*/ .dirs={
  //      		   "FullSim_SMS-T2tt_mStop-425_mLSP-325", "FullSim_SMS-T2tt_mStop-500_mLSP-325",
  //      		   "FullSim_SMS-T2tt_mStop-850_mLSP-100" 
  //      		 } });
  signal_selected.push_back(signal_all[0]);
  for (int i=0; i<4; ++i) signal_fastsim.push_back(signal_all[i]);
  for (int i=0; i<3; ++i) signal_gluino .push_back(signal_all[i]);
  for (int i=3; i<4; ++i) signal_stop   .push_back(signal_all[i]);
  T5ttcc.push_back(signal_all[0]);
  T5tttt.push_back(signal_all[1]);
  T1tttt.push_back(signal_all[2]);
  T2tt  .push_back(signal_all[3]);

  //"T5ttttDeg (M_{#tilde{g}}=1TeV)","1",/*Black*/
  //"T1tttt (M_{#tilde{g}}=1.5TeV)", "862",/*Azure*/
  //"T1tttt (M_{#tilde{g}}=1.5TeV, M_{#tilde{#chi}^{0}}=100GeV)", "841",/*Teal*/     
  //"T1tttt (M_{#tilde{g}}=1.2TeV, M_{#tilde{#chi}^{0}}=800GeV)", "843",/*DarkTeal*/ 
  //"T5ttttDeg (M_{#tilde{g}}=1TeV, 2,3-body)",                   "12", /*DarkGrey*/ 
  //"T2ttDeg (M_{#tilde{t}}=350GeV)",                             "434",/*Cyan*/     

  // Sample postfixes
  if (debug) std::cout<<"Analysis::define_histo_options: ok5"<<std::endl;
  static const PostfixOptions all_samples_opt=get_pf_opts_({data_all, bkg_all, signal_all}, sample);
  sh.AddNewPostfix("AllSamples", [] { return all_samples_opt.index; }, all_samples_opt.postfixes, all_samples_opt.legends, all_samples_opt.colors);

  static const PostfixOptions plot_samples_opt=get_pf_opts_({data_selected, signal_fastsim, bkg_selected}, sample);
  sh.AddNewPostfix("StackPlot", [&d] { 
		     // Select gluino/stop mass to give ~1k events with 40 fb^-1 (xsec~=0.025pb)
		     if (plot_samples_opt.index>=1 && plot_samples_opt.index<=3) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1; // T1ttt/T5ttcc/T5tttt
		     } else if (plot_samples_opt.index==4) {
		       if (d.evt.SUSY_Stop_Mass  != 850 || d.evt.SUSY_LSP_Mass != 100) return (size_t)-1; // T2tt
		     }
		     return plot_samples_opt.index; 
		   }, plot_samples_opt.postfixes, plot_samples_opt.legends, plot_samples_opt.colors);

  std::vector<Sample> background;
  std::vector<std::string> background_dirs;
  for (auto bkg : bkg_selected) for (auto dir : bkg.dirs) background_dirs.push_back(dir);
  background.push_back({ .postfix="Background", .legend="Background", .color="1", .dirs=background_dirs });
  static const PostfixOptions background_opt = get_pf_opts_({background}, sample);
  sh.AddNewPostfix("Background",  [] { return background_opt.index; }, background_opt.postfixes, background_opt.legends, background_opt.colors);

  static const PostfixOptions gluino_signalscans_opt = get_pf_opts_({signal_gluino}, sample);
  sh.AddNewPostfix("GluinoSignalScans",  [] { return gluino_signalscans_opt.index; }, gluino_signalscans_opt.postfixes, gluino_signalscans_opt.legends, gluino_signalscans_opt.colors);

  static const PostfixOptions stop_signalscans_opt = get_pf_opts_({signal_stop}, sample);
  sh.AddNewPostfix("StopSignalScans",  [] { return stop_signalscans_opt.index; }, stop_signalscans_opt.postfixes, stop_signalscans_opt.legends, stop_signalscans_opt.colors);

  static const PostfixOptions Bkg_T5ttcc_opt=get_pf_opts_({background, T5ttcc}, sample);
  static const PostfixOptions Bkg_T5tttt_opt=get_pf_opts_({background, T5tttt}, sample);
  static const PostfixOptions Bkg_T1tttt_opt=get_pf_opts_({background, T1tttt}, sample);
  static const PostfixOptions Bkg_T2tt_opt  =get_pf_opts_({background, T2tt},   sample);

  static const PostfixOptions T5ttcc_opt = get_pf_opts_({T5ttcc}, sample);
  sh.AddNewPostfix("T5ttcc",  [] { return T5ttcc_opt.index; }, T5ttcc_opt.postfixes, T5ttcc_opt.legends, T5ttcc_opt.colors);
  static const PostfixOptions T5tttt_opt = get_pf_opts_({T5tttt}, sample);
  sh.AddNewPostfix("T5tttt",  [] { return T5tttt_opt.index; }, T5tttt_opt.postfixes, T5tttt_opt.legends, T5tttt_opt.colors);
  static const PostfixOptions T1tttt_opt = get_pf_opts_({T1tttt}, sample);
  sh.AddNewPostfix("T1tttt",  [] { return T1tttt_opt.index; }, T1tttt_opt.postfixes, T1tttt_opt.legends, T1tttt_opt.colors);
  static const PostfixOptions T2tt_opt = get_pf_opts_({T2tt}, sample);
  sh.AddNewPostfix("T2tt",  [] { return T2tt_opt.index; }, T2tt_opt.postfixes, T2tt_opt.legends, T2tt_opt.colors);

  static const PostfixOptions background_signal_opt = get_pf_opts_({background, signal_selected}, sample);
  sh.AddNewPostfix("Background_Signal", [&d] { 
		     if (background_signal_opt.index==1) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1;
		     }
		     return background_signal_opt.index;
		   }, background_signal_opt.postfixes, background_signal_opt.legends, "633,601");

  static const PostfixOptions signals_opt = get_pf_opts_({signal_all}, sample);
  sh.AddNewPostfix("Signals",  [&d] { 
		     // Select gluino/stop mass to give ~1k events with 40 fb^-1
		     if (signals_opt.index<3) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1; // T1ttt/T5ttcc/T5tttt
		     } else if (signals_opt.index==3) {
		       if (d.evt.SUSY_Stop_Mass  != 850 || d.evt.SUSY_LSP_Mass != 100) return (size_t)-1; // T2tt - Same as FullSim point
		     }
		     return signals_opt.index; 
		   }, signals_opt.postfixes, signals_opt.legends, signals_opt.colors);

  static const PostfixOptions signals_background_opt = get_pf_opts_({signal_all, background}, sample);
  sh.AddNewPostfix("Signals_Background",  [&d] { 
		     // Select gluino/stop mass to give ~1k events with 40 fb^-1
		     if (signals_background_opt.index<3) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1; // T1ttt/T5ttcc/T5tttt
		     } else if (signals_background_opt.index==3) {
		       if (d.evt.SUSY_Stop_Mass  != 850 || d.evt.SUSY_LSP_Mass != 100) return (size_t)-1; // T2tt - Same as FullSim point
		     }
		     return signals_background_opt.index; 
		   }, signals_background_opt.postfixes, signals_background_opt.legends, signals_background_opt.colors);

  static const PostfixOptions ttbar_signal_opt = get_pf_opts_({ttbar_selected, signal_selected}, sample);
  sh.AddNewPostfix("TT_Signal",  [&d] { 
		     if (ttbar_signal_opt.index==0) return (size_t)0;
		     else if (ttbar_signal_opt.index==1) {
		       if (d.evt.SUSY_LSP_Mass == 300 && d.evt.SUSY_Gluino_Mass==1400) return (size_t)1;
		     }
		     return (size_t)-1;
		   }, "TTbar;T5ttcc_Mlsp300_Mglu1400", "t#bar{t};T5ttcc M_{#tilde{g}}=1.4TeV", "1,633");
  sh.AddNewPostfix("TT_SignalPoints",  [&d] { 
		     if (ttbar_signal_opt.index==0) return (size_t)0;
		     else if (ttbar_signal_opt.index==1) {
		       if (d.evt.SUSY_LSP_Mass == 300) {
			 if      (d.evt.SUSY_Gluino_Mass== 900) return (size_t)1;
			 else if (d.evt.SUSY_Gluino_Mass==1100) return (size_t)2;
			 else if (d.evt.SUSY_Gluino_Mass==1300) return (size_t)3;
			 else if (d.evt.SUSY_Gluino_Mass==1500) return (size_t)4;
			 else if (d.evt.SUSY_Gluino_Mass==1700) return (size_t)5;
		       }
		     }
		     return (size_t)-1;
		   }, "TTbar;T5ttcc_Mlsp300_Mglu[900to1700++200]", "t#bar{t};T5ttcc M_{#tilde{g}}=[0.9to1.7++0.2]TeV", "1,"+col5_green_to_red);

  if (debug) std::cout<<"Analysis::define_histo_options: ok6"<<std::endl;
  static const PostfixOptions mgluinopoints_opt = get_pf_opts_({signal_gluino}, sample);
  sh.AddNewPostfix("MGluinoPoints",  [&d] { 
		     if (mgluinopoints_opt.index==(size_t)-1) return (size_t)-1;
		     else {
		       if (d.evt.SUSY_LSP_Mass != 300) return (size_t)-1;
		       if      (d.evt.SUSY_Gluino_Mass== 900) return (size_t)0;
		       else if (d.evt.SUSY_Gluino_Mass==1100) return (size_t)1;
		       else if (d.evt.SUSY_Gluino_Mass==1300) return (size_t)2;
		       else if (d.evt.SUSY_Gluino_Mass==1500) return (size_t)3;
		       else if (d.evt.SUSY_Gluino_Mass==1700) return (size_t)4;
		       else return (size_t)-1;
		     }
		   }, "Mlsp300_Mglu[900to1700++200]", "M_{#tilde{#chi}^{0}}=300GeV, M_{#tilde{g}}=[0.9to1.7++0.2]TeV", col5_green_to_red);

  static const PostfixOptions mstoppoints_opt = get_pf_opts_({signal_stop}, sample);
  sh.AddNewPostfix("MStopPoints",  [&d] { 
		     if (mstoppoints_opt.index==(size_t)-1) return (size_t)-1;
		     else {
		       if (d.evt.SUSY_LSP_Mass != 100) return (size_t)-1;
		       if      (d.evt.SUSY_Stop_Mass== 600) return (size_t)0;
		       else if (d.evt.SUSY_Stop_Mass== 800) return (size_t)1;
		       else if (d.evt.SUSY_Stop_Mass==1000) return (size_t)2;
		       else if (d.evt.SUSY_Stop_Mass==1200) return (size_t)3;
		       else return (size_t)-1;
		     }
		   }, "Mlsp100_Mstop[600to1200++200]", "M_{#tilde{#chi}^{0}}=100GeV, M_{#tilde{t}}=[0.6to1.2++0.2]TeV", col4_cyan_to_red);

  if (debug) std::cout<<"Analysis::define_histo_options: ok7"<<std::endl;
  static const PostfixOptions data_mc_opt = get_pf_opts_({data_selected, background}, sample);
  sh.AddNewPostfix("Data_MC",  [] { return data_mc_opt.index; }, data_mc_opt.postfixes, data_mc_opt.legends, "1,633");

  static const PostfixOptions single_lep_opt = get_pf_opts_({single_ele, single_mu}, sample);
  sh.AddNewPostfix("SingleEle_SingleMu", [] { return single_lep_opt.index; }, single_lep_opt.postfixes, single_lep_opt.legends, "1,633");

  static const PostfixOptions triggers_opt = get_pf_opts_({data_selected, single_ele, met, single_mu, background}, sample);
  sh.AddNewPostfix("Triggers", [&d]()
		   {
		     if (triggers_opt.index==0) {
		       //bool Pass_any_PFHT = //(d.hlt.AK8PFJet450==1) ||
		       //  (d.hlt.PFHT400==1) || (d.hlt.PFHT475==1) || (d.hlt.PFHT600==1) || 
		       //  (d.hlt.PFHT650==1) || (d.hlt.PFHT800==1) ||
		       //  (d.hlt.PFHT550_4JetPt50==1) ||(d.hlt.PFHT650_4JetPt50==1) || (d.hlt.PFHT750_4JetPt50==1);
		       //// JetHT Dataset: Pass any low threshold HT Trigger
		       //if (Pass_any_PFHT) return (size_t)0;
		       return (size_t)0;
		     } else if (triggers_opt.index==1) { // SingleElectron
		       if ((d.hlt.Ele23_WPLoose_Gsf==1||d.hlt.Ele27_WPTight_Gsf==1)&&nEleTight==1) return (size_t)1;
		     } else if (triggers_opt.index==2) { // MET
		       if (d.hlt.PFMET120_PFMHT120_IDTight==1&&d.met.Pt[0]>200&&nLepVeto==0&&d.evt.NIsoTrk==0) return (size_t)1;
		     } else if (triggers_opt.index==3) { // SingleMuon
		       if ((d.hlt.IsoMu24==1||d.hlt.IsoTkMu24==1)&&nMuTight==1) return (size_t)2;
		     } else if (triggers_opt.index==4) { // Background MC
		       return (size_t)3;
		     }
		     return (size_t)-1; 
		   //}, "JetHT;Ele23or27;IsoMu24;MET120;MC", "JetHT (All events);Ele23/27, 1 Tight Ele;IsoMu24, 1 Tight muon;PFMET120_PFMHT120, MET>200, lep veto;Simulation", "1,417,601,618,633");
		   }, "JetHT;SingleEle_MET;SingleMu;MC", "JetHT (All events);SingleEle+MET;SingleMu;Simulation", "1,417,601,633");
  static const PostfixOptions trigger_opt = get_pf_opts_({single_ele, met}, sample);
  sh.AddNewPostfix("Trigger", [&d]()
		   {
		     if (trigger_opt.index==0) { // SingleElectron
		       if ((d.hlt.Ele23_WPLoose_Gsf==1||d.hlt.Ele27_WPTight_Gsf==1)&&nEleTight==1) return (size_t)0;
		     } else if (trigger_opt.index==1) { // MET
		       if (d.hlt.PFMET120_PFMHT120_IDTight==1&&d.met.Pt[0]>200&&nLepVeto==0&&d.evt.NIsoTrk==0) return (size_t)0;
		     }
		     return (size_t)-1; 
		   }, "SingleEle_MET", "SingleEle + MET", "1");
  //sh.AddNewPostfix("Triggers", [&d] { return triggers_opt.index; },
  //      	   "JetHT;SingleEle;SingleMu;MET;MC", "Data: JetHT;Data: SingleEle;Data: SingleMu;Data: MET;MC: t#bar{t}", "1,417,601,618,633");
  
  // Systematics Postfixes
  sh.AddNewPostfix("Syst", [&syst_index] { return syst_index; }, std::string(";Syst[0to")+std::to_string(syst_nSyst)+"]", std::string(";systematics [0to")+std::to_string(syst_nSyst)+"]", "1-999");
  if (syst_nSyst>998) utils::error("Error: Too large number of systematics, define more colors!");
  if (debug) std::cout<<"Analysis::define_histo_options: sample postfixes ok"<<std::endl;

  // Cut names
  std::map<std::string, std::string> legname;
  legname["3Jet"]        = "Njet#geq3";
  legname["MR_R2"]       = "MR, R^{2}";
  legname["HLT"]         = "HLT";
  legname["0Ele"]        = "ele veto";
  legname["0Mu"]         = "muon veto";
  legname["0IsoTrk"]     = "isol trk veto";
  legname["1b"]          = "Nb#geq1";
  legname["1W"]          = "NW#geq1";
  legname["mDPhi"]       = "#Delta#phi";
  legname["InvmDPhi"]    = "inv. #Delta#phi";
  legname["0b"]          = "b-tag veto";
  legname["1aW"]         = "NW(anti-tag)#geq1";
  legname["InvmDPhi0p3"] = "#Delta#phi<0.3";
  legname["1Lep"]        = "Nlep=1";
  legname["MT"]          = "m_{T}";
  legname["1mW"]         = "NW(mass-tag)#geq1";
  legname["MR_R2ll"]     = "MR, R^{2}";
  legname["2Lep"]        = "Nlep=2";
  legname["OppCharge"]   = "#sumq_{lep}=0";
  legname["mDPhill"]     = "#Delta#phi";
  legname["Mll"]         = "|m_{ll} - m_{Z}| < 10 GeV";
  legname["1Top"]        = "Ntop#geq1";
  std::map<char, std::string> regionname;
  regionname['S'] = "Signal region";
  regionname['s'] = "S' region";
  regionname['Q'] = "QCD enriched region";
  regionname['q'] = "Q' region";
  regionname['T'] = "Top enriched region";
  regionname['W'] = "W enriched region";
  regionname['Z'] = "Z enriched region";
  regionname['t'] = "Boosted top region";

  // Cut Postfixes
  sh.AddNewPostfix("BaselineCuts", [] { return 0; }, "BaselineCuts", "Baseline cuts", "1");
  all_cuts.push_back("BaselineCuts");
  for (const auto& region : analysis_cuts) {
    std::string cutflow_str = "";
    sh.AddNewPostfix(std::string(1,region.first), [this,region] { return apply_all_cuts(region.first) ? 0 : (size_t)-1; },
		     std::string(1,region.first), regionname[region.first], "1");
    for (size_t i=0, n=region.second.size(); i<n; ++i) {


      // Cuts in order 1-N: "PassNCuts[search region]"
      sh.AddNewPostfix(std::string(1,region.first)+"_"+std::to_string(i+1)+"Cuts", [this,i,region] { return apply_ncut(region.first, i) ? 0 : (size_t)-1; },
		       std::string(1,region.first)+"_"+std::to_string(i+1)+"Cuts", std::string(1,region.first)+" region, first "+std::to_string(i+1)+" cuts", "1");
      all_cuts.push_back(std::string(1,region.first)+"_"+std::to_string(i+1)+"Cuts");
      cutflow_str += region.second[i].name+std::string(1,region.first)+";";
      // N-1 Cuts: "[search region]_Excl[cut]"
      sh.AddNewPostfix(std::string(1,region.first)+"_Excl"+region.second[i].name, [this,i,region] { 
			 unsigned int mask = (1<<region.second.size())-1 - (1<<i); 
			 return ((cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
		       std::string(1,region.first)+"_Excl"+region.second[i].name, regionname[region.first]+", no "+legname[region.second[i].name]+" cut", "1");
      // N-2 Cuts: "[search region]_Excl[cut1][cut2]"
      for (size_t j=i+1, n=region.second.size(); j<n; ++j)
	sh.AddNewPostfix(std::string(1,region.first)+"_Excl"+region.second[i].name+region.second[j].name, [this,i,j,region] { 
			   unsigned int mask = (1<<region.second.size())-1 - (1<<i) - (1<<j); 
			   return ((cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
		         std::string(1,region.first)+"_Excl"+region.second[i].name+region.second[j].name, regionname[region.first]+", no "+legname[region.second[i].name]+", "+legname[region.second[j].name]+" cut", "1");
    }
    // Stackable Cut Histos: "CutFlow"
    sh.AddNewPostfix("CutFlow"+std::string(1,region.first), [this,region] { for (size_t i=0, n=region.second.size(); i<n; ++i) if (!region.second[i].func()) return i; return region.second.size(); }, 
		     cutflow_str+"PassAll"+std::string(1,region.first), cutflow_str+regionname[region.first], col10+col10);
  }
  sh.AddNewPostfix("TriggerPreSelection",  [this] { return apply_cuts('W', {W_3Jet, W_MR_R2})?0:(size_t)-1; }, "TriggerPreSelection", "Preselection", "1");
  sh.AddNewPostfix("TriggerPreSelPlus1mW", [this] { return apply_cuts('W', {W_3Jet, W_MR_R2, W_1mW})?0:(size_t)-1; }, "TriggerPreSelPlus1mW", "Preselection + 1mW", "1");
  //sh.AddNewPostfix("AllLepIsolated",          [this] { return allVetoLepIsolated;                             }, "NotAllLepIsoLated;AllLepIsolated", "Not all lepton isolated;All lepton isolated", "633,418");
  //sh.AddNewPostfix("LepInsideJet",            [this] { return isLepInsideJet;                                 }, "LeptonOutsideJet;LeptonInsideJet", "Lepton not inside jet;Lepton inside jet", "633,418");
  
  // Individual Cuts implemented as Postfixes
  // Triggers
  sh.AddNewPostfix("JetHT",      [this,&d] { 
		     if (sample.find("SingleElectron")!=std::string::npos) return (size_t)-1;
		     else if (sample.find("SingleMuon")!=std::string::npos) return (size_t)-1;
		     else if (sample.find("MET")!=std::string::npos) return (size_t)-1;
		     //else if (sample.find("JetHT")!=std::string::npos) return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1 ? 0 : (size_t)-1;
		     else if (sample.find("JetHT")!=std::string::npos) return (size_t)0;
		     else return (size_t)0;
		   }, "JetHT",      "", "1");
  sh.AddNewPostfix("Blind",      [this,&d] { 
		     if (sample.find("SingleElectron")!=std::string::npos) return (size_t)-1;
		     else if (sample.find("SingleMuon")!=std::string::npos) return (size_t)-1;
		     else if (sample.find("MET")!=std::string::npos) return (size_t)-1;
		     else if (sample.find("JetHT")!=std::string::npos) return (size_t)-1;
		     else return (size_t)0;
		   }, "BlindData", "", "1");
  sh.AddNewPostfix("PFHT475",         [&d] { if (d.hlt.PFHT475==-9999) return (size_t)-1; else return (size_t)d.hlt.PFHT475; }, "NoPassHLT_PFHT475;PassHLT_PFHT475", "Do not pass HLT_PFHT475;Pass HLT_PFHT475", "633;418");

  // AK4 Jet Postfixes
  sh.AddNewPostfix("Jets",    [&d] {  size_t i=itJet[d.jetsAK4.it];        return (i<4)?i:(size_t)-1; }, "Jet[1to5]",  "1st Jet;2nd Jet;3rd Jet;[4to5]th Jet", col5_red_to_green);
  sh.AddNewPostfix("BTags",   [&d] {  size_t i=itMediumBTag[d.jetsAK4.it]; return (i<4)?i:(size_t)-1; }, "BTag[1to5]", "1st b;2nd b;3rd b;[4to5]th b",         col5_red_to_green);



  // AK8 Jet Postfixes
  sh.AddNewPostfix("JetsAK8",  [&d] {  size_t i=itJetAK8[d.jetsAK8.it];    return (i<4)?i:(size_t)-1; }, "Jet[1to4]",     "1st Jet;2nd Jet;3rd Jet;4th Jet",                     col4_red_to_cyan);
  sh.AddNewPostfix("mWs",      [&d] {  size_t i=itWPreTag[d.jetsAK8.it];   return (i<4)?i:(size_t)-1; }, "mW[1to4]", "1st W-masstag;2nd W-masstag;3rd W-masstag;4th W-masstag", col4_red_to_cyan);
  sh.AddNewPostfix("aWs",      [&d] {  size_t i=itTightWAntiTag[d.jetsAK8.it];   return (i<4)?i:(size_t)-1; }, "aW[1to4]", "1st W-antitag;2nd W-antitag;3rd W-antitag;4th W-antitag", col4_red_to_cyan);
  sh.AddNewPostfix("Ws",       [&d] {  size_t i=itTightWTag[d.jetsAK8.it]; return (i<4)?i:(size_t)-1; }, "W[1to4]",       "1st W;2nd W;3rd W;4th W",                             col4_red_to_cyan);
  sh.AddNewPostfix("Jet1AK8Pt450",  [&d] {  return d.jetsAK8.Pt[iJetAK8[0]]>450 ? 0 : (size_t)-1; }, "Jet1AK8_Pt450",  "1st jet p_{T} (AK8) > 450", "1");
  sh.AddNewPostfix("Jet1AK8Pt500",  [&d] {  return d.jetsAK8.Pt[iJetAK8[0]]>500 ? 0 : (size_t)-1; }, "Jet1AK8_Pt500",  "1st jet p_{T} (AK8) > 500", "1");
  sh.AddNewPostfix("Jet1AK8Mass65", [&d] {  return softDropMassW[iJetAK8[0]]>65 ? 0 : (size_t)-1; }, "Jet1AK8_Mass65", "1st jet M_{SD} (AK8) > 65", "1");
  sh.AddNewPostfix("Tau21Tagged",   [&d] {  return tau21[d.jetsAK8.it]<W_TAU21_TIGHT_CUT; }, "Tau21AntiTag;Tau21Tag", "#tau_{2}/#tau_{1} anti-tagged;#tau_{2}/#tau_{1} tagged", "633,418");

  // Event
  sh.AddNewPostfix("RBins",          [&d] { return (size_t)((d.evt.R>=0.1)+(d.evt.R>=0.2)+(d.evt.R>=0.4)); }, "R0to0p1;R0p1to0p2;R0p2to0p4;R0p4", "0.0<R<0.1;0.1<R<0.2;0.2<R<0.4;R>=0.4", "1,4,418,633");
  sh.AddNewPostfix("OtherUnisoLep",  [] { return std::min(nLepVetoNoIso-nLepSelect,(unsigned int)1); }, "NoOtherUnisoLep;OtherUnisoLep", "0 other unisol. lepton;#geq1 other unisol. lepton", "418,633");
  sh.AddNewPostfix("OtherLooseLep",  [] { return std::min(nLepVeto     -nLepSelect,(unsigned int)1); }, "NoOtherLep;OtherLep",           "0 other loose lepton;#geq1 other loose lepton", "633,418");
  sh.AddNewPostfix("R2Bins",         [&d] { return (size_t)((d.evt.R2>=0.08)+(d.evt.R2>=0.12)+(d.evt.R2>=0.16)+(d.evt.R2>=0.24)+(d.evt.R2>=0.5)); }, 
		   "R2_0to0p08;R2_0p08to0p12;R2_0p12to0p16;R2_0p16to0p24;R2_0p24to0p5;R2_0p5", 
		   "R^{2}#in[0,0.08[;R^{2}#in[0.08,0.12[;R^{2}#in[0.12,0.16[;R^{2}#in[0.16,0.24[;R^{2}#in[0.24,0.5[;R^{2}#in[0.5,1[", col6_rainbow_dark);
  sh.AddNewPostfix("R2llBins",       [] { return (size_t)((R2_ll>=0.08)+(R2_ll>=0.12)+(R2_ll>=0.16)+(R2_ll>=0.24)+(R2_ll>=0.5)); }, 
		   "R2ll_0to0p08;R2ll_0p08to0p12;R2ll_0p12to0p16;R2ll_0p16to0p24;R2ll_0p24to0p5;R2ll_0p5", 
		   "R_{ll}^{2}#in[0,0.08[;R_{ll}^{2}#in[0.08,0.12[;R_{ll}^{2}#in[0.12,0.16[;R_{ll}^{2}#in[0.16,0.24[;R_{ll}^{2}#in[0.24,0.5[;R_{ll}^{2}#in[0.5,1[", col6_rainbow_dark);
  sh.AddNewPostfix("Ele_Muon",       [] {  return (size_t)(nEleSelect==1 ? 0 : nMuSelect==1 ? 1 : -1); }, "EleOnly;MuOnly", "1 ele;1 muon", "1,2");
  sh.AddNewPostfix("2Ele_2Muon",     [] {  return (size_t)(nEleLoose==2 ? 0 : nMuLoose==2 ? 1 : -1); }, "EleOnly;MuOnly", "2 ele;2 muon", "1,2");
  sh.AddNewPostfix("NJet35",         [] {  return (size_t)(nJet<3 ? -1 : nJet>5); }, "NJet3to5;NJet6", "3 #leq N_{jet} #leq 5 ;6 #leq N_{jet}", "1,2");
  if (debug) std::cout<<"Analysis::define_histo_options: postfixes ok"<<std::endl;

  // Weights
  sh.AddNewPostfix("NoPUWeight",     [&d] { return 0; }, "NoPUWeight",   "No pile-up reweighting", "1");
  sh.AddNewPostfix("NoTrigWeight",   [&d] { return 0; }, "NoTrigWeight", "No trigger weighting",   "1");
  sh.AddNewPostfix("NoEleSF",        [&d] { return 0; }, "NoEleSF",      "No ele SF",              "1");
  sh.AddNewPostfix("NoMuonSF",       [&d] { return 0; }, "NoMuonSF",     "No muon SF",             "1");
  sh.AddNewPostfix("NoBTagSF",       [&d] { return 0; }, "NoBTagSF",     "No b-tag SF",            "1");
  sh.AddNewPostfix("NoWTagSF",       [&d] { return 0; }, "NoWTagSF",     "No W-tag SF",            "1");

  // Bins
  std::vector<double> E   = {0, 100, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 5000, 10000};
  std::vector<double> Pt  = {0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 1400, 2000, 3000, 4000, 5000, 10000};
  std::vector<double> PtF = {0, 200, 300, 400, 600, 1000, 2000, 5000};
  std::vector<double> M   = {0, 10, 20, 30, 40, 50, 65, 75, 85, 95, 105, 120, 135, 150, 165, 180, 195, 210, 230, 260, 300, 500, 1000};
  std::vector<double> MW  = {65, 75, 85, 95, 105};
  std::vector<double> CSV = {0, 0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 };
  std::vector<double> MDP;
  for (double x=0.0; x< 1.8; x+=0.1) MDP.push_back(x);
  for (double x=1.8; x< 2.4; x+=0.2) MDP.push_back(x);
  for (double x=2.4; x<=3.2; x+=0.4) MDP.push_back(x);
  std::vector<double> NVTX(1,0);
  for (double x=6;  x<  40;  x+=2) NVTX.push_back(x);
  for (double x=40; x<=100;  x+=5) NVTX.push_back(x);
  std::vector<double> R   = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.7, 0.8, 1.0, 1.2, 2.0 };
  std::vector<double> MR  = {0, 600, 800, 1000, 1200, 1600, 2000, 4000, 10000};
  std::vector<double> MTR = {0, 100, 200, 300, 400, 600, 800, 1000, 1200, 1600, 2000, 4000};
  std::vector<double> MET = {0, 100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000};
  std::vector<double> R2  = {0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.5, 1.0, 5.0};
  std::vector<double> HT  = {0, 200, 300, 400, 500, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1200, 1500, 2000, 2500, 3000, 4000, 10000};
  std::vector<double> HTB = {400, 500, 600, 700, 750, 800, 850, 900, 950, 1000, 1500, 10000}; // 2D Trigger Eff
  std::vector<double> PtB = {200, 300, 400, 450, 500, 550, 600, 1000, 10000}; // 2D Trigger Eff

  // Bin Postfixes
  std::stringstream HT_pf, HT_leg;
  for (size_t i=0, n=HTB.size(); i<n-1; ++i) {
    HT_pf<<"HT"<<HTB[i]<<"to"<<HTB[i+1];
    HT_leg<<"H_{T} #subset ["<<HTB[i]<<","<<HTB[i+1]<<"[";
    if (i!=n-2) { HT_pf<<";"; HT_leg<<";"; }
  }
  sh.AddNewPostfix("HTBins", [HTB] { for (size_t i=0, n=HTB.size(); i<n-1; ++i) if (AK4_Ht>=HTB[i]&&AK4_Ht<HTB[i+1]) return i; return (size_t)-1; },
		   HT_pf.str(), HT_leg.str(), col12+col12);
  std::stringstream AK8Pt_pf, AK8Pt_leg;
  for (size_t i=0, n=PtB.size(); i<n-1; ++i) {
    AK8Pt_pf<<"Jet1AK8Pt"<<PtB[i]<<"to"<<PtB[i+1];
    AK8Pt_leg<<"AK8 jet1 p_{T} #subset ["<<PtB[i]<<","<<PtB[i+1]<<"[";
    if (i!=n-2) { AK8Pt_pf<<";"; AK8Pt_leg<<";"; }
  }
  sh.AddNewPostfix("Jet1AK8PtBins", [PtB,&d] { if (nJetAK8<1) return (size_t)-1;
		     for (size_t i=0, n=PtB.size(); i<n-1; ++i) if (d.jetsAK8.Pt[iJetAK8[0]]>=PtB[i]&&d.jetsAK8.Pt[iJetAK8[0]]<PtB[i+1]) return i; 
		     return (size_t)-1; }, AK8Pt_pf.str(), AK8Pt_leg.str(), col8+col8);


  // --------------------------------------------------------------------
  //                         Fill Parameters
  // --------------------------------------------------------------------

  sh.AddNewFillParam("Bin", { .nbin=1,   .bins={0,1}, .fill=[&d] { return 0; }, .axis_title="Bin"}); // For averages/counts

  // AK4 Jets
  /* Variables for ID
     sh.AddNewFillParam("JetPhotonE",           { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.PhotonEnergy[d.jetsAK4.it];              }, .axis_title="Jet Photon Energy (GeV)"});
     sh.AddNewFillParam("JetElectronE",         { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.ElectronEnergy[d.jetsAK4.it];            }, .axis_title="Jet Electron Energy (GeV)"});
     sh.AddNewFillParam("JetMuonE",             { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.MuonEnergy[d.jetsAK4.it];                }, .axis_title="Jet Muon Energy (GeV)"});
     sh.AddNewFillParam("JetChargedMuE",        { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.ChargeMuEnergy[d.jetsAK4.it];            }, .axis_title="Jet Charged Mu Energy (GeV)"});
     sh.AddNewFillParam("JetChargedEmE",        { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.chargedEmEnergy[d.jetsAK4.it];           }, .axis_title="Jet Charged Em Energy (GeV)"});
     sh.AddNewFillParam("JetChargedHadronE",    { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.chargedHadronEnergy[d.jetsAK4.it];       }, .axis_title="Jet Charged Hadron Energy (GeV)"});
     sh.AddNewFillParam("JetNeutralHadronE",    { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.neutralHadronEnergy[d.jetsAK4.it];       }, .axis_title="Jet Neutral Hadron Energy (GeV)"});
     sh.AddNewFillParam("JetNeutralEmE",        { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.neutralEmEnergy[d.jetsAK4.it];           }, .axis_title="Jet Neutral Em Energy (GeV)"});
     sh.AddNewFillParam("JetHFHadronE",         { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.HFHadronEnergy[d.jetsAK4.it];            }, .axis_title="Jet HF Hadron Energy (GeV)"});
     sh.AddNewFillParam("JetHFEME",             { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4.HFEMEnergy[d.jetsAK4.it];                }, .axis_title="Jet HF EM Energy (GeV)"});
     sh.AddNewFillParam("JetNumberOfDaughters", { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4.numberOfDaughters[d.jetsAK4.it];         }, .axis_title="Jet Number Of Daughters"});
     sh.AddNewFillParam("JetPhotonMult",        { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4.photonMultiplicity[d.jetsAK4.it];        }, .axis_title="Jet Photon Multiplicity"});
     sh.AddNewFillParam("JetElectronMult",      { .nbin=  10, .bins={     0,     10}, .fill=[&d] { return d.jetsAK4.electronMultiplicity[d.jetsAK4.it];      }, .axis_title="Jet Electron Multiplicity"});
     sh.AddNewFillParam("JetMuonMult",          { .nbin=  10, .bins={     0,     10}, .fill=[&d] { return d.jetsAK4.muonMultiplicity[d.jetsAK4.it];          }, .axis_title="Jet Muon Multiplicity"});
     sh.AddNewFillParam("JetNeutralMult",       { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4.neutralMultiplicity[d.jetsAK4.it];       }, .axis_title="Jet Neutral Multiplicity"});
     sh.AddNewFillParam("JetChargedMult",       { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4.chargedMultiplicity[d.jetsAK4.it];       }, .axis_title="Jet Charged Multiplicity"});
     sh.AddNewFillParam("JetChargedHadronMult", { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4.ChargedHadronMultiplicity[d.jetsAK4.it]; }, .axis_title="Jet Charged Hadron Multiplicity"});
     sh.AddNewFillParam("JetNeutralHadronMult", { .nbin=  40, .bins={     0,     40}, .fill=[&d] { return d.jetsAK4.neutralHadronMultiplicity[d.jetsAK4.it]; }, .axis_title="Jet Neutral Hadron Multiplicity"});
     sh.AddNewFillParam("JetHFHadronMult",      { .nbin=  80, .bins={     0,     80}, .fill=[&d] { return d.jetsAK4.HFHadronMultiplicity[d.jetsAK4.it];      }, .axis_title="Jet HF Hadron Multiplicity"});
     sh.AddNewFillParam("JetHFEMMult",          { .nbin=  50, .bins={     0,     50}, .fill=[&d] { return d.jetsAK4.HFEMMultiplicity[d.jetsAK4.it];          }, .axis_title="Jet HF EM Multiplicity"});
  */
  sh.AddNewFillParam("JetPtBins",            { .nbin=Pt.size()-1,   .bins=Pt,      .fill=[&d] { return d.jetsAK4.Pt[d.jetsAK4.it];           }, .axis_title="Jet p_{T} (GeV)", .def_range={200,2000} });
  sh.AddNewFillParam("JetPtFewBins",         { .nbin=PtF.size()-1,  .bins=PtF,     .fill=[&d] { return d.jetsAK4.Pt[d.jetsAK4.it];           }, .axis_title="Jet p_{T} (GeV)", .def_range={200,2000} });
  sh.AddNewFillParam("JetPtOneBin",          { .nbin=   1, .bins={     0,   5000}, .fill=[&d] { return d.jetsAK4.Pt[d.jetsAK4.it];           }, .axis_title="Jet p_{T} (GeV)"});
  sh.AddNewFillParam("JetPt",                { .nbin= 200, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK4.Pt[d.jetsAK4.it];           }, .axis_title="Jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("JetEta",               { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK4.Eta[d.jetsAK4.it];          }, .axis_title="Jet #eta",        .def_range={-2.4,2.4}});
  sh.AddNewFillParam("JetPhi",               { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK4.Phi[d.jetsAK4.it];          }, .axis_title="Jet #phi"});
  sh.AddNewFillParam("JetCSV",               { .nbin=  20, .bins={     0,   1.00}, .fill=[&d] { return std::min(d.jetsAK4.CSVv2[d.jetsAK4.it],float(0.999)); }, .axis_title="Jet CSV"});
  // BJets
  sh.AddNewFillParam("BJetPtBins",           { .nbin=Pt.size()-1,   .bins=Pt,      .fill=[&d] { return d.jetsAK4.Pt[d.jetsAK4.it];           }, .axis_title="B-jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("BJetPt",               { .nbin= 200, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK4.Pt[d.jetsAK4.it];           }, .axis_title="B-jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("BJetEta",              { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK4.Eta[d.jetsAK4.it];          }, .axis_title="B-jet #eta",        .def_range={-2.4,2.4}});
  sh.AddNewFillParam("BJetPhi",              { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK4.Phi[d.jetsAK4.it];          }, .axis_title="B-jet #phi"});
  sh.AddNewFillParam("BJetCSV",              { .nbin=  20, .bins={     0,   1.00}, .fill=[&d] { return std::min(d.jetsAK4.CSVv2[d.jetsAK4.it],float(0.999)); }, .axis_title="B-jet CSV"});

  // AK8 Jets
  sh.AddNewFillParam("JetAK8PtOneBin",       { .nbin=   1, .bins={   200,   5000}, .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="AK8 jet p_{T} (GeV)"});
  sh.AddNewFillParam("JetAK8PtFewBins",      { .nbin=PtF.size()-1, .bins=PtF,      .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="AK8 jet p_{T} (GeV)", .def_range={200,2000} });
  sh.AddNewFillParam("JetAK8PtBins",         { .nbin=Pt.size()-1 , .bins=Pt,       .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="AK8 jet p_{T} (GeV)", .def_range={200,2000} });
  sh.AddNewFillParam("JetAK8Pt",             { .nbin= 200, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="AK8 jet p_{T} (GeV)", .def_range={200,2000} });
  sh.AddNewFillParam("JetAK8Eta",            { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8.Eta[d.jetsAK8.it];          }, .axis_title="AK8 jet #eta",        .def_range={-2.4,2.4}});
  sh.AddNewFillParam("JetAK8Phi",            { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8.Phi[d.jetsAK8.it];          }, .axis_title="AK8 jet #phi"});
  sh.AddNewFillParam("JetAK8Mass",           { .nbin= 200, .bins={     0,   2000}, .fill=[&d] { return softDropMassW[d.jetsAK8.it];          }, .axis_title="AK8 jet soft-drop mass (GeV)", .def_range={0,400}});
  sh.AddNewFillParam("JetAK8MassTop",        { .nbin= 200, .bins={     0,   2000}, .fill=[&d] { return softDropMassTop[d.jetsAK8.it];        }, .axis_title="AK8 jet soft-drop mass (GeV)", .def_range={0,400}});
  /*
    sh.AddNewFillParam("JetAK8PrunedMass",     { .nbin= 400, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8.prunedMass[d.jetsAK8.it];   }, .axis_title="AK8 jet pruned mass (GeV)",    .def_range={0,150}});
    sh.AddNewFillParam("JetAK8FilteredMass",   { .nbin= 400, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8.filteredMass[d.jetsAK8.it]; }, .axis_title="AK8 jet filtered mass (GeV)",  .def_range={0,150}});
    sh.AddNewFillParam("JetAK8TrimmedMass",    { .nbin= 400, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8.trimmedMass[d.jetsAK8.it];  }, .axis_title="AK8 jet trimmed mass (GeV)",   .def_range={0,150}});
  */
#if VER == 0
  sh.AddNewFillParam("JetAK8Tau1",           { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return std::min(d.jetsAK8.tau1[d.jetsAK8.it], float(0.999));    }, .axis_title="AK8 jet #tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau2",           { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return std::min(d.jetsAK8.tau2[d.jetsAK8.it], float(0.999));    }, .axis_title="AK8 jet #tau_{2}"});
  sh.AddNewFillParam("JetAK8Tau3",           { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return std::min(d.jetsAK8.tau3[d.jetsAK8.it], float(0.999));    }, .axis_title="AK8 jet #tau_{3}"});
  sh.AddNewFillParam("MaxAK8SubjetCSV",      { .nbin=  20, .bins={     0,   1.00}, .fill=[&d] { return std::min(maxSubjetCSV[d.jetsAK8.it], float(0.999));      }, .axis_title="Max. AK8 subjet CSV"});
  sh.AddNewFillParam("MaxAK8SubJetCSVBins",  { .nbin=CSV.size()-1, .bins=CSV,      .fill=[&d] { return std::min(maxSubjetCSV[d.jetsAK8.it], float(0.999));      }, .axis_title="Max. AK8 subjet CSV"});
#else
  sh.AddNewFillParam("JetAK8Tau1",           { .nbin= 100, .bins={     0,      1}, .fill=[&d] { return std::min(d.jetsAK8.tau1Puppi[d.jetsAK8.it], float(0.999));         }, .axis_title="AK8 jet #tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau2",           { .nbin= 100, .bins={     0,      1}, .fill=[&d] { return std::min(d.jetsAK8.tau2Puppi[d.jetsAK8.it], float(0.999));         }, .axis_title="AK8 jet #tau_{2}"});
  sh.AddNewFillParam("JetAK8Tau3",           { .nbin= 100, .bins={     0,      1}, .fill=[&d] { return std::min(d.jetsAK8.tau3Puppi[d.jetsAK8.it], float(0.999));         }, .axis_title="AK8 jet #tau_{3}"});
  sh.AddNewFillParam("MaxAK8SubjetCSV",      { .nbin=  20, .bins={     0,   1.00}, .fill=[&d] { return std::min(d.jetsAK8.maxSubjetCSVv2[d.jetsAK8.it], float(0.999));    }, .axis_title="Max. AK8 subjet CSV"});
  sh.AddNewFillParam("MaxAK8SubJetCSVBins",  { .nbin=CSV.size()-1, .bins=CSV,      .fill=[&d] { return std::min(d.jetsAK8.maxSubjetCSVv2[d.jetsAK8.it], float(0.999));    }, .axis_title="Max. AK8 subjet CSV"});
#endif
  sh.AddNewFillParam("JetAK8Tau21",          { .nbin=  20, .bins={     0,      1}, .fill=[&d] { return std::min(tau21[d.jetsAK8.it], double(0.999));                       }, .axis_title="AK8 jet #tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau31",          { .nbin=  20, .bins={     0,      1}, .fill=[&d] { return std::min(tau31[d.jetsAK8.it], double(0.999));                       }, .axis_title="AK8 jet #tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau32",          { .nbin=  20, .bins={     0,      1}, .fill=[&d] { return std::min(tau32[d.jetsAK8.it], double(0.999));                       }, .axis_title="AK8 jet #tau_{3}/#tau_{2}"});
  // mWs
  sh.AddNewFillParam("mWPtBins",        { .nbin=Pt.size()-1 , .bins=Pt,       .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="Mass-tagged W p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("mWPt",            { .nbin= 200, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="Mass-tagged W p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("mWEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8.Eta[d.jetsAK8.it];          }, .axis_title="Mass-tagged W #eta",        .def_range={-2.4,2.4}});
  sh.AddNewFillParam("mWPhi",           { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8.Phi[d.jetsAK8.it];          }, .axis_title="Mass-tagged W #phi"});
  sh.AddNewFillParam("mWTau21",         { .nbin=  20, .bins={     0,      1}, .fill=[&d] { return tau21[d.jetsAK8.it];                  }, .axis_title="Mass-tagged W #tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("mWMass",          { .nbin=M.size()-1, .bins=M,          .fill=[&d] { return softDropMassW[d.jetsAK8.it];          }, .axis_title="Mass-tagged W M_{Soft-Drop} (GeV)"});
  // aWs
  sh.AddNewFillParam("aWPtBins",        { .nbin=Pt.size()-1 , .bins=Pt,       .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="Anti-tagged W p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("aWPt",            { .nbin= 200, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="Anti-tagged W p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("aWEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8.Eta[d.jetsAK8.it];          }, .axis_title="Anti-tagged W #eta",        .def_range={-2.4,2.4}});
  sh.AddNewFillParam("aWPhi",           { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8.Phi[d.jetsAK8.it];          }, .axis_title="Anti-tagged W #phi"});
  sh.AddNewFillParam("aWTau21",         { .nbin=  20, .bins={     0,      1}, .fill=[&d] { return tau21[d.jetsAK8.it];                  }, .axis_title="Anti-tagged W #tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("aWMass",          { .nbin=M.size()-1, .bins=M,          .fill=[&d] { return softDropMassW[d.jetsAK8.it];          }, .axis_title="Anti-tagged W M_{Soft-Drop} (GeV)"});
  // Ws
  sh.AddNewFillParam("WPtBins",         { .nbin=Pt.size()-1, .bins=Pt,        .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="Tagged W p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("WPt",             { .nbin= 200, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8.Pt[d.jetsAK8.it];           }, .axis_title="Tagged W p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("WEta",            { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8.Eta[d.jetsAK8.it];          }, .axis_title="Tagged W #eta",        .def_range={-2.4,2.4}});
  sh.AddNewFillParam("WPhi",            { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8.Phi[d.jetsAK8.it];          }, .axis_title="Tagged W #phi"});
  sh.AddNewFillParam("WTau21",          { .nbin=  20, .bins={     0,      1}, .fill=[&d] { return tau21[d.jetsAK8.it];                  }, .axis_title="Tagged W #tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("WMass",           { .nbin=M.size()-1, .bins=M,          .fill=[&d] { return softDropMassW[d.jetsAK8.it];          }, .axis_title="Tagged W M_{Soft-Drop} (GeV)"});

  // Leptons
  sh.AddNewFillParam("VetoElePt",       { .nbin= 200, .bins={     0,   1000}, .fill=[&d] { return d.ele.Pt[d.ele.it];                   }, .axis_title="Loose Electron p_{T} (GeV)", .def_range={0,500}});
  sh.AddNewFillParam("VetoEleEta",      { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.ele.Eta[d.ele.it];                  }, .axis_title="Loose Electron #eta (GeV)",  .def_range={-2.5,2.5}});
  sh.AddNewFillParam("VetoMuPt",        { .nbin= 200, .bins={     0,   1000}, .fill=[&d] { return d.mu.Pt[d.mu.it];                     }, .axis_title="Loose Muon p_{T} (GeV)",     .def_range={0,500}});
  sh.AddNewFillParam("VetoMuEta",       { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.mu.Eta[d.mu.it];                    }, .axis_title="Loose Muon #eta (GeV)",      .def_range={-2.4,2.4}});

  sh.AddNewFillParam("ElePt",           { .nbin= 200, .bins={     0,   1000}, .fill=[&d] { return d.ele.Pt[d.ele.it];                   }, .axis_title="Tight Electron p_{T} (GeV)", .def_range={0,250}});
  sh.AddNewFillParam("EleEta",          { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.ele.Eta[d.ele.it];                  }, .axis_title="Tight Electron #eta (GeV)",  .def_range={-2.5,2.5}});
  sh.AddNewFillParam("EleJetDR",        { .nbin=  60, .bins={     0,      6}, .fill=[&d] { return eleJetDR[d.ele.it];                   }, .axis_title="#DeltaR (ele, jet)",         .def_range={0,4}});
  sh.AddNewFillParam("EleJetPt",        { .nbin= 200, .bins={     0,   1000}, .fill=[&d] { return eleJetPt[d.ele.it];                   }, .axis_title="p_{T, nearest jet to ele}"});
  sh.AddNewFillParam("EleJetDPhi",      { .nbin=MDP.size()-1, .bins=MDP,      .fill=[&d] { return eleJetDPhi[d.ele.it];                 }, .axis_title="#Delta#phi (ele, jet)"});
  sh.AddNewFillParam("Ele1JetDPhi",     { .nbin=MDP.size()-1, .bins=MDP,      .fill=[&d] { return nEleSelect<1 ? -9999 : eleJetDPhi[iEleSelect[0]]; }, .axis_title="#Delta#phi (1st ele, jet)"});
  sh.AddNewFillParam("Ele2JetDPhi",     { .nbin=MDP.size()-1, .bins=MDP,      .fill=[&d] { return nEleSelect<2 ? -9999 : eleJetDPhi[iEleSelect[1]]; }, .axis_title="#Delta#phi (2nd ele, jet)"});

  sh.AddNewFillParam("MuPt",            { .nbin= 200, .bins={     0,   1000}, .fill=[&d] { return d.mu.Pt[d.mu.it];                     }, .axis_title="Tight Muon p_{T} (GeV)",     .def_range={0,500}});
  sh.AddNewFillParam("MuEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.mu.Eta[d.mu.it];                    }, .axis_title="Tight Muon #eta (GeV)",      .def_range={-2.4,2.4}});
  sh.AddNewFillParam("MuJetDR",         { .nbin=  60, .bins={     0,      6}, .fill=[&d] { return muJetDR[d.mu.it];                     }, .axis_title="#DeltaR (muon, jet)",        .def_range={0,4}});
  sh.AddNewFillParam("MuJetPt",         { .nbin= 200, .bins={     0,   1000}, .fill=[&d] { return muJetPt[d.mu.it];                     }, .axis_title="p_{T, nearest jet to muon}"});
  sh.AddNewFillParam("MuJetDPhi",       { .nbin=MDP.size()-1, .bins=MDP,      .fill=[&d] { return muJetDPhi[d.mu.it];                   }, .axis_title="#Delta#phi (muon, jet)"});
  sh.AddNewFillParam("Mu1JetDPhi",      { .nbin=MDP.size()-1, .bins=MDP,      .fill=[&d] { return nMuSelect<1 ? -9999 : muJetDPhi[iMuSelect[0]]; }, .axis_title="#Delta#phi (1st muon, jet)"});
  sh.AddNewFillParam("Mu2JetDPhi",      { .nbin=MDP.size()-1, .bins=MDP,      .fill=[&d] { return nMuSelect<2 ? -9999 : muJetDPhi[iMuSelect[1]]; }, .axis_title="#Delta#phi (2nd muon, jet)"});

  /*
    sh.AddNewFillParam("EleDRJet",             { .nbin=  60, .bins={    0,       6}, .fill=[&d] { return d.evt.EleDRJet[d.ele.it];                       }, .axis_title="#DeltaR (e, jet)"});
    sh.AddNewFillParam("EleRelPtJet",          { .nbin=  50, .bins={    0,     500}, .fill=[&d] { return d.evt.EleRelPtJet[d.ele.it];                    }, .axis_title="p_{T}^{rel} (e, jet) (GeV)"});
    sh.AddNewFillParam("MuDRJet",              { .nbin=  60, .bins={    0,       6}, .fill=[&d] { return d.evt.MuDRJet[d.mu.it];                         }, .axis_title="#DeltaR (#mu, jet)"});
    sh.AddNewFillParam("MuRelPtJet",           { .nbin=  50, .bins={    0,     500}, .fill=[&d] { return d.evt.MuRelPtJet[d.mu.it];                      }, .axis_title="p_{T}^{rel} (#mu, jet) (GeV)"});
  */



  // Event
  // Cuts
  // Object counts
  sh.AddNewFillParam("NVtx",                 { .nbin=NVTX.size()-1, .bins=NVTX,    .fill=[&d] { return d.evt.NGoodVtx;          }, .axis_title="N_{Vertices}",         .def_range={0,50}});
  sh.AddNewFillParam("NJet",                 { .nbin=  50, .bins={    0,      50}, .fill=[&d] { return nJet;                    }, .axis_title="N_{Jet}",              .def_range={2,20}});
  sh.AddNewFillParam("NJetAK8",              { .nbin=  10, .bins={    0,      10}, .fill=[&d] { return nJetAK8;                 }, .axis_title="N_{AK8 jet}",          .def_range={1,10}});
  sh.AddNewFillParam("NBTag",                { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nMediumBTag;             }, .axis_title="N_{b}",                .def_range={0,8}});
  sh.AddNewFillParam("NLooseBTag",           { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nLooseBTag;              }, .axis_title="N_{b, loose tag}",     .def_range={0,8}});
  sh.AddNewFillParam("NTightBTag",           { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nTightBTag;              }, .axis_title="N_{b, tight tag}",     .def_range={0,5}});
  sh.AddNewFillParam("NmW",                  { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nWPreTag;                }, .axis_title="N_{W, mass-tag}",       .def_range={0,5}});
  sh.AddNewFillParam("NaW",                  { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nTightWAntiTag;          }, .axis_title="N_{W, anti-tag}",       .def_range={0,5}});
  sh.AddNewFillParam("NW",                   { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nTightWTag;              }, .axis_title="N_{W}",                .def_range={0,5}});
  sh.AddNewFillParam("NLooseW",              { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nLooseWTag;              }, .axis_title="N_{W, loose tag}",     .def_range={0,5}});
  sh.AddNewFillParam("NHadTopTag",           { .nbin=   8, .bins={    0,       8}, .fill=[&d] { return nHadTopTag;              }, .axis_title="N_{top (had.)}",       .def_range={0,5}});
  sh.AddNewFillParam("NLepVeto",             { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nLepVeto;                }, .axis_title="N_{lepton, Veto}",     .def_range={0,5}});
  sh.AddNewFillParam("NEleVeto",             { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nEleVeto;                }, .axis_title="N_{ele, Veto}",        .def_range={0,5}});
  sh.AddNewFillParam("NMuVeto",              { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nMuVeto;                 }, .axis_title="N_{muon, Veto}",       .def_range={0,5}});
  sh.AddNewFillParam("NLepLoose",            { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nLepLoose;               }, .axis_title="N_{lepton, Loose}",     .def_range={0,5}});
  sh.AddNewFillParam("NEleLoose",            { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nEleLoose;               }, .axis_title="N_{ele, Loose}",        .def_range={0,5}});
  sh.AddNewFillParam("NMuLoose",             { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nMuLoose;                }, .axis_title="N_{muon, Loose}",       .def_range={0,5}});
  sh.AddNewFillParam("NIsoTrk",              { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return d.evt.NIsoTrk;           }, .axis_title="N_{iso trk}",          .def_range={0,5}});
  sh.AddNewFillParam("NLep",                 { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nLepSelect;              }, .axis_title="N_{lepton}",           .def_range={0,5}});
  sh.AddNewFillParam("NEle",                 { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nEleSelect;              }, .axis_title="N_{ele}",              .def_range={0,5}});
  sh.AddNewFillParam("NMu",                  { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nMuSelect;               }, .axis_title="N_{muon}",             .def_range={0,5}});
  // Razor
  sh.AddNewFillParam("R",                    { .nbin=  40, .bins={    0,     2.0}, .fill=[&d] { return d.evt.R;                 }, .axis_title="R",                    .def_range={0,1}});
  sh.AddNewFillParam("RFine",                { .nbin= 200, .bins={    0,     2.0}, .fill=[&d] { return d.evt.R;                 }, .axis_title="R",                    .def_range={0,1}});
  sh.AddNewFillParam("RBins",                { .nbin=R.size()-1, .bins=R,          .fill=[&d] { return d.evt.R;                 }, .axis_title="R",                    .def_range={0,1}});
  //sh.AddNewFillParam("MR",                   { .nbin= 100, .bins={    0,   10000}, .fill=[&d] { return d.evt.MR;                }, .axis_title="M_{R} (GeV)",          .def_range={0,2000}});
  //sh.AddNewFillParam("R2",                   { .nbin=  80, .bins={    0,     4.0}, .fill=[&d] { return d.evt.R2;                }, .axis_title="R^{2}",                .def_range={0,1}});
  sh.AddNewFillParam("MR",                   { .nbin=MR.size()-1, .bins=MR,        .fill=[&d] { return d.evt.MR;                }, .axis_title="M_{R} (GeV)",          .def_range={0,4000}});
  //sh.AddNewFillParam("MTR",                  { .nbin=  80, .bins={    0,    4000}, .fill=[&d] { return d.evt.MTR;               }, .axis_title="M_{T}^{R} (GeV)",      .def_range={0,2000}});
  sh.AddNewFillParam("MTR",                  { .nbin=MTR.size()-1, .bins=MTR,      .fill=[&d] { return d.evt.MTR;               }, .axis_title="M_{T}^{R} (GeV)",      .def_range={0,2000}});
  sh.AddNewFillParam("R2",                   { .nbin=R2.size()-1, .bins=R2,        .fill=[&d] { return d.evt.R2;                }, .axis_title="R^{2}",                .def_range={0,1}});
  sh.AddNewFillParam("MTRll",                { .nbin=MTR.size()-1, .bins=MTR,      .fill=[&d] { return MTR_ll;                  }, .axis_title="M_{T,ll}^{R} (GeV)",   .def_range={0,2000}});
  sh.AddNewFillParam("R2ll",                 { .nbin=R2.size()-1, .bins=R2,        .fill=[&d] { return R2_ll;                   }, .axis_title="R_{ll}^{2}",           .def_range={0,1}});
  // HT
  //sh.AddNewFillParam("HT",                   { .nbin=HT.size()-1, .bins=HT,        .fill=[&d] { return AK4_Ht;             }, .axis_title="H_{T} (GeV)",          .def_range={200, 4000}});
  sh.AddNewFillParam("HT",                   { .nbin= 100, .bins={    0,   10000}, .fill=[&d] { return AK4_Ht;             }, .axis_title="H_{T} (GeV)",          .def_range={400, 3000}});
  sh.AddNewFillParam("OnlineHT",             { .nbin= 100, .bins={    0,   10000}, .fill=[&d] { return AK4_HtOnline;       }, .axis_title="H_{T}^{HLT} (GeV)",    .def_range={400, 3000}});
  sh.AddNewFillParam("HTBins",               { .nbin=HTB.size()-1, .bins=HTB,      .fill=[&d] { return AK4_Ht;             }, .axis_title="H_{T} (GeV)",          .def_range={400, 1500}});
  sh.AddNewFillParam("GenHT",                { .nbin=HT.size()-1, .bins=HT,        .fill=[&d] { return d.evt.Gen_Ht;            }, .axis_title="H_{T}^{Gen} (GeV)",    .def_range={0, 2000}});
  sh.AddNewFillParam("AK8HT",                { .nbin=HT.size()-1, .bins=HT,        .fill=[&d] { return AK8_Ht;             }, .axis_title="H_{T}^{AK8} (GeV)",    .def_range={0, 2000}});
  // MET
  sh.AddNewFillParam("MET",                  { .nbin=MET.size()-1, .bins=MET,      .fill=[&d] { return d.met.Pt[0];             }, .axis_title="#slash{E}_{T} (GeV)", .def_range={0,2000}});
  sh.AddNewFillParam("METll",                { .nbin=MET.size()-1, .bins=MET,      .fill=[&d] { return MET_ll;                  }, .axis_title="#slash{E}_{T,ll} (GeV)", .def_range={0,2000}});
  sh.AddNewFillParam("Met",                  { .nbin=  80, .bins={    0,    4000}, .fill=[&d] { return d.met.Pt[0];             }, .axis_title="MET (GeV)", .def_range={0,2000}});
  // DPhi
  //sh.AddNewFillParam("MinDeltaPhi",          { .nbin=  64, .bins={    0,     3.2}, .fill=[]   { return minDeltaPhi;             }, .axis_title="#Delta#phi_{min}"});
  sh.AddNewFillParam("MinDeltaPhi",          { .nbin=MDP.size()-1, .bins=MDP,      .fill=[]   { return minDeltaPhi;             }, .axis_title="#Delta#phi_{min}"});
  sh.AddNewFillParam("MinDeltaPhill",        { .nbin=MDP.size()-1, .bins=MDP,      .fill=[]   { return minDeltaPhi_ll;          }, .axis_title="#Delta#phi_{min,ll}"});
  sh.AddNewFillParam("DeltaPhiLLMET",        { .nbin=MDP.size()-1, .bins=MDP,      .fill=[]   { return dPhi_ll_met;             }, .axis_title="#Delta#phi (ll, MET)"});
  sh.AddNewFillParam("DeltaPhiLLJet",        { .nbin=MDP.size()-1, .bins=MDP,      .fill=[]   { return dPhi_ll_jet;             }, .axis_title="#Delta#phi_{min} (ll, jet)"});
  sh.AddNewFillParam("DeltaRWb",             { .nbin=  60, .bins={    0,       6}, .fill=[]   { return minDeltaR_W_b;           }, .axis_title="#DeltaR_{min} (W, b)"});
  // MT/Mll
  sh.AddNewFillParam("MT",                   { .nbin= 100, .bins={    0,    2000}, .fill=[]   { return MT;                      }, .axis_title="m_{T} (GeV)",  .def_range={0,500}});
  sh.AddNewFillParam("Mll",                  { .nbin=  50, .bins={    0,     500}, .fill=[]   { return M_ll;                    }, .axis_title="m_{ll} (GeV)", .def_range={0,200}});
  // SUSY
  sh.AddNewFillParam("MGluino",              { .nbin= 121, .bins={-12.5, 3012.5 }, .fill=[&d] { return d.evt.SUSY_Gluino_Mass;  }, .axis_title="M_{#tilde{g}} (GeV)",        .def_range={550,2350}});
  sh.AddNewFillParam("MStop",                { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[&d] { return d.evt.SUSY_Stop_Mass;    }, .axis_title="M_{#tilde{s}} (GeV)",        .def_range={  0,1650}});
  sh.AddNewFillParam("MLSP",                 { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[&d] { return d.evt.SUSY_LSP_Mass;     }, .axis_title="M_{#tilde{#chi}^{0}} (GeV)", .def_range={  0,1650}});
  sh.AddNewFillParam("StopLSPMassDiff",      { .nbin= 400, .bins={0, 2000 },       .fill=[&d] { return d.evt.SUSY_Stop_Mass-d.evt.SUSY_LSP_Mass; }, .axis_title="M_{#tilde{s}}-M_{#tilde{#chi}^{0}} (GeV)"});
  // AK8 JetN
  sh.AddNewFillParam("Jet1AK8Mass",         { .nbin=M.size()-1, .bins=M,           .fill=[&d] { return (nJetAK8<1) ? -9999. : softDropMassW[iJetAK8[0]]; }, .axis_title="Leading AK8 jet M_{Soft-Drop} (GeV)",    .def_range={0, 300}});
  sh.AddNewFillParam("Jet2AK8Mass",         { .nbin=M.size()-1, .bins=M,           .fill=[&d] { return (nJetAK8<2) ? -9999. : softDropMassW[iJetAK8[1]]; }, .axis_title="Subleading AK8 jet M_{Soft-Drop} (GeV)", .def_range={0, 300}});
  sh.AddNewFillParam("Jet1AK8Pt",           { .nbin=  100, .bins={    0,   10000}, .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8.Pt[iJetAK8[0]];           }, .axis_title="Leading AK8 jet p_{T} (GeV)",    .def_range={200, 1000}});
  sh.AddNewFillParam("Jet1AK8PtBins",       { .nbin=PtB.size()-1, .bins=PtB,       .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8.Pt[iJetAK8[0]];           }, .axis_title="Leading AK8 jet p_{T} (GeV)",    .def_range={200, 1000}});
  sh.AddNewFillParam("Jet2AK8PtBins",       { .nbin=PtB.size()-1, .bins=PtB,       .fill=[&d] { return (nJetAK8<2) ? -9999. : d.jetsAK8.Pt[iJetAK8[1]];           }, .axis_title="Subleading AK8 jet p_{T} (GeV)", .def_range={200, 1000}});
  sh.AddNewFillParam("Jet1AK8Eta",          { .nbin=   80, .bins={   -4,       4}, .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8.Eta[iJetAK8[0]];          }, .axis_title="Leading AK8 jet #eta",    .def_range={-3, 3}});
  sh.AddNewFillParam("Jet2AK8Eta",          { .nbin=   80, .bins={   -4,       4}, .fill=[&d] { return (nJetAK8<2) ? -9999. : d.jetsAK8.Eta[iJetAK8[1]];          }, .axis_title="Subleading AK8 jet #eta", .def_range={-3, 3}});
  sh.AddNewFillParam("Jet1AK8Tau32",        { .nbin=   50, .bins={    0,       1}, .fill=[&d] { return (nJetAK8<1) ? -9999. : tau32[iJetAK8[0]];                       }, .axis_title="Leading AK8 jet #tau_{32}"});
  sh.AddNewFillParam("Jet2AK8Tau32",        { .nbin=   50, .bins={    0,       1}, .fill=[&d] { return (nJetAK8<2) ? -9999. : tau32[iJetAK8[1]];                       }, .axis_title="Subleading AK8 jet #tau_{32}"});
  sh.AddNewFillParam("mW1Mass",             { .nbin=MW.size()-1, .bins=MW,         .fill=[&d] { return (nWPreTag<1)? -9999. : softDropMassW[iWPreTag[0]]; }, .axis_title="Mass-tagged W M_{Soft-Drop} (GeV)"});
#if VER == 0
  sh.AddNewFillParam("Jet1AK8BTagCSV",      { .nbin=  101, .bins={    0,    1.01}, .fill=[&d] { return (nJetAK8<1) ? -9999. : maxSubjetCSV[iJetAK8[0]];            }, .axis_title="Leading AK8 jet - Max. Subjet CSV",    .def_range={0,1}});
  sh.AddNewFillParam("Jet2AK8BTagCSV",      { .nbin=  101, .bins={    0,    1.01}, .fill=[&d] { return (nJetAK8<2) ? -9999. : maxSubjetCSV[iJetAK8[1]];            }, .axis_title="Subleading AK8 jet - Max. Subjet CSV", .def_range={0,1}});
#else
  sh.AddNewFillParam("Jet1AK8BTagCSV",      { .nbin=  101, .bins={    0,    1.01}, .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8.maxSubjetCSVv2[iJetAK8[0]];     }, .axis_title="Leading AK8 jet - Max. Subjet CSV",    .def_range={0,1}});
  sh.AddNewFillParam("Jet2AK8BTagCSV",      { .nbin=  101, .bins={    0,    1.01}, .fill=[&d] { return (nJetAK8<2) ? -9999. : d.jetsAK8.maxSubjetCSVv2[iJetAK8[1]];     }, .axis_title="Subleading AK8 jet - Max. Subjet CSV", .def_range={0,1}});
#endif
  // gen particles
  sh.AddNewFillParam("GenWPt",              { .nbin=  200, .bins={    0,   10000}, .fill=[&d] { return d.gen.Pt[d.gen.it];  }, .axis_title="Gen-W p_{T} (GeV)",   .def_range={000, 2000}});
  sh.AddNewFillParam("GenWPtBins",          { .nbin=Pt.size()-1, .bins=Pt,         .fill=[&d] { return d.gen.Pt[d.gen.it];  }, .axis_title="Gen-W p_{T} (GeV)",   .def_range={000, 2000}});
  sh.AddNewFillParam("GenTopPt",            { .nbin=  200, .bins={    0,   10000}, .fill=[&d] { return d.gen.Pt[d.gen.it];  }, .axis_title="Gen-top p_{T} (GeV)", .def_range={000, 2000}});

  if (debug) std::cout<<"Analysis::define_histo_options: non-special fillparams ok"<<std::endl;

  // SPECIAL
  // Special Y/Z axis parameters:
  sh.AddSpecial({ .name="Counts",                         .name_plus_1d="Syst",                          .axis="Counts (Incl Syst Unc)",              .axis_plus_1d="Systematics variation index"});
  /*
    sh.AddSpecial({ .name="MergedTopFraction",              .name_plus_1d="IsGenTopMerged",                .axis="Fraction of Merged Tops",             .axis_plus_1d="Gen. b and W Merged in R<0.8 (bool)"});
    sh.AddSpecial({ .name="JetFindingEfficiency",           .name_plus_1d="HasJet",                        .axis="Jet finding Efficiency",              .axis_plus_1d="Found AK8 jet (bool)"});
    sh.AddSpecial({ .name="TopFindingEfficiency",           .name_plus_1d="HasHadTopTaggedJet",            .axis="Top finding Efficiency",              .axis_plus_1d="Has hadronic top-tagged jet (bool)"});
    sh.AddSpecial({ .name="TopTagEfficiency",               .name_plus_1d="JetIsHadTopTagged",             .axis="Top-tagging Efficiency",              .axis_plus_1d="Jet is hadronic top-tagged (bool)"});
    sh.AddSpecial({ .name="MisTagRate",                     .name_plus_1d="JetHasNoGenTop",                .axis="Mis-tag Rate",                        .axis_plus_1d="Jet is mis-matched (bool)"});
  */
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360",             .name_plus_1d="HLT_AK8PFJet360_TrimMass30",    .axis="#epsilon_{HLT_AK8PFJet360_TrimMass30}",                   .axis_plus_1d="HLT_AK8PFJet360_TrimMass30"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet450",             .name_plus_1d="HLT_Ak8PFJet450",               .axis="#epsilon_{HLT_AK8PFJet450}",                              .axis_plus_1d="HLT_AK8PFJet450"});
  sh.AddSpecial({ .name="HLTEff_AK8PFHT700_TrimMass50",   .name_plus_1d="HLT_AK8PFHT700_TrimMass50",     .axis="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}",          .axis_plus_1d="HLT_AK8PFHT700_TrimR0p1PT0p03Mass50"});
  sh.AddSpecial({ .name="HLTEff_PFHT750_4JetPt50",        .name_plus_1d="HLT_PFHT750_4JetPt50",          .axis="#epsilon_{HLT_PFHT750_4JetPt50}",                         .axis_plus_1d="HLT_PFHT750_4JetPt50"});
  sh.AddSpecial({ .name="HLTEff_PFHT800or900",                 .name_plus_1d="HLT_PFHT800or900",                   .axis="#epsilon_{HLT_PFHT800or900}",                                  .axis_plus_1d="HLT_PFHT800or900"});
  //sh.AddSpecial({ .name="HLTEff_AK8DiPFJet250_200",       .name_plus_1d="HLT_AK8DiPFJet250_200",         .axis="#epsilon_{HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20}", .axis_plus_1d="HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20"});
  //sh.AddSpecial({ .name="HLTEff_Rsq0p25",                 .name_plus_1d="HLT_Rsq0p25",                   .axis="#epsilon_{HLT_Rsq0p25}",                                  .axis_plus_1d="HLT_Rsq0p25"});
  //sh.AddSpecial({ .name="HLTEff_RsqMR270_Rsq0p09_MR200",  .name_plus_1d="HLT_RsqMR270_Rsq0p09_MR200",    .axis="#epsilon_{HLT_RsqMR270_Rsq0p09_MR200}",                   .axis_plus_1d="HLT_RsqMR270_Rsq0p09_MR200"});
  sh.AddSpecial({ .name="HLTEff_AK8PFHT700orPFHT800or900",     .name_plus_1d="HLT_AK8PFHT700_or_PFHT800or900",     .axis="#epsilon_{HLT_AK8PFHT700 OR HLT_PFHT800or900}",                .axis_plus_1d="HLT_AK8PFHT700 OR HLT_PFHT800or900"});
  sh.AddSpecial({ .name="HLTEff_PFJet450orPFHT800or900",       .name_plus_1d="HLT_PFJet450_or_PFHT800or900",       .axis="#epsilon_{HLT_PFJet450 OR HLT_PFHT800or900}",                  .axis_plus_1d="HLT_PFJet450 OR HLT_PFHT800or900"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet450orPFHT800or900",    .name_plus_1d="HLT_AK8PFJet450_or_PFHT800or900",    .axis="#epsilon_{HLT_AK8PFJet450 OR HLT_PFHT800or900}",               .axis_plus_1d="HLT_AK8PFJet450 OR HLT_PFHT800or900"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet450orAK8PFHT700", .name_plus_1d="HLT_AK8PFJet450_or_AK8PFHT700", .axis="#epsilon_{HLT_AK8PFJet450 OR HLT_AK8PFHT700}",            .axis_plus_1d="HLT_AK8PFJet450 OR HLT_AK8PFHT700"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360orPFHT800or900",    .name_plus_1d="HLT_AK8PFJet360_or_PFHT800or900",    .axis="#epsilon_{HLT_AK8PFJet360 OR HLT_PFHT800or900}",               .axis_plus_1d="HLT_AK8PFJet360 OR HLT_PFHT800or900"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360orAK8PFHT700", .name_plus_1d="HLT_AK8PFJet360_or_AK8PFHT700", .axis="#epsilon_{HLT_AK8PFJet360 OR HLT_AK8PFHT700}",            .axis_plus_1d="HLT_AK8PFJet360 OR HLT_AK8PFHT700"});
  sh.AddSpecial({ .name="WTaggingEfficiency",        .name_plus_1d="PassWTag",            .axis="W-tagging Efficiency",         .axis_plus_1d="Pass W Tag"});
  sh.AddSpecial({ .name="SignalSelectionEfficiency",    .name_plus_1d="PassSignalSelection",    .axis="Signal Selection Efficiency - W ana",    .axis_plus_1d="Pass Signal Selection - W"});
  sh.AddSpecial({ .name="TopSignalSelectionEfficiency", .name_plus_1d="PassTopSignalSelection", .axis="Signal Selection Efficiency - top ana",  .axis_plus_1d="Pass Signal Selection - top"});
  sh.AddSpecial({ .name="SignalSignificance_T5ttcc",    .name_plus_1d="Bkg_T5ttcc",          .axis="S/#sqrt{S+B} - T5ttcc",        .axis_plus_1d="Background, Signal - T5ttcc"});
  sh.AddSpecial({ .name="SignalSignificance_T5tttt",    .name_plus_1d="Bkg_T5tttt",          .axis="S/#sqrt{S+B} - T5tttt",        .axis_plus_1d="Background, Signal - T5tttt"});
  sh.AddSpecial({ .name="SignalSignificance_T1tttt",    .name_plus_1d="Bkg_T1tttt",          .axis="S/#sqrt{S+B} - T1tttt",        .axis_plus_1d="Background, Signal - T1tttt"});
  sh.AddSpecial({ .name="SignalSignificance_T2tt",      .name_plus_1d="Bkg_T2tt",            .axis="S/#sqrt{S+B} - T2tt",          .axis_plus_1d="Background, Signal - T2tt"  });

  sh.AddNewFillParam("Counts",                         { .nbin= 1+syst_nSyst, .bins={-0.5, syst_nSyst+0.5}, .fill=[&syst_index] { return syst_index; }, .axis_title="Counts (Incl Syst Unc)"});
  /*
    sh.AddNewFillParam("MergedTopFraction",              { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return top_Children_Within_Cone[d.jetsAK8.it]; }, .axis_title="Fraction of Merged Tops" });
    sh.AddNewFillParam("JetFindingEfficiency",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return has_Matched_Jet[d.gen.it];                   }, .axis_title="Jet finding Efficiency" });
    sh.AddNewFillParam("TopFindingEfficiency",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return has_Matched_Tagged_Jet[d.gen.it];            }, .axis_title="Top finding Efficiency" });
    sh.AddNewFillParam("TopTagEfficiency",               { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return passHadTopTag[d.jetsAK8.it];            }, .axis_title="Top-tagging Efficiency" });
  */
  sh.AddNewFillParam("HLTEff_AK8PFJet360",             { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet360_TrimMass30;                   }, .axis_title="#epsilon_{HLT_AK8PFJet360_TrimMass30}",                   .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet450",             { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet450;                              }, .axis_title="#epsilon_{HLT_AK8PFJet450}",                              .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFHT700_TrimMass50",   { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50;          }, .axis_title="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}",          .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_PFHT750_4JetPt50",        { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.PFHT750_4JetPt50;                         }, .axis_title="#epsilon_{HLT_PFHT750_4JetPt50}",                         .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_PFHT800or900",                 { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.PFHT800==1 || d.hlt.PFHT900==1;           }, .axis_title="#epsilon_{HLT_PFHT800or900}",                                  .def_range={0,1} });
  //sh.AddNewFillParam("HLTEff_AK8DiPFJet250_200",       { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8DiPFJet250_200_TrimMass30_BTagCSV_p20; }, .axis_title="#epsilon_{HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20}", .def_range={0,1} });
  //sh.AddNewFillParam("HLTEff_Rsq0p25",                 { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.Rsq0p25;                                  }, .axis_title="#epsilon_{HLT_Rsq0p25}",                                  .def_range={0,1} });
  //sh.AddNewFillParam("HLTEff_RsqMR270_Rsq0p09_MR200",  { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.RsqMR270_Rsq0p09_MR200;                   }, .axis_title="#epsilon_{HLT_RsqMR270_Rsq0p09_MR200}",                   .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFHT700orPFHT800or900",     { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1;                }, .axis_title="#epsilon_{HLT_AK8PFHT700 OR HLT_PFHT800or900}",  .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_PFJet450orPFHT800or900",       { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1;                                       }, .axis_title="#epsilon_{HLT_PFJet450 OR HLT_PFHT800or900}",    .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet450orPFHT800or900",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1;                                    }, .axis_title="#epsilon_{HLT_AK8PFJet450 OR HLT_PFHT800or900}", .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet450orAK8PFHT700", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet450==1 || d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1;                                    }, .axis_title="#epsilon_{HLT_AK8PFJet450 OR HLT_AK8PFHT700}", .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet360orPFHT800or900",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet360_TrimMass30==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1;                         }, .axis_title="#epsilon_{HLT_AK8PFJet360 OR HLT_PFHT800or900}",     .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet360orAK8PFHT700", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet360_TrimMass30==1 || d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; }, .axis_title="#epsilon_{HLT_AK8PFJet360 OR HLT_AK8PFHT700}",  .def_range={0,1} });
  sh.AddNewFillParam("WTaggingEfficiency",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return genHadWPassWTag[itGenHadW[d.gen.it]]; }, .axis_title="W-tagging Efficiency",         .def_range={0,1} });
  sh.AddNewFillParam("SignalSelectionEfficiency",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return pass_all_cuts['S'];   }, .axis_title="Signal Selection Efficiency - W ana"});
  sh.AddNewFillParam("TopSignalSelectionEfficiency", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return pass_all_cuts['t'];   }, .axis_title="Signal Selection Efficiency - top ana"});
  sh.AddNewFillParam("SignalSignificance_T5ttcc",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return Bkg_T5ttcc_opt.index; }, .axis_title="S/#sqrt{S+B} - T5ttcc", .def_range={0,10}});
  sh.AddNewFillParam("SignalSignificance_T5tttt",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return Bkg_T5tttt_opt.index; }, .axis_title="S/#sqrt{S+B} - T5tttt", .def_range={0,10}});
  sh.AddNewFillParam("SignalSignificance_T1tttt",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return Bkg_T1tttt_opt.index; }, .axis_title="S/#sqrt{S+B} - T1tttt", .def_range={0,10}});
  sh.AddNewFillParam("SignalSignificance_T2tt",      { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return Bkg_T2tt_opt.index;   }, .axis_title="S/#sqrt{S+B} - T2tt",   .def_range={0,10}});

  if (debug) std::cout<<"Analysis::define_histo_options: fillparams ok"<<std::endl;
}


//_______________________________________________________
//                 List of Histograms

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  /*
    RazorBoost cuts:
    - MET filters + 1 Vtx
    - 3 jets
    - jet1 pt > 200
    - MR > 800
    - R2 > 0.08
    - nele == 0
    - nmu == 0
    - ntautrkIso = 0
    - nbtag > 0
    - nWtag > 0
    - minDPhi > 0.5 (small jet cone size)

    RazorBoost plot order:
    
    - cleaning + HLT
    + Vertices (before and after reweighting) - num
    Objects selections:
    + AK4                        - num, pt
    + bs                         - num,
    + AK8                        - num, mass, mass vs pt, dau1 vs dau2 pt, dau1 vs dau2 mass
    + W-masstags                 - num, dau1 vs dau2 pt, dau1 vs dau2 mass, massdrop(/pt), yasym, mdr(/pt) vs yasym
    + Ws                         - num, numb vs num, DRAK4
    + MET/Razor                  - MET, MR, R2, R2 vs MR, R2 vs MET (Also add lep/mu to MET) + MinDPhi plot combinations
    + muons                      - num
    + electrons                  - num
    + gen                        - numtop, pt(top), DR(W,b), pt(W), TT - had/semilep/dilep, DR(q1,q2), many others

  */

  //__________________________________
  //        Define Smarthistos
  //__________________________________

  // Define histo types (for different object to loop on, and different cuts to apply)
  sh.AddHistoType("gen W",   "Gen-Ws");
  sh.AddHistoType("gen top", "Gen-tops");
  sh.AddHistoType("AK4",     "Jets");
  sh.AddHistoType("b",       "b-tagged jets");
  sh.AddHistoType("b loose", "Loose b-tagged jets");
  sh.AddHistoType("AK8",     "AK8 jets");
  sh.AddHistoType("mW",      "Mass-tagged Ws");
  sh.AddHistoType("aW",      "Anti-tagged Ws");
  sh.AddHistoType("W",       "Tagged Ws");
  sh.AddHistoType("ele",     "Electrons");
  sh.AddHistoType("ele veto","Veto electrons");
  sh.AddHistoType("mu",      "Muons");
  sh.AddHistoType("mu veto", "Veto muons");
  sh.AddHistoType("evt",     "Events");
  sh.AddHistoType("syst",    "Systematic variations");

  // Histo options
  std::string d = "HISTE1";
  std::string o_stk_d = "LogSumw2Stack5AddRatioTwoCol57AddIntApproval15";
  std::string o_stk_s = "LogSumw2Stack5AddRatioTwoCol57AddIntApproval45";
  std::string o_1or2d_d = "Sumw2Approval15";
  std::string o_1or2d_s = "Sumw2Approval45";
  std::string o_norm_d = "Sumw2NormApproval15";
  std::string o_norm_s = "Sumw2NormApproval45";
  std::vector<double> r_stk  = {0,0, 1.01e-2,1e6,  0.4,0.9};
  std::vector<double> r_stk2 = {0,0, 1.01e-2,1e8,  0.4,0.9};
  std::string Stack = "StackPlot";

  // -------------------------------------------------------------------------
  //                                   Trigger

  //for (const auto& cut : {"TriggerPreSelection","TriggerPreSelPlus1mW","S_Excl1W","S_ExclMR_R2","S","T","Q","W"}) {
  //for (const auto& cut : {"TriggerPreSelection","TriggerPreSelPlus1mW","S"}) {
  for (const auto& cut : {"TriggerPreSelection"}) {
    // (AK8)HT triggers
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_R2_vs_MR",              .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_mW1Mass_vs_AK8HT",      .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Mass_vs_AK8HT",  .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_mW1Mass_vs_HTBins",     .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Mass_vs_HTBins", .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_AK8HT",                 .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_HT",                    .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_OnlineHT",              .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Bin",                   .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Pt",             .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Mass",           .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_mW1Mass",               .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_HT",                                  .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_OnlineHT",                            .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_Bin",                                 .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_Jet1AK8Pt",                           .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_Jet1AK8Mass",                         .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_mW1Mass",                             .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_R2_vs_MR",                            .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_mW1Mass_vs_HTBins",                   .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_Jet1AK8Mass_vs_HTBins",               .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0,1} });
    //for (const auto& iso : {"OtherUnisoLep", "OtherLooseLep"})    
    //  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_HT",                                .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });

    //sh.AddHistos("evt", { .fill="HLTEff_PFHT750_4JetPt50_vs_R2_vs_MR",                   .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_PFHT750_4JetPt50_vs_HT",                         .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_PFHT750_4JetPt50_vs_Bin",                        .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    
    // AK8/B2G triggers
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_R2_vs_MR",                        .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8PtBins_vs_Jet1AK8Mass",    .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8PtBins_vs_mW1Mass",        .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_mW1Mass",                         .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt",                       .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt",                       .pfs={"Triggers",cut,"Jet1AK8Pt450"},  .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt",                       .pfs={"Triggers",cut,"Jet1AK8Pt500"},  .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Mass",                     .pfs={"Triggers",cut},                  .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Mass",                     .pfs={"Triggers",cut,"Jet1AK8Mass65"}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_HT",                              .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_OnlineHT",                        .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,"Jet1AK8Pt450"},  .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,"Jet1AK8Pt500"},  .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,"Jet1AK8Mass65"}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,"Jet1AK8Mass65","Jet1AK8Pt450"}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,"Jet1AK8Mass65","Jet1AK8Pt500"}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_R2_vs_MR",                        .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_mW1Mass",                         .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_HT",                              .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_OnlineHT",                        .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_Bin",                             .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    
    //sh.AddHistos("evt", { .fill="HLTEff_AK8DiPFJet250_200_vs_R2_vs_MR",                  .pfs={"Triggers",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8DiPFJet250_200_vs_mW1Mass",                   .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8DiPFJet250_200_vs_Bin",                       .pfs={"Triggers",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,0, 0.45,0.45} });

    // Suggestion: AK8PFJet450 OR PFHT800or900
    //for (auto trigger_comb : std::vector<std::string>({"HLTEff_AK8PFHT700orPFHT800or900", "HLTEff_AK8PFJet450orPFHT800or900", "HLTEff_AK8PFJet450orAK8PFHT700", "HLTEff_AK8PFJet360orPFHT800or900","HLTEff_AK8PFJet360orAK8PFHT700"})) {
    for (auto trigger_comb : std::vector<std::string>({"HLTEff_AK8PFJet450orPFHT800or900", "HLTEff_AK8PFJet360orPFHT800or900"})) {
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Bin",                     .pfs={                "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_R2_vs_MR",                .pfs={                "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8PtBins_vs_HTBins", .pfs={                "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass_vs_HTBins",   .pfs={                "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_mW1Mass_vs_HTBins",       .pfs={                "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_HT",                      .pfs={                "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,2000, 0,1, 0.44,0.53} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_OnlineHT",                .pfs={                "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,2000, 0,1, 0.44,0.53} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass",             .pfs={                "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt",               .pfs={                "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_mW1Mass",                 .pfs={                "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_HT",                      .pfs={"Jet1AK8PtBins","Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt",               .pfs={"HTBins"       ,"Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Bin",                     .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_R2_vs_MR",                .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8PtBins_vs_HTBins", .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass_vs_HTBins",   .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_mW1Mass_vs_HTBins",       .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_HT",                      .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_OnlineHT",                .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt",               .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass",             .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_mW1Mass",                 .pfs={"NJet35",       "Trigger",cut}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_HT",                      .pfs={"Jet1AK8PtBins","Trigger",cut,"NJet35"}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt",               .pfs={"HTBins"       ,"Trigger",cut,"NJet35"}, .cuts={}, .draw="PE1",  .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
    }
    
    sh.AddHistos("evt", { .fill="HT",              .pfs={"Jet1AK8PtBins","Trigger",cut},          .cuts={}, .draw="HIST", .opt=o_1or2d_d+"TwoCol66", .ranges={0,0, 0,0, 0.35,0.9} });
    sh.AddHistos("evt", { .fill="HT",              .pfs={"Jet1AK8PtBins","Trigger",cut,"NJet35"}, .cuts={}, .draw="HIST", .opt=o_1or2d_d+"TwoCol66", .ranges={0,0, 0,0, 0.35,0.9} });

    sh.AddHistos("evt", { .fill="Jet1AK8Pt",       .pfs={"HTBins",       "Trigger",cut},          .cuts={}, .draw="HIST", .opt=o_1or2d_d+"TwoCol44", .ranges={0,0, 0,0, 0.35,0.9} });
    sh.AddHistos("evt", { .fill="Jet1AK8Pt",       .pfs={"HTBins",       "Trigger",cut,"NJet35"}, .cuts={}, .draw="HIST", .opt=o_1or2d_d+"TwoCol44", .ranges={0,0, 0,0, 0.35,0.9} });

    sh.AddHistos("evt", { .fill="Jet1AK8Pt_vs_HT", .pfs={"Triggers",cut},          .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0} });
    sh.AddHistos("evt", { .fill="Jet1AK8Pt_vs_HT", .pfs={"Triggers",cut,"NJet35"}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,0, 0,0} });
  }


  // -------------------------------------------------------------------------
  //                              Selected AK4 Jets

  for (auto region : {'S', 's', 'T','W','Q', 'q', 'Z', 't'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut = std::string(1,region);
    std::vector<std::string> showdata = {"JetHT"};
    if (region=='S'||region=='t') showdata.push_back("Blind");
    for (auto data : showdata ) {
      std::string opt = (data=="Blind") ? o_stk_s : o_stk_d;
      sh.AddHistos("AK4",  { .fill="JetPtBins",           .pfs={Stack,data,cut},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("AK4",  { .fill="JetPt",               .pfs={Stack,data,cut},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("AK4",  { .fill="JetEta",              .pfs={Stack,data,cut},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("AK4",  { .fill="JetPhi",              .pfs={Stack,data,cut},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("AK4",  { .fill="JetCSV",              .pfs={Stack,data,cut},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    }
  }

  // -------------------------------------------------------------------------
  //                                  Leptons
  
  // Veto Leptons
  for (auto region : {'S', 's', 'Q', 'q', 'T', 'W', 't'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut = std::string(1,region);
    std::vector<std::string> showdata = {"JetHT"};
    if (region=='S'||region=='t') showdata.push_back("Blind");
    if (region=='T'||region=='W') for (auto data : showdata)  {
      sh.AddHistos("evt",      { .fill="NEle",                 .pfs={Stack,data,cut+"_Excl1LepMT"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("ele veto", { .fill="VetoElePt",            .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("ele veto", { .fill="VetoEleEta",           .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("evt",      { .fill="NMu",                  .pfs={Stack,data,cut+"_Excl1LepMT"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("mu veto",  { .fill="VetoMuPt",             .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("mu veto",  { .fill="VetoMuEta",            .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("evt",      { .fill="NLep",                 .pfs={Stack,data,cut+"_Excl1LepMT"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    } else for (auto data : showdata) {
      std::string opt = (data=="Blind") ? o_stk_s : o_stk_d;
      sh.AddHistos("evt",      { .fill="NEleVeto",             .pfs={Stack,data,cut+"_Excl0Ele0IsoTrk"}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("ele veto", { .fill="VetoElePt",            .pfs={Stack,data,cut+"_Excl0Ele0IsoTrk"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("ele veto", { .fill="VetoEleEta",           .pfs={Stack,data,cut+"_Excl0Ele0IsoTrk"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("evt",      { .fill="NMuVeto",              .pfs={Stack,data,cut+"_Excl0Mu0IsoTrk"},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("mu veto",  { .fill="VetoMuPt",             .pfs={Stack,data,cut+"_Excl0Mu0IsoTrk"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
      sh.AddHistos("mu veto",  { .fill="VetoMuEta",            .pfs={Stack,data,cut+"_Excl0Mu0IsoTrk"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    }
  }

  // Selected Leptons
  for (auto region : {'Z'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut = std::string(1,region);
    sh.AddHistos("evt",  { .fill="NEle",                       .pfs={Stack,"JetHT",cut+"_ExclMR_R2ll2Lep"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NMu",                        .pfs={Stack,"JetHT",cut+"_ExclMR_R2ll2Lep"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NLep",                       .pfs={Stack,"JetHT",cut+"_ExclMR_R2ll2Lep"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});

    sh.AddHistos("ele",  { .fill="ElePt",                      .pfs={Stack,"JetHT",cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("ele",  { .fill="EleEta",                     .pfs={Stack,"JetHT",cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("mu",   { .fill="MuPt",                       .pfs={Stack,"JetHT",cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("mu",   { .fill="MuEta",                      .pfs={Stack,"JetHT",cut},                    .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});

    // DPhi debug
    sh.AddHistos("ele",  { .fill="EleJetPt",                   .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("ele",  { .fill="EleJetDR",                   .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("ele",  { .fill="EleJetDPhi",                 .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("ele",  { .fill="ElePt_vs_EleJetPt",          .pfs={    "Data_MC",cut+"_ExclmDPhill"},     .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("ele",  { .fill="EleJetDR_vs_EleJetPt",       .pfs={    "Data_MC",cut+"_ExclmDPhill"},     .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Ele1JetDPhi",                .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Ele2JetDPhi",                .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Ele2JetDPhi_vs_Ele1JetDPhi", .pfs={    "Data_MC",cut+"_ExclmDPhill"},     .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});

    sh.AddHistos("mu",   { .fill="MuJetPt",                    .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("mu",   { .fill="MuJetDR",                    .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("mu",   { .fill="MuJetDPhi",                  .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("mu",   { .fill="MuPt_vs_MuJetPt",            .pfs={    "Data_MC",cut+"_ExclmDPhill"},     .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("mu",   { .fill="MuJetDR_vs_MuJetPt",         .pfs={    "Data_MC",cut+"_ExclmDPhill"},     .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Mu1JetDPhi",                 .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Mu2JetDPhi",                 .pfs={Stack,"JetHT",cut+"_ExclmDPhill"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Mu2JetDPhi_vs_Mu1JetDPhi",   .pfs={    "Data_MC",cut+"_ExclmDPhill"},     .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
  }

  // ----------------------------------------------------------------------------------------------
  //                                          W ANALYSIS
  //-----------------------------------------------------------------------------------------------

  o_stk_d = "LogSumw2Stack5AddRatioTwoCol57AddIntApproval16";
  o_stk_s = "LogSumw2Stack5AddRatioTwoCol57AddIntApproval46";
  o_1or2d_d = "Sumw2Approval16";
  o_1or2d_s = "Sumw2Approval46";
  o_norm_d = "Sumw2NormApproval16";
  o_norm_s = "Sumw2NormApproval46";

  // -------------------------------------------------------------------------
  //                                     bs

  // Selected b-tags
  for (auto region : {'S', 's', 'T', 'Z', 't'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut1 = std::string(1,region);
    std::string cut2 = cut1;
    if (region=='S'||region=='s'||region=='T') cut2 += "_Excl1b";
    std::vector<std::string> showdata = {"JetHT"};
    if (region=='S'||region=='t') showdata.push_back("Blind");
    for (auto data : showdata ) {
      std::string opt = (data=="Blind") ? o_stk_s : o_stk_d;
      sh.AddHistos("b",    { .fill="BJetPtBins",          .pfs={Stack,data,cut1}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("b",    { .fill="BJetPt",              .pfs={Stack,data,cut1}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("b",    { .fill="BJetEta",             .pfs={Stack,data,cut1}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("b",    { .fill="BJetPhi",             .pfs={Stack,data,cut1}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("b",    { .fill="BJetCSV",             .pfs={Stack,data,cut1}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      sh.AddHistos("evt",  { .fill="NBTag",               .pfs={Stack,data,cut1}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      if (region!='Z'&&region!='t')
	sh.AddHistos("evt",  { .fill="NBTag",               .pfs={Stack,data,cut2}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    }
  }

  // Veto b-tags
  for (auto region : {'Q', 'q', 'W', 'Z'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut1 = std::string(1,region);
    std::string cut2 = cut1;
    if (region=='Q'||region=='q'||region=='W') cut2 += "_Excl0b";
    sh.AddHistos("evt",  { .fill="NLooseBTag",          .pfs={Stack,"JetHT",cut1}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    if (region!='Z')
      sh.AddHistos("evt",  { .fill="NLooseBTag",          .pfs={Stack,"JetHT",cut2}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  }

  // -------------------------------------------------------------------------
  //                                 AK8/Ws Jets

  for (auto region : {'S', 's', 'T','W','Q', 'q', 'Z', 't'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut1 = std::string(1,region);
    std::string cut2 = cut1;
    if      (region=='S'||region=='s'||region=='T') cut2 += "_Excl1W";
    else if (region=='Q'||region=='q') cut2 += "_Excl1aW";
    else if (region=='W'||region=='Z') cut2 += "_Excl1mW";
    else if (region=='t') cut2 += "_Excl1Top";
    std::vector<std::string> showdata = {"JetHT"};
    if (region=='S'||region=='t') showdata.push_back("Blind");
    for (auto cut : { cut1, cut2 }) {
      for (auto data : showdata ) {
	std::string opt = (data=="Blind") ? o_stk_s : o_stk_d;
	// Njet
	sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
	sh.AddHistos("AK8",  { .fill="JetAK8PtBins",       .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Pt",           .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Eta",          .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Phi",          .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Mass",         .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Mass",         .pfs={Stack,data,cut,"Tau21Tagged"}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Tau21",        .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="JetAK8Tau32",        .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("AK8",  { .fill="MaxAK8SubjetCSV",    .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});

	// Mass-tagged Ws
	sh.AddHistos("evt",  { .fill="NmW",                .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
	sh.AddHistos("mW",   { .fill="mWPtBins",           .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("mW",   { .fill="mWPt",               .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("mW",   { .fill="mWEta",              .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("mW",   { .fill="mWPhi",              .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("mW",   { .fill="mWTau21",            .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("mW",   { .fill="mWMass",             .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});

	// Anti-tagged Ws
	sh.AddHistos("evt",  { .fill="NaW",                .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
	sh.AddHistos("aW",   { .fill="aWPtBins",           .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("aW",   { .fill="aWPt",               .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("aW",   { .fill="aWEta",              .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("aW",   { .fill="aWPhi",              .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("aW",   { .fill="aWTau21",            .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("aW",   { .fill="aWMass",             .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});

	// Tagged Ws
	sh.AddHistos("W",    { .fill="WPtBins",            .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("W",    { .fill="WPt",                .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("W",    { .fill="WEta",               .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("W",    { .fill="WPhi",               .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("W",    { .fill="WTau21",             .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("W",    { .fill="WMass",              .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk2});
	sh.AddHistos("evt",  { .fill="NW",                 .pfs={Stack,data,cut},               .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      }
      sh.AddHistos("evt",  { .fill="NW",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={0,0, 0,1, 0.32,0.90}});
      sh.AddHistos("evt",  { .fill="NW",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={0,0, 0,1, 0.32,0.90}});
    }
  }

  // -------------------------------------------------------------------------
  //                              W GenInfo

  sh.SetHistoWeights({ [this] { return sf_weight['S']; } });
  sh.AddHistos("gen W",   { .fill="GenWPt",                                       .pfs={"TT_SignalPoints"},   .cuts={}, .draw=d,    .opt=o_1or2d_s+"Norm",.ranges={0,2000, 0,0.2, 0.6,0.9}});
  sh.AddHistos("gen W",   { .fill="WTaggingEfficiency_vs_GenWPtBins",             .pfs={"TT_Signal"},         .cuts={}, .draw="PE1",.opt=o_1or2d_s,       .ranges={0,2000, 0,1.0}});

  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MGluino",    .pfs={"T5ttcc"},          .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MGluino",    .pfs={"T5tttt"},          .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MGluino",    .pfs={"T1tttt"},          .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MStop",      .pfs={"T2tt"},            .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MGluino",    .pfs={"T5ttcc","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MGluino",    .pfs={"T5tttt","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MGluino",    .pfs={"T1tttt","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSelectionEfficiency_vs_MLSP_vs_MStop",      .pfs={"T2tt"  ,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});

  sh.AddHistos("evt",   { .fill="SignalSignificance_T5ttcc_vs_MLSP_vs_MGluino",    .pfs={"S"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T5tttt_vs_MLSP_vs_MGluino",    .pfs={"S"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T1tttt_vs_MLSP_vs_MGluino",    .pfs={"S"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T2tt_vs_MLSP_vs_MStop",        .pfs={"S"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T5ttcc_vs_MLSP_vs_MGluino",    .pfs={"S","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T5tttt_vs_MLSP_vs_MGluino",    .pfs={"S","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T1tttt_vs_MLSP_vs_MGluino",    .pfs={"S","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T2tt_vs_MLSP_vs_MStop",        .pfs={"S","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});

  // -------------------------------------------------------------------------
  //                           Signal Region: S and S'
  
  sh.SetHistoWeights({ [this] { return sf_weight['S']; } });

  for (const auto& cut : {"S_ExclmDPhi", "S", "s"}) {
    std::string data = std::string(cut)=="S" ? "Blind" :"JetHT";
    std::string opt  = std::string(cut)=="S" ? o_stk_s : o_stk_d;
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut,"R2Bins"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut,"R2Bins","NJet35"},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    if (std::string(cut) != "S" && std::string(cut) != "S_ExclMR_R2" && std::string(cut) != "s")
      sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={ "Data_MC",cut},                    .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={ "Data_MC",cut,"NJet35"},           .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    // MGlunio/MStop plots
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR",                .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR",                .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MGluinoPoints","GluinoSignalScans","R2Bins",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MStopPoints",  "StopSignalScans"  ,"R2Bins",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={"T5ttcc",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,1700, 0,1400}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={"T5tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={800,2300, 0,1600}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={"T1tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,2300, 0,1600}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MStop",      .pfs={"T2tt"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={150,1200, 0, 650}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",          .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",          .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",          .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",          .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR_vs_MR",          .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR_vs_MR",          .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    //sh.AddHistos("evt",  { .fill="NJet",               .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJet",               .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    // Avg plots
    //sh.AddHistos("evt",  { .fill="AvgMet_vs_MLSP_vs_MGluino",             .pfs={"T5ttcc",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,1700, 0,1400}});
    //sh.AddHistos("evt",  { .fill="AvgMet_vs_MLSP_vs_MGluino",             .pfs={"T5tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={800,2300, 0,1600}});
    //sh.AddHistos("evt",  { .fill="AvgMet_vs_MLSP_vs_MGluino",             .pfs={"T1tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,2300, 0,1600}});
    //sh.AddHistos("evt",  { .fill="AvgMet_vs_MLSP_vs_MStop",               .pfs={"T2tt"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={150,1200, 0, 650}});
    //sh.AddHistos("evt",  { .fill="AvgStopLSPMassDiff_vs_MLSP_vs_MGluino", .pfs={"T5ttcc",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,1700, 0,1400}});
    //sh.AddHistos("evt",  { .fill="AvgStopLSPMassDiff_vs_MLSP_vs_MGluino", .pfs={"T5tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={800,2300, 0,1600}});
    //sh.AddHistos("evt",  { .fill="AvgStopLSPMassDiff_vs_MLSP_vs_MStop",   .pfs={"T2tt"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={150,1200, 0, 650}});
  }

  // Unskimmed plots
  sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"Blind","S_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"Blind","S_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"Blind","S_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"Blind","S_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","S_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","S_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","S_ExclMR_R2","R2Bins"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","S_ExclMR_R2","R2Bins","NJet35"},  .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"Blind","S_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"Blind","S_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"Blind","S_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"Blind","S_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background","S_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background","S_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background","S_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background","S_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background","S_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background","S_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background","S_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background","S_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});

  // N-1 Cut plots
  // S
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","S_Excl3Jet"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","S_Excl0Ele"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","S_Excl0Mu"},               .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","S_Excl0IsoTrk"},           .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","S_Excl1b"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"Blind","S_Excl1W"},                .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","S_Excl1W"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","S_Excl1W"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","S"},                       .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","S_ExclmDPhi"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","S_Excl0Ele",   "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","S_Excl0Mu",    "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","S_Excl0IsoTrk","NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","S_Excl1b",     "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","S_Excl1W",     "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","S_Excl1W",     "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","S",            "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","S_ExclmDPhi",  "NJet35"},  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // S'
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","s_Excl3Jet"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","s_Excl0Ele"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","s_Excl0Mu"},               .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","s_Excl0IsoTrk"},           .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","s_Excl1b"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","s_Excl1W"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","s_Excl1W"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","s"},                       .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","s_ExclInvmDPhi"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","s_Excl0Ele",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","s_Excl0Mu",     "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","s_Excl0IsoTrk", "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","s_Excl1b",      "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","s_Excl1W",      "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","s_Excl1W",      "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","s",             "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","s_ExclInvmDPhi","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});

  // N-1 Weights
  // Number of vertices
  sh.AddHistos("evt",  { .fill="NVtx",        .pfs={Stack,"JetHT","S_3Cuts"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NVtx",        .pfs={Stack,"JetHT","S_6Cuts"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NVtx",        .pfs={Stack,"JetHT","S"},       .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // Same plots with no pile-up reweighting
  sh.SetHistoWeights({ [this] { return w_nm1['S'][1]; } });
  sh.AddHistos("evt",  { .fill="NVtx",        .pfs={Stack,"JetHT","S_3Cuts","NoPUWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NVtx",        .pfs={Stack,"JetHT","S_6Cuts","NoPUWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NVtx",        .pfs={Stack,"JetHT","S"      ,"NoPUWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // No trigger efficiency
  sh.SetHistoWeights({ [this] { return w_nm1['S'][5]; } });
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","S_Excl3JetHLT",   "NoTrigWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",          .pfs={Stack,"JetHT","S_ExclHLT",       "NoTrigWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",         .pfs={Stack,"JetHT","S_ExclHLT",       "NoTrigWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","S_ExclHLT",       "NoTrigWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","S_ExclHLT",       "NoTrigWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","S_ExclHLT",       "NoTrigWeight"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",    .pfs={"Signals_Background","S_ExclHLT","NoTrigWeight"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
  // No scale factors
  sh.SetHistoWeights({ [this] { return w_nm1['S'][6]; } });
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","S_Excl0Ele","NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","S_Excl3Jet","NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",          .pfs={Stack,"JetHT","S",         "NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",         .pfs={Stack,"JetHT","S",         "NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","S",         "NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","S",         "NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","S",         "NoEleSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",    .pfs={"Signals_Background","S",  "NoEleSF"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
  sh.SetHistoWeights({ [this] { return w_nm1['S'][7]; } });
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","S_Excl0Mu", "NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","S_Excl3Jet","NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",          .pfs={Stack,"JetHT","S",         "NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",         .pfs={Stack,"JetHT","S",         "NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","S",         "NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","S",         "NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","S",         "NoMuonSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",    .pfs={"Signals_Background","S",  "NoMuonSF"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
  sh.SetHistoWeights({ [this] { return w_nm1['S'][8]; } });
  sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","S_Excl1b",  "NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","S_Excl3Jet","NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",          .pfs={Stack,"JetHT","S",         "NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",         .pfs={Stack,"JetHT","S",         "NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","S",         "NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","S",         "NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","S",         "NoBTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",    .pfs={"Signals_Background","S",  "NoBTagSF"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
  sh.SetHistoWeights({ [this] { return w_nm1['S'][9]; } });
  sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","S_Excl1W",  "NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","S_Excl3Jet","NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",          .pfs={Stack,"JetHT","S",         "NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",         .pfs={Stack,"JetHT","S",         "NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","S",         "NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","S",         "NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","S",         "NoWTagSF"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",    .pfs={"Signals_Background","S",  "NoWTagSF"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});

  
  // -------------------------------------------------------------------------
  //                       QCD enriched Region: Q and Q'

  sh.SetHistoWeights({ [this] { return sf_weight['Q']; } });

  for (const auto& cut : 
//    {   "Q_Excl3Jet", "Q_Excl0Ele", "Q_Excl0Mu", "Q_Excl0IsoTrk",
//	"Q_Excl0b",   "Q_Excl1aW",  "Q_ExclInvmDPhi0p3", "Q",
//	"Q_1Cuts", "Q_2Cuts", "Q_3Cuts", "Q_4Cuts",
//	"Q_5Cuts", "Q_6Cuts", "Q_7Cuts", "Q_8Cuts" }) {
    {   "Q_Excl0b",   "Q_Excl1aW",  "Q_ExclInvmDPhi0p3", "Q", "q" }) {
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut},               .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
  }
  
  // N-1 plots
  // Q
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","Q_Excl3Jet"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","Q_ExclMR_R2"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","Q_ExclMR_R2","R2Bins"},       .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","Q_ExclMR_R2"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","Q_ExclMR_R2"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","Q_Excl0Ele"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","Q_Excl0Mu"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","Q_Excl0IsoTrk"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","Q_Excl0b"},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NaW",         .pfs={Stack,"JetHT","Q_Excl1aW"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","Q_Excl1aW"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","Q_ExclInvmDPhi0p3"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","Q_Excl0Ele",       "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","Q_Excl0Mu",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","Q_Excl0IsoTrk",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","Q_Excl0b",         "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NaW",         .pfs={Stack,"JetHT","Q_Excl1aW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","Q_Excl1aW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","Q_ExclInvmDPhi0p3","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // Q'
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","q_Excl3Jet"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","q_ExclMR_R2"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","q_ExclMR_R2","R2Bins"},       .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","q_ExclMR_R2"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","q_ExclMR_R2"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","q_Excl0Ele"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","q_Excl0Mu"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","q_Excl0IsoTrk"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","q_Excl0b"},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NaW",         .pfs={Stack,"JetHT","q_Excl1aW"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","q_Excl1aW"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","q_ExclmDPhi"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","q_Excl0Ele",       "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","q_Excl0Mu",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","q_Excl0IsoTrk",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","q_Excl0b",         "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NaW",         .pfs={Stack,"JetHT","q_Excl1aW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","q_Excl1aW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  
  // -------------------------------------------------------------------------
  //                          Top enriched Region: T

  sh.SetHistoWeights({ [this] { return sf_weight['T']; } });

  for (const auto& cut : {"T"}) {
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut},               .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    // Ele or Muon
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"Ele_Muon"},    .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
  }

  // N-1 plots
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","T_Excl3Jet"},               .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","T_ExclMR_R2"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","T_ExclMR_R2","R2Bins"},     .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","T_ExclMR_R2"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","T_ExclMR_R2"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NEle",        .pfs={Stack,"JetHT","T_Excl1LepMT"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NMu",         .pfs={Stack,"JetHT","T_Excl1LepMT"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLep",        .pfs={Stack,"JetHT","T_Excl1LepMT"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","T_Excl1b"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","T_Excl1W"},                .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","T"},                        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","T_ExclmDPhi"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","T_ExclmDPhiMT"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","T_ExclMT"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","T_ExclmDPhiMT"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEle",        .pfs={Stack,"JetHT","T_Excl1LepMT", "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMu",         .pfs={Stack,"JetHT","T_Excl1LepMT", "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLep",        .pfs={Stack,"JetHT","T_Excl1LepMT", "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","T_Excl1b",     "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","T_Excl1W",     "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","T_Excl1W",     "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","T",            "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","T_ExclmDPhi",  "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","T_ExclmDPhiMT","NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","T_ExclMT",     "NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","T_ExclmDPhiMT","NJet35"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // Ele or Muon
  sh.AddHistos("evt",  { .fill="NEle",        .pfs={Stack,"JetHT","T_Excl1LepMT", "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMu",         .pfs={Stack,"JetHT","T_Excl1LepMT", "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLep",        .pfs={Stack,"JetHT","T_Excl1LepMT", "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NBTag",       .pfs={Stack,"JetHT","T_Excl1b",     "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",          .pfs={Stack,"JetHT","T_Excl1W",     "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","T_Excl1W",     "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaRWb",    .pfs={Stack,"JetHT","T",            "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","T_ExclmDPhi",  "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","T_ExclmDPhiMT","Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","T_ExclMT",     "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","T_ExclmDPhiMT","Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});



  // -------------------------------------------------------------------------
  //                          W enriched Region: W


  sh.SetHistoWeights({ [this] { return sf_weight['W']; } });

  for (const auto& cut : {"W"}) {
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut},               .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"log",.ranges={}});
    // Ele or Muon
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"Ele_Muon"},        .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"Ele_Muon"},    .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
  }

  // N-1 plots
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","W_Excl3Jet"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","W_ExclMR_R2"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","W_ExclMR_R2","R2Bins"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","W_ExclMR_R2"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","W_ExclMR_R2"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NEle",        .pfs={Stack,"JetHT","W_Excl1LepMT"},           .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NMu",         .pfs={Stack,"JetHT","W_Excl1LepMT"},           .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLep",        .pfs={Stack,"JetHT","W_Excl1LepMT"},           .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","W_Excl0b"},               .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NmW",          .pfs={Stack,"JetHT","W_Excl1W"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","W_ExclmDPhi"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","W_ExclmDPhiMT"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","W_ExclMT"},               .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","W_ExclmDPhiMT"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEle",        .pfs={Stack,"JetHT","W_Excl1LepMT", "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMu",         .pfs={Stack,"JetHT","W_Excl1LepMT", "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLep",        .pfs={Stack,"JetHT","W_Excl1LepMT", "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","W_Excl0b",     "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NmW",         .pfs={Stack,"JetHT","W_Excl1mW",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","W_Excl1mW",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","W_ExclmDPhi",  "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","W_ExclmDPhiMT","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","W_ExclMT",     "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","W_ExclmDPhiMT","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // Ele or Muon
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","W_Excl3Jet",   "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",          .pfs={Stack,"JetHT","W_ExclMR_R2",  "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",         .pfs={Stack,"JetHT","W_ExclMR_R2",  "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",          .pfs={Stack,"JetHT","W_ExclMR_R2",  "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLooseBTag",  .pfs={Stack,"JetHT","W_Excl0b",     "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NmW",         .pfs={Stack,"JetHT","W_Excl1mW",    "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","W_Excl1mW",    "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","W_ExclmDPhi",  "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","W_ExclmDPhiMT","Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","W_ExclMT",     "Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MT",          .pfs={Stack,"JetHT","W_ExclmDPhiMT","Ele_Muon"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});


  // -------------------------------------------------------------------------
  //                          Z enriched Region: Z

  sh.SetHistoWeights({ [this] { return sf_weight['Z']; } });

  for (const auto& cut : {"Z"}) {
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="METll",                   .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="METll",                   .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"R2Bins","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTRll",                   .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTRll",                   .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2ll",                    .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2ll",                    .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut},                   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"NJet35"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut},               .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_METll",             .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_METll",             .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2ll_vs_METll",           .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2ll_vs_METll",           .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2ll_vs_MR",              .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2ll_vs_MR",              .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut},            .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"},   .cuts={},.draw="COLZ",.opt=o_1or2d_s+"Log",.ranges={}});
    // Ele or Muon
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="METll",                   .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTRll",                   .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2ll",                    .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,"JetHT",cut,"2Ele_2Muon"},      .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={"Data_MC","JetHT",cut,"2Ele_2Muon"},  .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
  }

  // N-1 plots
  sh.AddHistos("evt",  { .fill="NJet",          .pfs={Stack,"JetHT","Z_Excl3Jet"},                 .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",            .pfs={Stack,"JetHT","Z_ExclMR_R2ll"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",            .pfs={Stack,"JetHT","Z_ExclMR_R2ll","R2llBins"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTRll",         .pfs={Stack,"JetHT","Z_ExclMR_R2ll"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2ll",          .pfs={Stack,"JetHT","Z_ExclMR_R2ll"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NEleLoose",     .pfs={Stack,"JetHT","Z_ExclMR_R2ll2Lep"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NMuLoose",      .pfs={Stack,"JetHT","Z_ExclMR_R2ll2Lep"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLepLoose",     .pfs={Stack,"JetHT","Z_ExclMR_R2ll2Lep"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NLooseBTag",    .pfs={Stack,"JetHT","Z"},                          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NBTag",         .pfs={Stack,"JetHT","Z"},                          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NW",            .pfs={Stack,"JetHT","Z_Excl1mW"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NmW",           .pfs={Stack,"JetHT","Z_Excl1mW"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",            .pfs={Stack,"JetHT","Z_ExclmDPhill"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NmW",           .pfs={Stack,"JetHT","Z_ExclmDPhill"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi",   .pfs={Stack,"JetHT","Z_ExclmDPhill"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhill", .pfs={Stack,"JetHT","Z_ExclmDPhill"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaPhiLLMET", .pfs={Stack,"JetHT","Z"},                          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaPhiLLMET", .pfs={Stack,"JetHT","Z_ExclmDPhill"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaPhiLLJet", .pfs={Stack,"JetHT","Z"},                          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="DeltaPhiLLJet", .pfs={Stack,"JetHT","Z_ExclmDPhill"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="Mll",           .pfs={Stack,"JetHT","Z_ExclMll"},                  .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="Mll",           .pfs={Stack,"JetHT","Z_ExclmDPhillMll"},           .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleLoose",     .pfs={Stack,"JetHT","Z_ExclMR_R2ll2Lep","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuLoose",      .pfs={Stack,"JetHT","Z_ExclMR_R2ll2Lep","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLepLoose",     .pfs={Stack,"JetHT","Z_ExclMR_R2ll2Lep","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLooseBTag",    .pfs={Stack,"JetHT","Z",                "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NBTag",         .pfs={Stack,"JetHT","Z",                "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",            .pfs={Stack,"JetHT","Z_Excl1mW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NmW",           .pfs={Stack,"JetHT","Z_Excl1mW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",       .pfs={Stack,"JetHT","Z_Excl1mW",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi",   .pfs={Stack,"JetHT","Z_ExclmDPhill",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhill", .pfs={Stack,"JetHT","Z_ExclmDPhill",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="Mll",           .pfs={Stack,"JetHT","Z_ExclMll",        "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  // Ele or Muon
  sh.AddHistos("evt",  { .fill="NJet",          .pfs={Stack,"JetHT","Z_Excl3Jet",   "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",            .pfs={Stack,"JetHT","Z_ExclMR_R2ll","2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTRll",         .pfs={Stack,"JetHT","Z_ExclMR_R2ll","2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2ll",          .pfs={Stack,"JetHT","Z_ExclMR_R2ll","2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NLooseBTag",    .pfs={Stack,"JetHT","Z",            "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NBTag",         .pfs={Stack,"JetHT","Z",            "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NW",            .pfs={Stack,"JetHT","Z_Excl1mW",    "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NmW",           .pfs={Stack,"JetHT","Z_Excl1mW",    "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("mW",   { .fill="mWTau21",       .pfs={Stack,"JetHT","Z_Excl1mW",    "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi",   .pfs={Stack,"JetHT","Z_ExclmDPhill","2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhill", .pfs={Stack,"JetHT","Z_ExclmDPhill","2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="Mll",           .pfs={Stack,"JetHT","Z_ExclMll",    "2Ele_2Muon"},   .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});

  // ----------------------------------------------------------------------------------------------
  //                                        TOP ANALYSIS
  //-----------------------------------------------------------------------------------------------

  o_stk_d = "LogSumw2Stack5AddRatioTwoCol57AddIntApproval17";
  o_stk_s = "LogSumw2Stack5AddRatioTwoCol57AddIntApproval47";
  o_1or2d_d = "Sumw2Approval17";
  o_1or2d_s = "Sumw2Approval47";
  o_norm_d = "Sumw2NormApproval17";
  o_norm_s = "Sumw2NormApproval47";

  // -------------------------------------------------------------------------
  //                           Selected tops
  
  for (auto region : {'S', 's', 'T','W','Q', 'q', 'Z', 't'}) {
    sh.SetHistoWeights({ [this,region] { return sf_weight[region]; } });
    std::string cut1 = std::string(1,region);
    std::string cut2 = cut1;
    if      (region=='S'||region=='s'||region=='T') cut2 += "_Excl1b1W";
    else if (region=='Q'||region=='q') cut2 += "_Excl0b1aW";
    else if (region=='W') cut2 += "_Excl0b1mW";
    else if (region=='Z') cut2 += "_Excl1mW";
    else if (region=='t') cut2 += "_Excl1Top";
    std::vector<std::string> showdata = {"JetHT"};
    if (region=='S'||region=='t') showdata.push_back("Blind");
    for (auto cut : { cut1, cut2 }) {
      for (auto data : showdata ) {
	std::string opt = (data=="Blind") ? o_stk_s : o_stk_d;
	sh.AddHistos("evt",  { .fill="NHadTopTag",         .pfs={Stack,data,cut},                          .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
      }
      sh.AddHistos("evt",  { .fill="NHadTopTag",         .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={0,0, 0,1, 0.32,0.90}});
      sh.AddHistos("evt",  { .fill="NHadTopTag",         .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={0,0, 0,1, 0.32,0.90}});
    }
  }

  // -------------------------------------------------------------------------
  //                             Top GenInfo

  sh.SetHistoWeights({ [this] { return sf_weight['t']; } });
  sh.AddHistos("gen top", { .fill="GenTopPt",                                      .pfs={"TT_SignalPoints"},   .cuts={}, .draw=d,    .opt=o_1or2d_s+"Norm",.ranges={0,2000, 0,0.2, 0.6,0.9}});

  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MGluino", .pfs={"T5ttcc"},          .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MGluino", .pfs={"T5tttt"},          .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MGluino", .pfs={"T1tttt"},          .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MStop",   .pfs={"T2tt"},            .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MGluino", .pfs={"T5ttcc","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MGluino", .pfs={"T5tttt","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MGluino", .pfs={"T1tttt","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="TopSignalSelectionEfficiency_vs_MLSP_vs_MStop",   .pfs={"T2tt"  ,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});

  sh.AddHistos("evt",   { .fill="SignalSignificance_T5ttcc_vs_MLSP_vs_MGluino",    .pfs={"t"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T5tttt_vs_MLSP_vs_MGluino",    .pfs={"t"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T1tttt_vs_MLSP_vs_MGluino",    .pfs={"t"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T2tt_vs_MLSP_vs_MStop",        .pfs={"t"},               .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T5ttcc_vs_MLSP_vs_MGluino",    .pfs={"t","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,1700, 0,1400, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T5tttt_vs_MLSP_vs_MGluino",    .pfs={"t","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={800,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T1tttt_vs_MLSP_vs_MGluino",    .pfs={"t","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={600,2300, 0,1600, 0,0, 0.02,0.95}});
  sh.AddHistos("evt",   { .fill="SignalSignificance_T2tt_vs_MLSP_vs_MStop",        .pfs={"t","NJet35"},      .cuts={},.draw="COLZ",.opt=o_1or2d_s, .ranges={150,1200, 0, 650, 0,0, 0.02,0.95}});

  // -------------------------------------------------------------------------
  //                  Top-tag Signal Region: t
  
  sh.SetHistoWeights({ [this] { return sf_weight['t']; } });

  for (const auto& cut : {"t"}) {
    std::string data = std::string(cut)=="t" ? "Blind" :"JetHT";
    std::string opt  = std::string(cut)=="t" ? o_stk_s : o_stk_d;
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut,"R2Bins"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,data,cut,"R2Bins","NJet35"},  .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background",cut,"NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="NJet",                    .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    if (std::string(cut) != "t")
      sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NJetAK8",                 .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="HTBins",                  .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,data,cut},                    .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins",           .pfs={Stack,data,cut,"NJet35"},           .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={ "Data_MC",cut},                    .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    sh.AddHistos("evt",  { .fill="Jet1AK8PtBins_vs_HTBins", .pfs={ "Data_MC",cut,"NJet35"},           .cuts={},.draw="COLZ",.opt=o_1or2d_d,.ranges={}});
    // MGlunio/MStop plots
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR",                .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR",                .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MGluinoPoints","GluinoSignalScans",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MStopPoints",  "StopSignalScans"  ,cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MGluinoPoints","GluinoSignalScans","R2Bins",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MStopPoints",  "StopSignalScans"  ,"R2Bins",cut}, .cuts={},.draw=d,.opt=o_norm_s,.ranges={}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={"T5ttcc",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,1700, 0,1400}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={"T5tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={800,2300, 0,1600}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={"T1tttt",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={600,2300, 0,1600}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MStop",      .pfs={"T2tt"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_s,.ranges={150,1200, 0, 650}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",          .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR_vs_MET",          .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",          .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MET",          .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR_vs_MR",          .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    sh.AddHistos("evt",  { .fill="MTR_vs_MR",          .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
    //sh.AddHistos("evt",  { .fill="NJet",               .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJet",               .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={"GluinoSignalScans","MGluinoPoints",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
    //sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={"StopSignalScans",  "MStopPoints"  ,cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_stk});
  }

  // Unskimmed plots
  sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"Blind","t_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="HT",                      .pfs={Stack,"Blind","t_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"Blind","t_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MET",                     .pfs={Stack,"Blind","t_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","t_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","t_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","t_ExclMR_R2","R2Bins"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR",                      .pfs={Stack,"Blind","t_ExclMR_R2","R2Bins","NJet35"},  .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"Blind","t_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MTR",                     .pfs={Stack,"Blind","t_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"Blind","t_ExclMR_R2"},                    .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="R2",                      .pfs={Stack,"Blind","t_ExclMR_R2","NJet35"},           .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background","t_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="MR_vs_MET",               .pfs={"Signals_Background","t_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background","t_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MET",               .pfs={"Signals_Background","t_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background","t_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="R2_vs_MR",                .pfs={"Signals_Background","t_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background","t_ExclMR_R2"},          .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});
  sh.AddHistos("evt",  { .fill="HT_vs_MR",                .pfs={"Signals_Background","t_ExclMR_R2","NJet35"}, .cuts={},.draw="COLZ",.opt=o_1or2d_d+"Log",.ranges={}});

  // N-1 Cut plots
  sh.AddHistos("evt",  { .fill="NJet",        .pfs={Stack,"JetHT","t_Excl3Jet"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","t_Excl0Ele"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","t_Excl0Mu"},              .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","t_Excl0IsoTrk"},          .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NHadTopTag",  .pfs={Stack,"Blind","t_Excl1Top"},             .cuts={},.draw=d,.opt=o_stk_s,.ranges=r_stk});
  //sh.AddHistos("evt",  { .fill="NHadTopTag",  .pfs={Stack,"JetHT","t_Excl1Top"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","t_Excl1Top"},             .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","t_ExclmDPhi"},            .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NEleVeto",    .pfs={Stack,"JetHT","t_Excl0Ele",   "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NMuVeto",     .pfs={Stack,"JetHT","t_Excl0Mu",    "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NIsoTrk",     .pfs={Stack,"JetHT","t_Excl0IsoTrk","NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NHadTopTag",  .pfs={Stack,"JetHT","t_Excl1Top",   "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  //sh.AddHistos("mW",   { .fill="mWTau21",     .pfs={Stack,"JetHT","t_Excl1Top",   "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="MinDeltaPhi", .pfs={Stack,"JetHT","t_ExclmDPhi",  "NJet35"}, .cuts={},.draw=d,.opt=o_stk_d,.ranges=r_stk});

  /*
  // --------------------------------------------------------------------------
  //                                   AK8 Jets

  // --------------------------------------
  //            Jet kinematics
  //         (pt, eta, phi, mass)

  bool Systematics_Only = false;

  if (!Systematics_Only) {

  // Data-MC Comparison - background composition
  for (const auto& cut : all_cuts) {
  // Jet quantities
  sh.AddHistos("evt",        { .fill="NJet",            .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt",        { .fill="NHadTopTag",      .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"JetHT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPhi",          .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPhi",          .pfs={Stack,cut,"JetHT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetEta",          .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetEta",          .pfs={Stack,cut,"JetHT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetTau32",        .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetTau32",        .pfs={Stack,cut,"JetHT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetSoftDropMass", .pfs={Stack,cut,"JetHT"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetSoftDropMass", .pfs={Stack,cut,"JetHT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  // Subjet BTags
  sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,"NSubjetBTag",cut,"JetHT"},                       .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,"NSubjetBTag",cut,"JetHT","PtOrder"},       .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"JetHT","JetPassSubjetBTag"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"JetHT","JetPassSubjetBTag","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetBTagCSV",      .pfs={Stack,cut,"JetHT"},                                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("AK8",      { .fill="JetBTagCSV",      .pfs={Stack,cut,"JetHT","PtOrder"},                     .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
      
  // Cutflow
  sh.AddHistos("cutflow", { .fill="Cutflow", .pfs={Stack,cut,"JetHT"}, .cuts={"PassAnaSelection"}, .draw=d, .opt=o_stk_d, .ranges={0,0, 0,0} });
  sh.AddHistos("cutflow", { .fill="Cutflow", .pfs={"Background", cut,"JetHT"}, .cuts={"PassAnaSelection"}, .draw=d, .opt="Sumw2Log",     .ranges={0,0, 0,0} });
  }
    
  // N-1 plots
  sh.AddHistos("evt", { .fill="Jet1Eta",     .pfs={Stack,"AllCuts_ExclJet1Eta",   "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet2Eta",     .pfs={Stack,"AllCuts_ExclJet2Eta",   "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet1Pt",      .pfs={Stack,"AllCuts_ExclJet1Pt",    "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet2Pt",      .pfs={Stack,"AllCuts_ExclJet2Pt",    "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet1Mass",    .pfs={Stack,"AllCuts_ExclJet1Mass",  "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet2Mass",    .pfs={Stack,"AllCuts_ExclJet2Mass",  "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet1BTagCSV", .pfs={Stack,"AllCuts_ExclJet1BTag",  "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet2BTagCSV", .pfs={Stack,"AllCuts_ExclJet2BTag",  "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet1Tau32",   .pfs={Stack,"AllCuts_ExclJet1Tau32", "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="Jet2Tau32",   .pfs={Stack,"AllCuts_ExclJet2Tau32", "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="DPhi",        .pfs={Stack,"AllCuts_ExclDeltaPhi",  "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="R",           .pfs={Stack,"AllCuts_ExclR",         "JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
    
  // --------------------------------------------------------------------------
  //                                 Trigger

  sh.AddHistos("evt", { .fill="GenHt",      .pfs={"AllSamples"},                           .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="AK4Ht", .pfs={"AllSamples"},                           .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="AK8Ht", .pfs={"AllSamples"},                           .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="GenHt",      .pfs={"AllSamples","JetHT"},              .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="AK4Ht", .pfs={"AllSamples","JetHT"},              .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="AK8Ht", .pfs={"AllSamples","JetHT"},              .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="GenHt",      .pfs={"AllSamples","Pass12Cuts","JetHT"},  .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });

  // Event composition - under different cuts
  for (const auto& cut : all_cuts) {
  sh.AddHistos("evt", { .fill="SumPt",      .pfs={Stack,cut},              .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="SumPt",      .pfs={Stack,cut,"JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="AK8Ht", .pfs={Stack,cut},              .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="AK8Ht", .pfs={Stack,cut,"JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="AK4Ht", .pfs={Stack,cut},              .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt", { .fill="AK4Ht", .pfs={Stack,cut,"JetHT"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });

  sh.AddHistos("evt", { .fill="MetPt",       .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="MinDeltaPhi", .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="MR",      .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MTR",     .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="R",       .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt", { .fill="DPhi",    .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt", { .fill="DEta",    .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt", { .fill="DR",      .pfs={Stack,cut,"JetHT"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  }
  sh.AddHistos("evt",   { .fill="SumPt",      .pfs={Stack,"Pass12Cuts","JetHT","PassDPhiCut"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt",   { .fill="AK8Ht", .pfs={Stack,"Pass12Cuts","JetHT","PassDPhiCut"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("evt",   { .fill="AK4Ht", .pfs={Stack,"Pass12Cuts","JetHT","PassDPhiCut"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });

  // N-1 plot for all cuts
    

  // --------------------------------------
  //               Efficiencies

  // 1D plots - sumpt - cutflow
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass11Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass11Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass11Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });

  // Check other triggers
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass5Cuts"},  .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass7Cuts"},  .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass9Cuts"},  .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass11Cuts"}, .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass12Cuts"}, .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass12Cuts","PassDPhiCut"}, .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPt",          .pfs={"Triggers","Pass5Cuts"},  .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPt",          .pfs={"Triggers","Pass7Cuts"},  .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPt",          .pfs={"Triggers","Pass9Cuts"},  .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPt",          .pfs={"Triggers","Pass11Cuts"}, .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPt",          .pfs={"Triggers","Pass12Cuts"}, .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPt",          .pfs={"Triggers","Pass12Cuts","PassDPhiCut"}, .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });

  // 1 Bin
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0.95,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0.95,1} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPtOneBin",          .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });
  sh.AddHistos("evt", { .fill="HLTEff_PFHT800or900_vs_SumPtOneBin",          .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0, 0,1, 0.45,0.45} });

  // 2D plots
  // sumpt vs mass1
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass1Cuts"},                .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass3Cuts"},                .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts"},                .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass7Cuts"},                .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass9Cuts"},                .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass11Cuts"},               .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });
  sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass12Cuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={0,200, 400,1200, 0,1} });

  
  // Main variables (Shape, area, Data-MC agreement)
  // Hadronic top selection
  //   // No Cut
  //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
  //   sh.AddHistos("AK8", { .fill="JetTau2",         .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetTau3",         .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} }); // Note
  //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note
  //   // Apply 1 Cut
  //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack,"JetMassCut"},  .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack,"JetTau32Cut"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack,"JetPtCut"},    .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack,"JetMassCut"},  .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack,"JetTau32Cut"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack,"JetPtCut"},    .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  //   // Apply 2 Cuts (N-1)
  //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack,"JetMassCut","JetTau32Cut"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack,"JetPtCut","JetMassCut"},    .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });  // Note
  //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack,"JetTau32Cut","JetPtCut"},   .cuts={},  .draw=d, .opt=o_stk_d, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note

  // 2D
  // No Cut
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",     .pfs={"AllSamples"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetPt",        .pfs={"AllSamples"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("AK8", { .fill="JetPt_vs_JetTau3",                .pfs={"AllSamples"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });
  // Apply 1 Cut
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",     .pfs={"AllSamples","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetPt",        .pfs={"AllSamples","JetTau32Cut"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("AK8", { .fill="JetPt_vs_JetTau32",               .pfs={"AllSamples","JetMassCut"},  .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });
  // Apply 2 Cuts
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",     .pfs={"AllSamples","SubjetBTagCut","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetPt",        .pfs={"AllSamples","SubjetBTagCut","JetTau32Cut"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("AK8", { .fill="JetPt_vs_JetTau32",               .pfs={"AllSamples","SubjetBTagCut","JetMassCut"},  .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });

  // Same plots, but use Gen Particle Truth
  sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} }); // No Cut
  sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau1",                               .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau2",                               .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau3",                               .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau21",                              .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau31",                              .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("AK8", { .fill="JetNSubJets",                           .pfs={"JetGenTruth","Signals_Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals_Background","JetMassCut"},  .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} }); // 1 Cut
  sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals_Background","JetTau32Cut"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals_Background","JetPtCut"},    .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals_Background","JetMassCut"},  .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals_Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals_Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals_Background","JetTau32Cut"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals_Background","JetPtCut"},    .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals_Background","JetMassCut","JetTau32Cut"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} }); // 2 Cuts
  sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals_Background","JetPtCut","JetMassCut"},    .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals_Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals_Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals_Background","JetTau32Cut","JetPtCut"},   .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",           .pfs={"JetGenTruth","Signals_Background"},             .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 2D No Cut
  sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",           .pfs={"JetGenTruth","Signals_Background","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 1 Cut

  // Fraction of merged sub-jets
  sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"TT_Signal"},                                   .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"TT_Signal","JetMatchedGenTopType"},            .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"TT_Signal"},                                   .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"TT_Signal","JetMatchedGenTopType"},            .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  
  // Top Tag/Finding Efficiency
  sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"TT_Signal"},              .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });
  sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"TT_Signal","TopSizeCut"}, .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });
  sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"TT_Signal"},              .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });
  sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"TT_Signal","TopSizeCut"}, .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });

  // --------------------------------------------------------------------------
  //                              Gen particles

  // Jet Finding Efficiency
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"TT_Signal"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });  // Note
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"TT_Signal"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"TT_Signal","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"TT_Signal","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"TT_Signal"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"TT_Signal"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"TT_Signal","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"TT_Signal","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });

  // --------------------------------------------------------------------------
  //                                   HLT

  // // Trigger Efficiencies
  //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPt",       .pfs={"Data_MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={800,2000, 0.8,1.05, 0.15,0.95} });
  //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPtOneBin", .pfs={"Data_MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0,      0.8,1.05, 0.15,0.95} });
  //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtBins",        .pfs={"Data_MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,2000,   0.8,1.05, 0.15,0.95} });
  //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtOneBin",      .pfs={"Data_MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0,      0.8,1.05, 0.15,0.95} });

  // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPt",       .pfs={"Background_Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={800,2000, 0.8,1.05, 0.15,0.95} });
  // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPtOneBin", .pfs={"Background_Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0,      0.8,1.05, 0.15,0.95} });
  // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtBins",        .pfs={"Background_Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,2000,   0.8,1.05, 0.15,0.95} });
  // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtOneBin",      .pfs={"Background_Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt=o_1or2d_d, .ranges={0,0,      0.8,1.05, 0.15,0.95} });

  // --------------------------------------------------------------------------
  //                                MET/Razor

  // MET, Razor Variables
  sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
  sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts","RBands"},           .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts","RBands","DPhiBands"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("baseline events", { .fill="MR",               .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("baseline events", { .fill="MR",               .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("baseline events", { .fill="MR",               .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("baseline events", { .fill="MTR",              .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("baseline events", { .fill="MTR",              .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("baseline events", { .fill="MTR",              .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1e-3,1e6} });
  
  // 2D Correlation plots
  sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples","DPhiBands"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background","DPhiBands"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgR_vs_DPhi",        .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="AvgTTHadR_vs_DPhi",   .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMR_vs_DPhi",       .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="AvgTTHadMR_vs_DPhi",  .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_DPhi",      .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","DPhiBands","AllSamples"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Tau32Cuts","AllSamples"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgR_vs_DPhi",        .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="AvgTTHadR_vs_DPhi",   .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMR_vs_DPhi",       .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events", { .fill="AvgTTHadMR_vs_DPhi",  .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_DPhi",      .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","DPhiBands","Background"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Tau32Cuts","Background"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  
  // Signal selection (Apply loose masstag selection)
  sh.AddHistos("baseline events",   { .fill="R",                 .pfs={Stack},                     .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={Stack},                     .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("baseline events",   { .fill="R",                 .pfs={Stack,"Tau32Cuts"},           .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={Stack,"Tau32Cuts"},           .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("baseline events",   { .fill="R",                 .pfs={Stack,"Tau32Cuts","DPhiBands"}, .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={Stack,"Tau32Cuts","RBands"},  .cuts={}, .draw=d, .opt=o_stk_d, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  
  // 3D Plots to get best signal cuts (Maximize Smin) --> input for B2GAnalyzer
  sh.AddHistos("baseline events",   { .fill="MR_vs_DPhiFine_vs_RFine",           .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events",   { .fill="MR_vs_DPhiFine_vs_TTHadRFine",      .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events",   { .fill="TTHadMR_vs_DPhiFine_vs_RFine",      .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("baseline events",   { .fill="HtAll_vs_DPhiFine_vs_RFine",        .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  
  // R plots - Distributions for All samples, All cut combinations, Ratios
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  // ABCD regions - DPhi<DPHI_CUT, define regions by: R and Ntop
  sh.AddHistos("baseline events",   { .fill="R",         .pfs={"ABCD,Signals","DPhiBands"},      .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="R",         .pfs={"ABCD","DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="R",         .pfs={"ABCD","DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("baseline events",   { .fill="R",         .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  //sh.AddHistos("baseline events",   { .fill="R",         .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
  sh.AddHistos("baseline events",   { .fill="RFine",     .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RFine",     .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","Tau32Cuts","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  
  // DPhi plots - Distributions for All samples, All cut combinations, Ratios
  sh.AddHistos("baseline events",   { .fill="DPhi",      .pfs={"ABCD,Signals","DPhiBands"},      .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands", "AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands", "Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",  "AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",  "Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  // Alternative ABCD regions - Ntop==2, define regions by: R and DPhi
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands",   "Tau32Cuts","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands",   "Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",    "Tau32Cuts","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",    "Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","RBins",    "AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","RBins",    "Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });

  // ABCD method + systematics
  //sh.AddHistos("baseline events syst",   { .fill="Counts_vs_DPhiBins",  .pfs={"RBins",    "Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  //sh.AddHistos("baseline events syst",   { .fill="Counts_vs_R",         .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  //sh.AddHistos("baseline events syst",   { .fill="Counts_vs_R",         .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2

  // Signal specific plots
  sh.AddHistos("baseline events",   { .fill="MLSP_vs_MGluino", .pfs={"DPhiBands","Tau32Cuts","RBands","Signals"}, .cuts={}, .draw="COLZ", .opt=o_1or2d_d, .ranges={} });
  } // End all plots

  // --------------------------------------------------------------------------
  //                             Systematics

  // Data-MC Comparison - background composition
  std::vector<std::string> cuts = { "Pass5Cuts", "Pass11Cuts", "Pass12Cuts", "Pass13Cuts", "Pass15Cuts", "Pass16Cuts" };
  for (auto cut : cuts) { 
  sh.AddHistos("all events syst",   { .fill="Counts_vs_NJet",            .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("all events syst",   { .fill="Counts_vs_NHadTopTag",      .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetPt",           .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetPt",           .pfs={Stack,cut,"JetHT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetPhi",          .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetPhi",          .pfs={Stack,cut,"JetHT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetEta",          .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetEta",          .pfs={Stack,cut,"JetHT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetTau32",        .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetTau32",        .pfs={Stack,cut,"JetHT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetSoftDropMass", .pfs={Stack,cut,"JetHT"},                 .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  sh.AddHistos("jetsAK8 syst", { .fill="Counts_vs_JetSoftDropMass", .pfs={Stack,cut,"JetHT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk_d, .ranges=r_stk });
  }

  // Background estimation
  //sh.AddHistos("baseline events",      { .fill="HtAll_vs_DPhiFine_vs_RFine",      .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  sh.AddHistos("baseline events",      { .fill="AK8Ht_vs_DPhiFine_vs_RFine", .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  sh.AddHistos("baseline events",      { .fill="DPhiBins",           .pfs={"RBins",     "Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_DPhiBins", .pfs={"RBins",     "Tau32Cuts","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  sh.AddHistos("baseline events",      { .fill="R",                  .pfs={"Tau32Cuts", "DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_R",        .pfs={"Tau32Cuts", "DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  sh.AddHistos("baseline events",      { .fill="R",                  .pfs={"Tau32Cuts", "DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_R",        .pfs={"Tau32Cuts", "DPhiBands","Background"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_RFine",    .pfs={"Tau32Cuts", "DPhiBands","AllSamples"}, .cuts={}, .draw=d, .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // python ABCD

  */
}

//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_analysis_histos(DataStruct& d, const unsigned int& syst_index, const double& weight)
{
  if (syst_index == 0) { // Default (no systematic variation)
    //__________________________________
    //         Fill Smarthistos
    //__________________________________
    while(d.jetsAK4.Loop()) if (passLooseJet     [d.jetsAK4.it]) sh.Fill("AK4");
    while(d.jetsAK8.Loop()) if (passLooseJetAK8  [d.jetsAK8.it]) sh.Fill("AK8");
    while(d.jetsAK4.Loop()) if (passMediumBTag   [d.jetsAK4.it]) sh.Fill("b");
    while(d.jetsAK4.Loop()) if (passLooseBTag    [d.jetsAK4.it]) sh.Fill("b loose");
    while(d.jetsAK8.Loop()) if (passWPreTag      [d.jetsAK8.it]) sh.Fill("mW");
    while(d.jetsAK8.Loop()) if (passTightWAntiTag[d.jetsAK8.it]) sh.Fill("aW");
    while(d.jetsAK8.Loop()) if (passTightWTag    [d.jetsAK8.it]) sh.Fill("W");
    while(d.ele.Loop())     if (passEleSelect    [d.ele.it])     sh.Fill("ele");
    while(d.ele.Loop())     if (passEleVeto      [d.ele.it])     sh.Fill("ele veto");
    while(d.mu.Loop())      if (passMuSelect     [d.mu.it])      sh.Fill("mu");
    while(d.mu.Loop())      if (passMuVeto       [d.mu.it])      sh.Fill("mu veto");
    while(d.gen.Loop())     if (passGenHadW      [d.gen.it])     sh.Fill("gen W");
    while(d.gen.Loop())     if (passGenHadTop    [d.gen.it])     sh.Fill("gen top");
    sh.Fill("evt");
  }
  sh.Fill("syst");
}

void
Analysis::load_analysis_histos(std::string inputfile)
{
  sh.Add(inputfile.c_str());
}

void
Analysis::save_analysis_histos(bool draw=0)
{
  if (draw) sh.DrawPlots();
  sh.Write();
}
