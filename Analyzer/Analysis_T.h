#ifndef VER
#define VER 2
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
  int NJetAK8 = 0;
  while(data.jetsAK8.Loop()) {
    size_t i = data.jetsAK8.it;
    // pt cut intentionally removed to accept all jets for systematics
    if ( data.jetsAK8.looseJetID[i] == 1 &&
         std::abs(data.jetsAK8.Eta[i])  <  JET_AK8_ETA_CUT ) {
      NJetAK8++;
    }
  }
  if (!(NJetAK8>=1)) return 0;
  if (!(data.evt.R2>=0.04)) return 0;

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


enum SCuts { S_3Jet, S_MRR2, S_HLT, S_0Ele, S_0Mu, S_0IsoTrk, S_1Top, S_dPhiRazor };
enum QCuts { Q_3Jet, Q_MRR2, Q_HLT, Q_0Ele, Q_0Mu, Q_0IsoTrk, Q_0b, Q_1aTop, Q_dPhiRazor };
enum TCuts { T_3Jet, T_MRR2, T_HLT, T_1Lep, W_0b,                  T_1Top, T_dPhiRazor, T_MT};
enum WCuts { W_3Jet, W_MRR2, W_HLT, W_1Lep,                   W_1mTop, W_dPhiRazor, W_MT};



void
Analysis::define_selections(const DataStruct& d)
{
  analysis_cuts.clear();

// Define here cuts that are common in all Signal/Control regions
  // MET Filters, etc. are already applied in AnalysisBase.h, See baseline_cuts


  // cut0: signal mass region
bool isT2tt = TString(sample).Contains("T2tt");
  if(isT2tt){ 
  baseline_cuts.push_back({ .name="signal_mass_selection",   .func = [&d]{ 
            //return d.evt.SUSY_Gluino_Mass==mGluino[num] && d.evt.SUSY_LSP_Mass==mLSP[num];
            //return d.evt.SUSY_Gluino_Mass==2000 && d.evt.SUSY_LSP_Mass==300;
            return d.evt.SUSY_Stop_Mass==850 && d.evt.SUSY_LSP_Mass==100;
	    } });}
bool isT5tttt = TString(sample).Contains("T5tttt");
  if(isT5tttt){ 
  baseline_cuts.push_back({ .name="signal_mass_selection",   .func = [&d]{ 
            //return d.evt.SUSY_Gluino_Mass==mGluino[num] && d.evt.SUSY_LSP_Mass==mLSP[num];
            return d.evt.SUSY_Gluino_Mass==1400 && d.evt.SUSY_LSP_Mass==300;
            //return d.evt.SUSY_Stop_Mass==800 && d.evt.SUSY_LSP_Mass==100;
	    } });}
bool isT1tttt = TString(sample).Contains("T1tttt");
  if(isT1tttt){ 
  baseline_cuts.push_back({ .name="signal_mass_selection",   .func = [&d]{ 
            return d.evt.SUSY_Gluino_Mass==1400 && d.evt.SUSY_LSP_Mass==300;
	    } });}
bool isT5ttcc = TString(sample).Contains("T5ttcc");
  if(isT5ttcc){ 
  baseline_cuts.push_back({ .name="signal_mass_selection",   .func = [&d]{ 
            //return d.evt.SUSY_Gluino_Mass==mGluino[num] && d.evt.SUSY_LSP_Mass==mLSP[num];
            return d.evt.SUSY_Gluino_Mass==1400 && d.evt.SUSY_LSP_Mass==300;
            //return d.evt.SUSY_Stop_Mass==800 && d.evt.SUSY_LSP_Mass==100;
} });}


 // P: Preselection region
  analysis_cuts['P'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['P'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['P'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['P'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['P'].push_back({ .name="HLT", .func = [this,&d] { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});

  // Q: QCD enriched control sample
 analysis_cuts['Q'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['Q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['Q'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['Q'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['Q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['Q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['Q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['Q'].push_back({ .name="1aTop",      .func = []    { return nHadTop0BAntiTag>=1;              }});
//analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                  }});
  analysis_cuts['Q'].push_back({ .name="DPhiHemi",   .func = []        { return dPhiRazor>=2.8;                  }});

   // Q': Dphi Control region of QCD enriched sample
  analysis_cuts['q'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['q'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['q'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['q'].push_back({ .name="1aTop",      .func = []    { return nHadTop0BAntiTag>=1;              }});
//analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['q'].push_back({ .name="DPhiHemi",   .func = []        { return dPhiRazor<2.8;                 }});

  // T: Top enriched control sample
  analysis_cuts['T'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['T'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['T'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['T'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['T'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['T'].push_back({ .name="1Top",       .func = []    { return nHadTopTag>=1;                    }});
//analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['T'].push_back({ .name="DPhiHemi",      .func = []    { return dPhiRazor<2.8;                    }});
  analysis_cuts['T'].push_back({ .name="MT", .func = [] { return MT_vetolep<100; }});

  // W: W enriched control sample
 analysis_cuts['W'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['W'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['W'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['W'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['W'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['W'].push_back({ .name="1mTop",      .func = []    { return nHadTop0BMassTag>=1;              }});
//analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['W'].push_back({ .name="DPhiHemi", .func = [] { return dPhiRazor<2.8; }});
  analysis_cuts['W'].push_back({ .name="MT",         .func = []        { return MT_vetolep>=30 && MT_vetolep<100; }});

 // S: Signal region
  analysis_cuts['S'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['S'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['S'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['S'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['S'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['S'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['S'].push_back({ .name="1Top",       .func = []    { return nHadTopTag>=1;                    }});
//analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['S'].push_back({ .name="DPhiHemi",   .func = []        { return dPhiRazor<2.8;                 }});

  // S': DPhi Control region of Signal region
  analysis_cuts['s'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['s'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['s'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['s'].push_back({ .name="R2",         .func = [&d]  { return d.evt.R2>=0.08;                   }});
  analysis_cuts['s'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['s'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['s'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['s'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['s'].push_back({ .name="1Top",       .func = []    { return nHadTopTag>=1;                    }});
//analysis_cuts['s'].push_back({ .name="InvmDPhi",   .func = []    { return minDeltaPhi<0.5;                  }});
  analysis_cuts['s'].push_back({ .name="DPhiHemi",   .func = []        { return dPhiRazor>=2.8;                  }});

  // Z: Z->ll enriched control region
  analysis_cuts['Z'].push_back({ .name="1JetAK8", .func = [] { return nJetAK8>=1; }});
  analysis_cuts['Z'].push_back({ .name="3Jet",       .func = []        { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['Z'].push_back({ .name="MRR2ll",     .func = [&d]      { return d.evt.MR>=800 && R2_ll>=0.08;     }});
  analysis_cuts['Z'].push_back({ .name="HLT",        .func = [this,&d] { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Z'].push_back({ .name="2Lep",       .func = []        { return (nEleSelect==2&&nMuVeto==0)||(nMuSelect==2&&nEleVeto==0); }});
  analysis_cuts['Z'].push_back({ .name="OppCharge",  .func = [&d]      { 
				   if (nEleSelect==2) return (d.ele.Charge[iEleSelect[0]] + d.ele.Charge[iEleSelect[1]])==0;
				   else if (nMuSelect==2) return (d.mu.Charge[iMuSelect[0]] + d.mu.Charge[iMuSelect[1]])==0;
				   return false;
				 }});
  analysis_cuts['Z'].push_back({ .name="nHTop0BMassTag", .func = []    { return nHadTop0BMassTag>=1;                     }});
  analysis_cuts['Z'].push_back({ .name="DPhiHemill",     .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['Z'].push_back({ .name="Mll",            .func = []    { return std::abs(M_ll-91.2)<10; }});

  // L: W->lv enriched control region for Top Analysis
  analysis_cuts['L'].push_back({ .name="3Jet",       .func = []        { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['L'].push_back({ .name="HLT",        .func = [this,&d] { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['L'].push_back({ .name="0b",         .func = []        { return nLooseBTag==0;                    }});
  analysis_cuts['L'].push_back({ .name="1mTop",      .func = []        { return nHadTop0BMassTag>=1;              }});
  analysis_cuts['L'].push_back({ .name="1Lep",       .func = []        { return nLepSelect==1; }}); //Lepton's pt are defined here. Check the AnalysisBase.h
  analysis_cuts['L'].push_back({ .name="Razor",     .func = [&d]      { return d.evt.MR>=800 && R2_1l>=0.08;     }}); //Determine Razor cuts for this specific region. Only difference is R2 for this specific region, because only razor variable which use MET is R2. We add lepton pt to MET so R2 must be redefined.
  analysis_cuts['L'].push_back({ .name="DPhiHemill",     .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['L'].push_back({ .name="MT",         .func = []        { return MT>=30 && MT<100; }});   

  
 // G: Photon enriched sample
  analysis_cuts['G'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['G'].push_back({ .name="1Pho",       .func = []    { return nPhotonSelect==1;                 }});
  analysis_cuts['G'].push_back({ .name="3Jet",       .func = []    { return nJetNoPho>=3;                     }});
  analysis_cuts['G'].push_back({ .name="MR",         .func = [&d]  { return MR_pho>=800;                      }});
  analysis_cuts['G'].push_back({ .name="R2",         .func = [&d]  { return R2_pho>=0.08;                     }});
  analysis_cuts['G'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['G'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['G'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['G'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['G'].push_back({ .name="1MTop",      .func = []    { return nHadTopMassTag>=1;                }});
//analysis_cuts['G'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi_pho>=0.5;             }});
  analysis_cuts['G'].push_back({ .name="mDPhi", .func = [] { return dPhiRazorNoPho<2.8; }});
  
   // G-1: Photon enriched sample (G without boosted object)
  analysis_cuts['g'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }});
  analysis_cuts['g'].push_back({ .name="1Pho",       .func = []    { return nPhotonSelect==1;                 }});
  analysis_cuts['g'].push_back({ .name="3Jet",       .func = []    { return nJetNoPho>=3;                     }});
  analysis_cuts['g'].push_back({ .name="MR",         .func = [&d]  { return MR_pho>=800;                      }});
  analysis_cuts['g'].push_back({ .name="R2",         .func = [&d]  { return R2_pho>=0.08;                     }});
  analysis_cuts['g'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['g'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['g'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['g'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['g'].push_back({ .name="mDPhi", .func = [] { return dPhiRazorNoPho<2.8; }});



//---------------------- -------------for W analisis--------------------------------------------------------------------------------------------------------------

  // A: S Signal region
  analysis_cuts['A'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['A'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['A'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['A'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['A'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['A'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['A'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['A'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['A'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['A'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});

  // a': s' DPhi Control region of Signal region
  analysis_cuts['a'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['a'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['a'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['a'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['a'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['a'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['a'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['a'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['a'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['s'].push_back({ .name="InvmDPhi",   .func = []    { return minDeltaPhi<0.5;                  }});
  analysis_cuts['a'].push_back({ .name="InvmDPhi", .func = [] { return dPhiRazor>=2.8; }});

   // B: Q QCD enriched control sample
  analysis_cuts['B'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['B'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['B'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['B'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['B'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['B'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['B'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['B'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['B'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;                }});
  //analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                  }});
  analysis_cuts['B'].push_back({ .name="InvmDPhi0p3",.func = []    { return dPhiRazor>=2.8;                  }});

  // b': q' Dphi Control region of QCD enriched sample
  analysis_cuts['b'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['b'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['b'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['b'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['b'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['b'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['b'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['b'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['b'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;                }});
  //analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
//analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.3 && minDeltaPhi<0.5;         }});

   // C: T Top enriched control sample
  analysis_cuts['C'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['C'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['C'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['C'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['C'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['C'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['C'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['C'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  analysis_cuts['C'].push_back({ .name="MT",         .func = []    { return MT_vetolep<100;                   }});

  // D: W enriched control sample
  analysis_cuts['D'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['D'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['D'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['D'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['D'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['D'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['D'].push_back({ .name="1mW",        .func = []    { return nWMassTag>=1;                     }});
  //analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['D'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  analysis_cuts['D'].push_back({ .name="MT", .func = [] { return MT_vetolep>=30 && MT_vetolep<100; }});

  // E: Z->ll enriched control sample
  analysis_cuts['E'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['E'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['E'].push_back({ .name="MR_R2ll",    .func = [&d]  { return d.evt.MR>=800 && R2_ll>=0.08;     }});
  analysis_cuts['E'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['E'].push_back({ .name="2Lep",       .func = []    { return (nEleSelect==2&&nMuVeto==0)||(nMuSelect==2&&nEleVeto==0); }});
  analysis_cuts['E'].push_back({ .name="OppCharge",  .func = [&d]  { 
				   if (nEleSelect==2) return (d.ele.Charge[iEleSelect[0]] + d.ele.Charge[iEleSelect[1]])==0;
				   else if (nMuSelect==2) return (d.mu.Charge[iMuSelect[0]] + d.mu.Charge[iMuSelect[1]])==0;
				   return false;
				 }});
  analysis_cuts['E'].push_back({ .name="1mW",        .func = []    { return nWMassTag>=1;                     }});
  //analysis_cuts['Z'].push_back({ .name="mDPhill",    .func = []    { return minDeltaPhi_ll>=0.5;              }});
  analysis_cuts['E'].push_back({ .name="mDPhill",    .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['E'].push_back({ .name="Mll", .func = [] { return std::abs(M_ll-91.2)<10; }});
  

  //M: W->lv enriched control region for W Analysis
  analysis_cuts['M'].push_back({ .name="3Jet",       .func = []        { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['M'].push_back({ .name="HLT",        .func = [this,&d] { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['M'].push_back({ .name="0b",         .func = []        { return nLooseBTag==0;                    }});
  analysis_cuts['M'].push_back({ .name="1mW",        .func = []        { return nWMassTag>=1; }}); 
  analysis_cuts['M'].push_back({ .name="1Lep",       .func = []        { return nLepSelect==1; }}); //Lepton's pt are defined here. Check the AnalysisBase.h
  analysis_cuts['M'].push_back({ .name="Razor" ,      .func = [&d]      { return d.evt.MR>=800 && R2_1l>=0.08;  }}); //Determine Razor cuts for this specific region. Only difference is R2 for this specific region, because only razor variable which use MET is R2. We add lepton pt to MET so R2 must be redefined.
  analysis_cuts['M'].push_back({ .name="DPhiHemill",     .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['M'].push_back({ .name="MT",         .func = []        { return MT>=30 && MT<100; }});

  

  //m: + jets for W Analysis
  analysis_cuts['m'].push_back({ .name="4Jet",       .func = []        { return nJet>=4;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['m'].push_back({ .name="HLT",        .func = [this,&d] { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['m'].push_back({ .name="0b",         .func = []        { return nLooseBTag==0;                    }});
  analysis_cuts['m'].push_back({ .name="1mW",        .func = []        { return nWMassTag>=1; }}); 
  analysis_cuts['m'].push_back({ .name="1Lep",       .func = []        { return nLepSelect==1; }}); //Lepton's pt are defined here. Check the AnalysisBase.h
  analysis_cuts['m'].push_back({ .name="Razor" ,      .func = [&d]      { return d.evt.MR>=800 && R2_1l>=0.08;  }}); //Determine Razor cuts for this specific region. Only difference is R2 for this specific region, because only razor variable which use MET is R2. We add lepton pt to MET so R2 must be redefined.
  analysis_cuts['m'].push_back({ .name="DPhiHemill",     .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['m'].push_back({ .name="MT",         .func = []        { return MT>=30 && MT<100; }});

 
}//End of define_selection




void
Analysis::apply_scale_factors(DataStruct& data, const unsigned int& s, const std::vector<std::vector<double> >& nSigmaSFs)
{
  bool isFastSim = TString(sample).Contains("FastSim");
  size_t i = 0;

  // Don't forget to specify the total number of sigmas you use in settings_*.h !

  // Electron SFs (4 sigmas - reco, id, iso, fastsim)
  double sf_ele_veto, sf_ele_loose, sf_ele_medium;
  std::tie(sf_ele_veto, sf_ele_loose, sf_ele_medium) = calc_ele_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s], isFastSim);
  i+=4;

  // Muon SFs (3 sigmas - tracking, fullsim, fastsim)
  double sf_muon_veto, sf_muon_loose, sf_muon_medium;
  std::tie(sf_muon_veto, sf_muon_loose, sf_muon_medium) =  calc_muon_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], isFastSim);
  i+=3;

  // b tagging SFs (1 sigma)
  std::pair<double, double> sf_btag = calc_b_tagging_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], isFastSim);
  double sf_btag_loose = sf_btag.first, sf_btag_medium = sf_btag.second;
i+=2;

  // W tagging SF  (3 sigma - fullsim, fastsim, mistag)
 double sf_w = calc_w_tagging_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s], isFastSim);
  // fake W tagging SFs (no varriation of systematics) 
  double sf_fake_mW = calc_fake_w_mass_tagging_sf(data, nSigmaSFs[i][s]);
  double sf_fake_aW = calc_fake_w_anti_tagging_sf(data, nSigmaSFs[i][s]);
i+=5;
 
  // top tagging SF (3 sigma - fullsim, fastsim, mistag)
 double sf_top = calc_top_tagging_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s], isFastSim);
  // fake top tagging SFs (3 sigmas)
  double sf_fake_0b_mTop = calc_fake_top_0b_mass_tagging_sf(data, nSigmaSFs[i][s]);
  double sf_fake_MTop = calc_fake_top_mass_tagging_sf(data, nSigmaSFs[i][s]);
  double sf_fake_aTop = calc_fake_top_anti_tagging_sf(data, nSigmaSFs[i][s]);

  i+=6;

  // Select scale factors to use
  for (auto& sf : scale_factors) sf.second.clear();



  scale_factors['S'].push_back(sf_ele_veto);
  scale_factors['S'].push_back(sf_muon_veto);
  scale_factors['S'].push_back(sf_top);
  

  scale_factors['s'] = scale_factors['S'];

  scale_factors['Q'].push_back(sf_ele_veto);
  scale_factors['Q'].push_back(sf_muon_veto);
  scale_factors['Q'].push_back(sf_fake_aTop);  // We only have SF for the tag
  scale_factors['Q'].push_back(sf_btag_loose);

  scale_factors['q'] = scale_factors['Q'];

  scale_factors['T'].push_back(sf_ele_veto);
  scale_factors['T'].push_back(sf_muon_veto);
  scale_factors['T'].push_back(sf_top);
  

  scale_factors['W'].push_back(sf_ele_veto);
  scale_factors['W'].push_back(sf_muon_veto);
  scale_factors['W'].push_back(sf_fake_0b_mTop); // We only have SF for the tag
  scale_factors['W'].push_back(sf_btag_loose);

  scale_factors['Z'].push_back(sf_ele_medium);
  scale_factors['Z'].push_back(sf_muon_medium);
  scale_factors['Z'].push_back(sf_fake_MTop);


 /* scale_factors['L'].push_back(sf_ele_veto);
  scale_factors['L'].push_back(sf_muon_veto);
  scale_factors['L'].push_back(sf_fake_MTop);
  scale_factors['L'].push_back(sf_btag_loose); */


  scale_factors['L'] = scale_factors['W'];

  scale_factors['G'].push_back(sf_ele_veto);
  scale_factors['G'].push_back(sf_muon_veto);
  scale_factors['G'].push_back(sf_fake_MTop); 
 
   // G-1 region (G without boosted object)
  scale_factors['g'].push_back(sf_ele_veto);
  scale_factors['g'].push_back(sf_muon_veto);

// ----------------for W analiysis---------------------------------------

  scale_factors['A'].push_back(sf_ele_veto);
  scale_factors['A'].push_back(sf_muon_veto);
  scale_factors['A'].push_back(sf_btag_medium);
  scale_factors['A'].push_back(sf_w);

  scale_factors['a'] = scale_factors['A']; 

  scale_factors['B'].push_back(sf_ele_veto);
  scale_factors['B'].push_back(sf_muon_veto);
  scale_factors['B'].push_back(sf_btag_loose);
  scale_factors['B'].push_back(sf_fake_aW);

  scale_factors['b'] = scale_factors['B'];

  scale_factors['C'].push_back(sf_ele_veto);
  scale_factors['C'].push_back(sf_muon_veto);
  scale_factors['C'].push_back(sf_btag_medium);
  scale_factors['C'].push_back(sf_w);

  scale_factors['D'].push_back(sf_ele_veto);
  scale_factors['D'].push_back(sf_muon_veto);
  scale_factors['D'].push_back(sf_btag_loose);
  scale_factors['D'].push_back(sf_fake_mW);
 
  
  scale_factors['M'] = scale_factors['D'];

  scale_factors['E'].push_back(sf_ele_medium);
  scale_factors['E'].push_back(sf_muon_medium);
  scale_factors['E'].push_back(sf_fake_mW);

  
  scale_factors['m'].push_back(sf_ele_medium);
  scale_factors['m'].push_back(sf_muon_medium);
  scale_factors['m'].push_back(sf_btag_loose);
  scale_factors['m'].push_back(sf_fake_mW);


  ////////////////////////////////////////////////////////////////
 ///        NEW SCALE FACTORS for NEW C.R. REGIONS (L & M)     //
////////////////////////////////////////////////////////////////



}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& d) {
 return apply_all_cuts('S');
 return 0;
}


//_______________________________________________________
//                 List of Histograms

TH2D* h_S_MRR2;
TH2D* h_S_MRR2_nj3;
TH2D* h_S_MRR2_nj35;
TH2D* h_S_MRR2_nj5;

TH2D* h_s_MRR2;
TH2D* h_s_MRR2_nj3;
TH2D* h_s_MRR2_nj35;
TH2D* h_s_MRR2_nj5;

TH2D* h_Q_MRR2;
TH2D* h_Q_MRR2_nj3;
TH2D* h_Q_MRR2_nj35;
TH2D* h_Q_MRR2_nj5;

TH2D* h_q_MRR2;
TH2D* h_q_MRR2_nj3;
TH2D* h_q_MRR2_nj35;
TH2D* h_q_MRR2_nj5;

TH2D* h_W_MRR2;
TH2D* h_W_MRR2_nj3;
TH2D* h_W_MRR2_nj35;
TH2D* h_W_MRR2_nj5;

TH2D* h_T_MRR2;
TH2D* h_T_MRR2_nj3;
TH2D* h_T_MRR2_nj35;
TH2D* h_T_MRR2_nj5;

TH2D* h_Z_MRR2;
TH2D* h_Z_MRR2_nj3;
TH2D* h_Z_MRR2_nj35;
TH2D* h_Z_MRR2_nj5;

TH2D* h_L_MRR2;
TH2D* h_M_MRR2;
TH2D* h_m_MRR2;
TH2D* h_G_MRR2; 

TH1D* h_MR_Z;
TH1D* h_R2_Z;
TH1D* h_MR_T;
TH1D* h_R2_T;
TH1D* h_MR_W;
TH1D* h_R2_W;
TH1D* h_MR_Q;
TH1D* h_R2_Q;
TH1D* h_MR_q;
TH1D* h_R2_q;
TH1D* h_MR_S;
TH1D* h_R2_S;
TH1D* h_MR_s;
TH1D* h_R2_s;
TH1D* h_MR_L;
TH1D* h_R2_L;
TH1D* h_MR_M;
TH1D* h_R2_M;
TH1D* h_MR_m;
TH1D* h_R2_m;
TH1D* h_MR_G;
TH1D* h_R2_G; 

TH1D* h_njet_S;
TH1D* h_njet_Q;
TH1D* h_njet_T;  // number of jets 
TH1D* h_njet_W;
TH1D* h_njet_Z;
TH1D* h_njet_L;
TH1D* h_njet_M;
TH1D* h_njet_m;
TH1D* h_njet_G;

TH1D* h_S_HT;
TH1D* h_Q_HT;
TH1D* h_T_HT;
TH1D* h_W_HT;  // HT 
TH1D* h_s_HT;
TH1D* h_q_HT;
TH1D* h_Z_HT;

TH1D* h_S_MET;
TH1D* h_Q_MET;
TH1D* h_T_MET; // MET 
TH1D* h_W_MET;
TH1D* h_s_MET;
TH1D* h_q_MET;
TH1D* h_Z_MET;

TH1D* h_jet1_AK8_pt_Q;
TH1D* h_jet2_AK8_pt_Q;
TH1D* h_jet3_AK8_pt_Q;
TH1D* h_jet1_AK8_pt_q;
TH1D* h_jet2_AK8_pt_q;
TH1D* h_jet3_AK8_pt_q;
TH1D* h_jet1_AK8_pt_S;
TH1D* h_jet2_AK8_pt_S;
TH1D* h_jet3_AK8_pt_S;   // PT 
TH1D* h_jet1_AK8_pt_s;
TH1D* h_jet2_AK8_pt_s;
TH1D* h_jet3_AK8_pt_s;
TH1D* h_jet1_AK8_pt_T;
TH1D* h_jet2_AK8_pt_T;
TH1D* h_jet3_AK8_pt_T;
TH1D* h_jet1_AK8_pt_W;
TH1D* h_jet2_AK8_pt_W;
TH1D* h_jet3_AK8_pt_W;
TH1D* h_jet1_AK8_pt_Z;
TH1D* h_jet2_AK8_pt_Z;
TH1D* h_jet3_AK8_pt_Z;

TH1D* h_jet1_AK8_Eta_Q;
TH1D* h_jet2_AK8_Eta_Q;
TH1D* h_jet3_AK8_Eta_Q;
TH1D* h_jet1_AK8_Eta_q;
TH1D* h_jet2_AK8_Eta_q;
TH1D* h_jet3_AK8_Eta_q;
TH1D* h_jet1_AK8_Eta_S;
TH1D* h_jet2_AK8_Eta_S;
TH1D* h_jet3_AK8_Eta_S;   // ETA 
TH1D* h_jet1_AK8_Eta_s;
TH1D* h_jet2_AK8_Eta_s;
TH1D* h_jet3_AK8_Eta_s;
TH1D* h_jet1_AK8_Eta_T;
TH1D* h_jet2_AK8_Eta_T;
TH1D* h_jet3_AK8_Eta_T;
TH1D* h_jet1_AK8_Eta_W;
TH1D* h_jet2_AK8_Eta_W;
TH1D* h_jet3_AK8_Eta_W;
TH1D* h_jet1_AK8_Eta_Z;
TH1D* h_jet2_AK8_Eta_Z;
TH1D* h_jet3_AK8_Eta_Z;



TH1D* h_DPhi_Q;
TH1D* h_DPhi_q;
TH1D* h_DPhi_S;
TH1D* h_DPhi_s;
TH1D* h_DPhi_T;
TH1D* h_DPhi_W;
TH1D* h_DPhi_Z;

TH1D* h_MTR_Q;
TH1D* h_MTR_q;
TH1D* h_MTR_S;
TH1D* h_MTR_s;
TH1D* h_MTR_T;
TH1D* h_MTR_W;
TH1D* h_MTR_Z;


TH1D* h_MT_Q;
TH1D* h_MT_q;
TH1D* h_MT_S;
TH1D* h_MT_s;
TH1D* h_MT_T;
TH1D* h_MT_W;
TH1D* h_MT_Z;



TH1D* h_AK8Jet1Pt_Top_S;
TH1D* h_AK8Jet1Pt_mTop_W;
TH1D* h_AK8Jet1Pt_aTop_Q;
TH1D* h_AK8Jet1Pt_Top_T;
TH1D* h_AK8Jet1Pt_mTop_Z;

TH1D* h_AK8Jet1Eta_Top_S;
TH1D* h_AK8Jet1Eta_mTop_W;
TH1D* h_AK8Jet1Eta_aTop_Q;
TH1D* h_AK8Jet1Eta_Top_T;
TH1D* h_AK8Jet1Eta_mTop_Z;

std::vector<TH2D*> vh_R2_MR_S;
std::vector<TH2D*> vh_R2_MR_Q;
std::vector<TH2D*> vh_R2_MR_T;
std::vector<TH2D*> vh_R2_MR_W;
std::vector<TH2D*> vh_R2_MR_Z;
std::vector<TH2D*> vh_R2_MR_G;
std::vector<TH2D*> vh_R2_MR_L;
std::vector<TH2D*> vh_R2_MR_q;
std::vector<TH2D*> vh_R2_MR_s;

//FOR W analysis
std::vector<TH2D*> vh_R2_MR_A;
std::vector<TH2D*> vh_R2_MR_B;
std::vector<TH2D*> vh_R2_MR_C;
std::vector<TH2D*> vh_R2_MR_D;
std::vector<TH2D*> vh_R2_MR_E;
//std::vector<TH2D*> vh_R2_MR_G;
std::vector<TH2D*> vh_R2_MR_M;
std::vector<TH2D*> vh_R2_MR_b;
std::vector<TH2D*> vh_R2_MR_a;


   


//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
 // h_gen_toppt      = new TH1D("h_gen_toppt", "h_gen_toppt", 50, 0, 1000);
  
  //int nbn_HT = 19;
  //int nbn_MET = 10;
 
 
  //int nbn_AK8J1pt = 15;
  int nbn_AK8j1toppt = 16;
  
  int nbn_MTR = 8;

  int nbn_MR = 7;
  int nbn_R2 = 7;

  Double_t HT_bins[19] = {0, 200, 300, 400, 500, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500, 2000, 3000};
 
  //Double_t bn_HT_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,900.,1000.,2500.};
 // Double_t* bn_HT = 0;
  //bn_HT = utils::getVariableBinEdges(nbn_HT+1,bn_HT_tmp);

  /*Double_t bn_MET_tmp[] = {0., 100., 200., 300., 400., 500., 600., 800., 1000., 1500., 2000.};
  Double_t* bn_MET = 0;
  bn_MET = utils::getVariableBinEdges(nbn_MET+1,bn_MET_tmp); */

  Double_t bn_MTR_tmp[] = {200., 300., 400., 600., 800., 1000., 1200., 1600., 2000.};
  Double_t* bn_MTR = 0;
  bn_MTR = utils::getVariableBinEdges(nbn_MTR+1,bn_MTR_tmp); 

  Double_t bn_MR_tmp[] = {0.,600.,800.,1000.,1200.,1600.,2000.,4000.};
  Double_t* bn_MR = 0;
  bn_MR = utils::getVariableBinEdges(nbn_MR+1,bn_MR_tmp);
  Double_t bn_R2_tmp[] = {0.,0.04,0.08,0.12,0.16,0.24,0.5,1.};
  Double_t* bn_R2 = 0;
  bn_R2 = utils:: getVariableBinEdges(nbn_R2+1,bn_R2_tmp); 

  

 // Double_t bn_AK8J1pt_tmp[] = {0.,200.,220.,240,260,280,300.,320.,340,360,380,400.,450,500.,700,1000.};
 // Double_t* bn_AK8J1pt = 0;
 // bn_AK8J1pt = utils::getVariableBinEdges(nbn_AK8J1pt+1,bn_AK8J1pt_tmp);

  Double_t bn_AK8j1toppt_tmp[] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 1000, 1400, 2000, 3000, 4000, 5000};
  Double_t* bn_AK8j1toppt = 0;
  bn_AK8j1toppt = utils::getVariableBinEdges(nbn_AK8j1toppt+1,bn_AK8j1toppt_tmp);

    //S Region; S, S34, Sge5-->  MRvsR2...3 2D plots for Closure Test
    h_S_MRR2_nj3     = new TH2D("S_MRR2_nj3", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=3 && nJet<5
    h_S_MRR2_nj35    = new TH2D("S_MRR2_nj35", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=5
    h_S_MRR2_nj5     = new TH2D("S_MRR2_nj5", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_S_MRR2 = new TH2D("S_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    
    
    //Q Region; Q, Q34, Qge5-->  MRvsR2 ... 3 2D plots for Closure Test
    h_Q_MRR2_nj3     = new TH2D("Q_MRR2_nj3",";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=3 && nJet<5
    h_Q_MRR2_nj35    = new TH2D("Q_MRR2_nj35", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=5
    h_Q_MRR2_nj5     = new TH2D("Q_MRR2_nj5", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_Q_MRR2 = new TH2D("Q_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    
    //W Region; W, W34, Wge5-->  MRvsR2 ...3 2D plots for Closure Test
    h_W_MRR2_nj3     = new TH2D("W_MRR2_nj3",";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=3 && nJet<5
    h_W_MRR2_nj35    = new TH2D("W_MRR2_nj35", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=5
    h_W_MRR2_nj5     = new TH2D("W_MRR2_nj5", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_W_MRR2 = new TH2D("W_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    
    //T Region; T, T34, Tge5-->  MRvsR2 ...3 2D plots for Closure Test
    h_T_MRR2_nj3     = new TH2D("T_MRR2_nj3",";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=3 && nJet<5
    h_T_MRR2_nj35     = new TH2D("T_MRR2_nj35", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=5
    h_T_MRR2_nj5     = new TH2D("T_MRR2_nj5", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_T_MRR2 = new TH2D("T_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    
    //q Region; q, q35, qge5--> MRvsR2 ...3 2D plots for Closure Test
    h_q_MRR2_nj3     = new TH2D("q_MRR2_nj3",";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=3 && nJet<5
    h_q_MRR2_nj35    = new TH2D("q_MRR2_nj35",";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=5
    h_q_MRR2_nj5     = new TH2D("q_MRR2_nj5",";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_q_MRR2 = new TH2D("q_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    
    //s Region; s3, s35, sge5--> MRvsR2 ...3 2D plots for Closure Test
    h_s_MRR2_nj3     = new TH2D("s_MRR2_nj3", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );

    h_s_MRR2_nj35    = new TH2D("s_MRR2_nj35", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );

    h_s_MRR2_nj5     = new TH2D("s_MRR2_nj5", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_s_MRR2 = new TH2D("s_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);

   //Z Region; Z, S34, Sge5-->  MRvsR2...3 2D plots for Closure Test
    h_Z_MRR2_nj3     = new TH2D("Z_MRR2_nj3", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=3 && nJet<5
    h_Z_MRR2_nj35    = new TH2D("Z_MRR2_nj35", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    //nJet>=5
    h_Z_MRR2_nj5     = new TH2D("Z_MRR2_nj5", ";MR;R2", nbn_MR, bn_MR, nbn_R2, bn_R2 );
    h_Z_MRR2 = new TH2D("Z_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);

    h_L_MRR2 = new TH2D("L_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    h_M_MRR2 = new TH2D("M_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2);
    h_m_MRR2 = new TH2D("m_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2); 
    h_G_MRR2 = new TH2D("G_MRR2", ";MR;R2",nbn_MR,bn_MR,nbn_R2,bn_R2); 

   h_MR_Z = new TH1D("MR_Z", ";Z_MR", nbn_MR, bn_MR);
   h_R2_Z = new TH1D("R2_Z", ";Z_R2", nbn_R2,bn_R2);
   h_MR_T = new TH1D("MR_T", ";T_MR", nbn_MR, bn_MR);
   h_R2_T = new TH1D("R2_T", ";T_R2", nbn_R2,bn_R2);
   h_MR_W = new TH1D("MR_W", ";W_MR", nbn_MR, bn_MR);
   h_R2_W = new TH1D("R2_W", ";W_R2", nbn_R2,bn_R2);
   h_MR_Q = new TH1D("MR_Q", ";Q_MR", nbn_MR, bn_MR);
   h_R2_Q = new TH1D("R2_Q", ";Q_R2", nbn_R2,bn_R2);
   h_MR_q = new TH1D("MR_q", ";q_MR", nbn_MR, bn_MR);
   h_R2_q = new TH1D("R2_q", ";q_R2", nbn_R2,bn_R2);
   h_MR_S = new TH1D("MR_S", ";S_MR", nbn_MR, bn_MR);
   h_R2_S = new TH1D("R2_S", ";S_R2", nbn_R2,bn_R2);
   h_MR_s = new TH1D("MR_s", ";s_MR", nbn_MR, bn_MR);
   h_R2_s = new TH1D("R2_s", ";s_R2", nbn_R2,bn_R2);
   h_MR_L = new TH1D("MR_L", ";L_MR", nbn_MR, bn_MR);
   h_R2_L = new TH1D("R2_L", ";L_R2", nbn_R2,bn_R2);
   h_MR_M = new TH1D("MR_M", ";M_MR", nbn_MR, bn_MR);
   h_R2_M = new TH1D("R2_M", ";M_R2", nbn_R2,bn_R2);
   h_MR_m = new TH1D("MR_m", ";m_MR", nbn_MR, bn_MR);
   h_R2_m = new TH1D("R2_m", ";m_R2", nbn_R2,bn_R2);
   h_MR_G = new TH1D("MR_G", ";G_MR", nbn_MR, bn_MR);
   h_R2_G = new TH1D("R2_G", ";G_R2", nbn_R2,bn_R2); 
    
    

   h_njet_S = new TH1D("njet_S",         ";N_{jet}",                20, 0,  20);
   h_njet_Q = new TH1D("njet_Q",         ";N_{jet}",                20, 0,  20);
   h_njet_T = new TH1D("njet_T",         ";N_{jet}",                20, 0,  20);
   h_njet_W = new TH1D("njet_W",         ";N_{jet}",                20, 0,  20);
   h_njet_Z = new TH1D("njet_Z",         ";N_{jet}",                20, 0,  20);
   h_njet_L = new TH1D("njet_L",         ";N_{jet}",                20, 0,  20);
   h_njet_M = new TH1D("njet_M",         ";N_{jet}",                20, 0,  20);
   h_njet_m = new TH1D("njet_m",         ";N_{jet}",                20, 0,  20);
   h_njet_G = new TH1D("njet_G",         ";N_{jet}",                20, 0,  20);

   h_S_HT           = new TH1D("S_HT",";S_H_{T}", 15, HT_bins);
   h_Q_HT           = new TH1D("Q_HT",";Q_H_{T}", 15, HT_bins);
   h_T_HT           = new TH1D("T_HT",";T_H_{T}", 15, HT_bins);
   h_W_HT           = new TH1D("W_HT",";W_H_{T}", 15, HT_bins);
   h_s_HT           = new TH1D("s_HT",";s_H_{T}", 15, HT_bins);
   h_q_HT           = new TH1D("q_HT",";q_H_{T}", 15, HT_bins);
   h_Z_HT           = new TH1D("Z_HT",";Z_H_{T}", 15, HT_bins);

   h_S_MET = new TH1D("S_MET", ";S_#slash{E}_{T} (GeV)",  400,0,2000);
   h_Q_MET = new TH1D("Q_MET", ";Q_#slash{E}_{T} (GeV)",  400,0,2000);
   h_T_MET = new TH1D("T_MET", ";T_#slash{E}_{T} (GeV)",  400,0,2000);
   h_W_MET = new TH1D("W_MET", ";W_#slash{E}_{T} (GeV)",  400,0,2000);
   h_s_MET = new TH1D("s_MET", ";s_#slash{E}_{T} (GeV)",  400,0,2000);
   h_q_MET = new TH1D("q_MET", ";q_#slash{E}_{T} (GeV)",  400,0,2000);
   h_Z_MET = new TH1D("Z_MET", ";Z_#slash{E}_{T} (GeV)",  400,0,2000); 
  
   h_jet1_AK8_pt_Q = new TH1D("jet1_AK8_pt_Q",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_Q = new TH1D("jet2_AK8_pt_Q",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_Q = new TH1D("jet3_AK8_pt_Q",      ";AK8 jet p_{T, jet3}", 200, 0,2000);
   h_jet1_AK8_pt_q = new TH1D("jet1_AK8_pt_q",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_q = new TH1D("jet2_AK8_pt_q",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_q = new TH1D("jet3_AK8_pt_q",      ";AK8 jet p_{T, jet3}", 200, 0,2000);
   h_jet1_AK8_pt_S = new TH1D("jet1_AK8_pt_S",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_S = new TH1D("jet2_AK8_pt_S",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_S = new TH1D("jet3_AK8_pt_S",      ";AK8 jet p_{T, jet3}", 200, 0,2000);
   h_jet1_AK8_pt_s = new TH1D("jet1_AK8_pt_s",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_s = new TH1D("jet2_AK8_pt_s",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_s = new TH1D("jet3_AK8_pt_s",      ";AK8 jet p_{T, jet3}", 200, 0,2000);
   h_jet1_AK8_pt_T = new TH1D("jet1_AK8_pt_T",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_T = new TH1D("jet2_AK8_pt_T",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_T = new TH1D("jet3_AK8_pt_T",      ";AK8 jet p_{T, jet3}", 200, 0,2000);
   h_jet1_AK8_pt_W = new TH1D("jet1_AK8_pt_W",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_W = new TH1D("jet2_AK8_pt_W",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_W = new TH1D("jet3_AK8_pt_W",      ";AK8 jet p_{T, jet3}", 200, 0,2000);
   h_jet1_AK8_pt_Z = new TH1D("jet1_AK8_pt_Z",      ";AK8 jet p_{T, jet1}", 200, 0,2000);
   h_jet2_AK8_pt_Z = new TH1D("jet2_AK8_pt_Z",      ";AK8 jet p_{T, jet2}", 200, 0,2000);
   h_jet3_AK8_pt_Z = new TH1D("jet3_AK8_pt_Z",      ";AK8 jet p_{T, jet3}", 200, 0,2000);

   h_jet1_AK8_Eta_Q = new TH1D("jet1_AK8_Eta_Q",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_Q = new TH1D("jet2_AK8_Eta_Q",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_Q = new TH1D("jet3_AK8_Eta_Q",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet1_AK8_Eta_q = new TH1D("jet1_AK8_Eta_q",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_q = new TH1D("jet2_AK8_Eta_q",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_q = new TH1D("jet3_AK8_Eta_q",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet1_AK8_Eta_S = new TH1D("jet1_AK8_Eta_S",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_S = new TH1D("jet2_AK8_Eta_S",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_S = new TH1D("jet3_AK8_Eta_S",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet1_AK8_Eta_s = new TH1D("jet1_AK8_Eta_s",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_s = new TH1D("jet2_AK8_Eta_s",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_s = new TH1D("jet3_AK8_Eta_s",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet1_AK8_Eta_T = new TH1D("jet1_AK8_Eta_T",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_T = new TH1D("jet2_AK8_Eta_T",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_T = new TH1D("jet3_AK8_Eta_T",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet1_AK8_Eta_W = new TH1D("jet1_AK8_Eta_W",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_W = new TH1D("jet2_AK8_Eta_W",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_W = new TH1D("jet3_AK8_Eta_W",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet1_AK8_Eta_Z = new TH1D("jet1_AK8_Eta_Z",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet2_AK8_Eta_Z = new TH1D("jet2_AK8_Eta_Z",      ";AK8 jet #eta", 48,-2.4,2.4);
   h_jet3_AK8_Eta_Z = new TH1D("jet3_AK8_Eta_Z",      ";AK8 jet #eta", 48,-2.4,2.4);


   h_DPhi_Q = new TH1D("DPhi_Q",  ";#Delta#phi_{megajets}", 32,0,3.2);
   h_DPhi_q = new TH1D("DPhi_q",  ";#Delta#phi_{megajets}", 32,0,3.2);
   h_DPhi_S = new TH1D("DPhi_S",  ";#Delta#phi_{megajets}", 32,0,3.2);
   h_DPhi_s = new TH1D("DPhi_s",  ";#Delta#phi_{megajets}", 32,0,3.2);
   h_DPhi_T = new TH1D("DPhi_T",  ";#Delta#phi_{megajets}", 32,0,3.2);
   h_DPhi_W = new TH1D("DPhi_W",  ";#Delta#phi_{megajets}", 32,0,3.2);
   h_DPhi_Z = new TH1D("DPhi_Z",  ";#Delta#phi_{megajets}", 32,0,3.2);

   h_MTR_Q = new TH1D("MTR_Q",  ";Q_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
   h_MTR_q = new TH1D("MTR_q",  ";q_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
   h_MTR_S = new TH1D("MTR_S",  ";S_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
   h_MTR_s = new TH1D("MTR_s",  ";s_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
   h_MTR_T = new TH1D("MTR_T",  ";T_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
   h_MTR_W = new TH1D("MTR_W",  ";W_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
   h_MTR_Z = new TH1D("MTR_Z",  ";Z_M_{T}^{R} (GeV)", nbn_MTR, bn_MTR);
 
   h_MT_Q = new TH1D("MT_Q",  ";Q_m_{T} (GeV)", 100,0,2000);
   h_MT_q = new TH1D("MT_q",  ";q_m_{T} (GeV)", 100,0,2000);
   h_MT_S = new TH1D("MT_S",  ";S_m_{T} (GeV)", 100,0,2000);
   h_MT_s = new TH1D("MT_s",  ";s_m_{T} (GeV)", 100,0,2000);
   h_MT_T = new TH1D("MT_T",  ";T_m_{T} (GeV)", 100,0,2000);
   h_MT_W = new TH1D("MT_W",  ";W_m_{T} (GeV)", 100,0,2000);
   h_MT_Z = new TH1D("MT_Z",  ";Z_m_{T} (GeV)", 100,0,2000); 


   

   h_AK8Jet1Pt_Top_S = new TH1D("AK8Jet1Pt_Top_S",";S_AK8 Top jet p_{T} [GeV]", nbn_AK8j1toppt,bn_AK8j1toppt);
   h_AK8Jet1Pt_mTop_W = new TH1D("AK8Jet1Pt_mTop_W",";W_AK8 Top jet p_{T} [GeV]", nbn_AK8j1toppt,bn_AK8j1toppt);
   h_AK8Jet1Pt_aTop_Q = new TH1D("AK8Jet1Pt_aTop_Q",";Q_AK8 Top jet p_{T} [GeV]", nbn_AK8j1toppt,bn_AK8j1toppt);
   h_AK8Jet1Pt_Top_T = new TH1D("AK8Jet1Pt_Top_T",";T_AK8 Top jet p_{T} [GeV]", nbn_AK8j1toppt,bn_AK8j1toppt);
   h_AK8Jet1Pt_mTop_Z = new TH1D("AK8Jet1Pt_mTop_Z",";Z_AK8 Top jet p_{T} [GeV]", nbn_AK8j1toppt,bn_AK8j1toppt);


  h_AK8Jet1Eta_Top_S = new TH1D("AK8Jet1Eta_Top_S",";S_AK8 Top jet #eta", 48,-2.4,2.4);
  h_AK8Jet1Eta_mTop_W = new TH1D("AK8Jet1Eta_mTop_W",";W_AK8 Top jet #eta", 48,-2.4,2.4);
  h_AK8Jet1Eta_aTop_Q = new TH1D("AK8Jet1Eta_aTop_Q",";Q_AK8 Top jet #eta", 48,-2.4,2.4);
  h_AK8Jet1Eta_Top_T = new TH1D("AK8Jet1Eta_Top_T",";T_AK8 Top jet #eta", 48,-2.4,2.4);
  h_AK8Jet1Eta_mTop_Z = new TH1D("AK8Jet1Eta_mTop_Z",";Z_AK8 Top jet #eta", 48,-2.4,2.4); 
   
// for systematic uncertainty

  for (unsigned int i=0; i<=syst_nSyst; ++i) {
  std::stringstream histoname, title;
  title<<"Systematic variation #="<<i;
  
  histoname<<"R2_MR_S_syst"<<i;
  vh_R2_MR_S.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_Q_syst"<<i;
  vh_R2_MR_Q.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str(""); 

  histoname<<"R2_MR_T_syst"<<i;
  vh_R2_MR_T.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_W_syst"<<i;
  vh_R2_MR_W.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_Z_syst"<<i;
  vh_R2_MR_Z.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_G_syst"<<i;
  vh_R2_MR_G.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_s_syst"<<i;
  vh_R2_MR_s.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_q_syst"<<i;
  vh_R2_MR_q.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_L_syst"<<i;
  vh_R2_MR_L.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");


// for W analysis

  histoname<<"R2_MR_A_syst"<<i;
  vh_R2_MR_A.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_B_syst"<<i;
  vh_R2_MR_B.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str(""); 

  histoname<<"R2_MR_C_syst"<<i;
  vh_R2_MR_C.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_D_syst"<<i;
  vh_R2_MR_D.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_E_syst"<<i;
  vh_R2_MR_E.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_a_syst"<<i;
  vh_R2_MR_a.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_b_syst"<<i;
  vh_R2_MR_b.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

  histoname<<"R2_MR_M_syst"<<i;
  vh_R2_MR_M.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR;R2").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
  histoname.str("");

}
}



//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_analysis_histos(DataStruct& data, const unsigned int& syst_index, const double& weight)
{
    double w = weight;

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
  //Switch on settings.varySystematics to be effective

   

      w = sf_weight['S'];

      if (apply_all_cuts('S')) {
          
          if(nJet>=3){
              h_S_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>=3 && nJet<5){
              h_S_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>5){
              h_S_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
          }
      }

   if (apply_all_cuts('S')) h_njet_S->Fill(nJet,w); // number of jets for  S region 
   if (apply_all_cuts('S')) h_S_HT -> Fill(AK8_Ht, w); // HT 
   if (apply_all_cuts('S')) h_S_MET->Fill(data.met.Pt.at(0),w);
   if (apply_all_cuts('S')) {
      h_jet1_AK8_pt_S->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_S->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_S->Fill(data.jetsAK8.Pt[iJet[2]], w);
     }

   if (apply_all_cuts('S')) {
      h_jet1_AK8_Eta_S->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_S->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_S->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }
   
   if (apply_all_cuts_except('S',"DPhiHemi" )) h_DPhi_S->Fill(dPhiRazor, w);

   if (apply_all_cuts('S')) h_MTR_S->Fill(data.evt.MTR, w);
   if (apply_all_cuts('S')) h_S_MRR2->Fill(data.evt.MR, data.evt.R2, w);
   if (apply_all_cuts('S')) h_AK8Jet1Pt_Top_S->Fill(data.jetsAK8.Pt[iHadTopTag[0]], w);
   if (apply_all_cuts('S')) h_AK8Jet1Eta_Top_S->Fill(data.jetsAK8.Eta[iHadTopTag[0]], w);
   if (apply_all_cuts('S')) h_MR_S->Fill(data.evt.MR, w);
   if (apply_all_cuts('S')) h_R2_S->Fill(data.evt.R2, w);
   if (apply_all_cuts('S')) h_MT_S->Fill(MT, w);
//S' Region 2D plots       
      w = sf_weight['s'];
      if (apply_all_cuts('s')){
          
          if(nJet>=3){
              h_s_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>=3 && nJet<5){
              h_s_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          
          if(nJet>5){
              h_s_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
      }
  if (apply_all_cuts('s')) h_s_HT -> Fill(AK8_Ht, w); // HT 
  if (apply_all_cuts('s')) h_s_MET->Fill(data.met.Pt.at(0),w);
  if (apply_all_cuts('s')) {
      h_jet1_AK8_pt_s->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_s->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_s->Fill(data.jetsAK8.Pt[iJet[2]], w);
      }
  
   if (apply_all_cuts('s')) {
      h_jet1_AK8_Eta_s->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_s->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_s->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }
   
   if (apply_all_cuts_except('s',"DPhiHemi" )) h_DPhi_s->Fill(dPhiRazor, w);
   if (apply_all_cuts('s')) h_MTR_s->Fill(data.evt.MTR, w); 
   if (apply_all_cuts('s')) h_s_MRR2->Fill(data.evt.MR, data.evt.R2, w); 
   if (apply_all_cuts('s')) h_MR_s->Fill(data.evt.MR, w);
   if (apply_all_cuts('s')) h_R2_s->Fill(data.evt.R2, w);
   if (apply_all_cuts('s')) h_MT_s->Fill(MT, w);
//W region 2D plots      
      w = sf_weight['W'];
      if (apply_all_cuts('W')){
          
          if(nJet>=3){
              h_W_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>=3 && nJet<5){
              h_W_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          
          if(nJet>5){
              h_W_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
      }
 
  if (apply_all_cuts('W')) h_njet_W->Fill(nJet,w);  // number of jets for W region 
  if (apply_all_cuts('W')) h_W_HT -> Fill(AK8_Ht, w); // HT 
  if (apply_all_cuts('W')) h_W_MET->Fill(data.met.Pt.at(0),w);
  if (apply_all_cuts('W')) {
      h_jet1_AK8_pt_W->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_W->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_W->Fill(data.jetsAK8.Pt[iJet[2]], w);
     }
   if (apply_all_cuts('W')) {
      h_jet1_AK8_Eta_W->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_W->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_W->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }
  
  if (apply_all_cuts_except('W',"DPhiHemi" )) h_DPhi_W->Fill(dPhiRazor, w);
  if (apply_all_cuts('W')) h_MTR_W->Fill(data.evt.MTR, w);
  if (apply_all_cuts('W')) h_W_MRR2->Fill(data.evt.MR, data.evt.R2, w);
  if (apply_all_cuts('W')) h_AK8Jet1Pt_mTop_W->Fill(data.jetsAK8.Pt[iHadTop0BMassTag[0]], w);
  if (apply_all_cuts('W')) h_AK8Jet1Eta_mTop_W->Fill(data.jetsAK8.Eta[iHadTop0BMassTag[0]], w);
  if (apply_all_cuts('W')) h_MR_W->Fill(data.evt.MR, w);
  if (apply_all_cuts('W')) h_R2_W->Fill(data.evt.R2, w);
  if (apply_all_cuts_except('W',"MT" )) h_MT_W->Fill(MT, w);  



//Q Region 2D plots      
      w = sf_weight['Q'];
      if (apply_all_cuts('Q')){
          
          if(nJet>=3){
              h_Q_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>=3 && nJet<5){
              h_Q_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          
          if(nJet>5){
              h_Q_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
      }

   if (apply_all_cuts('Q')) h_njet_Q->Fill(nJet,w); // number of jets for Q region 
   if (apply_all_cuts('Q')) h_Q_HT -> Fill(AK8_Ht, w); // HT 
   if (apply_all_cuts('Q')) h_Q_MET->Fill(data.met.Pt.at(0),w);
   if (apply_all_cuts('Q')){
      h_jet1_AK8_pt_Q->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_Q->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_Q->Fill(data.jetsAK8.Pt[iJet[2]], w);                      
   }
 
   if (apply_all_cuts('Q')) {
      h_jet1_AK8_Eta_Q->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_Q->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_Q->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }

  if (apply_all_cuts_except('Q',"DPhiHemi" )) h_DPhi_Q->Fill(dPhiRazor, w);
  if (apply_all_cuts('Q')) h_MTR_Q->Fill(data.evt.MTR, w);
  if (apply_all_cuts('Q')) h_Q_MRR2->Fill(data.evt.MR, data.evt.R2, w);
  if (apply_all_cuts('Q')) h_AK8Jet1Pt_aTop_Q->Fill(data.jetsAK8.Pt[iHadTop0BAntiTag[0]], w);
  if (apply_all_cuts('Q')) h_AK8Jet1Eta_aTop_Q->Fill(data.jetsAK8.Eta[iHadTop0BAntiTag[0]], w); 
  if (apply_all_cuts('Q')) h_MR_Q->Fill(data.evt.MR, w);
  if (apply_all_cuts('Q')) h_R2_Q->Fill(data.evt.R2, w);
  if (apply_all_cuts('Q')) h_MT_Q->Fill(MT, w);

//Q' Region 2D plots      
      w = sf_weight['q'];
      if (apply_all_cuts('q')){
          
          if(nJet>=3){
              h_q_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>=3 && nJet<5){
              h_q_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          
          if(nJet>5){
              h_q_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
      }

   if (apply_all_cuts('q')) h_q_HT -> Fill(AK8_Ht, w); // HT  
   if (apply_all_cuts('q')) h_q_MET->Fill(data.met.Pt.at(0),w);
   if (apply_all_cuts('q')) {
      h_jet1_AK8_pt_q->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_q->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_q->Fill(data.jetsAK8.Pt[iJet[2]], w);
      }
   if (apply_all_cuts('q')) {
      h_jet1_AK8_Eta_q->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_q->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_q->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }

   if (apply_all_cuts_except('q',"DPhiHemi" )) h_DPhi_q->Fill(dPhiRazor, w);
   if (apply_all_cuts('q')) h_MTR_q->Fill(data.evt.MTR, w); 
   if (apply_all_cuts('q')) h_q_MRR2->Fill(data.evt.MR, data.evt.R2, w);
   if (apply_all_cuts('q')) h_MR_q->Fill(data.evt.MR, w);
   if (apply_all_cuts('q')) h_R2_q->Fill(data.evt.R2, w);
   if (apply_all_cuts('q')) h_MT_q->Fill(MT, w); 
//T region 2D plots      
      w = sf_weight['T'];
      if (apply_all_cuts('T')){
          
          if(nJet>=3){
              h_T_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          if(nJet>=3 && nJet<5){
              h_T_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          
          if(nJet>5){
              h_T_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
      }
  if (apply_all_cuts('T')) h_njet_T->Fill(nJet,w);  // number of jets for W region data.met.Phi[0] 
  if (apply_all_cuts('T')) h_T_HT -> Fill(AK8_Ht, w); // HT
  if (apply_all_cuts('T')) h_T_MET->Fill(data.met.Pt.at(0),w); 
  if (apply_all_cuts('T')) {
      h_jet1_AK8_pt_T->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_T->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_T->Fill(data.jetsAK8.Pt[iJet[2]], w);
      }
   if (apply_all_cuts('T')) {
      h_jet1_AK8_Eta_T->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_T->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_T->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }
 
  if (apply_all_cuts_except('T',"DPhiHemi" )) h_DPhi_T->Fill(dPhiRazor, w);
  //if (apply_all_cuts('T')) h_DPhi_T->Fill(dPhiRazor, w);
  if (apply_all_cuts('T')) h_MTR_T->Fill(data.evt.MTR, w);
  if (apply_all_cuts('T')) h_T_MRR2->Fill(data.evt.MR, data.evt.R2, w);
  if (apply_all_cuts('T')) h_AK8Jet1Pt_Top_T->Fill(data.jetsAK8.Pt[iHadTopTag[0]], w);
  if (apply_all_cuts('T')) h_AK8Jet1Eta_Top_T->Fill(data.jetsAK8.Eta[iHadTopTag[0]], w);
  if (apply_all_cuts('T')) h_MR_T->Fill(data.evt.MR, w);
  if (apply_all_cuts('T')) h_R2_T->Fill(data.evt.R2, w);
  if (apply_all_cuts_except('T',"MT" )) h_MT_T->Fill(MT, w); 


// Z enriched region

   w = sf_weight['Z'];

   if (apply_all_cuts('Z')){
          
          if(nJet>=3){
              h_Z_MRR2_nj3 -> Fill(data.evt.MR, data.evt.R2, w);
          }
          
          if(nJet>=3 && nJet<5){
              h_Z_MRR2_nj35 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
          
          if(nJet>5){
              h_Z_MRR2_nj5 -> Fill(data.evt.MR, data.evt.R2, w);
              
          }
      }
  

   if (apply_all_cuts('Z')) h_njet_Z->Fill(nJet,w);
   if (apply_all_cuts('Z')) h_Z_HT -> Fill(AK8_Ht, w); // HT
   if (apply_all_cuts('Z')) h_Z_MET->Fill(data.met.Pt.at(0),w);
   if (apply_all_cuts('Z')) {
      h_jet1_AK8_pt_Z->Fill(data.jetsAK8.Pt[iJet[0]], w);
      h_jet2_AK8_pt_Z->Fill(data.jetsAK8.Pt[iJet[1]], w);
      h_jet3_AK8_pt_Z->Fill(data.jetsAK8.Pt[iJet[2]], w);
      } 
   if (apply_all_cuts('Z')) {
      h_jet1_AK8_Eta_Z->Fill(data.jetsAK8.Eta[iJet[0]], w);
      h_jet2_AK8_Eta_Z->Fill(data.jetsAK8.Eta[iJet[1]], w);
      h_jet3_AK8_Eta_Z->Fill(data.jetsAK8.Eta[iJet[2]], w);
     }
   if (apply_all_cuts_except('Z',"DPhiHemill" )) h_DPhi_Z->Fill(dPhiRazor, w);
   if (apply_all_cuts('Z')) h_MR_Z->Fill(data.evt.MR, w);
   if (apply_all_cuts('Z')) h_R2_Z->Fill(R2_ll, w);
   if (apply_all_cuts('Z')) h_Z_MRR2-> Fill(data.evt.MR, R2_ll, w);
   if (apply_all_cuts('Z'))h_MTR_Z->Fill(data.evt.MTR, w);
   if (apply_all_cuts('Z')) h_AK8Jet1Pt_mTop_Z->Fill(data.jetsAK8.Pt[iHadTop0BMassTag[0]], w);
   if (apply_all_cuts('Z')) h_AK8Jet1Eta_mTop_Z->Fill(data.jetsAK8.Eta[iHadTop0BMassTag[0]], w);
   if (apply_all_cuts_except('Z',"Mll" )) h_MT_Z->Fill(MT, w); 



// L enriched region

  
 w = sf_weight['L'];

   if (apply_all_cuts('L')) h_MR_L->Fill(data.evt.MR, w);
   if (apply_all_cuts('L')) h_R2_L->Fill(R2_1l, w);
   if (apply_all_cuts('L')) h_L_MRR2-> Fill(data.evt.MR, R2_1l, w);
   if (apply_all_cuts('L')) h_njet_L->Fill(nJet,w);

// M enriched region

   w = sf_weight['M'];

   if (apply_all_cuts('M')) h_MR_M->Fill(data.evt.MR, w);
   if (apply_all_cuts('M')) h_R2_M->Fill(R2_1l, w);
   if (apply_all_cuts('M')) h_M_MRR2-> Fill(data.evt.MR, R2_1l, w);
   if (apply_all_cuts('M')) h_njet_M->Fill(nJet, w);


// m region 4 jets   for W analysis

   w = sf_weight['m'];  

  if (apply_all_cuts('m')) h_m_MRR2-> Fill(data.evt.MR, R2_1l, w);
  if (apply_all_cuts('m')) h_MR_m->Fill(data.evt.MR, w);
  if (apply_all_cuts('m')) h_R2_m->Fill(R2_1l, w);
  if (apply_all_cuts('m')) h_njet_m->Fill(nJet, w);

   w = sf_weight['G'];  

  if (apply_all_cuts('G')) h_G_MRR2-> Fill(MR_pho, R2_pho, w);
  if (apply_all_cuts('G')) h_MR_G->Fill(MR_pho, w);
  if (apply_all_cuts('G')) h_R2_G->Fill(R2_pho, w);
  if (apply_all_cuts('G')) h_njet_G->Fill(nJet, w);
 
}

// systematic uncertainty fill 

 w = sf_weight['S'];
 if (apply_all_cuts('S')) vh_R2_MR_S[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['Q'];
 if (apply_all_cuts('Q')) vh_R2_MR_Q[syst_index]->Fill(data.evt.MR, data.evt.R2, w);


 w = sf_weight['T'];
 if (apply_all_cuts('T')) vh_R2_MR_T[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['W'];
 if (apply_all_cuts('W')) vh_R2_MR_W[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['Z'];
 if (apply_all_cuts('Z')) vh_R2_MR_Z[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['G'];
 if (apply_all_cuts('G')) vh_R2_MR_G[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['s'];
 if (apply_all_cuts('s')) vh_R2_MR_s[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['q'];
 if (apply_all_cuts('q')) vh_R2_MR_q[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

  w = sf_weight['L'];
 if (apply_all_cuts('L')) vh_R2_MR_L[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
 

// for W nalysisis systematic uncertainty fill


  w = sf_weight['A'];
 if (apply_all_cuts('A')) vh_R2_MR_A[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['B'];
 if (apply_all_cuts('B')) vh_R2_MR_B[syst_index]->Fill(data.evt.MR, data.evt.R2, w);


 w = sf_weight['C'];
 if (apply_all_cuts('C')) vh_R2_MR_C[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['D'];
 if (apply_all_cuts('D')) vh_R2_MR_D[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['E'];
 if (apply_all_cuts('E')) vh_R2_MR_E[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['a'];
 if (apply_all_cuts('a')) vh_R2_MR_a[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

 w = sf_weight['b'];
 if (apply_all_cuts('b')) vh_R2_MR_b[syst_index]->Fill(data.evt.MR, data.evt.R2, w);

  w = sf_weight['M'];
 if (apply_all_cuts('M')) vh_R2_MR_M[syst_index]->Fill(data.evt.MR, data.evt.R2, w);




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
