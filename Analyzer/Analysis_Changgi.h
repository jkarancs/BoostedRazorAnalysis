#ifndef VER
#define VER 2
#endif

#include "TLorentzVector.h"
#include "TMath.h"
#include "common/AnalysisBase.h"

//_______________________________________________________
//                  Calculate variables

void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{
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
  return 1;
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

  //baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  //baseline_cuts.push_back({ .name="Skim_R2",         .func = [&d]  { return d.evt.R2>=0.04;                  }}); // New skim cut introduced in 2017 february

  // S: Signal region
  analysis_cuts['S'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['S'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['S'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['S'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['S'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['S'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['S'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['S'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});

  // S': DPhi Control region of Signal region
  analysis_cuts['s'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['s'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['s'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['s'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['s'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['s'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['s'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['s'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['s'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['s'].push_back({ .name="InvmDPhi",   .func = []    { return minDeltaPhi<0.5;                  }});
  analysis_cuts['s'].push_back({ .name="InvmDPhi",   .func = []    { return dPhiRazor>=2.8;                  }});

  // W boson tag fake rate
  analysis_cuts['Y'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['Y'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  //analysis_cuts['Y'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['Y'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Y'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['Y'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['Y'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['Y'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['Y'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['Y'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                  }});
  analysis_cuts['Y'].push_back({ .name="InvmDPhi0p3",.func = []    { return dPhiRazor>=2.8;                  }});

  // Top quark anti tag fake rate
  analysis_cuts['K'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['K'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  //analysis_cutsK'I'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['K'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['K'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['K'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['K'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['K'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['K'].push_back({ .name="1a0bTop",    .func = []    { return nHadTop0BAntiTag>=1;              }});
  analysis_cuts['K'].push_back({ .name="InvmDPhi0p3",.func = []    { return dPhiRazor>=2.8;                  }});

  // k: Photon enriched sample
  analysis_cuts['k'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['k'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  //analysis_cuts['k'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  //analysis_cuts['k'].push_back({ .name="R2",         .func = [&d]  { return R2_pho>=0.08;                     }});
  analysis_cuts['k'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1         || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['k'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['k'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['k'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['k'].push_back({ .name="1Pho",       .func = []    { return nPhotonSelect==1;                 }});
  analysis_cuts['k'].push_back({ .name="1a0bTop",    .func = []    { return nHadTop0BAntiTag>=1;              }});
  analysis_cuts['k'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                    }});

  // Gen Level Information
  bool isFFsim = TString(sample).Contains("TTJets_madgraphMLM");
  if(isFFsim){ 
    analysis_cuts['g'].push_back({ .name="genW",     .func = []    { return npreWTag > 0;                     }});
    analysis_cuts['g'].push_back({ .name="mWTag",    .func = []    { return nmWTag > 0;                       }});
    analysis_cuts['g'].push_back({ .name="WTag",     .func = []    { return nWTag > 0;                        }});
    analysis_cuts['B'].push_back({ .name="genTop",   .func = []    { return npreTopTag > 0;                   }});
    analysis_cuts['B'].push_back({ .name="mTopTag",  .func = []    { return nmTopTag > 0;                     }});
    analysis_cuts['B'].push_back({ .name="TopTag",   .func = []    { return nTopTag > 0;                      }});
  }
  
  // Q: QCD enriched control sample
  analysis_cuts['Q'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['Q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['Q'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['Q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['Q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['Q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['Q'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;                }});
  //analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return minDeltaPhi<0.3;                  }});
  analysis_cuts['Q'].push_back({ .name="InvmDPhi0p3",.func = []    { return dPhiRazor>=2.8;                  }});

  // Q': Dphi Control region of QCD enriched sample
  analysis_cuts['q'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['q'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['q'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['q'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['q'].push_back({ .name="1aW",        .func = []    { return nTightWAntiTag>=1;                }});
  //analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  //analysis_cuts['q'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.3 && minDeltaPhi<0.5;         }});

  // T: Top enriched control sample
  analysis_cuts['T'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['T'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['T'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['T'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['T'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                   }});
  analysis_cuts['T'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                    }});
  //analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['T'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  analysis_cuts['T'].push_back({ .name="MT",         .func = []    { return MT_vetolep<100;                   }});

  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['W'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['W'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['W'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['W'].push_back({ .name="1mW",        .func = []    { return nWMassTag>=1;                     }});
  //analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['W'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  analysis_cuts['W'].push_back({ .name="MT",         .func = []    { return MT_vetolep>=30 && MT_vetolep<100; }});

  // Z: Z->ll enriched control sample
  analysis_cuts['Z'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['Z'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['Z'].push_back({ .name="MR_R2ll",    .func = [&d]  { return d.evt.MR>=800 && R2_ll>=0.08;     }});
  analysis_cuts['Z'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['Z'].push_back({ .name="2Lep",       .func = []    { return (nEleSelect==2&&nMuVeto==0)||(nMuSelect==2&&nEleVeto==0); }});
  analysis_cuts['Z'].push_back({ .name="OppCharge",  .func = [&d]  { 
				   if (nEleSelect==2) return (d.ele.Charge[iEleSelect[0]] + d.ele.Charge[iEleSelect[1]])==0;
				   else if (nMuSelect==2) return (d.mu.Charge[iMuSelect[0]] + d.mu.Charge[iMuSelect[1]])==0;
				   return false;
				 }});
  analysis_cuts['Z'].push_back({ .name="1mW",        .func = []    { return nWMassTag>=1;                     }});
  //analysis_cuts['Z'].push_back({ .name="mDPhill",    .func = []    { return minDeltaPhi_ll>=0.5;              }});
  analysis_cuts['Z'].push_back({ .name="mDPhill",    .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['Z'].push_back({ .name="Mll",        .func = []    { return std::abs(M_ll-91.2)<10;           }});

  // G: Photon enriched sample
  analysis_cuts['G'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['G'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                       }}); // Similar         to pt>200, one AK8 jet has pt>200
  analysis_cuts['G'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }});
  analysis_cuts['G'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  analysis_cuts['G'].push_back({ .name="R2",         .func = [&d]  { return R2_pho>=0.08;                     }});
  analysis_cuts['G'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1         || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['G'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['G'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['G'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['G'].push_back({ .name="1Pho",       .func = []    { return nPhotonSelect==1;                 }});
  //analysis_cuts['G'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['G'].push_back({ .name="1mW",        .func = []    { return nWMassTag>=1;                     }});
  //analysis_cuts['G'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi_pho>=0.5;             }});
  analysis_cuts['G'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                    }});

  // z: Photon enriched sample
  analysis_cuts['z'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['z'].push_back({ .name="1JetAK8",    .func = []    { return nJetAK8>=1;                    }}); // Similar         to pt>200, one AK8 jet has pt>200
  analysis_cuts['z'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                       }});
  //analysis_cuts['z'].push_back({ .name="MR",         .func = [&d]  { return d.evt.MR>=800;                    }});
  //analysis_cuts['z'].push_back({ .name="R2",         .func = [&d]  { return R2_pho>=0.08;                     }});
  //analysis_cuts['z'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.Photon165_HE10==1 : 1; }});
  analysis_cuts['z'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1         || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['z'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                   }});
  analysis_cuts['z'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                    }});
  analysis_cuts['z'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;              }});
  //analysis_cuts['z'].push_back({ .name="1Pho",       .func = []    { return nPhotonSelect==1;              }});
  //analysis_cuts['z'].push_back({ .name="1FakePho",       .func = []    { return nPhotonFake==1;         }});
  //analysis_cuts['z'].push_back({ .name="1mW",        .func = []    { return nWMassTag>=1;                     }});
  //analysis_cuts['G'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  //analysis_cuts['z'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                    }});

  // t: Boosted Top Signal region
  analysis_cuts['A'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['A'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['A'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['A'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['A'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['A'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['A'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['A'].push_back({ .name="1Top",       .func = []    { return nHadTopTag>=1;                    }});
  //analysis_cuts['A'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['A'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});

  // a: Boosted Top' Signal region
  analysis_cuts['a'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['a'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['a'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['a'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['a'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['a'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['a'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['a'].push_back({ .name="1Top",       .func = []    { return nHadTopTag>=1;                    }});
  //analysis_cuts['a'].push_back({ .name="InvmDPhi",      .func = []    { return minDeltaPhi<0.5;                 }});
  analysis_cuts['a'].push_back({ .name="InvmDPhi",      .func = []    { return dPhiRazor>=2.8;                 }});

  // t: Boosted Top Q region
  analysis_cuts['X'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['X'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['X'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['X'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['X'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['X'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['X'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['X'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['X'].push_back({ .name="1aTop",       .func = []   { return nHadTop0BAntiTag>=1;                    }});
  //analysis_cuts['X'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi<0.3;                 }});
  analysis_cuts['X'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor>=2.8;                 }});

  // t: Boosted Top Q' region
  analysis_cuts['x'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['x'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['x'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['x'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['x'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                      }});
  analysis_cuts['x'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                       }});
  analysis_cuts['x'].push_back({ .name="0IsoTrk",    .func = [&d]  { return d.evt.NIsoTrk==0;                 }});
  analysis_cuts['x'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['x'].push_back({ .name="1aTop",       .func = []   { return nHadTop0BAntiTag>=1;                    }});
  //analysis_cuts['x'].push_back({ .name="InvmDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['x'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});

  // t: Boosted Top T region
  analysis_cuts['V'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['V'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['V'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['V'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['V'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['V'].push_back({ .name="1Top",       .func = []   { return nHadTopTag>=1;                    }});
  //analysis_cuts['V'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['V'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  analysis_cuts['V'].push_back({ .name="MT",         .func = []    { return     MT_vetolep<100;               }});

  // t: Boosted Top W region
  analysis_cuts['C'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['C'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['C'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;  }});
  analysis_cuts['C'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['C'].push_back({ .name="1Lep",       .func = []    { return nLepVeto==1;                      }});
  analysis_cuts['C'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                    }});
  analysis_cuts['C'].push_back({ .name="1mTop",       .func = []   { return nHadTop0BMassTag>=1;                    }});
  //analysis_cuts['C'].push_back({ .name="mDPhi",      .func = []    { return minDeltaPhi>=0.5;                 }});
  analysis_cuts['C'].push_back({ .name="mDPhi",      .func = []    { return dPhiRazor<2.8;                 }});
  analysis_cuts['C'].push_back({ .name="MT",         .func = []    { return MT_vetolep>=30 && MT_vetolep<100; }});

  // Z: Z->ll enriched control sample
  analysis_cuts['F'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['F'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['F'].push_back({ .name="MR_R2ll",    .func = [&d]  { return d.evt.MR>=800 && R2_ll>=0.08;     }});
  analysis_cuts['F'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});
  analysis_cuts['F'].push_back({ .name="2Lep",       .func = []    { return (nEleSelect==2&&nMuVeto==0)||(nMuSelect==2&&nEleVeto==0); }});
  analysis_cuts['F'].push_back({ .name="OppCharge",  .func = [&d]  { 
				   if (nEleSelect==2) return (d.ele.Charge[iEleSelect[0]] + d.ele.Charge[iEleSelect[1]])==0;
				   else if (nMuSelect==2) return (d.mu.Charge[iMuSelect[0]] + d.mu.Charge[iMuSelect[1]])==0;
				   return false;
				 }});
  analysis_cuts['F'].push_back({ .name="1mTop",        .func = []    { return nHadTopMassTag>=1;                     }});
  analysis_cuts['F'].push_back({ .name="mDPhill",    .func = []    { return dPhiRazor<2.8;              }});
  analysis_cuts['F'].push_back({ .name="Mll",        .func = []    { return std::abs(M_ll-91.2)<10;           }});

  // preselection
  analysis_cuts['P'].push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>200
  analysis_cuts['P'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                          }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['P'].push_back({ .name="MR_R2",    .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08;     }});
  analysis_cuts['P'].push_back({ .name="HLT",   .func = [this,&d]  { return isData ? d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1 || d.hlt.PFHT900==1 : 1; }});

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
  double sf_w = calc_w_tagging_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], isFastSim);
  i+=2;

  // fake W tagging SFs (no varriation of systematics) 
  double sf_fake_mW = calc_fake_w_mass_tagging_sf(data);
  double sf_fake_aW = calc_fake_w_anti_tagging_sf(data);

  double sf_fake_mTop = calc_fake_top_mass_tagging_sf(data);
  double sf_fake_aTop = calc_fake_top_anti_tagging_sf(data);

  // b tagging SFs (1 sigma)
  std::pair<double, double> sf_btag = calc_b_tagging_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], isFastSim);
  double sf_btag_loose = sf_btag.first, sf_btag_medium = sf_btag.second;
  i+=2;

  // top tagging SF (1 sigma)
  double sf_top = calc_top_tagging_sf(data, nSigmaSFs[i][s], nSigmaSFs[i+1][s], isFastSim);
  i+=2;

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
  scale_factors['Q'].push_back(sf_fake_aW);

  scale_factors['q'] = scale_factors['Q'];

  scale_factors['T'].push_back(sf_ele_medium);
  scale_factors['T'].push_back(sf_muon_medium);
  scale_factors['T'].push_back(sf_btag_medium);
  scale_factors['T'].push_back(sf_w);

  scale_factors['W'].push_back(sf_ele_medium);
  scale_factors['W'].push_back(sf_muon_medium);
  scale_factors['W'].push_back(sf_btag_loose);
  scale_factors['W'].push_back(sf_fake_mW);

  scale_factors['Z'].push_back(sf_ele_medium);
  scale_factors['Z'].push_back(sf_muon_medium);
  scale_factors['Z'].push_back(sf_fake_mW);

  scale_factors['G'].push_back(sf_ele_medium);
  scale_factors['G'].push_back(sf_muon_medium);
  scale_factors['G'].push_back(sf_fake_mW);

  scale_factors['z'].push_back(sf_ele_medium);
  scale_factors['z'].push_back(sf_muon_medium);
  //scale_factors['z'].push_back(sf_fake_mW);

  scale_factors['Y'].push_back(sf_ele_veto);
  scale_factors['Y'].push_back(sf_muon_veto);
  scale_factors['Y'].push_back(sf_btag_loose);

  scale_factors['K'].push_back(sf_ele_veto);
  scale_factors['K'].push_back(sf_muon_veto);
  scale_factors['K'].push_back(sf_btag_loose);
  scale_factors['k'].push_back(sf_ele_veto);
  scale_factors['k'].push_back(sf_muon_veto);

  // Top analysis
  scale_factors['A'].push_back(sf_ele_veto);
  scale_factors['A'].push_back(sf_muon_veto);
  scale_factors['A'].push_back(sf_top);
  // Top analysis
  scale_factors['a'].push_back(sf_ele_veto);
  scale_factors['a'].push_back(sf_muon_veto);
  scale_factors['a'].push_back(sf_top);
  // Top analysis
  scale_factors['X'].push_back(sf_ele_veto);
  scale_factors['X'].push_back(sf_muon_veto);
  scale_factors['X'].push_back(sf_btag_loose);
  scale_factors['X'].push_back(sf_fake_aTop);
  // Top analysis
  scale_factors['x'].push_back(sf_ele_veto);
  scale_factors['x'].push_back(sf_muon_veto);
  scale_factors['x'].push_back(sf_btag_loose);
  scale_factors['x'].push_back(sf_fake_aTop);
  // Top analysis
  scale_factors['C'].push_back(sf_ele_medium);
  scale_factors['C'].push_back(sf_muon_medium);
  scale_factors['C'].push_back(sf_btag_loose);
  scale_factors['C'].push_back(sf_fake_mTop);
  // Top analysis
  scale_factors['V'].push_back(sf_ele_medium);
  scale_factors['V'].push_back(sf_muon_medium);
  scale_factors['V'].push_back(sf_top);
  // Top analysis
  scale_factors['F'].push_back(sf_ele_medium);
  scale_factors['F'].push_back(sf_muon_medium);
  scale_factors['F'].push_back(sf_fake_mTop);

/*
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
*/
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& data) {
  return apply_all_cuts('S');
  //return 0;
}

//_______________________________________________________
//                 List of Histograms

TH1D* h_njet_S;
TH1D* h_njet_Q;
TH1D* h_njet_T;
TH1D* h_njet_W;
TH1D* h_njet_S_MRR2;
TH1D* h_njet_S_3Jet;
TH1D* h_njet_Q_3Jet;
TH1D* h_njet_W_3Jet;
TH1D* h_njet_T_3Jet;
TH1D* h_njetAK8_S;
TH1D* h_njetAK8_Q;
TH1D* h_njetAK8_W;
TH1D* h_njetAK8_T;

TH1D* h_njet;
TH1D* h_nb;
TH1D* h_nw;
TH1D* h_ht_gen;
TH1D* h_ht_AK4;
TH1D* h_ht_AK8;

TH1D* h_ht_AK4_S;
TH1D* h_ht_AK8_S;
TH1D* h_jet1_pt_S;
TH1D* h_jet2_pt_S;
TH1D* h_jet3_pt_S;
TH1D* h_MR_S;
TH1D* h_MTR_S;
TH1D* h_R_S;
TH1D* h_R2_S;
TH1D* h_tau21_S;
TH1D* h_MET_S;
TH1D* h_AK8_jet1_pt_S;
TH1D* h_AK8_jet1_eta_S;
TH1D* h_MET_phi_S;
TH1D* h_AK8_jet1_pt_S_nj35;
TH1D* h_AK8_jet1_eta_S_nj35;
TH1D* h_MET_phi_S_nj35;
TH1D* h_AK8_jet1_pt_S_nj6;
TH1D* h_AK8_jet1_eta_S_nj6;
TH1D* h_MET_phi_S_nj6;

TH1D* h_ht_AK4_S_nj35;
TH1D* h_ht_AK8_S_nj35;
TH1D* h_jet1_pt_S_nj35;
TH1D* h_jet2_pt_S_nj35;
TH1D* h_jet3_pt_S_nj35;
TH1D* h_MR_S_nj35;
TH1D* h_MTR_S_nj35;
TH1D* h_R_S_nj35;
TH1D* h_R2_S_nj35;
TH1D* h_tau21_S_nj35;
TH1D* h_MET_S_nj35;

TH1D* h_ht_AK4_S_nj6;
TH1D* h_ht_AK8_S_nj6;
TH1D* h_jet1_pt_S_nj6;
TH1D* h_jet2_pt_S_nj6;
TH1D* h_jet3_pt_S_nj6;
TH1D* h_MR_S_nj6;
TH1D* h_MTR_S_nj6;
TH1D* h_R_S_nj6;
TH1D* h_R2_S_nj6;
TH1D* h_tau21_S_nj6;
TH1D* h_MET_S_nj6;

TH1D* h_ht_AK4_s;
TH1D* h_ht_AK8_s;
TH1D* h_jet1_pt_s;
TH1D* h_jet2_pt_s;
TH1D* h_jet3_pt_s;
TH1D* h_MR_s;
TH1D* h_MTR_s;
TH1D* h_R_s;
TH1D* h_R2_s;
TH1D* h_tau21_s;
TH1D* h_MET_s;
TH1D* h_AK8_jet1_pt_s;
TH1D* h_AK8_jet1_eta_s;
TH1D* h_MET_phi_s;
TH1D* h_AK8_jet1_pt_s_nj35;
TH1D* h_AK8_jet1_eta_s_nj35;
TH1D* h_MET_phi_s_nj35;
TH1D* h_AK8_jet1_pt_s_nj6;
TH1D* h_AK8_jet1_eta_s_nj6;
TH1D* h_MET_phi_s_nj6;

TH1D* h_ht_AK4_s_nj35;
TH1D* h_ht_AK8_s_nj35;
TH1D* h_jet1_pt_s_nj35;
TH1D* h_jet2_pt_s_nj35;
TH1D* h_jet3_pt_s_nj35;
TH1D* h_MR_s_nj35;
TH1D* h_MTR_s_nj35;
TH1D* h_R_s_nj35;
TH1D* h_R2_s_nj35;
TH1D* h_tau21_s_nj35;
TH1D* h_MET_s_nj35;

TH1D* h_ht_AK4_s_nj6;
TH1D* h_ht_AK8_s_nj6;
TH1D* h_jet1_pt_s_nj6;
TH1D* h_jet2_pt_s_nj6;
TH1D* h_jet3_pt_s_nj6;
TH1D* h_MR_s_nj6;
TH1D* h_MTR_s_nj6;
TH1D* h_R_s_nj6;
TH1D* h_R2_s_nj6;
TH1D* h_tau21_s_nj6;
TH1D* h_MET_s_nj6;

TH1D* h_ht_AK4_G;
TH1D* h_ht_AK8_G;
TH1D* h_jet1_pt_G;
TH1D* h_jet2_pt_G;
TH1D* h_jet3_pt_G;
TH1D* h_MR_G;
TH1D* h_MTR_G;
TH1D* h_R_G;
TH1D* h_R2_G;
TH1D* h_tau21_G;
TH1D* h_MET_G;
TH1D* h_AK8_jet1_pt_G;
TH1D* h_AK8_jet1_eta_G;
TH1D* h_MET_phi_G;
TH1D* h_AK8_jet1_pt_G_nj35;
TH1D* h_AK8_jet1_eta_G_nj35;
TH1D* h_MET_phi_G_nj35;
TH1D* h_AK8_jet1_pt_G_nj6;
TH1D* h_AK8_jet1_eta_G_nj6;
TH1D* h_MET_phi_G_nj6;

TH1D* h_ht_AK4_G_nj35;
TH1D* h_ht_AK8_G_nj35;
TH1D* h_jet1_pt_G_nj35;
TH1D* h_jet2_pt_G_nj35;
TH1D* h_jet3_pt_G_nj35;
TH1D* h_MR_G_nj35;
TH1D* h_MTR_G_nj35;
TH1D* h_R_G_nj35;
TH1D* h_R2_G_nj35;
TH1D* h_tau21_G_nj35;
TH1D* h_MET_G_nj35;

TH1D* h_ht_AK4_G_nj6;
TH1D* h_ht_AK8_G_nj6;
TH1D* h_jet1_pt_G_nj6;
TH1D* h_jet2_pt_G_nj6;
TH1D* h_jet3_pt_G_nj6;
TH1D* h_MR_G_nj6;
TH1D* h_MTR_G_nj6;
TH1D* h_R_G_nj6;
TH1D* h_R2_G_nj6;
TH1D* h_tau21_G_nj6;
TH1D* h_MET_G_nj6;

TH1D* h_ht_AK4_Z;
TH1D* h_ht_AK8_Z;
TH1D* h_jet1_pt_Z;
TH1D* h_jet2_pt_Z;
TH1D* h_jet3_pt_Z;
TH1D* h_MR_Z;
TH1D* h_MTR_Z;
TH1D* h_R_Z;
TH1D* h_R2_Z;
TH1D* h_tau21_Z;
TH1D* h_MET_Z;
TH1D* h_AK8_jet1_pt_Z;
TH1D* h_AK8_jet1_eta_Z;
TH1D* h_MET_phi_Z;
TH1D* h_AK8_jet1_pt_Z_nj35;
TH1D* h_AK8_jet1_eta_Z_nj35;
TH1D* h_MET_phi_Z_nj35;
TH1D* h_AK8_jet1_pt_Z_nj6;
TH1D* h_AK8_jet1_eta_Z_nj6;
TH1D* h_MET_phi_Z_nj6;

TH1D* h_ht_AK4_Z_nj35;
TH1D* h_ht_AK8_Z_nj35;
TH1D* h_jet1_pt_Z_nj35;
TH1D* h_jet2_pt_Z_nj35;
TH1D* h_jet3_pt_Z_nj35;
TH1D* h_MR_Z_nj35;
TH1D* h_MTR_Z_nj35;
TH1D* h_R_Z_nj35;
TH1D* h_R2_Z_nj35;
TH1D* h_tau21_Z_nj35;
TH1D* h_MET_Z_nj35;

TH1D* h_ht_AK4_Z_nj6;
TH1D* h_ht_AK8_Z_nj6;
TH1D* h_jet1_pt_Z_nj6;
TH1D* h_jet2_pt_Z_nj6;
TH1D* h_jet3_pt_Z_nj6;
TH1D* h_MR_Z_nj6;
TH1D* h_MTR_Z_nj6;
TH1D* h_R_Z_nj6;
TH1D* h_R2_Z_nj6;
TH1D* h_tau21_Z_nj6;
TH1D* h_MET_Z_nj6;

TH1D* h_ht_AK4_Q;
TH1D* h_ht_AK8_Q;
TH1D* h_jet1_pt_Q;
TH1D* h_jet2_pt_Q;
TH1D* h_jet3_pt_Q;
TH1D* h_MR_Q;
TH1D* h_MTR_Q;
TH1D* h_R_Q;
TH1D* h_R2_Q;
TH1D* h_tau21_Q;
TH1D* h_MET_Q;
TH1D* h_AK8_jet1_pt_Q;
TH1D* h_AK8_jet1_eta_Q;
TH1D* h_MET_phi_Q;
TH1D* h_AK8_jet1_pt_Q_nj35;
TH1D* h_AK8_jet1_eta_Q_nj35;
TH1D* h_MET_phi_Q_nj35;
TH1D* h_AK8_jet1_pt_Q_nj6;
TH1D* h_AK8_jet1_eta_Q_nj6;
TH1D* h_MET_phi_Q_nj6;

TH1D* h_ht_AK4_nodPhi_Q;
TH1D* h_ht_AK8_nodPhi_Q;
TH1D* h_jet1_pt_nodPhi_Q;
TH1D* h_jet2_pt_nodPhi_Q;
TH1D* h_jet3_pt_nodPhi_Q;
TH1D* h_MR_nodPhi_Q;
TH1D* h_MTR_nodPhi_Q;
TH1D* h_R_nodPhi_Q;
TH1D* h_R2_nodPhi_Q;
TH1D* h_tau21_nodPhi_Q;
TH1D* h_MET_nodPhi_Q;
TH1D* h_dPhiRazor_nodPhi_Q;

TH1D* h_ht_AK4_notau21_Q;
TH1D* h_ht_AK8_notau21_Q;
TH1D* h_jet1_pt_notau21_Q;
TH1D* h_jet2_pt_notau21_Q;
TH1D* h_jet3_pt_notau21_Q;
TH1D* h_MR_notau21_Q;
TH1D* h_MTR_notau21_Q;
TH1D* h_R_notau21_Q;
TH1D* h_R2_notau21_Q;
TH1D* h_tau21_notau21_Q;
TH1D* h_MET_notau21_Q;

TH1D* h_ht_AK4_nonb_Q;
TH1D* h_ht_AK8_nonb_Q;
TH1D* h_jet1_pt_nonb_Q;
TH1D* h_jet2_pt_nonb_Q;
TH1D* h_jet3_pt_nonb_Q;
TH1D* h_MR_nonb_Q;
TH1D* h_MTR_nonb_Q;
TH1D* h_R_nonb_Q;
TH1D* h_R2_nonb_Q;
TH1D* h_tau21_nonb_Q;
TH1D* h_MET_nonb_Q;

TH1D* h_ht_AK4_Q_nj35;
TH1D* h_ht_AK8_Q_nj35;
TH1D* h_jet1_pt_Q_nj35;
TH1D* h_jet2_pt_Q_nj35;
TH1D* h_jet3_pt_Q_nj35;
TH1D* h_MR_Q_nj35;
TH1D* h_MTR_Q_nj35;
TH1D* h_R_Q_nj35;
TH1D* h_R2_Q_nj35;
TH1D* h_tau21_Q_nj35;
TH1D* h_MET_Q_nj35;

TH1D* h_ht_AK4_nodPhi_Q_nj35;
TH1D* h_ht_AK8_nodPhi_Q_nj35;
TH1D* h_jet1_pt_nodPhi_Q_nj35;
TH1D* h_jet2_pt_nodPhi_Q_nj35;
TH1D* h_jet3_pt_nodPhi_Q_nj35;
TH1D* h_MR_nodPhi_Q_nj35;
TH1D* h_MTR_nodPhi_Q_nj35;
TH1D* h_R_nodPhi_Q_nj35;
TH1D* h_R2_nodPhi_Q_nj35;
TH1D* h_tau21_nodPhi_Q_nj35;
TH1D* h_MET_nodPhi_Q_nj35;
TH1D* h_dPhiRazor_nodPhi_Q_nj35;

TH1D* h_ht_AK4_notau21_Q_nj35;
TH1D* h_ht_AK8_notau21_Q_nj35;
TH1D* h_jet1_pt_notau21_Q_nj35;
TH1D* h_jet2_pt_notau21_Q_nj35;
TH1D* h_jet3_pt_notau21_Q_nj35;
TH1D* h_MR_notau21_Q_nj35;
TH1D* h_MTR_notau21_Q_nj35;
TH1D* h_R_notau21_Q_nj35;
TH1D* h_R2_notau21_Q_nj35;
TH1D* h_tau21_notau21_Q_nj35;
TH1D* h_MET_notau21_Q_nj35;

TH1D* h_ht_AK4_nonb_Q_nj35;
TH1D* h_ht_AK8_nonb_Q_nj35;
TH1D* h_jet1_pt_nonb_Q_nj35;
TH1D* h_jet2_pt_nonb_Q_nj35;
TH1D* h_jet3_pt_nonb_Q_nj35;
TH1D* h_MR_nonb_Q_nj35;
TH1D* h_MTR_nonb_Q_nj35;
TH1D* h_R_nonb_Q_nj35;
TH1D* h_R2_nonb_Q_nj35;
TH1D* h_tau21_nonb_Q_nj35;
TH1D* h_MET_nonb_Q_nj35;

TH1D* h_ht_AK4_Q_nj6;
TH1D* h_ht_AK8_Q_nj6;
TH1D* h_jet1_pt_Q_nj6;
TH1D* h_jet2_pt_Q_nj6;
TH1D* h_jet3_pt_Q_nj6;
TH1D* h_MR_Q_nj6;
TH1D* h_MTR_Q_nj6;
TH1D* h_R_Q_nj6;
TH1D* h_R2_Q_nj6;
TH1D* h_tau21_Q_nj6;
TH1D* h_MET_Q_nj6;

TH1D* h_ht_AK4_nodPhi_Q_nj6;
TH1D* h_ht_AK8_nodPhi_Q_nj6;
TH1D* h_jet1_pt_nodPhi_Q_nj6;
TH1D* h_jet2_pt_nodPhi_Q_nj6;
TH1D* h_jet3_pt_nodPhi_Q_nj6;
TH1D* h_MR_nodPhi_Q_nj6;
TH1D* h_MTR_nodPhi_Q_nj6;
TH1D* h_R_nodPhi_Q_nj6;
TH1D* h_R2_nodPhi_Q_nj6;
TH1D* h_tau21_nodPhi_Q_nj6;
TH1D* h_MET_nodPhi_Q_nj6;
TH1D* h_dPhiRazor_nodPhi_Q_nj6;

TH1D* h_ht_AK4_notau21_Q_nj6;
TH1D* h_ht_AK8_notau21_Q_nj6;
TH1D* h_jet1_pt_notau21_Q_nj6;
TH1D* h_jet2_pt_notau21_Q_nj6;
TH1D* h_jet3_pt_notau21_Q_nj6;
TH1D* h_MR_notau21_Q_nj6;
TH1D* h_MTR_notau21_Q_nj6;
TH1D* h_R_notau21_Q_nj6;
TH1D* h_R2_notau21_Q_nj6;
TH1D* h_tau21_notau21_Q_nj6;
TH1D* h_MET_notau21_Q_nj6;

TH1D* h_ht_AK4_nonb_Q_nj6;
TH1D* h_ht_AK8_nonb_Q_nj6;
TH1D* h_jet1_pt_nonb_Q_nj6;
TH1D* h_jet2_pt_nonb_Q_nj6;
TH1D* h_jet3_pt_nonb_Q_nj6;
TH1D* h_MR_nonb_Q_nj6;
TH1D* h_MTR_nonb_Q_nj6;
TH1D* h_R_nonb_Q_nj6;
TH1D* h_R2_nonb_Q_nj6;
TH1D* h_tau21_nonb_Q_nj6;
TH1D* h_MET_nonb_Q_nj6;

TH1D* h_ht_AK4_T;
TH1D* h_ht_AK8_T;
TH1D* h_jet1_pt_T;
TH1D* h_jet2_pt_T;
TH1D* h_jet3_pt_T;
TH1D* h_MR_T;
TH1D* h_MTR_T;
TH1D* h_R_T;
TH1D* h_R2_T;
TH1D* h_tau21_T;
TH1D* h_MET_T;
TH1D* h_AK8_jet1_pt_T;
TH1D* h_AK8_jet1_eta_T;
TH1D* h_MET_phi_T;
TH1D* h_AK8_jet1_pt_T_nj35;
TH1D* h_AK8_jet1_eta_T_nj35;
TH1D* h_MET_phi_T_nj35;
TH1D* h_AK8_jet1_pt_T_nj6;
TH1D* h_AK8_jet1_eta_T_nj6;
TH1D* h_MET_phi_T_nj6;

TH1D* h_ht_AK4_noMT_T;
TH1D* h_ht_AK8_noMT_T;
TH1D* h_jet1_pt_noMT_T;
TH1D* h_jet2_pt_noMT_T;
TH1D* h_jet3_pt_noMT_T;
TH1D* h_MR_noMT_T;
TH1D* h_MTR_noMT_T;
TH1D* h_R_noMT_T;
TH1D* h_R2_noMT_T;
TH1D* h_tau21_noMT_T;
TH1D* h_MET_noMT_T;

TH1D* h_ht_AK4_T_nj35;
TH1D* h_ht_AK8_T_nj35;
TH1D* h_jet1_pt_T_nj35;
TH1D* h_jet2_pt_T_nj35;
TH1D* h_jet3_pt_T_nj35;
TH1D* h_MR_T_nj35;
TH1D* h_MTR_T_nj35;
TH1D* h_R_T_nj35;
TH1D* h_R2_T_nj35;
TH1D* h_tau21_T_nj35;
TH1D* h_MET_T_nj35;

TH1D* h_ht_AK4_noMT_T_nj35;
TH1D* h_ht_AK8_noMT_T_nj35;
TH1D* h_jet1_pt_noMT_T_nj35;
TH1D* h_jet2_pt_noMT_T_nj35;
TH1D* h_jet3_pt_noMT_T_nj35;
TH1D* h_MR_noMT_T_nj35;
TH1D* h_MTR_noMT_T_nj35;
TH1D* h_R_noMT_T_nj35;
TH1D* h_R2_noMT_T_nj35;
TH1D* h_tau21_noMT_T_nj35;
TH1D* h_MET_noMT_T_nj35;

TH1D* h_ht_AK4_T_nj6;
TH1D* h_ht_AK8_T_nj6;
TH1D* h_jet1_pt_T_nj6;
TH1D* h_jet2_pt_T_nj6;
TH1D* h_jet3_pt_T_nj6;
TH1D* h_MR_T_nj6;
TH1D* h_MTR_T_nj6;
TH1D* h_R_T_nj6;
TH1D* h_R2_T_nj6;
TH1D* h_tau21_T_nj6;
TH1D* h_MET_T_nj6;

TH1D* h_ht_AK4_noMT_T_nj6;
TH1D* h_ht_AK8_noMT_T_nj6;
TH1D* h_jet1_pt_noMT_T_nj6;
TH1D* h_jet2_pt_noMT_T_nj6;
TH1D* h_jet3_pt_noMT_T_nj6;
TH1D* h_MR_noMT_T_nj6;
TH1D* h_MTR_noMT_T_nj6;
TH1D* h_R_noMT_T_nj6;
TH1D* h_R2_noMT_T_nj6;
TH1D* h_tau21_noMT_T_nj6;
TH1D* h_MET_noMT_T_nj6;

TH1D* h_ht_AK4_W;
TH1D* h_ht_AK8_W;
TH1D* h_jet1_pt_W;
TH1D* h_jet2_pt_W;
TH1D* h_jet3_pt_W;
TH1D* h_MR_W;
TH1D* h_MTR_W;
TH1D* h_R_W;
TH1D* h_R2_W;
TH1D* h_tau21_W;
TH1D* h_MET_W;
TH1D* h_AK8_jet1_pt_W;
TH1D* h_AK8_jet1_eta_W;
TH1D* h_MET_phi_W;
TH1D* h_AK8_jet1_pt_W_nj35;
TH1D* h_AK8_jet1_eta_W_nj35;
TH1D* h_MET_phi_W_nj35;
TH1D* h_AK8_jet1_pt_W_nj6;
TH1D* h_AK8_jet1_eta_W_nj6;
TH1D* h_MET_phi_W_nj6;

TH1D* h_ht_AK4_noMT_W;
TH1D* h_ht_AK8_noMT_W;
TH1D* h_jet1_pt_noMT_W;
TH1D* h_jet2_pt_noMT_W;
TH1D* h_jet3_pt_noMT_W;
TH1D* h_MR_noMT_W;
TH1D* h_MTR_noMT_W;
TH1D* h_R_noMT_W;
TH1D* h_R2_noMT_W;
TH1D* h_tau21_noMT_W;
TH1D* h_MET_noMT_W;

TH1D* h_ht_AK4_nonb_W;
TH1D* h_ht_AK8_nonb_W;
TH1D* h_jet1_pt_nonb_W;
TH1D* h_jet2_pt_nonb_W;
TH1D* h_jet3_pt_nonb_W;
TH1D* h_MR_nonb_W;
TH1D* h_MTR_nonb_W;
TH1D* h_R_nonb_W;
TH1D* h_R2_nonb_W;
TH1D* h_tau21_nonb_W;
TH1D* h_MET_nonb_W;

TH1D* h_ht_AK4_W_nj35;
TH1D* h_ht_AK8_W_nj35;
TH1D* h_jet1_pt_W_nj35;
TH1D* h_jet2_pt_W_nj35;
TH1D* h_jet3_pt_W_nj35;
TH1D* h_MR_W_nj35;
TH1D* h_MTR_W_nj35;
TH1D* h_R_W_nj35;
TH1D* h_R2_W_nj35;
TH1D* h_tau21_W_nj35;
TH1D* h_MET_W_nj35;

TH1D* h_ht_AK4_noMT_W_nj35;
TH1D* h_ht_AK8_noMT_W_nj35;
TH1D* h_jet1_pt_noMT_W_nj35;
TH1D* h_jet2_pt_noMT_W_nj35;
TH1D* h_jet3_pt_noMT_W_nj35;
TH1D* h_MR_noMT_W_nj35;
TH1D* h_MTR_noMT_W_nj35;
TH1D* h_R_noMT_W_nj35;
TH1D* h_R2_noMT_W_nj35;
TH1D* h_tau21_noMT_W_nj35;
TH1D* h_MET_noMT_W_nj35;

TH1D* h_ht_AK4_nonb_W_nj35;
TH1D* h_ht_AK8_nonb_W_nj35;
TH1D* h_jet1_pt_nonb_W_nj35;
TH1D* h_jet2_pt_nonb_W_nj35;
TH1D* h_jet3_pt_nonb_W_nj35;
TH1D* h_MR_nonb_W_nj35;
TH1D* h_MTR_nonb_W_nj35;
TH1D* h_R_nonb_W_nj35;
TH1D* h_R2_nonb_W_nj35;
TH1D* h_tau21_nonb_W_nj35;
TH1D* h_MET_nonb_W_nj35;

TH1D* h_ht_AK4_W_nj6;
TH1D* h_ht_AK8_W_nj6;
TH1D* h_jet1_pt_W_nj6;
TH1D* h_jet2_pt_W_nj6;
TH1D* h_jet3_pt_W_nj6;
TH1D* h_MR_W_nj6;
TH1D* h_MTR_W_nj6;
TH1D* h_R_W_nj6;
TH1D* h_R2_W_nj6;
TH1D* h_tau21_W_nj6;
TH1D* h_MET_W_nj6;

TH1D* h_ht_AK4_noMT_W_nj6;
TH1D* h_ht_AK8_noMT_W_nj6;
TH1D* h_jet1_pt_noMT_W_nj6;
TH1D* h_jet2_pt_noMT_W_nj6;
TH1D* h_jet3_pt_noMT_W_nj6;
TH1D* h_MR_noMT_W_nj6;
TH1D* h_MTR_noMT_W_nj6;
TH1D* h_R_noMT_W_nj6;
TH1D* h_R2_noMT_W_nj6;
TH1D* h_tau21_noMT_W_nj6;
TH1D* h_MET_noMT_W_nj6;

TH1D* h_ht_AK4_nonb_W_nj6;
TH1D* h_ht_AK8_nonb_W_nj6;
TH1D* h_jet1_pt_nonb_W_nj6;
TH1D* h_jet2_pt_nonb_W_nj6;
TH1D* h_jet3_pt_nonb_W_nj6;
TH1D* h_MR_nonb_W_nj6;
TH1D* h_MTR_nonb_W_nj6;
TH1D* h_R_nonb_W_nj6;
TH1D* h_R2_nonb_W_nj6;
TH1D* h_tau21_nonb_W_nj6;
TH1D* h_MET_nonb_W_nj6;

TH1D* h_ht_AK4_q;
TH1D* h_ht_AK8_q;
TH1D* h_jet1_pt_q;
TH1D* h_jet2_pt_q;
TH1D* h_jet3_pt_q;
TH1D* h_MR_q;
TH1D* h_MTR_q;
TH1D* h_R_q;
TH1D* h_R2_q;
TH1D* h_tau21_q;
TH1D* h_MET_q;
TH1D* h_AK8_jet1_pt_q;
TH1D* h_AK8_jet1_eta_q;
TH1D* h_MET_phi_q;
TH1D* h_AK8_jet1_pt_q_nj35;
TH1D* h_AK8_jet1_eta_q_nj35;
TH1D* h_MET_phi_q_nj35;
TH1D* h_AK8_jet1_pt_q_nj6;
TH1D* h_AK8_jet1_eta_q_nj6;
TH1D* h_MET_phi_q_nj6;

TH1D* h_ht_AK4_q_nj35;
TH1D* h_ht_AK8_q_nj35;
TH1D* h_jet1_pt_q_nj35;
TH1D* h_jet2_pt_q_nj35;
TH1D* h_jet3_pt_q_nj35;
TH1D* h_MR_q_nj35;
TH1D* h_MTR_q_nj35;
TH1D* h_R_q_nj35;
TH1D* h_R2_q_nj35;
TH1D* h_tau21_q_nj35;
TH1D* h_MET_q_nj35;

TH1D* h_ht_AK4_q_nj6;
TH1D* h_ht_AK8_q_nj6;
TH1D* h_jet1_pt_q_nj6;
TH1D* h_jet2_pt_q_nj6;
TH1D* h_jet3_pt_q_nj6;
TH1D* h_MR_q_nj6;
TH1D* h_MTR_q_nj6;
TH1D* h_R_q_nj6;
TH1D* h_R2_q_nj6;
TH1D* h_tau21_q_nj6;
TH1D* h_MET_q_nj6;


TH1D* h_softDropMass;
TH1D* h_StopMass;
TH1D* h_GluinoMass;
TH1D* h_LSPMass;
TH2D* h_GluinoLSPMass;

TH1D* h_ht_AK4_topana_S;
TH1D* h_ht_AK8_topana_S;
TH1D* h_jet1_pt_topana_S;
TH1D* h_jet2_pt_topana_S;
TH1D* h_jet3_pt_topana_S;
TH1D* h_MR_topana_S;
TH1D* h_MTR_topana_S;
TH1D* h_R_topana_S;
TH1D* h_R2_topana_S;
TH1D* h_tau21_topana_S;
TH1D* h_MET_topana_S;
TH1D* h_AK8_jet1_pt_topana_S;
TH1D* h_AK8_jet1_eta_topana_S;
TH1D* h_MET_phi_topana_S;

TH1D* h_ht_AK4_topana_s;
TH1D* h_ht_AK8_topana_s;
TH1D* h_jet1_pt_topana_s;
TH1D* h_jet2_pt_topana_s;
TH1D* h_jet3_pt_topana_s;
TH1D* h_MR_topana_s;
TH1D* h_MTR_topana_s;
TH1D* h_R_topana_s;
TH1D* h_R2_topana_s;
TH1D* h_tau21_topana_s;
TH1D* h_MET_topana_s;
TH1D* h_AK8_jet1_pt_topana_s;
TH1D* h_AK8_jet1_eta_topana_s;
TH1D* h_MET_phi_topana_s;

TH1D* h_ht_AK4_topana_Q;
TH1D* h_ht_AK8_topana_Q;
TH1D* h_jet1_pt_topana_Q;
TH1D* h_jet2_pt_topana_Q;
TH1D* h_jet3_pt_topana_Q;
TH1D* h_MR_topana_Q;
TH1D* h_MTR_topana_Q;
TH1D* h_R_topana_Q;
TH1D* h_R2_topana_Q;
TH1D* h_tau21_topana_Q;
TH1D* h_MET_topana_Q;
TH1D* h_AK8_jet1_pt_topana_Q;
TH1D* h_AK8_jet1_eta_topana_Q;
TH1D* h_MET_phi_topana_Q;

TH1D* h_ht_AK4_topana_q;
TH1D* h_ht_AK8_topana_q;
TH1D* h_jet1_pt_topana_q;
TH1D* h_jet2_pt_topana_q;
TH1D* h_jet3_pt_topana_q;
TH1D* h_MR_topana_q;
TH1D* h_MTR_topana_q;
TH1D* h_R_topana_q;
TH1D* h_R2_topana_q;
TH1D* h_tau21_topana_q;
TH1D* h_MET_topana_q;
TH1D* h_AK8_jet1_pt_topana_q;
TH1D* h_AK8_jet1_eta_topana_q;
TH1D* h_MET_phi_topana_q;

TH1D* h_ht_AK4_topana_T;
TH1D* h_ht_AK8_topana_T;
TH1D* h_jet1_pt_topana_T;
TH1D* h_jet2_pt_topana_T;
TH1D* h_jet3_pt_topana_T;
TH1D* h_MR_topana_T;
TH1D* h_MTR_topana_T;
TH1D* h_R_topana_T;
TH1D* h_R2_topana_T;
TH1D* h_tau21_topana_T;
TH1D* h_MET_topana_T;
TH1D* h_AK8_jet1_pt_topana_T;
TH1D* h_AK8_jet1_eta_topana_T;
TH1D* h_MET_phi_topana_T;

TH1D* h_ht_AK4_topana_W;
TH1D* h_ht_AK8_topana_W;
TH1D* h_jet1_pt_topana_W;
TH1D* h_jet2_pt_topana_W;
TH1D* h_jet3_pt_topana_W;
TH1D* h_MR_topana_W;
TH1D* h_MTR_topana_W;
TH1D* h_R_topana_W;
TH1D* h_R2_topana_W;
TH1D* h_tau21_topana_W;
TH1D* h_MET_topana_W;
TH1D* h_AK8_jet1_pt_topana_W;
TH1D* h_AK8_jet1_eta_topana_W;
TH1D* h_MET_phi_topana_W;

TH1D* h_ht_AK4_topana_Z;
TH1D* h_ht_AK8_topana_Z;
TH1D* h_jet1_pt_topana_Z;
TH1D* h_jet2_pt_topana_Z;
TH1D* h_jet3_pt_topana_Z;
TH1D* h_MR_topana_Z;
TH1D* h_MTR_topana_Z;
TH1D* h_R_topana_Z;
TH1D* h_R2_topana_Z;
TH1D* h_tau21_topana_Z;
TH1D* h_MET_topana_Z;
TH1D* h_AK8_jet1_pt_topana_Z;
TH1D* h_AK8_jet1_eta_topana_Z;
TH1D* h_MET_phi_topana_Z;

TH2D* h_R2_MR;

TH2D* h_R2_MR_topana_S;
TH2D* h_R2_MR_topana_s;
TH2D* h_R2_MR_topana_Q;
TH2D* h_R2_MR_topana_q;
TH2D* h_R2_MR_topana_T;
TH2D* h_R2_MR_topana_W;
TH2D* h_R2_MR_topana_Z;

TH2D* h_R2_MR_s;
TH2D* h_MR_MET_s;
TH2D* h_R2_MET_s;
TH2D* h_AK8Jet1pT_MET_s;
TH2D* h_R2_MR_q;
TH2D* h_MR_MET_q;
TH2D* h_R2_MET_q;
TH2D* h_AK8Jet1pT_MET_q;
TH2D* h_R2_MR_S;
TH2D* h_MR_MET_S;
TH2D* h_R2_MET_S;
TH2D* h_AK8Jet1pT_MET_S;
TH2D* h_R2_MR_Q;
TH2D* h_MR_MET_Q;
TH2D* h_R2_MET_Q;
TH2D* h_AK8Jet1pT_MET_Q;
TH2D* h_R2_MR_T;
TH2D* h_MR_MET_T;
TH2D* h_R2_MET_T;
TH2D* h_AK8Jet1pT_MET_T;
TH2D* h_R2_MR_W;
TH2D* h_MR_MET_W;
TH2D* h_R2_MET_W;
TH2D* h_AK8Jet1pT_MET_W;
TH2D* h_R2_MR_Z;
TH2D* h_MR_MET_Z;
TH2D* h_R2_MET_Z;
TH2D* h_AK8Jet1pT_MET_Z;
TH2D* h_R2_MR_G;
TH2D* h_MR_MET_G;
TH2D* h_R2_MET_G;
TH2D* h_AK8Jet1pT_MET_G;

TH2D* h_R2_MR_s_nj35;
TH2D* h_R2_MR_s_nj6;
TH2D* h_R2_MR_q_nj35;
TH2D* h_R2_MR_q_nj6;
TH2D* h_R2_MR_S_nj35;
TH2D* h_R2_MR_S_nj6;
TH2D* h_R2_MR_Q_nj35;
TH2D* h_R2_MR_Q_nj6;
TH2D* h_R2_MR_T_nj35;
TH2D* h_R2_MR_T_nj6;
TH2D* h_R2_MR_W_nj35;
TH2D* h_R2_MR_W_nj6;
TH2D* h_R2_MR_Z_nj35;
TH2D* h_R2_MR_Z_nj6;
TH2D* h_R2_MR_G_nj35;
TH2D* h_R2_MR_G_nj6;

TH1D *h_HT_pre;
TH1D *h_HT_pre_pass;

TH2D *h_HT_j1pt_pre;
TH2D *h_HT_j1pt_pre_pass;

TH1D* h_MR_0mW_EB_DPP;
TH1D* h_R2_0mW_EB_DPP;
TH2D* h_R2_MR_0mW_EB_DPP;
TH1D* h_MR_0mW_EB_FPP;
TH1D* h_R2_0mW_EB_FPP;
TH2D* h_R2_MR_0mW_EB_FPP;
TH1D* h_MR_0mW_EB_FP;
TH1D* h_R2_0mW_EB_FP;
TH2D* h_R2_MR_0mW_EB_FP;
TH1D* h_MR_1mW_EB_DPP;
TH1D* h_R2_1mW_EB_DPP;
TH2D* h_R2_MR_1mW_EB_DPP;
TH1D* h_MR_1mW_EB_FPP;
TH1D* h_R2_1mW_EB_FPP;
TH2D* h_R2_MR_1mW_EB_FPP;
TH1D* h_MR_1mW_EB_FP;
TH1D* h_R2_1mW_EB_FP;
TH2D* h_R2_MR_1mW_EB_FP;
TH1D* h_MR_0mW_EE_DPP;
TH1D* h_R2_0mW_EE_DPP;
TH2D* h_R2_MR_0mW_EE_DPP;
TH1D* h_MR_0mW_EE_FPP;
TH1D* h_R2_0mW_EE_FPP;
TH2D* h_R2_MR_0mW_EE_FPP;
TH1D* h_MR_0mW_EE_FP;
TH1D* h_R2_0mW_EE_FP;
TH2D* h_R2_MR_0mW_EE_FP;
TH1D* h_MR_1mW_EE_DPP;
TH1D* h_R2_1mW_EE_DPP;
TH2D* h_R2_MR_1mW_EE_DPP;
TH1D* h_MR_1mW_EE_FPP;
TH1D* h_R2_1mW_EE_FPP;
TH2D* h_R2_MR_1mW_EE_FPP;
TH1D* h_MR_1mW_EE_FP;
TH1D* h_R2_1mW_EE_FP;
TH2D* h_R2_MR_1mW_EE_FP;

TH1D* h_MR_zinv_0mW_EB;
TH1D* h_R2_zinv_0mW_EB;
TH2D* h_R2_MR_zinv_0mW_EB;
TH1D* h_MR_zinv_0mW_EE;
TH1D* h_R2_zinv_0mW_EE;
TH2D* h_R2_MR_zinv_0mW_EE;
TH1D* h_MR_zinv_1mW_EB;
TH1D* h_R2_zinv_1mW_EB;
TH2D* h_R2_MR_zinv_1mW_EB;
TH1D* h_MR_zinv_1mW_EE;
TH1D* h_R2_zinv_1mW_EE;
TH2D* h_R2_MR_zinv_1mW_EE;
TH1D* h_nmW_zinv_EB;
TH1D* h_MR_zinv_EB;
TH1D* h_R2_zinv_EB;
TH2D* h_R2_MR_zinv_EB;
TH1D* h_nmW_zinv_EE;
TH1D* h_MR_zinv_EE;
TH1D* h_R2_zinv_EE;
TH2D* h_R2_MR_zinv_EE;

TH1D* h_ChargedIso_MR_900_R2_010_SR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_010_SR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_010_SR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_010_SR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_010_SR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_014_SR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_014_SR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_014_SR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_014_SR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_014_SR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_020_SR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_020_SR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_020_SR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_020_SR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_020_SR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_037_SR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_037_SR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_037_SR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_037_SR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_037_SR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_075_SR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_075_SR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_075_SR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_075_SR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_075_SR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_010_CR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_010_CR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_010_CR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_010_CR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_010_CR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_014_CR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_014_CR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_014_CR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_014_CR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_014_CR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_020_CR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_020_CR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_020_CR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_020_CR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_020_CR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_037_CR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_037_CR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_037_CR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_037_CR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_037_CR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_075_CR_mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_075_CR_mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_075_CR_mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_075_CR_mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_075_CR_mW_EB;
TH1D* h_ChargedIso_MR_900_R2_010_SR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_010_SR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_010_SR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_010_SR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_010_SR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_014_SR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_014_SR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_014_SR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_014_SR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_014_SR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_020_SR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_020_SR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_020_SR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_020_SR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_020_SR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_037_SR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_037_SR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_037_SR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_037_SR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_037_SR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_075_SR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_075_SR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_075_SR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_075_SR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_075_SR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_010_CR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_010_CR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_010_CR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_010_CR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_010_CR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_014_CR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_014_CR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_014_CR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_014_CR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_014_CR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_020_CR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_020_CR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_020_CR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_020_CR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_020_CR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_037_CR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_037_CR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_037_CR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_037_CR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_037_CR_mW_EE;
TH1D* h_ChargedIso_MR_900_R2_075_CR_mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_075_CR_mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_075_CR_mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_075_CR_mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_075_CR_mW_EE;

TH1D* h_ChargedIso_MR_900_R2_010_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_010_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_010_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_010_SR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_010_SR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_014_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_014_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_014_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_014_SR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_014_SR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_020_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_020_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_020_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_020_SR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_020_SR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_037_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_037_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_037_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_037_SR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_037_SR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_075_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_075_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_075_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_075_SR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_075_SR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_010_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_010_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_010_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_010_CR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_010_CR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_014_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_014_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_014_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_014_CR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_014_CR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_020_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_020_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_020_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_020_CR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_020_CR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_037_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_037_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_037_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_037_CR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_037_CR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_075_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_075_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_075_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_075_CR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_075_CR_0mW_EB;
TH1D* h_ChargedIso_MR_900_R2_010_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_010_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_010_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_010_SR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_010_SR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_014_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_014_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_014_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_014_SR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_014_SR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_020_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_020_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_020_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_020_SR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_020_SR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_037_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_037_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_037_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_037_SR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_037_SR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_075_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_075_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_075_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_075_SR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_075_SR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_010_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_010_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_010_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_010_CR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_010_CR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_014_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_014_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_014_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_014_CR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_014_CR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_020_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_020_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_020_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_020_CR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_020_CR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_037_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_037_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_037_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_037_CR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_037_CR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_075_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_R2_075_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_R2_075_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_R2_075_CR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_R2_075_CR_1mW_EB;
TH1D* h_ChargedIso_MR_900_R2_010_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_010_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_010_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_010_SR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_010_SR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_014_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_014_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_014_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_014_SR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_014_SR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_020_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_020_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_020_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_020_SR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_020_SR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_037_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_037_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_037_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_037_SR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_037_SR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_075_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_075_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_075_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_075_SR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_075_SR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_010_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_010_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_010_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_010_CR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_010_CR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_014_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_014_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_014_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_014_CR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_014_CR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_020_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_020_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_020_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_020_CR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_020_CR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_037_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_037_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_037_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_037_CR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_037_CR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_075_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_075_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_075_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_075_CR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_075_CR_0mW_EE;
TH1D* h_ChargedIso_MR_900_R2_010_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_010_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_010_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_010_SR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_010_SR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_014_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_014_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_014_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_014_SR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_014_SR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_020_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_020_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_020_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_020_SR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_020_SR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_037_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_037_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_037_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_037_SR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_037_SR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_075_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_075_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_075_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_075_SR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_075_SR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_010_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_010_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_010_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_010_CR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_010_CR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_014_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_014_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_014_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_014_CR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_014_CR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_020_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_020_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_020_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_020_CR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_020_CR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_037_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_037_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_037_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_037_CR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_037_CR_1mW_EE;
TH1D* h_ChargedIso_MR_900_R2_075_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_R2_075_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_R2_075_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_R2_075_CR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_R2_075_CR_1mW_EE;

TH1D* h_ChargedIso_MR_900_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_SR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_SR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_SR_0mW_EB;
TH1D* h_ChargedIso_R2_010_SR_0mW_EB;
TH1D* h_ChargedIso_R2_014_SR_0mW_EB;
TH1D* h_ChargedIso_R2_020_SR_0mW_EB;
TH1D* h_ChargedIso_R2_037_SR_0mW_EB;
TH1D* h_ChargedIso_R2_075_SR_0mW_EB;
TH1D* h_ChargedIso_MR_900_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1100_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1400_CR_0mW_EB;
TH1D* h_ChargedIso_MR_1800_CR_0mW_EB;
TH1D* h_ChargedIso_MR_3000_CR_0mW_EB;
TH1D* h_ChargedIso_R2_010_CR_0mW_EB;
TH1D* h_ChargedIso_R2_014_CR_0mW_EB;
TH1D* h_ChargedIso_R2_020_CR_0mW_EB;
TH1D* h_ChargedIso_R2_037_CR_0mW_EB;
TH1D* h_ChargedIso_R2_075_CR_0mW_EB;
TH1D* h_ChargedIso_MR_900_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_SR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_SR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_SR_1mW_EB;
TH1D* h_ChargedIso_R2_010_SR_1mW_EB;
TH1D* h_ChargedIso_R2_014_SR_1mW_EB;
TH1D* h_ChargedIso_R2_020_SR_1mW_EB;
TH1D* h_ChargedIso_R2_037_SR_1mW_EB;
TH1D* h_ChargedIso_R2_075_SR_1mW_EB;
TH1D* h_ChargedIso_MR_900_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1100_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1400_CR_1mW_EB;
TH1D* h_ChargedIso_MR_1800_CR_1mW_EB;
TH1D* h_ChargedIso_MR_3000_CR_1mW_EB;
TH1D* h_ChargedIso_R2_010_CR_1mW_EB;
TH1D* h_ChargedIso_R2_014_CR_1mW_EB;
TH1D* h_ChargedIso_R2_020_CR_1mW_EB;
TH1D* h_ChargedIso_R2_037_CR_1mW_EB;
TH1D* h_ChargedIso_R2_075_CR_1mW_EB;
TH1D* h_ChargedIso_MR_900_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_SR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_SR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_SR_0mW_EE;
TH1D* h_ChargedIso_R2_010_SR_0mW_EE;
TH1D* h_ChargedIso_R2_014_SR_0mW_EE;
TH1D* h_ChargedIso_R2_020_SR_0mW_EE;
TH1D* h_ChargedIso_R2_037_SR_0mW_EE;
TH1D* h_ChargedIso_R2_075_SR_0mW_EE;
TH1D* h_ChargedIso_MR_900_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1100_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1400_CR_0mW_EE;
TH1D* h_ChargedIso_MR_1800_CR_0mW_EE;
TH1D* h_ChargedIso_MR_3000_CR_0mW_EE;
TH1D* h_ChargedIso_R2_010_CR_0mW_EE;
TH1D* h_ChargedIso_R2_014_CR_0mW_EE;
TH1D* h_ChargedIso_R2_020_CR_0mW_EE;
TH1D* h_ChargedIso_R2_037_CR_0mW_EE;
TH1D* h_ChargedIso_R2_075_CR_0mW_EE;
TH1D* h_ChargedIso_MR_900_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_SR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_SR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_SR_1mW_EE;
TH1D* h_ChargedIso_R2_010_SR_1mW_EE;
TH1D* h_ChargedIso_R2_014_SR_1mW_EE;
TH1D* h_ChargedIso_R2_020_SR_1mW_EE;
TH1D* h_ChargedIso_R2_037_SR_1mW_EE;
TH1D* h_ChargedIso_R2_075_SR_1mW_EE;
TH1D* h_ChargedIso_MR_900_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1100_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1400_CR_1mW_EE;
TH1D* h_ChargedIso_MR_1800_CR_1mW_EE;
TH1D* h_ChargedIso_MR_3000_CR_1mW_EE;
TH1D* h_ChargedIso_R2_010_CR_1mW_EE;
TH1D* h_ChargedIso_R2_014_CR_1mW_EE;
TH1D* h_ChargedIso_R2_020_CR_1mW_EE;
TH1D* h_ChargedIso_R2_037_CR_1mW_EE;
TH1D* h_ChargedIso_R2_075_CR_1mW_EE;

TH1D* h_ChargedIso_MR_900_SR_1mTop_EB;
TH1D* h_ChargedIso_MR_1100_SR_1mTop_EB;
TH1D* h_ChargedIso_MR_1400_SR_1mTop_EB;
TH1D* h_ChargedIso_MR_1800_SR_1mTop_EB;
TH1D* h_ChargedIso_MR_3000_SR_1mTop_EB;
TH1D* h_ChargedIso_R2_010_SR_1mTop_EB;
TH1D* h_ChargedIso_R2_014_SR_1mTop_EB;
TH1D* h_ChargedIso_R2_020_SR_1mTop_EB;
TH1D* h_ChargedIso_R2_037_SR_1mTop_EB;
TH1D* h_ChargedIso_R2_075_SR_1mTop_EB;
TH1D* h_ChargedIso_MR_900_CR_1mTop_EB;
TH1D* h_ChargedIso_MR_1100_CR_1mTop_EB;
TH1D* h_ChargedIso_MR_1400_CR_1mTop_EB;
TH1D* h_ChargedIso_MR_1800_CR_1mTop_EB;
TH1D* h_ChargedIso_MR_3000_CR_1mTop_EB;
TH1D* h_ChargedIso_R2_010_CR_1mTop_EB;
TH1D* h_ChargedIso_R2_014_CR_1mTop_EB;
TH1D* h_ChargedIso_R2_020_CR_1mTop_EB;
TH1D* h_ChargedIso_R2_037_CR_1mTop_EB;
TH1D* h_ChargedIso_R2_075_CR_1mTop_EB;
TH1D* h_ChargedIso_MR_900_SR_1mTop_EE;
TH1D* h_ChargedIso_MR_1100_SR_1mTop_EE;
TH1D* h_ChargedIso_MR_1400_SR_1mTop_EE;
TH1D* h_ChargedIso_MR_1800_SR_1mTop_EE;
TH1D* h_ChargedIso_MR_3000_SR_1mTop_EE;
TH1D* h_ChargedIso_R2_010_SR_1mTop_EE;
TH1D* h_ChargedIso_R2_014_SR_1mTop_EE;
TH1D* h_ChargedIso_R2_020_SR_1mTop_EE;
TH1D* h_ChargedIso_R2_037_SR_1mTop_EE;
TH1D* h_ChargedIso_R2_075_SR_1mTop_EE;
TH1D* h_ChargedIso_MR_900_CR_1mTop_EE;
TH1D* h_ChargedIso_MR_1100_CR_1mTop_EE;
TH1D* h_ChargedIso_MR_1400_CR_1mTop_EE;
TH1D* h_ChargedIso_MR_1800_CR_1mTop_EE;
TH1D* h_ChargedIso_MR_3000_CR_1mTop_EE;
TH1D* h_ChargedIso_R2_010_CR_1mTop_EE;
TH1D* h_ChargedIso_R2_014_CR_1mTop_EE;
TH1D* h_ChargedIso_R2_020_CR_1mTop_EE;
TH1D* h_ChargedIso_R2_037_CR_1mTop_EE;
TH1D* h_ChargedIso_R2_075_CR_1mTop_EE;

TH1D* h_AK8Jet1Pt_W_fakerate;
TH1D* h_AK8Jet1Pt_no_W_fakerate;
TH1D* h_AK8Jet1Pt_mW_fakerate;
TH1D* h_AK8Jet1Pt_no_mW_fakerate;
TH1D* h_AK8Jet1Pt_m0bW_fakerate;
TH1D* h_AK8Jet1Pt_no_m0bW_fakerate;
TH1D* h_AK8Jet1Pt_aW_fakerate;
TH1D* h_AK8Jet1Pt_no_aW_fakerate;

TH1D* h_AK8Jet1Pt_Top_fakerate;
TH1D* h_AK8Jet1Pt_no_Top_fakerate;
TH1D* h_AK8Jet1Pt_mTop_fakerate;
TH1D* h_AK8Jet1Pt_no_mTop_fakerate;
TH1D* h_AK8Jet1Pt_m0bTop_fakerate;
TH1D* h_AK8Jet1Pt_no_m0bTop_fakerate;
TH1D* h_AK8Jet1Pt_aTop_fakerate;
TH1D* h_AK8Jet1Pt_no_aTop_fakerate;

TH2D* h_AK8Jet1Pt_Eta_W_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_W_fakerate;
TH2D* h_AK8Jet1Pt_Eta_mW_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_mW_fakerate;
TH2D* h_AK8Jet1Pt_Eta_m0bW_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_m0bW_fakerate;
TH2D* h_AK8Jet1Pt_Eta_aW_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_aW_fakerate;

TH2D* h_AK8Jet1Pt_Eta_Top_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_Top_fakerate;
TH2D* h_AK8Jet1Pt_Eta_mTop_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_mTop_fakerate;
TH2D* h_AK8Jet1Pt_Eta_m0bTop_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_m0bTop_fakerate;
TH2D* h_AK8Jet1Pt_Eta_aTop_fakerate;
TH2D* h_AK8Jet1Pt_Eta_no_aTop_fakerate;

TH1D* h_tau1_GenW_W;
TH1D* h_tau1_GenW_no_W;
TH1D* h_tau1_GenTop_Top;
TH1D* h_tau1_GenTop_no_Top;

TH1D* h_tau2_GenW_W;
TH1D* h_tau2_GenW_no_W;
TH1D* h_tau2_GenTop_Top;
TH1D* h_tau2_GenTop_no_Top;

TH1D* h_tau3_GenW_W;
TH1D* h_tau3_GenW_no_W;
TH1D* h_tau3_GenTop_Top;
TH1D* h_tau3_GenTop_no_Top;

TH1D* h_tau21_GenW_W;
TH1D* h_tau21_GenW_no_W;
TH1D* h_tau21_GenTop_Top;
TH1D* h_tau21_GenTop_no_Top;

TH1D* h_tau32_GenW_W;
TH1D* h_tau32_GenW_no_W;
TH1D* h_tau32_GenTop_Top;
TH1D* h_tau32_GenTop_no_Top;

TH1D* h_SubjetBTag_GenTop_Top;
TH1D* h_SubjetBTag_GenTop_no_Top;

TH1D* h_AK8Jet1Pt_GenW_W;
TH1D* h_AK8Jet1Pt_GenW_no_W;
TH1D* h_AK8Jet1Pt_GenTop_Top;
TH1D* h_AK8Jet1Pt_GenTop_no_Top;

TH2D* h_AK8Jet1Pt_eta_GenW_W;
TH2D* h_AK8Jet1Pt_eta_GenW_no_W;
TH2D* h_AK8Jet1Pt_eta_GenTop_Top;
TH2D* h_AK8Jet1Pt_eta_GenTop_no_Top;

std::vector<TH1D*> vh_ht_AK4_S;
std::vector<TH1D*> vh_jet1_pt_S;
std::vector<TH1D*> vh_jet2_pt_S;
std::vector<TH1D*> vh_jet3_pt_S;
std::vector<TH1D*> vh_MR_S;
std::vector<TH1D*> vh_MTR_S;
std::vector<TH1D*> vh_R2_S;
std::vector<TH1D*> vh_tau21_S;
std::vector<TH1D*> vh_MET_S;
std::vector<TH1D*> vh_AK8_jet1_pt_S;
std::vector<TH1D*> vh_AK8_jet1_eta_S;
std::vector<TH1D*> vh_MET_phi_S;
std::vector<TH1D*> vh_AK8_jet1_pt_S_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_S_nj35;
std::vector<TH1D*> vh_MET_phi_S_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_S_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_S_nj6;
std::vector<TH1D*> vh_MET_phi_S_nj6;

std::vector<TH1D*> vh_ht_AK4_S_nj35;
std::vector<TH1D*> vh_jet1_pt_S_nj35;
std::vector<TH1D*> vh_jet2_pt_S_nj35;
std::vector<TH1D*> vh_jet3_pt_S_nj35;
std::vector<TH1D*> vh_MR_S_nj35;
std::vector<TH1D*> vh_MTR_S_nj35;
std::vector<TH1D*> vh_R2_S_nj35;
std::vector<TH1D*> vh_tau21_S_nj35;
std::vector<TH1D*> vh_MET_S_nj35;

std::vector<TH1D*> vh_ht_AK4_S_nj6;
std::vector<TH1D*> vh_jet1_pt_S_nj6;
std::vector<TH1D*> vh_jet2_pt_S_nj6;
std::vector<TH1D*> vh_jet3_pt_S_nj6;
std::vector<TH1D*> vh_MR_S_nj6;
std::vector<TH1D*> vh_MTR_S_nj6;
std::vector<TH1D*> vh_R2_S_nj6;
std::vector<TH1D*> vh_tau21_S_nj6;
std::vector<TH1D*> vh_MET_S_nj6;

std::vector<TH1D*> vh_ht_AK4_s;
std::vector<TH1D*> vh_jet1_pt_s;
std::vector<TH1D*> vh_jet2_pt_s;
std::vector<TH1D*> vh_jet3_pt_s;
std::vector<TH1D*> vh_MR_s;
std::vector<TH1D*> vh_MTR_s;
std::vector<TH1D*> vh_R2_s;
std::vector<TH1D*> vh_tau21_s;
std::vector<TH1D*> vh_MET_s;
std::vector<TH1D*> vh_AK8_jet1_pt_s;
std::vector<TH1D*> vh_AK8_jet1_eta_s;
std::vector<TH1D*> vh_MET_phi_s;
std::vector<TH1D*> vh_AK8_jet1_pt_s_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_s_nj35;
std::vector<TH1D*> vh_MET_phi_s_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_s_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_s_nj6;
std::vector<TH1D*> vh_MET_phi_s_nj6;

std::vector<TH1D*> vh_ht_AK4_s_nj35;
std::vector<TH1D*> vh_jet1_pt_s_nj35;
std::vector<TH1D*> vh_jet2_pt_s_nj35;
std::vector<TH1D*> vh_jet3_pt_s_nj35;
std::vector<TH1D*> vh_MR_s_nj35;
std::vector<TH1D*> vh_MTR_s_nj35;
std::vector<TH1D*> vh_R2_s_nj35;
std::vector<TH1D*> vh_tau21_s_nj35;
std::vector<TH1D*> vh_MET_s_nj35;

std::vector<TH1D*> vh_ht_AK4_s_nj6;
std::vector<TH1D*> vh_jet1_pt_s_nj6;
std::vector<TH1D*> vh_jet2_pt_s_nj6;
std::vector<TH1D*> vh_jet3_pt_s_nj6;
std::vector<TH1D*> vh_MR_s_nj6;
std::vector<TH1D*> vh_MTR_s_nj6;
std::vector<TH1D*> vh_R2_s_nj6;
std::vector<TH1D*> vh_tau21_s_nj6;
std::vector<TH1D*> vh_MET_s_nj6;

std::vector<TH1D*> vh_ht_AK4_G;
std::vector<TH1D*> vh_jet1_pt_G;
std::vector<TH1D*> vh_jet2_pt_G;
std::vector<TH1D*> vh_jet3_pt_G;
std::vector<TH1D*> vh_MR_G;
std::vector<TH1D*> vh_MTR_G;
std::vector<TH1D*> vh_R2_G;
std::vector<TH1D*> vh_tau21_G;
std::vector<TH1D*> vh_MET_G;
std::vector<TH1D*> vh_AK8_jet1_pt_G;
std::vector<TH1D*> vh_AK8_jet1_eta_G;
std::vector<TH1D*> vh_MET_phi_G;
std::vector<TH1D*> vh_AK8_jet1_pt_G_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_G_nj35;
std::vector<TH1D*> vh_MET_phi_G_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_G_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_G_nj6;
std::vector<TH1D*> vh_MET_phi_G_nj6;

std::vector<TH1D*> vh_ht_AK4_G_nj35;
std::vector<TH1D*> vh_jet1_pt_G_nj35;
std::vector<TH1D*> vh_jet2_pt_G_nj35;
std::vector<TH1D*> vh_jet3_pt_G_nj35;
std::vector<TH1D*> vh_MR_G_nj35;
std::vector<TH1D*> vh_MTR_G_nj35;
std::vector<TH1D*> vh_R2_G_nj35;
std::vector<TH1D*> vh_tau21_G_nj35;
std::vector<TH1D*> vh_MET_G_nj35;

std::vector<TH1D*> vh_ht_AK4_G_nj6;
std::vector<TH1D*> vh_jet1_pt_G_nj6;
std::vector<TH1D*> vh_jet2_pt_G_nj6;
std::vector<TH1D*> vh_jet3_pt_G_nj6;
std::vector<TH1D*> vh_MR_G_nj6;
std::vector<TH1D*> vh_MTR_G_nj6;
std::vector<TH1D*> vh_R2_G_nj6;
std::vector<TH1D*> vh_tau21_G_nj6;
std::vector<TH1D*> vh_MET_G_nj6;

std::vector<TH1D*> vh_ht_AK4_Z;
std::vector<TH1D*> vh_jet1_pt_Z;
std::vector<TH1D*> vh_jet2_pt_Z;
std::vector<TH1D*> vh_jet3_pt_Z;
std::vector<TH1D*> vh_MR_Z;
std::vector<TH1D*> vh_MTR_Z;
std::vector<TH1D*> vh_R2_Z;
std::vector<TH1D*> vh_tau21_Z;
std::vector<TH1D*> vh_MET_Z;
std::vector<TH1D*> vh_AK8_jet1_pt_Z;
std::vector<TH1D*> vh_AK8_jet1_eta_Z;
std::vector<TH1D*> vh_MET_phi_Z;
std::vector<TH1D*> vh_AK8_jet1_pt_Z_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_Z_nj35;
std::vector<TH1D*> vh_MET_phi_Z_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_Z_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_Z_nj6;
std::vector<TH1D*> vh_MET_phi_Z_nj6;

std::vector<TH1D*> vh_ht_AK4_Z_nj35;
std::vector<TH1D*> vh_jet1_pt_Z_nj35;
std::vector<TH1D*> vh_jet2_pt_Z_nj35;
std::vector<TH1D*> vh_jet3_pt_Z_nj35;
std::vector<TH1D*> vh_MR_Z_nj35;
std::vector<TH1D*> vh_MTR_Z_nj35;
std::vector<TH1D*> vh_R2_Z_nj35;
std::vector<TH1D*> vh_tau21_Z_nj35;
std::vector<TH1D*> vh_MET_Z_nj35;

std::vector<TH1D*> vh_ht_AK4_Z_nj6;
std::vector<TH1D*> vh_jet1_pt_Z_nj6;
std::vector<TH1D*> vh_jet2_pt_Z_nj6;
std::vector<TH1D*> vh_jet3_pt_Z_nj6;
std::vector<TH1D*> vh_MR_Z_nj6;
std::vector<TH1D*> vh_MTR_Z_nj6;
std::vector<TH1D*> vh_R2_Z_nj6;
std::vector<TH1D*> vh_tau21_Z_nj6;
std::vector<TH1D*> vh_MET_Z_nj6;

std::vector<TH1D*> vh_ht_AK4_Q;
std::vector<TH1D*> vh_jet1_pt_Q;
std::vector<TH1D*> vh_jet2_pt_Q;
std::vector<TH1D*> vh_jet3_pt_Q;
std::vector<TH1D*> vh_MR_Q;
std::vector<TH1D*> vh_MTR_Q;
std::vector<TH1D*> vh_R2_Q;
std::vector<TH1D*> vh_tau21_Q;
std::vector<TH1D*> vh_MET_Q;
std::vector<TH1D*> vh_AK8_jet1_pt_Q;
std::vector<TH1D*> vh_AK8_jet1_eta_Q;
std::vector<TH1D*> vh_MET_phi_Q;
std::vector<TH1D*> vh_AK8_jet1_pt_Q_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_Q_nj35;
std::vector<TH1D*> vh_MET_phi_Q_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_Q_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_Q_nj6;
std::vector<TH1D*> vh_MET_phi_Q_nj6;

std::vector<TH1D*> vh_ht_AK4_Q_nj35;
std::vector<TH1D*> vh_jet1_pt_Q_nj35;
std::vector<TH1D*> vh_jet2_pt_Q_nj35;
std::vector<TH1D*> vh_jet3_pt_Q_nj35;
std::vector<TH1D*> vh_MR_Q_nj35;
std::vector<TH1D*> vh_MTR_Q_nj35;
std::vector<TH1D*> vh_R2_Q_nj35;
std::vector<TH1D*> vh_tau21_Q_nj35;
std::vector<TH1D*> vh_MET_Q_nj35;

std::vector<TH1D*> vh_ht_AK4_Q_nj6;
std::vector<TH1D*> vh_jet1_pt_Q_nj6;
std::vector<TH1D*> vh_jet2_pt_Q_nj6;
std::vector<TH1D*> vh_jet3_pt_Q_nj6;
std::vector<TH1D*> vh_MR_Q_nj6;
std::vector<TH1D*> vh_MTR_Q_nj6;
std::vector<TH1D*> vh_R2_Q_nj6;
std::vector<TH1D*> vh_tau21_Q_nj6;
std::vector<TH1D*> vh_MET_Q_nj6;

std::vector<TH1D*> vh_ht_AK4_T;
std::vector<TH1D*> vh_jet1_pt_T;
std::vector<TH1D*> vh_jet2_pt_T;
std::vector<TH1D*> vh_jet3_pt_T;
std::vector<TH1D*> vh_MR_T;
std::vector<TH1D*> vh_MTR_T;
std::vector<TH1D*> vh_R2_T;
std::vector<TH1D*> vh_tau21_T;
std::vector<TH1D*> vh_MET_T;
std::vector<TH1D*> vh_AK8_jet1_pt_T;
std::vector<TH1D*> vh_AK8_jet1_eta_T;
std::vector<TH1D*> vh_MET_phi_T;
std::vector<TH1D*> vh_AK8_jet1_pt_T_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_T_nj35;
std::vector<TH1D*> vh_MET_phi_T_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_T_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_T_nj6;
std::vector<TH1D*> vh_MET_phi_T_nj6;

std::vector<TH1D*> vh_ht_AK4_T_nj35;
std::vector<TH1D*> vh_jet1_pt_T_nj35;
std::vector<TH1D*> vh_jet2_pt_T_nj35;
std::vector<TH1D*> vh_jet3_pt_T_nj35;
std::vector<TH1D*> vh_MR_T_nj35;
std::vector<TH1D*> vh_MTR_T_nj35;
std::vector<TH1D*> vh_R2_T_nj35;
std::vector<TH1D*> vh_tau21_T_nj35;
std::vector<TH1D*> vh_MET_T_nj35;

std::vector<TH1D*> vh_ht_AK4_T_nj6;
std::vector<TH1D*> vh_jet1_pt_T_nj6;
std::vector<TH1D*> vh_jet2_pt_T_nj6;
std::vector<TH1D*> vh_jet3_pt_T_nj6;
std::vector<TH1D*> vh_MR_T_nj6;
std::vector<TH1D*> vh_MTR_T_nj6;
std::vector<TH1D*> vh_R2_T_nj6;
std::vector<TH1D*> vh_tau21_T_nj6;
std::vector<TH1D*> vh_MET_T_nj6;

std::vector<TH1D*> vh_ht_AK4_W;
std::vector<TH1D*> vh_jet1_pt_W;
std::vector<TH1D*> vh_jet2_pt_W;
std::vector<TH1D*> vh_jet3_pt_W;
std::vector<TH1D*> vh_MR_W;
std::vector<TH1D*> vh_MTR_W;
std::vector<TH1D*> vh_R2_W;
std::vector<TH1D*> vh_tau21_W;
std::vector<TH1D*> vh_MET_W;
std::vector<TH1D*> vh_AK8_jet1_pt_W;
std::vector<TH1D*> vh_AK8_jet1_eta_W;
std::vector<TH1D*> vh_MET_phi_W;
std::vector<TH1D*> vh_AK8_jet1_pt_W_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_W_nj35;
std::vector<TH1D*> vh_MET_phi_W_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_W_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_W_nj6;
std::vector<TH1D*> vh_MET_phi_W_nj6;

std::vector<TH1D*> vh_ht_AK4_W_nj35;
std::vector<TH1D*> vh_jet1_pt_W_nj35;
std::vector<TH1D*> vh_jet2_pt_W_nj35;
std::vector<TH1D*> vh_jet3_pt_W_nj35;
std::vector<TH1D*> vh_MR_W_nj35;
std::vector<TH1D*> vh_MTR_W_nj35;
std::vector<TH1D*> vh_R2_W_nj35;
std::vector<TH1D*> vh_tau21_W_nj35;
std::vector<TH1D*> vh_MET_W_nj35;

std::vector<TH1D*> vh_ht_AK4_W_nj6;
std::vector<TH1D*> vh_jet1_pt_W_nj6;
std::vector<TH1D*> vh_jet2_pt_W_nj6;
std::vector<TH1D*> vh_jet3_pt_W_nj6;
std::vector<TH1D*> vh_MR_W_nj6;
std::vector<TH1D*> vh_MTR_W_nj6;
std::vector<TH1D*> vh_R2_W_nj6;
std::vector<TH1D*> vh_tau21_W_nj6;
std::vector<TH1D*> vh_MET_W_nj6;

std::vector<TH1D*> vh_ht_AK4_q;
std::vector<TH1D*> vh_jet1_pt_q;
std::vector<TH1D*> vh_jet2_pt_q;
std::vector<TH1D*> vh_jet3_pt_q;
std::vector<TH1D*> vh_MR_q;
std::vector<TH1D*> vh_MTR_q;
std::vector<TH1D*> vh_R2_q;
std::vector<TH1D*> vh_tau21_q;
std::vector<TH1D*> vh_MET_q;
std::vector<TH1D*> vh_AK8_jet1_pt_q;
std::vector<TH1D*> vh_AK8_jet1_eta_q;
std::vector<TH1D*> vh_MET_phi_q;
std::vector<TH1D*> vh_AK8_jet1_pt_q_nj35;
std::vector<TH1D*> vh_AK8_jet1_eta_q_nj35;
std::vector<TH1D*> vh_MET_phi_q_nj35;
std::vector<TH1D*> vh_AK8_jet1_pt_q_nj6;
std::vector<TH1D*> vh_AK8_jet1_eta_q_nj6;
std::vector<TH1D*> vh_MET_phi_q_nj6;

std::vector<TH1D*> vh_ht_AK4_q_nj35;
std::vector<TH1D*> vh_jet1_pt_q_nj35;
std::vector<TH1D*> vh_jet2_pt_q_nj35;
std::vector<TH1D*> vh_jet3_pt_q_nj35;
std::vector<TH1D*> vh_MR_q_nj35;
std::vector<TH1D*> vh_MTR_q_nj35;
std::vector<TH1D*> vh_R2_q_nj35;
std::vector<TH1D*> vh_tau21_q_nj35;
std::vector<TH1D*> vh_MET_q_nj35;

std::vector<TH1D*> vh_ht_AK4_q_nj6;
std::vector<TH1D*> vh_jet1_pt_q_nj6;
std::vector<TH1D*> vh_jet2_pt_q_nj6;
std::vector<TH1D*> vh_jet3_pt_q_nj6;
std::vector<TH1D*> vh_MR_q_nj6;
std::vector<TH1D*> vh_MTR_q_nj6;
std::vector<TH1D*> vh_R2_q_nj6;
std::vector<TH1D*> vh_tau21_q_nj6;
std::vector<TH1D*> vh_MET_q_nj6;

std::vector<TH2D*> vh_R2_MR_s;
std::vector<TH2D*> vh_R2_MR_q;
std::vector<TH2D*> vh_R2_MR_S;
std::vector<TH2D*> vh_R2_MR_Q;
std::vector<TH2D*> vh_R2_MR_T;
std::vector<TH2D*> vh_R2_MR_W;
std::vector<TH2D*> vh_R2_MR_Z;
std::vector<TH2D*> vh_R2_MR_G;

std::vector<TH2D*> vh_R2_MR_s_nj35;
std::vector<TH2D*> vh_R2_MR_s_nj6;
std::vector<TH2D*> vh_R2_MR_q_nj35;
std::vector<TH2D*> vh_R2_MR_q_nj6;
std::vector<TH2D*> vh_R2_MR_S_nj35;
std::vector<TH2D*> vh_R2_MR_S_nj6;
std::vector<TH2D*> vh_R2_MR_Q_nj35;
std::vector<TH2D*> vh_R2_MR_Q_nj6;
std::vector<TH2D*> vh_R2_MR_T_nj35;
std::vector<TH2D*> vh_R2_MR_T_nj6;
std::vector<TH2D*> vh_R2_MR_W_nj35;
std::vector<TH2D*> vh_R2_MR_W_nj6;
std::vector<TH2D*> vh_R2_MR_Z_nj35;
std::vector<TH2D*> vh_R2_MR_Z_nj6;
std::vector<TH2D*> vh_R2_MR_G_nj35;
std::vector<TH2D*> vh_R2_MR_G_nj6;


//std::vector<TH1D*> vh_jet1_pt;

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{

  //double htbn = 20;
  //double htmn = 0;
 // double htmx = 1500;

  //double j1ptbn = 20;
  //double j1ptmn = 0;
  //double j1ptmx = 1000;

  //double MRbn = 7;
  //double MRmn = 0;
  //double MRmx = 4000;
  //double Rbn = 7;
  //double Rmn = 0;
  //double Rmx = 1;

  // Variable binning
  int nbn_HT = 19;
  int nbnHT = 10;
  int nbnj1pt = 7;
  //int nbn_j1pt = 12;
  //int nbn_AK8j1pt = 13;
  int nbn_AK8J1pt = 15;
  int nbn_MR = 7;
  int nbn_R = 6;
  int nbn_R2 = 7;
  int nbn_eta = 2;
  //Double_t  HT_bins[20]  = {0, 200, 300, 400, 500, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500, 2000, 4000, 10000};
  Double_t bn_HTtmp[] = {400.,500.,600.,700.,750.,800.,850.,900.,950.,1000.,1500.};
  Double_t* bnHT = 0;
  bnHT = utils::getVariableBinEdges(nbnHT+1,bn_HTtmp);
  Double_t bn_j1pttmp[] = {200.,300.,400.,450.,500.,550.,600.,1000.};
  Double_t* bnj1pt = 0;
  bnj1pt = utils::getVariableBinEdges(nbnj1pt+1,bn_j1pttmp);

  Double_t bn_HT_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,900.,1000.,2500.};
  Double_t* bn_HT = 0;
  bn_HT = utils::getVariableBinEdges(nbn_HT+1,bn_HT_tmp);
  //Double_t bn_j1pt_tmp[] = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,700.,1000.};
  //Double_t* bn_j1pt = 0;
 // bn_j1pt = utils::getVariableBinEdges(nbn_j1pt+1,bn_j1pt_tmp);
  //Double_t bn_AK8j1pt_tmp[] = {0.,200.,220.,240,260,280,300.,320.,340,360,380,400.,500.,1000.};
  //Double_t* bn_AK8j1pt = 0;
  //bn_AK8j1pt = utils::getVariableBinEdges(nbn_AK8j1pt+1,bn_AK8j1pt_tmp);
  //Double_t bn_AK8J1pt_tmp[] = {200.,250,350.,1000.};
  Double_t bn_AK8J1pt_tmp[] = {0.,200.,220.,240,260,280,300.,320.,340,360,380,400.,450,500.,700,1000.};
  Double_t* bn_AK8J1pt = 0;
  bn_AK8J1pt = utils::getVariableBinEdges(nbn_AK8J1pt+1,bn_AK8J1pt_tmp);
  Double_t bn_MR_tmp[] = {0.,600.,800.,1000.,1200.,1600.,2000.,4000.};
  Double_t* bn_MR = 0;
  bn_MR = utils::getVariableBinEdges(nbn_MR+1,bn_MR_tmp);
  Double_t bn_R_tmp[] = {0.,0.08,.12,.16,0.24,0.5,1.};
  Double_t* bn_R = 0;
  bn_R = utils::getVariableBinEdges(nbn_R+1,bn_R_tmp);
  Double_t bn_R2_tmp[] = {0.,0.04,0.08,0.12,0.16,0.24,0.5,1.};
  Double_t* bn_R2 = 0;
  bn_R2 = utils::getVariableBinEdges(nbn_R2+1,bn_R2_tmp);
  Double_t bn_eta_tmp[] = {0.,1.5,2.5};
  Double_t* bn_eta = 0;
  bn_eta = utils::getVariableBinEdges(nbn_eta+1,bn_eta_tmp);

  h_njet_S = new TH1D("njet_S",         ";N_{jet}",                20, 0,  20);
  h_njet_Q = new TH1D("njet_Q",         ";N_{jet}",                20, 0,  20);
  h_njet_T = new TH1D("njet_T",         ";N_{jet}",                20, 0,  20);
  h_njet_W = new TH1D("njet_W",         ";N_{jet}",                20, 0,  20);
  h_njet_S_MRR2 = new TH1D("njet_S_MRR2",         ";N_{jet}",                20, 0,  20);
  h_njet_S_3Jet = new TH1D("njet_S_3Jet",         ";N_{jet}",                20, 0,  20);
  h_njet_Q_3Jet = new TH1D("njet_Q_3Jet",         ";N_{jet}",                20, 0,  20);
  h_njet_W_3Jet = new TH1D("njet_W_3Jet",         ";N_{jet}",                20, 0,  20);
  h_njet_T_3Jet = new TH1D("njet_T_3Jet",         ";N_{jet}",                20, 0,  20);
  h_njetAK8_S = new TH1D("njetAK8_S",         ";N_{jet}",                20, 0,  20);
  h_njetAK8_Q = new TH1D("njetAK8_Q",         ";N_{jet}",                20, 0,  20);
  h_njetAK8_W = new TH1D("njetAK8_W",         ";N_{jet}",                20, 0,  20);
  h_njetAK8_T = new TH1D("njetAK8_T",         ";N_{jet}",                20, 0,  20);

  h_njet         = new TH1D("njet",         ";N_{jet}",                20, 0,  20);
  h_nw           = new TH1D("nw",           ";N_{W tag}",              20, 0,  20);
  h_nb           = new TH1D("nb",           ";N_{b tag}",              20, 0,  20);
  h_ht_gen       = new TH1D("ht_gen",       ";H_{T}^{gen}",            nbn_HT, bn_HT);
  h_ht_AK4  = new TH1D("ht_AK4",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8  = new TH1D("ht_AK8",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);

  h_ht_AK4_S = new TH1D("ht_AK4_S",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_S = new TH1D("ht_AK8_S",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_S = new TH1D("jet1_pt_S",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_S = new TH1D("jet2_pt_S",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_S = new TH1D("jet3_pt_S",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_S = new TH1D("MR_S",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_S = new TH1D("MTR_S",  ";MTR_{AK4}",        200, 0,2000);
  h_R_S = new TH1D("R_S",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_S = new TH1D("R2_S",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_S = new TH1D("tau21_S", ";tau21", 200,0,1);
  h_MET_S = new TH1D("MET_S", ";MET", 400,0,2000);
  h_AK8_jet1_pt_S = new TH1D("AK8_jet1_pt_S",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_S = new TH1D("AK8_jet1_eta_S",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_S = new TH1D("MET_phi_S", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_S_nj35 = new TH1D("AK8_jet1_pt_S_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_S_nj35 = new TH1D("AK8_jet1_eta_S_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_S_nj35 = new TH1D("MET_phi_S_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_S_nj6 = new TH1D("AK8_jet1_pt_S_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_S_nj6 = new TH1D("AK8_jet1_eta_S_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_S_nj6 = new TH1D("MET_phi_S_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_S_nj35 = new TH1D("ht_AK4_S_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_S_nj35 = new TH1D("ht_AK8_S_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_S_nj35 = new TH1D("jet1_pt_S_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_S_nj35 = new TH1D("jet2_pt_S_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_S_nj35 = new TH1D("jet3_pt_S_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_S_nj35 = new TH1D("MR_S_nj35",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_S_nj35 = new TH1D("MTR_S_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_S_nj35 = new TH1D("R_S_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_S_nj35 = new TH1D("R2_S_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_S_nj35 = new TH1D("tau21_S_nj35", ";tau21", 200,0,1);
  h_MET_S_nj35 = new TH1D("MET_S_nj35", ";MET", 400,0,2000);

  h_ht_AK4_S_nj6 = new TH1D("ht_AK4_S_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_S_nj6 = new TH1D("ht_AK8_S_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_S_nj6 = new TH1D("jet1_pt_S_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_S_nj6 = new TH1D("jet2_pt_S_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_S_nj6 = new TH1D("jet3_pt_S_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_S_nj6 = new TH1D("MR_S_nj6",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_S_nj6 = new TH1D("MTR_S_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_S_nj6 = new TH1D("R_S_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_S_nj6 = new TH1D("R2_S_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_S_nj6 = new TH1D("tau21_S_nj6", ";tau21", 200,0,1);
  h_MET_S_nj6 = new TH1D("MET_S_nj6", ";MET", 400,0,2000);

  h_ht_AK4_s = new TH1D("ht_AK4_s",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_s = new TH1D("ht_AK8_s",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_s = new TH1D("jet1_pt_s",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_s = new TH1D("jet2_pt_s",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_s = new TH1D("jet3_pt_s",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_s = new TH1D("MR_s",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_s = new TH1D("MTR_s",  ";MTR_{AK4}",        200, 0,2000);
  h_R_s = new TH1D("R_s",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_s = new TH1D("R2_s",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_s = new TH1D("tau21_s", ";tau21", 200,0,1);
  h_MET_s = new TH1D("MET_s", ";MET", 400,0,2000);
  h_AK8_jet1_pt_s = new TH1D("AK8_jet1_pt_s",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_s = new TH1D("AK8_jet1_eta_s",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_s = new TH1D("MET_phi_s", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_s_nj35 = new TH1D("AK8_jet1_pt_s_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_s_nj35 = new TH1D("AK8_jet1_eta_s_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_s_nj35 = new TH1D("MET_phi_s_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_s_nj6 = new TH1D("AK8_jet1_pt_s_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_s_nj6 = new TH1D("AK8_jet1_eta_s_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_s_nj6 = new TH1D("MET_phi_s_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_s_nj35 = new TH1D("ht_AK4_s_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_s_nj35 = new TH1D("ht_AK8_s_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_s_nj35 = new TH1D("jet1_pt_s_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_s_nj35 = new TH1D("jet2_pt_s_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_s_nj35 = new TH1D("jet3_pt_s_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_s_nj35 = new TH1D("MR_s_nj35",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_s_nj35 = new TH1D("MTR_s_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_s_nj35 = new TH1D("R_s_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_s_nj35 = new TH1D("R2_s_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_s_nj35 = new TH1D("tau21_s_nj35", ";tau21", 200,0,1);
  h_MET_s_nj35 = new TH1D("MET_s_nj35", ";MET", 400,0,2000);

  h_ht_AK4_s_nj6 = new TH1D("ht_AK4_s_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_s_nj6 = new TH1D("ht_AK8_s_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_s_nj6 = new TH1D("jet1_pt_s_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_s_nj6 = new TH1D("jet2_pt_s_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_s_nj6 = new TH1D("jet3_pt_s_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_s_nj6 = new TH1D("MR_s_nj6",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_s_nj6 = new TH1D("MTR_s_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_s_nj6 = new TH1D("R_s_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_s_nj6 = new TH1D("R2_s_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_s_nj6 = new TH1D("tau21_s_nj6", ";tau21", 200,0,1);
  h_MET_s_nj6 = new TH1D("MET_s_nj6", ";MET", 400,0,2000);

  h_ht_AK4_Z = new TH1D("ht_AK4_Z",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_Z = new TH1D("ht_AK8_Z",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_Z = new TH1D("jet1_pt_Z",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Z = new TH1D("jet2_pt_Z",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Z = new TH1D("jet3_pt_Z",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Z = new TH1D("MR_Z",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_Z = new TH1D("MTR_Z",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Z = new TH1D("R_Z",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_Z = new TH1D("R2_Z",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_Z = new TH1D("tau21_Z", ";tau21", 200,0,1);
  h_MET_Z = new TH1D("MET_Z", ";MET", 400,0,2000);
  h_AK8_jet1_pt_Z = new TH1D("AK8_jet1_pt_Z",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_Z = new TH1D("AK8_jet1_eta_Z",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_Z = new TH1D("MET_phi_Z", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_Z_nj35 = new TH1D("AK8_jet1_pt_Z_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_Z_nj35 = new TH1D("AK8_jet1_eta_Z_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_Z_nj35 = new TH1D("MET_phi_Z_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_Z_nj6 = new TH1D("AK8_jet1_pt_Z_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_Z_nj6 = new TH1D("AK8_jet1_eta_Z_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_Z_nj6 = new TH1D("MET_phi_Z_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_Z_nj35 = new TH1D("ht_AK4_Z_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_Z_nj35 = new TH1D("ht_AK8_Z_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_Z_nj35 = new TH1D("jet1_pt_Z_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Z_nj35 = new TH1D("jet2_pt_Z_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Z_nj35 = new TH1D("jet3_pt_Z_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Z_nj35 = new TH1D("MR_Z_nj35",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_Z_nj35 = new TH1D("MTR_Z_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Z_nj35 = new TH1D("R_Z_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_Z_nj35 = new TH1D("R2_Z_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_Z_nj35 = new TH1D("tau21_Z_nj35", ";tau21", 200,0,1);
  h_MET_Z_nj35 = new TH1D("MET_Z_nj35", ";MET", 400,0,2000);

  h_ht_AK4_Z_nj6 = new TH1D("ht_AK4_Z_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_Z_nj6 = new TH1D("ht_AK8_Z_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_Z_nj6 = new TH1D("jet1_pt_Z_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Z_nj6 = new TH1D("jet2_pt_Z_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Z_nj6 = new TH1D("jet3_pt_Z_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Z_nj6 = new TH1D("MR_Z_nj6",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_Z_nj6 = new TH1D("MTR_Z_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Z_nj6 = new TH1D("R_Z_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_Z_nj6 = new TH1D("R2_Z_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_Z_nj6 = new TH1D("tau21_Z_nj6", ";tau21", 200,0,1);
  h_MET_Z_nj6 = new TH1D("MET_Z_nj6", ";MET", 400,0,2000);

  h_ht_AK4_G = new TH1D("ht_AK4_G",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_G = new TH1D("ht_AK8_G",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_G = new TH1D("jet1_pt_G",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_G = new TH1D("jet2_pt_G",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_G = new TH1D("jet3_pt_G",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_G = new TH1D("MR_G",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_G = new TH1D("MTR_G",  ";MTR_{AK4}",        200, 0,2000);
  h_R_G = new TH1D("R_G",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_G = new TH1D("R2_G",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_G = new TH1D("tau21_G", ";tau21", 200,0,1);
  h_MET_G = new TH1D("MET_G", ";MET", 400,0,2000);
  h_AK8_jet1_pt_G = new TH1D("AK8_jet1_pt_G",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_G = new TH1D("AK8_jet1_eta_G",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_G = new TH1D("MET_phi_G", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_G_nj35 = new TH1D("AK8_jet1_pt_G_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_G_nj35 = new TH1D("AK8_jet1_eta_G_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_G_nj35 = new TH1D("MET_phi_G_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_G_nj6 = new TH1D("AK8_jet1_pt_G_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_G_nj6 = new TH1D("AK8_jet1_eta_G_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_G_nj6 = new TH1D("MET_phi_G_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_G_nj35 = new TH1D("ht_AK4_G_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_G_nj35 = new TH1D("ht_AK8_G_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_G_nj35 = new TH1D("jet1_pt_G_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_G_nj35 = new TH1D("jet2_pt_G_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_G_nj35 = new TH1D("jet3_pt_G_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_G_nj35 = new TH1D("MR_G_nj35",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_G_nj35 = new TH1D("MTR_G_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_G_nj35 = new TH1D("R_G_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_G_nj35 = new TH1D("R2_G_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_G_nj35 = new TH1D("tau21_G_nj35", ";tau21", 200,0,1);
  h_MET_G_nj35 = new TH1D("MET_G_nj35", ";MET", 400,0,2000);

  h_ht_AK4_G_nj6 = new TH1D("ht_AK4_G_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_G_nj6 = new TH1D("ht_AK8_G_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_G_nj6 = new TH1D("jet1_pt_G_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_G_nj6 = new TH1D("jet2_pt_G_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_G_nj6 = new TH1D("jet3_pt_G_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_G_nj6 = new TH1D("MR_G_nj6",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_G_nj6 = new TH1D("MTR_G_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_G_nj6 = new TH1D("R_G_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_G_nj6 = new TH1D("R2_G_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_G_nj6 = new TH1D("tau21_G_nj6", ";tau21", 200,0,1);
  h_MET_G_nj6 = new TH1D("MET_G_nj6", ";MET", 400,0,2000);

  //h_MR_S_nj6 = new TH1D("MR_S_nj6",";MR_{AK4}",nbn_MR,bn_MR);
  //h_MR_S_nj35 = new TH1D("MR_S_nj35",";MR_{AK4}",nbn_MR,bn_MR);
  //h_R2_S_nj6 = new TH1D("R2_S_nj6",";R2_{AK4}",nbn_R2,bn_R2);
  //h_R2_S_nj35 = new TH1D("R2_S_nj35",";R2_{AK4}",nbn_R2,bn_R2);

  //h_MR_Q_nj6 = new TH1D("MR_Q_nj6",";MR_{AK4}",nbn_MR,bn_MR);
  //h_MR_Q_nj35 = new TH1D("MR_Q_nj35",";MR_{AK4}",nbn_MR,bn_MR);
  //h_R2_Q_nj35 = new TH1D("R2_Q_nj35",";R2_{AK4}",nbn_R2,bn_R2);
  //h_R2_Q_nj6 = new TH1D("R2_Q_nj6",";R2_{AK4}",nbn_R2,bn_R2);

  //h_MR_T_nj6 = new TH1D("MR_T_nj6",";MR_{AK4}",nbn_MR,bn_MR);
  //h_MR_T_nj35 = new TH1D("MR_T_nj35",";MR_{AK4}",nbn_MR,bn_MR);
  //h_R2_T_nj35 = new TH1D("R2_T_nj35",";R2_{AK4}",nbn_R2,bn_R2);
  //h_R2_T_nj6 = new TH1D("R2_T_nj6",";R2_{AK4}",nbn_R2,bn_R2);

  //h_MR_W_nj6 = new TH1D("MR_W_nj6",";MR_{AK4}",nbn_MR,bn_MR);
  //h_MR_W_nj35 = new TH1D("MR_W_nj35",";MR_{AK4}",nbn_MR,bn_MR);
  //h_R2_W_nj35 = new TH1D("R2_W_nj35",";R2_{AK4}",nbn_R2,bn_R2);
  //h_R2_W_nj6 = new TH1D("R2_W_nj6",";R2_{AK4}",nbn_R2,bn_R2);

  h_ht_AK4_Q = new TH1D("ht_AK4_Q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_Q = new TH1D("ht_AK8_Q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_Q = new TH1D("jet1_pt_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Q = new TH1D("jet2_pt_Q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Q = new TH1D("jet3_pt_Q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Q = new TH1D("MR_Q",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_Q = new TH1D("MTR_Q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Q = new TH1D("R_Q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_Q = new TH1D("R2_Q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_Q = new TH1D("tau21_Q", ";tau21", 200,0,1);
  h_MET_Q = new TH1D("MET_Q", ";MET", 400,0,2000);
  h_AK8_jet1_pt_Q = new TH1D("AK8_jet1_pt_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_Q = new TH1D("AK8_jet1_eta_Q",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_Q = new TH1D("MET_phi_Q", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_Q_nj35 = new TH1D("AK8_jet1_pt_Q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_Q_nj35 = new TH1D("AK8_jet1_eta_Q_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_Q_nj35 = new TH1D("MET_phi_Q_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_Q_nj6 = new TH1D("AK8_jet1_pt_Q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_Q_nj6 = new TH1D("AK8_jet1_eta_Q_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_Q_nj6 = new TH1D("MET_phi_Q_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_nodPhi_Q = new TH1D("ht_AK4_nodPhi_Q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nodPhi_Q = new TH1D("ht_AK8_nodPhi_Q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nodPhi_Q = new TH1D("jet1_pt_nodPhi_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nodPhi_Q = new TH1D("jet2_pt_nodPhi_Q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nodPhi_Q = new TH1D("jet3_pt_nodPhi_Q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nodPhi_Q = new TH1D("MR_nodPhi_Q",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nodPhi_Q = new TH1D("MTR_nodPhi_Q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nodPhi_Q = new TH1D("R_nodPhi_Q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nodPhi_Q = new TH1D("R2_nodPhi_Q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nodPhi_Q = new TH1D("tau21_nodPhi_Q", ";tau21", 200,0,1);
  h_MET_nodPhi_Q = new TH1D("MET_nodPhi_Q", ";MET", 400,0,2000);
  h_dPhiRazor_nodPhi_Q = new TH1D("dPhiRazor_nodPhi_Q", ";#Delta#phi_{megajets}", 64,0,3.2);

  h_ht_AK4_notau21_Q = new TH1D("ht_AK4_notau21_Q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_notau21_Q = new TH1D("ht_AK8_notau21_Q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_notau21_Q = new TH1D("jet1_pt_notau21_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_notau21_Q = new TH1D("jet2_pt_notau21_Q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_notau21_Q = new TH1D("jet3_pt_notau21_Q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_notau21_Q = new TH1D("MR_notau21_Q",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_notau21_Q = new TH1D("MTR_notau21_Q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_notau21_Q = new TH1D("R_notau21_Q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_notau21_Q = new TH1D("R2_notau21_Q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_notau21_Q = new TH1D("tau21_notau21_Q", ";tau21", 200,0,1);
  h_MET_notau21_Q = new TH1D("MET_notau21_Q", ";MET", 400,0,2000);

  h_ht_AK4_nonb_Q = new TH1D("ht_AK4_nonb_Q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nonb_Q = new TH1D("ht_AK8_nonb_Q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nonb_Q = new TH1D("jet1_pt_nonb_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nonb_Q = new TH1D("jet2_pt_nonb_Q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nonb_Q = new TH1D("jet3_pt_nonb_Q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nonb_Q = new TH1D("MR_nonb_Q",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nonb_Q = new TH1D("MTR_nonb_Q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nonb_Q = new TH1D("R_nonb_Q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nonb_Q = new TH1D("R2_nonb_Q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nonb_Q = new TH1D("tau21_nonb_Q", ";tau21", 200,0,1);
  h_MET_nonb_Q = new TH1D("MET_nonb_Q", ";MET", 400,0,2000);

  h_ht_AK4_Q_nj35 = new TH1D("ht_AK4_Q_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_Q_nj35 = new TH1D("ht_AK8_Q_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_Q_nj35 = new TH1D("jet1_pt_Q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Q_nj35 = new TH1D("jet2_pt_Q_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Q_nj35 = new TH1D("jet3_pt_Q_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Q_nj35 = new TH1D("MR_Q_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_Q_nj35 = new TH1D("MTR_Q_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Q_nj35 = new TH1D("R_Q_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_Q_nj35 = new TH1D("R2_Q_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_Q_nj35 = new TH1D("tau21_Q_nj35", ";tau21", 200,0,1);
  h_MET_Q_nj35 = new TH1D("MET_Q_nj35", ";MET", 400,0,2000);

  h_ht_AK4_nodPhi_Q_nj35 = new TH1D("ht_AK4_nodPhi_Q_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nodPhi_Q_nj35 = new TH1D("ht_AK8_nodPhi_Q_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nodPhi_Q_nj35 = new TH1D("jet1_pt_nodPhi_Q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nodPhi_Q_nj35 = new TH1D("jet2_pt_nodPhi_Q_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nodPhi_Q_nj35 = new TH1D("jet3_pt_nodPhi_Q_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nodPhi_Q_nj35 = new TH1D("MR_nodPhi_Q_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nodPhi_Q_nj35 = new TH1D("MTR_nodPhi_Q_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nodPhi_Q_nj35 = new TH1D("R_nodPhi_Q_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nodPhi_Q_nj35 = new TH1D("R2_nodPhi_Q_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nodPhi_Q_nj35 = new TH1D("tau21_nodPhi_Q_nj35", ";tau21", 200,0,1);
  h_MET_nodPhi_Q_nj35 = new TH1D("MET_nodPhi_Q_nj35", ";MET", 400,0,2000);
  h_dPhiRazor_nodPhi_Q_nj35 = new TH1D("dPhiRazor_nodPhi_Q_nj35", ";#Delta#phi_{megajets}", 64,0,3.2);

  h_ht_AK4_notau21_Q_nj35 = new TH1D("ht_AK4_notau21_Q_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_notau21_Q_nj35 = new TH1D("ht_AK8_notau21_Q_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_notau21_Q_nj35 = new TH1D("jet1_pt_notau21_Q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_notau21_Q_nj35 = new TH1D("jet2_pt_notau21_Q_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_notau21_Q_nj35 = new TH1D("jet3_pt_notau21_Q_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_notau21_Q_nj35 = new TH1D("MR_notau21_Q_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_notau21_Q_nj35 = new TH1D("MTR_notau21_Q_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_notau21_Q_nj35 = new TH1D("R_notau21_Q_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_notau21_Q_nj35 = new TH1D("R2_notau21_Q_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_notau21_Q_nj35 = new TH1D("tau21_notau21_Q_nj35", ";tau21", 200,0,1);
  h_MET_notau21_Q_nj35 = new TH1D("MET_notau21_Q_nj35", ";MET", 400,0,2000);

  h_ht_AK4_nonb_Q_nj35 = new TH1D("ht_AK4_nonb_Q_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nonb_Q_nj35 = new TH1D("ht_AK8_nonb_Q_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nonb_Q_nj35 = new TH1D("jet1_pt_nonb_Q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nonb_Q_nj35 = new TH1D("jet2_pt_nonb_Q_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nonb_Q_nj35 = new TH1D("jet3_pt_nonb_Q_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nonb_Q_nj35 = new TH1D("MR_nonb_Q_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nonb_Q_nj35 = new TH1D("MTR_nonb_Q_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nonb_Q_nj35 = new TH1D("R_nonb_Q_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nonb_Q_nj35 = new TH1D("R2_nonb_Q_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nonb_Q_nj35 = new TH1D("tau21_nonb_Q_nj35", ";tau21", 200,0,1);
  h_MET_nonb_Q_nj35 = new TH1D("MET_nonb_Q_nj35", ";MET", 400,0,2000);

  h_ht_AK4_Q_nj6 = new TH1D("ht_AK4_Q_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_Q_nj6 = new TH1D("ht_AK8_Q_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_Q_nj6 = new TH1D("jet1_pt_Q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_Q_nj6 = new TH1D("jet2_pt_Q_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_Q_nj6 = new TH1D("jet3_pt_Q_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_Q_nj6 = new TH1D("MR_Q_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_Q_nj6 = new TH1D("MTR_Q_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_Q_nj6 = new TH1D("R_Q_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_Q_nj6 = new TH1D("R2_Q_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_Q_nj6 = new TH1D("tau21_Q_nj6", ";tau21", 200,0,1);
  h_MET_Q_nj6 = new TH1D("MET_Q_nj6", ";MET", 400,0,2000);

  h_ht_AK4_nodPhi_Q_nj6 = new TH1D("ht_AK4_nodPhi_Q_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nodPhi_Q_nj6 = new TH1D("ht_AK8_nodPhi_Q_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nodPhi_Q_nj6 = new TH1D("jet1_pt_nodPhi_Q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nodPhi_Q_nj6 = new TH1D("jet2_pt_nodPhi_Q_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nodPhi_Q_nj6 = new TH1D("jet3_pt_nodPhi_Q_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nodPhi_Q_nj6 = new TH1D("MR_nodPhi_Q_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nodPhi_Q_nj6 = new TH1D("MTR_nodPhi_Q_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nodPhi_Q_nj6 = new TH1D("R_nodPhi_Q_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nodPhi_Q_nj6 = new TH1D("R2_nodPhi_Q_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nodPhi_Q_nj6 = new TH1D("tau21_nodPhi_Q_nj6", ";tau21", 200,0,1);
  h_MET_nodPhi_Q_nj6 = new TH1D("MET_nodPhi_Q_nj6", ";MET", 400,0,2000);
  h_dPhiRazor_nodPhi_Q_nj6 = new TH1D("dPhiRazor_nodPhi_Q_nj6", ";#Delta#phi_{megajets}", 64,0,3.2);

  h_ht_AK4_notau21_Q_nj6 = new TH1D("ht_AK4_notau21_Q_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_notau21_Q_nj6 = new TH1D("ht_AK8_notau21_Q_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_notau21_Q_nj6 = new TH1D("jet1_pt_notau21_Q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_notau21_Q_nj6 = new TH1D("jet2_pt_notau21_Q_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_notau21_Q_nj6 = new TH1D("jet3_pt_notau21_Q_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_notau21_Q_nj6 = new TH1D("MR_notau21_Q_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_notau21_Q_nj6 = new TH1D("MTR_notau21_Q_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_notau21_Q_nj6 = new TH1D("R_notau21_Q_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_notau21_Q_nj6 = new TH1D("R2_notau21_Q_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_notau21_Q_nj6 = new TH1D("tau21_notau21_Q_nj6", ";tau21", 200,0,1);
  h_MET_notau21_Q_nj6 = new TH1D("MET_notau21_Q_nj6", ";MET", 400,0,2000);

  h_ht_AK4_nonb_Q_nj6 = new TH1D("ht_AK4_nonb_Q_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nonb_Q_nj6 = new TH1D("ht_AK8_nonb_Q_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nonb_Q_nj6 = new TH1D("jet1_pt_nonb_Q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nonb_Q_nj6 = new TH1D("jet2_pt_nonb_Q_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nonb_Q_nj6 = new TH1D("jet3_pt_nonb_Q_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nonb_Q_nj6 = new TH1D("MR_nonb_Q_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nonb_Q_nj6 = new TH1D("MTR_nonb_Q_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nonb_Q_nj6 = new TH1D("R_nonb_Q_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nonb_Q_nj6 = new TH1D("R2_nonb_Q_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nonb_Q_nj6 = new TH1D("tau21_nonb_Q_nj6", ";tau21", 200,0,1);
  h_MET_nonb_Q_nj6 = new TH1D("MET_nonb_Q_nj6", ";MET", 400,0,2000);

  h_ht_AK4_T = new TH1D("ht_AK4_T",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_T = new TH1D("ht_AK8_T",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_T = new TH1D("jet1_pt_T",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_T = new TH1D("jet2_pt_T",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_T = new TH1D("jet3_pt_T",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_T = new TH1D("MR_T",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_T = new TH1D("MTR_T",  ";MTR_{AK4}",        200, 0,2000);
  h_R_T = new TH1D("R_T",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_T = new TH1D("R2_T",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_T = new TH1D("tau21_T", ";tau21", 200,0,1);
  h_MET_T = new TH1D("MET_T", ";MET", 400,0,2000);
  h_AK8_jet1_pt_T = new TH1D("AK8_jet1_pt_T",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_T = new TH1D("AK8_jet1_eta_T",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_T = new TH1D("MET_phi_T", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_T_nj35 = new TH1D("AK8_jet1_pt_T_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_T_nj35 = new TH1D("AK8_jet1_eta_T_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_T_nj35 = new TH1D("MET_phi_T_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_T_nj6 = new TH1D("AK8_jet1_pt_T_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_T_nj6 = new TH1D("AK8_jet1_eta_T_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_T_nj6 = new TH1D("MET_phi_T_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_noMT_T = new TH1D("ht_AK4_noMT_T",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_noMT_T = new TH1D("ht_AK8_noMT_T",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_noMT_T = new TH1D("jet1_pt_noMT_T",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_noMT_T = new TH1D("jet2_pt_noMT_T",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_noMT_T = new TH1D("jet3_pt_noMT_T",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_noMT_T = new TH1D("MR_noMT_T",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_noMT_T = new TH1D("MTR_noMT_T",  ";MTR_{AK4}",        200, 0,2000);
  h_R_noMT_T = new TH1D("R_noMT_T",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_noMT_T = new TH1D("R2_noMT_T",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_noMT_T = new TH1D("tau21_noMT_T", ";tau21", 200,0,1);
  h_MET_noMT_T = new TH1D("MET_noMT_T", ";MET", 400,0,2000);

  h_ht_AK4_T_nj35 = new TH1D("ht_AK4_T_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_T_nj35 = new TH1D("ht_AK8_T_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_T_nj35 = new TH1D("jet1_pt_T_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_T_nj35 = new TH1D("jet2_pt_T_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_T_nj35 = new TH1D("jet3_pt_T_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_T_nj35 = new TH1D("MR_T_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_T_nj35 = new TH1D("MTR_T_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_T_nj35 = new TH1D("R_T_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_T_nj35 = new TH1D("R2_T_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_T_nj35 = new TH1D("tau21_T_nj35", ";tau21", 200,0,1);
  h_MET_T_nj35 = new TH1D("MET_T_nj35", ";MET", 400,0,2000);

  h_ht_AK4_noMT_T_nj35 = new TH1D("ht_AK4_noMT_T_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_noMT_T_nj35 = new TH1D("ht_AK8_noMT_T_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_noMT_T_nj35 = new TH1D("jet1_pt_noMT_T_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_noMT_T_nj35 = new TH1D("jet2_pt_noMT_T_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_noMT_T_nj35 = new TH1D("jet3_pt_noMT_T_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_noMT_T_nj35 = new TH1D("MR_noMT_T_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_noMT_T_nj35 = new TH1D("MTR_noMT_T_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_noMT_T_nj35 = new TH1D("R_noMT_T_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_noMT_T_nj35 = new TH1D("R2_noMT_T_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_noMT_T_nj35 = new TH1D("tau21_noMT_T_nj35", ";tau21", 200,0,1);
  h_MET_noMT_T_nj35 = new TH1D("MET_noMT_T_nj35", ";MET", 400,0,2000);

  h_ht_AK4_T_nj6 = new TH1D("ht_AK4_T_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_T_nj6 = new TH1D("ht_AK8_T_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_T_nj6 = new TH1D("jet1_pt_T_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_T_nj6 = new TH1D("jet2_pt_T_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_T_nj6 = new TH1D("jet3_pt_T_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_T_nj6 = new TH1D("MR_T_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_T_nj6 = new TH1D("MTR_T_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_T_nj6 = new TH1D("R_T_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_T_nj6 = new TH1D("R2_T_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_T_nj6 = new TH1D("tau21_T_nj6", ";tau21", 200,0,1);
  h_MET_T_nj6 = new TH1D("MET_T_nj6", ";MET", 400,0,2000);

  h_ht_AK4_noMT_T_nj6 = new TH1D("ht_AK4_noMT_T_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_noMT_T_nj6 = new TH1D("ht_AK8_noMT_T_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_noMT_T_nj6 = new TH1D("jet1_pt_noMT_T_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_noMT_T_nj6 = new TH1D("jet2_pt_noMT_T_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_noMT_T_nj6 = new TH1D("jet3_pt_noMT_T_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_noMT_T_nj6 = new TH1D("MR_noMT_T_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_noMT_T_nj6 = new TH1D("MTR_noMT_T_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_noMT_T_nj6 = new TH1D("R_noMT_T_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_noMT_T_nj6 = new TH1D("R2_noMT_T_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_noMT_T_nj6 = new TH1D("tau21_noMT_T_nj6", ";tau21", 200,0,1);
  h_MET_noMT_T_nj6 = new TH1D("MET_noMT_T_nj6", ";MET", 400,0,2000);

  h_ht_AK4_W = new TH1D("ht_AK4_W",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_W = new TH1D("ht_AK8_W",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_W = new TH1D("jet1_pt_W",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_W = new TH1D("jet2_pt_W",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_W = new TH1D("jet3_pt_W",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_W = new TH1D("MR_W",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_W = new TH1D("MTR_W",  ";MTR_{AK4}",        200, 0,2000);
  h_R_W = new TH1D("R_W",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_W = new TH1D("R2_W",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_W = new TH1D("tau21_W", ";tau21", 200,0,1);
  h_MET_W = new TH1D("MET_W", ";MET", 400,0,2000);
  h_AK8_jet1_pt_W = new TH1D("AK8_jet1_pt_W",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_W = new TH1D("AK8_jet1_eta_W",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_W = new TH1D("MET_phi_W", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_W_nj35 = new TH1D("AK8_jet1_pt_W_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_W_nj35 = new TH1D("AK8_jet1_eta_W_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_W_nj35 = new TH1D("MET_phi_W_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_W_nj6 = new TH1D("AK8_jet1_pt_W_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_W_nj6 = new TH1D("AK8_jet1_eta_W_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_W_nj6 = new TH1D("MET_phi_W_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_noMT_W = new TH1D("ht_AK4_noMT_W",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_noMT_W = new TH1D("ht_AK8_noMT_W",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_noMT_W = new TH1D("jet1_pt_noMT_W",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_noMT_W = new TH1D("jet2_pt_noMT_W",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_noMT_W = new TH1D("jet3_pt_noMT_W",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_noMT_W = new TH1D("MR_noMT_W",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_noMT_W = new TH1D("MTR_noMT_W",  ";MTR_{AK4}",        200, 0,2000);
  h_R_noMT_W = new TH1D("R_noMT_W",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_noMT_W = new TH1D("R2_noMT_W",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_noMT_W = new TH1D("tau21_noMT_W", ";tau21", 200,0,1);
  h_MET_noMT_W = new TH1D("MET_noMT_W", ";MET", 400,0,2000);

  h_ht_AK4_nonb_W = new TH1D("ht_AK4_nonb_W",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nonb_W = new TH1D("ht_AK8_nonb_W",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nonb_W = new TH1D("jet1_pt_nonb_W",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nonb_W = new TH1D("jet2_pt_nonb_W",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nonb_W = new TH1D("jet3_pt_nonb_W",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nonb_W = new TH1D("MR_nonb_W",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nonb_W = new TH1D("MTR_nonb_W",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nonb_W = new TH1D("R_nonb_W",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nonb_W = new TH1D("R2_nonb_W",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nonb_W = new TH1D("tau21_nonb_W", ";tau21", 200,0,1);
  h_MET_nonb_W = new TH1D("MET_nonb_W", ";MET", 400,0,2000);

  h_ht_AK4_W_nj35 = new TH1D("ht_AK4_W_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_W_nj35 = new TH1D("ht_AK8_W_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_W_nj35 = new TH1D("jet1_pt_W_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_W_nj35 = new TH1D("jet2_pt_W_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_W_nj35 = new TH1D("jet3_pt_W_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_W_nj35 = new TH1D("MR_W_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_W_nj35 = new TH1D("MTR_W_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_W_nj35 = new TH1D("R_W_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_W_nj35 = new TH1D("R2_W_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_W_nj35 = new TH1D("tau21_W_nj35", ";tau21", 200,0,1);
  h_MET_W_nj35 = new TH1D("MET_W_nj35", ";MET", 400,0,2000);

  h_ht_AK4_noMT_W_nj35 = new TH1D("ht_AK4_noMT_W_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_noMT_W_nj35 = new TH1D("ht_AK8_noMT_W_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_noMT_W_nj35 = new TH1D("jet1_pt_noMT_W_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_noMT_W_nj35 = new TH1D("jet2_pt_noMT_W_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_noMT_W_nj35 = new TH1D("jet3_pt_noMT_W_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_noMT_W_nj35 = new TH1D("MR_noMT_W_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_noMT_W_nj35 = new TH1D("MTR_noMT_W_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_noMT_W_nj35 = new TH1D("R_noMT_W_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_noMT_W_nj35 = new TH1D("R2_noMT_W_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_noMT_W_nj35 = new TH1D("tau21_noMT_W_nj35", ";tau21", 200,0,1);
  h_MET_noMT_W_nj35 = new TH1D("MET_noMT_W_nj35", ";MET", 400,0,2000);

  h_ht_AK4_nonb_W_nj35 = new TH1D("ht_AK4_nonb_W_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nonb_W_nj35 = new TH1D("ht_AK8_nonb_W_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nonb_W_nj35 = new TH1D("jet1_pt_nonb_W_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nonb_W_nj35 = new TH1D("jet2_pt_nonb_W_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nonb_W_nj35 = new TH1D("jet3_pt_nonb_W_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nonb_W_nj35 = new TH1D("MR_nonb_W_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nonb_W_nj35 = new TH1D("MTR_nonb_W_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nonb_W_nj35 = new TH1D("R_nonb_W_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nonb_W_nj35 = new TH1D("R2_nonb_W_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nonb_W_nj35 = new TH1D("tau21_nonb_W_nj35", ";tau21", 200,0,1);
  h_MET_nonb_W_nj35 = new TH1D("MET_nonb_W_nj35", ";MET", 400,0,2000);

  h_ht_AK4_W_nj6 = new TH1D("ht_AK4_W_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_W_nj6 = new TH1D("ht_AK8_W_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_W_nj6 = new TH1D("jet1_pt_W_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_W_nj6 = new TH1D("jet2_pt_W_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_W_nj6 = new TH1D("jet3_pt_W_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_W_nj6 = new TH1D("MR_W_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_W_nj6 = new TH1D("MTR_W_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_W_nj6 = new TH1D("R_W_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_W_nj6 = new TH1D("R2_W_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_W_nj6 = new TH1D("tau21_W_nj6", ";tau21", 200,0,1);
  h_MET_W_nj6 = new TH1D("MET_W_nj6", ";MET", 400,0,2000);

  h_ht_AK4_noMT_W_nj6 = new TH1D("ht_AK4_noMT_W_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_noMT_W_nj6 = new TH1D("ht_AK8_noMT_W_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_noMT_W_nj6 = new TH1D("jet1_pt_noMT_W_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_noMT_W_nj6 = new TH1D("jet2_pt_noMT_W_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_noMT_W_nj6 = new TH1D("jet3_pt_noMT_W_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_noMT_W_nj6 = new TH1D("MR_noMT_W_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_noMT_W_nj6 = new TH1D("MTR_noMT_W_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_noMT_W_nj6 = new TH1D("R_noMT_W_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_noMT_W_nj6 = new TH1D("R2_noMT_W_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_noMT_W_nj6 = new TH1D("tau21_noMT_W_nj6", ";tau21", 200,0,1);
  h_MET_noMT_W_nj6 = new TH1D("MET_noMT_W_nj6", ";MET", 400,0,2000);

  h_ht_AK4_nonb_W_nj6 = new TH1D("ht_AK4_nonb_W_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_nonb_W_nj6 = new TH1D("ht_AK8_nonb_W_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_nonb_W_nj6 = new TH1D("jet1_pt_nonb_W_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_nonb_W_nj6 = new TH1D("jet2_pt_nonb_W_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_nonb_W_nj6 = new TH1D("jet3_pt_nonb_W_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_nonb_W_nj6 = new TH1D("MR_nonb_W_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_nonb_W_nj6 = new TH1D("MTR_nonb_W_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_nonb_W_nj6 = new TH1D("R_nonb_W_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_nonb_W_nj6 = new TH1D("R2_nonb_W_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_nonb_W_nj6 = new TH1D("tau21_nonb_W_nj6", ";tau21", 200,0,1);
  h_MET_nonb_W_nj6 = new TH1D("MET_nonb_W_nj6", ";MET", 400,0,2000);

  h_ht_AK4_q = new TH1D("ht_AK4_q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_q = new TH1D("ht_AK8_q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_q = new TH1D("jet1_pt_q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_q = new TH1D("jet2_pt_q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_q = new TH1D("jet3_pt_q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_q = new TH1D("MR_q",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_q = new TH1D("MTR_q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_q = new TH1D("R_q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_q = new TH1D("R2_q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_q = new TH1D("tau21_q", ";tau21", 200,0,1);
  h_MET_q = new TH1D("MET_q", ";MET", 400,0,2000);
  h_AK8_jet1_pt_q = new TH1D("AK8_jet1_pt_q",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_q = new TH1D("AK8_jet1_eta_q",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_q = new TH1D("MET_phi_q", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_q_nj35 = new TH1D("AK8_jet1_pt_q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_q_nj35 = new TH1D("AK8_jet1_eta_q_nj35",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_q_nj35 = new TH1D("MET_phi_q_nj35", ";#phi_{MET}", 640, -3.2,3.2);
  h_AK8_jet1_pt_q_nj6 = new TH1D("AK8_jet1_pt_q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_q_nj6 = new TH1D("AK8_jet1_eta_q_nj6",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_q_nj6 = new TH1D("MET_phi_q_nj6", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_q_nj35 = new TH1D("ht_AK4_q_nj35",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_q_nj35 = new TH1D("ht_AK8_q_nj35",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_q_nj35 = new TH1D("jet1_pt_q_nj35",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_q_nj35 = new TH1D("jet2_pt_q_nj35",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_q_nj35 = new TH1D("jet3_pt_q_nj35",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_q_nj35 = new TH1D("MR_q_nj35",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_q_nj35 = new TH1D("MTR_q_nj35",  ";MTR_{AK4}",        200, 0,2000);
  h_R_q_nj35 = new TH1D("R_q_nj35",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_q_nj35 = new TH1D("R2_q_nj35",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_q_nj35 = new TH1D("tau21_q_nj35", ";tau21", 200,0,1);
  h_MET_q_nj35 = new TH1D("MET_q_nj35", ";MET", 400,0,2000);

  h_ht_AK4_q_nj6 = new TH1D("ht_AK4_q_nj6",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_q_nj6 = new TH1D("ht_AK8_q_nj6",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_q_nj6 = new TH1D("jet1_pt_q_nj6",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_q_nj6 = new TH1D("jet2_pt_q_nj6",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_q_nj6 = new TH1D("jet3_pt_q_nj6",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_q_nj6 = new TH1D("MR_q_nj6",   ";MR_{AK4}",         nbn_MR,bn_MR);
  h_MTR_q_nj6 = new TH1D("MTR_q_nj6",  ";MTR_{AK4}",        200, 0,2000);
  h_R_q_nj6 = new TH1D("R_q_nj6",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_q_nj6 = new TH1D("R2_q_nj6",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_q_nj6 = new TH1D("tau21_q_nj6", ";tau21", 200,0,1);
  h_MET_q_nj6 = new TH1D("MET_q_nj6", ";MET", 400,0,2000);


  h_softDropMass = new TH1D("softDropMass", "", 100,0,500);
  h_GluinoLSPMass = new TH2D("GluinoLSPMass","",34,600,2300,64,0,1600);
  h_StopMass = new TH1D("StopMass", "", 92,0,2300);
  h_GluinoMass = new TH1D("GluinoMass", "", 34,600,2300);
  h_LSPMass = new TH1D("LSPMass", "", 64,0,1600);

  h_R2_MR = new TH2D("R2_MR", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);

  h_R2_MR_s = new TH2D("R2_MR_s", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_s = new TH2D("R2_MET_s", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_s = new TH2D("MR_MET_s", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_s = new TH2D("AK8Jet1pT_MET_s", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_q = new TH2D("R2_MR_q", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_q = new TH2D("R2_MET_q", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_q = new TH2D("MR_MET_q", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_q = new TH2D("AK8Jet1pT_MET_q", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_S = new TH2D("R2_MR_S", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_S = new TH2D("R2_MET_S", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_S = new TH2D("MR_MET_S", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_S = new TH2D("AK8Jet1pT_MET_S", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_Q = new TH2D("R2_MR_Q", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_Q = new TH2D("R2_MET_Q", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_Q = new TH2D("MR_MET_Q", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_Q = new TH2D("AK8Jet1pT_MET_Q", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_T = new TH2D("R2_MR_T", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_T = new TH2D("R2_MET_T", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_T = new TH2D("MR_MET_T", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_T = new TH2D("AK8Jet1pT_MET_T", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_W = new TH2D("R2_MR_W", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_W = new TH2D("R2_MET_W", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_W = new TH2D("MR_MET_W", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_W = new TH2D("AK8Jet1pT_MET_W", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_Z = new TH2D("R2_MR_Z", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_Z = new TH2D("R2_MET_Z", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_Z = new TH2D("MR_MET_Z", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_Z = new TH2D("AK8Jet1pT_MET_Z", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_G = new TH2D("R2_MR_G", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MET_G = new TH2D("R2_MET_G", ";MET_{AK4};R2_{AK4}",40,0,2000,nbn_R2,bn_R2);
  h_MR_MET_G = new TH2D("MR_MET_G", ";MET_{AK4};MR_{AK4}",40,0,2000,nbn_MR,bn_MR);
  h_AK8Jet1pT_MET_G = new TH2D("AK8Jet1pT_MET_G", ";MET_{AK4};p_{T,AK8Jet1}",40,0,2000,40,0,2000);
  h_R2_MR_s_nj35 = new TH2D("R2_MR_s_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_s_nj6 = new TH2D("R2_MR_s_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_q_nj35 = new TH2D("R2_MR_q_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_q_nj6 = new TH2D("R2_MR_q_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_S_nj35 = new TH2D("R2_MR_S_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_S_nj6 = new TH2D("R2_MR_S_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_Q_nj35 = new TH2D("R2_MR_Q_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_Q_nj6 = new TH2D("R2_MR_Q_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_T_nj35 = new TH2D("R2_MR_T_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_T_nj6 = new TH2D("R2_MR_T_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_W_nj35 = new TH2D("R2_MR_W_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_W_nj6 = new TH2D("R2_MR_W_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_Z_nj35 = new TH2D("R2_MR_Z_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_Z_nj6 = new TH2D("R2_MR_Z_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_G_nj35 = new TH2D("R2_MR_G_nj35", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_R2_MR_G_nj6 = new TH2D("R2_MR_G_nj6", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  
  h_ht_AK4_topana_S = new TH1D("ht_AK4_topana_S",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_S = new TH1D("ht_AK8_topana_S",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_S = new TH1D("jet1_pt_topana_S",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_S = new TH1D("jet2_pt_topana_S",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_S = new TH1D("jet3_pt_topana_S",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_S = new TH1D("MR_topana_S",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_S = new TH1D("MTR_topana_S",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_S = new TH1D("R_topana_S",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_S = new TH1D("R2_topana_S",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_S = new TH1D("tau21_topana_S", ";tau21", 200,0,1);
  h_MET_topana_S = new TH1D("MET_topana_S", ";MET", 400,0,2000);
  h_R2_MR_topana_S = new TH2D("R2_MR_topana_S", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_S = new TH1D("AK8_jet1_pt_topana_S",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_S = new TH1D("AK8_jet1_eta_topana_S",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_S = new TH1D("MET_phi_topana_S", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_topana_s = new TH1D("ht_AK4_topana_s",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_s = new TH1D("ht_AK8_topana_s",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_s = new TH1D("jet1_pt_topana_s",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_s = new TH1D("jet2_pt_topana_s",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_s = new TH1D("jet3_pt_topana_s",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_s = new TH1D("MR_topana_s",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_s = new TH1D("MTR_topana_s",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_s = new TH1D("R_topana_s",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_s = new TH1D("R2_topana_s",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_s = new TH1D("tau21_topana_s", ";tau21", 200,0,1);
  h_MET_topana_s = new TH1D("MET_topana_s", ";MET", 400,0,2000);
  h_R2_MR_topana_s = new TH2D("R2_MR_topana_s", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_s = new TH1D("AK8_jet1_pt_topana_s",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_s = new TH1D("AK8_jet1_eta_topana_s",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_s = new TH1D("MET_phi_topana_s", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_topana_Q = new TH1D("ht_AK4_topana_Q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_Q = new TH1D("ht_AK8_topana_Q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_Q = new TH1D("jet1_pt_topana_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_Q = new TH1D("jet2_pt_topana_Q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_Q = new TH1D("jet3_pt_topana_Q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_Q = new TH1D("MR_topana_Q",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_Q = new TH1D("MTR_topana_Q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_Q = new TH1D("R_topana_Q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_Q = new TH1D("R2_topana_Q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_Q = new TH1D("tau21_topana_Q", ";tau21", 200,0,1);
  h_MET_topana_Q = new TH1D("MET_topana_Q", ";MET", 400,0,2000);
  h_R2_MR_topana_Q = new TH2D("R2_MR_topana_Q", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_Q = new TH1D("AK8_jet1_pt_topana_Q",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_Q = new TH1D("AK8_jet1_eta_topana_Q",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_Q = new TH1D("MET_phi_topana_Q", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_topana_q = new TH1D("ht_AK4_topana_q",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_q = new TH1D("ht_AK8_topana_q",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_q = new TH1D("jet1_pt_topana_q",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_q = new TH1D("jet2_pt_topana_q",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_q = new TH1D("jet3_pt_topana_q",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_q = new TH1D("MR_topana_q",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_q = new TH1D("MTR_topana_q",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_q = new TH1D("R_topana_q",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_q = new TH1D("R2_topana_q",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_q = new TH1D("tau21_topana_q", ";tau21", 200,0,1);
  h_MET_topana_q = new TH1D("MET_topana_q", ";MET", 400,0,2000);
  h_R2_MR_topana_q = new TH2D("R2_MR_topana_q", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_q = new TH1D("AK8_jet1_pt_topana_q",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_q = new TH1D("AK8_jet1_eta_topana_q",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_q = new TH1D("MET_phi_topana_q", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_topana_T = new TH1D("ht_AK4_topana_T",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_T = new TH1D("ht_AK8_topana_T",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_T = new TH1D("jet1_pt_topana_T",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_T = new TH1D("jet2_pt_topana_T",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_T = new TH1D("jet3_pt_topana_T",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_T = new TH1D("MR_topana_T",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_T = new TH1D("MTR_topana_T",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_T = new TH1D("R_topana_T",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_T = new TH1D("R2_topana_T",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_T = new TH1D("tau21_topana_T", ";tau21", 200,0,1);
  h_MET_topana_T = new TH1D("MET_topana_T", ";MET", 400,0,2000);
  h_R2_MR_topana_T = new TH2D("R2_MR_topana_T", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_T = new TH1D("AK8_jet1_pt_topana_T",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_T = new TH1D("AK8_jet1_eta_topana_T",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_T = new TH1D("MET_phi_topana_T", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_topana_W = new TH1D("ht_AK4_topana_W",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_W = new TH1D("ht_AK8_topana_W",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_W = new TH1D("jet1_pt_topana_W",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_W = new TH1D("jet2_pt_topana_W",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_W = new TH1D("jet3_pt_topana_W",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_W = new TH1D("MR_topana_W",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_W = new TH1D("MTR_topana_W",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_W = new TH1D("R_topana_W",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_W = new TH1D("R2_topana_W",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_W = new TH1D("tau21_topana_W", ";tau21", 200,0,1);
  h_MET_topana_W = new TH1D("MET_topana_W", ";MET", 400,0,2000);
  h_R2_MR_topana_W = new TH2D("R2_MR_topana_W", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_W = new TH1D("AK8_jet1_pt_topana_W",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_W = new TH1D("AK8_jet1_eta_topana_W",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_W = new TH1D("MET_phi_topana_W", ";#phi_{MET}", 640, -3.2,3.2);

  h_ht_AK4_topana_Z = new TH1D("ht_AK4_topana_Z",  ";H_{T}",                  nbn_HT, bn_HT);
  h_ht_AK8_topana_Z = new TH1D("ht_AK8_topana_Z",  ";H_{T}^{AK8}",            nbn_HT, bn_HT);
  h_jet1_pt_topana_Z = new TH1D("jet1_pt_topana_Z",      ";p_{T, jet1}",            200, 0,2000);
  h_jet2_pt_topana_Z = new TH1D("jet2_pt_topana_Z",      ";p_{T, jet2}",            200, 0,2000);
  h_jet3_pt_topana_Z = new TH1D("jet3_pt_topana_Z",      ";p_{T, jet3}",            200, 0,2000);
  h_MR_topana_Z = new TH1D("MR_topana_Z",   ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_MTR_topana_Z = new TH1D("MTR_topana_Z",  ";MTR_{AK4}",        200, 0,2000);
  h_R_topana_Z = new TH1D("R_topana_Z",    ";R_{AK4}",          nbn_R,bn_R);
  h_R2_topana_Z = new TH1D("R2_topana_Z",   ";R2_{AK4}",         nbn_R2,bn_R2);
  h_tau21_topana_Z = new TH1D("tau21_topana_Z", ";tau21", 200,0,1);
  h_MET_topana_Z = new TH1D("MET_topana_Z", ";MET", 400,0,2000);
  h_R2_MR_topana_Z = new TH2D("R2_MR_topana_Z", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_AK8_jet1_pt_topana_Z = new TH1D("AK8_jet1_pt_topana_Z",      ";p_{T, jet1}",            200, 0,2000);
  h_AK8_jet1_eta_topana_Z = new TH1D("AK8_jet1_eta_topana_Z",      ";#eta_{jet1}",            640, -3.2,3.2);
  h_MET_phi_topana_Z = new TH1D("MET_phi_topana_Z", ";#phi_{MET}", 640, -3.2,3.2);

  h_MR_0mW_EB_DPP = new TH1D("MR_0mW_EB_DPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_0mW_EB_DPP = new TH1D("R2_0mW_EB_DPP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_0mW_EB_DPP = new TH2D("R2_MR_0mW_EB_DPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_0mW_EB_FPP = new TH1D("MR_0mW_EB_FPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_0mW_EB_FPP = new TH1D("R2_0mW_EB_FPP", ";R2_{AK4}",         nbn_R2,bn_R2);   
  h_R2_MR_0mW_EB_FPP = new TH2D("R2_MR_0mW_EB_FPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_0mW_EB_FP = new TH1D("MR_0mW_EB_FP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_0mW_EB_FP = new TH1D("R2_0mW_EB_FP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_0mW_EB_FP  = new TH2D("R2_MR_0mW_EB_FP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_1mW_EB_DPP = new TH1D("MR_1mW_EB_DPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_1mW_EB_DPP = new TH1D("R2_1mW_EB_DPP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_1mW_EB_DPP = new TH2D("R2_MR_1mW_EB_DPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_1mW_EB_FPP = new TH1D("MR_1mW_EB_FPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_1mW_EB_FPP = new TH1D("R2_1mW_EB_FPP", ";R2_{AK4}",         nbn_R2,bn_R2);   
  h_R2_MR_1mW_EB_FPP = new TH2D("R2_MR_1mW_EB_FPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_1mW_EB_FP = new TH1D("MR_1mW_EB_FP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_1mW_EB_FP = new TH1D("R2_1mW_EB_FP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_1mW_EB_FP  = new TH2D("R2_MR_1mW_EB_FP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_0mW_EE_DPP = new TH1D("MR_0mW_EE_DPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_0mW_EE_DPP = new TH1D("R2_0mW_EE_DPP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_0mW_EE_DPP = new TH2D("R2_MR_0mW_EE_DPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_0mW_EE_FPP = new TH1D("MR_0mW_EE_FPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_0mW_EE_FPP = new TH1D("R2_0mW_EE_FPP", ";R2_{AK4}",         nbn_R2,bn_R2);   
  h_R2_MR_0mW_EE_FPP = new TH2D("R2_MR_0mW_EE_FPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_0mW_EE_FP = new TH1D("MR_0mW_EE_FP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_0mW_EE_FP = new TH1D("R2_0mW_EE_FP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_0mW_EE_FP  = new TH2D("R2_MR_0mW_EE_FP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_1mW_EE_DPP = new TH1D("MR_1mW_EE_DPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_1mW_EE_DPP = new TH1D("R2_1mW_EE_DPP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_1mW_EE_DPP = new TH2D("R2_MR_1mW_EE_DPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_1mW_EE_FPP = new TH1D("MR_1mW_EE_FPP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_1mW_EE_FPP = new TH1D("R2_1mW_EE_FPP", ";R2_{AK4}",         nbn_R2,bn_R2);   
  h_R2_MR_1mW_EE_FPP = new TH2D("R2_MR_1mW_EE_FPP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_1mW_EE_FP = new TH1D("MR_1mW_EE_FP", ";MR_{AK4}",     nbn_MR,    bn_MR);
  h_R2_1mW_EE_FP = new TH1D("R2_1mW_EE_FP", ";R2_{AK4}",         nbn_R2,bn_R2);
  h_R2_MR_1mW_EE_FP  = new TH2D("R2_MR_1mW_EE_FP", ";MR_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);

  h_MR_zinv_0mW_EB = new TH1D("MR_zinv_0mW_EB",";MR_{AK4}",nbn_MR,bn_MR);
  h_R2_zinv_0mW_EB = new TH1D("R2_zinv_0mW_EB",";R2_{AK4}",nbn_R2,bn_R2);
  h_R2_MR_zinv_0mW_EB  = new TH2D("R2_MR_zinv_0mW_EB", ";MR_zinv_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_zinv_0mW_EE = new TH1D("MR_zinv_0mW_EE",";MR_{AK4}",nbn_MR,bn_MR);
  h_R2_zinv_0mW_EE = new TH1D("R2_zinv_0mW_EE",";R2_{AK4}",nbn_R2,bn_R2);
  h_R2_MR_zinv_0mW_EE  = new TH2D("R2_MR_zinv_0mW_EE", ";MR_zinv_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_zinv_1mW_EB = new TH1D("MR_zinv_1mW_EB",";MR_{AK4}",nbn_MR,bn_MR);
  h_R2_zinv_1mW_EB = new TH1D("R2_zinv_1mW_EB",";R2_{AK4}",nbn_R2,bn_R2);
  h_R2_MR_zinv_1mW_EB  = new TH2D("R2_MR_zinv_1mW_EB", ";MR_zinv_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_MR_zinv_1mW_EE = new TH1D("MR_zinv_1mW_EE",";MR_{AK4}",nbn_MR,bn_MR);
  h_R2_zinv_1mW_EE = new TH1D("R2_zinv_1mW_EE",";R2_{AK4}",nbn_R2,bn_R2);
  h_R2_MR_zinv_1mW_EE  = new TH2D("R2_MR_zinv_1mW_EE", ";MR_zinv_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_nmW_zinv_EB = new TH1D("nmW_zinv_EB",";nmW",10,0,10);
  h_MR_zinv_EB = new TH1D("MR_zinv_EB",";MR_{AK4}",nbn_MR,bn_MR);
  h_R2_zinv_EB = new TH1D("R2_zinv_EB",";R2_{AK4}",nbn_R2,bn_R2);
  h_R2_MR_zinv_EB  = new TH2D("R2_MR_zinv_EB", ";MR_zinv_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);
  h_nmW_zinv_EE = new TH1D("nmW_zinv_EE",";nmW",10,0,10);
  h_MR_zinv_EE = new TH1D("MR_zinv_EE",";MR_{AK4}",nbn_MR,bn_MR);
  h_R2_zinv_EE = new TH1D("R2_zinv_EE",";R2_{AK4}",nbn_R2,bn_R2);
  h_R2_MR_zinv_EE  = new TH2D("R2_MR_zinv_EE", ";MR_zinv_{AK4};R2_{AK4}",nbn_MR,bn_MR,nbn_R2,bn_R2);

  h_ChargedIso_MR_900_R2_010_SR_mW_EB = new TH1D("ChargedIso_MR_900_R2_010_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_SR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_010_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_SR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_010_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_SR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_010_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_SR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_010_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_SR_mW_EB = new TH1D("ChargedIso_MR_900_R2_014_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_SR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_014_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_SR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_014_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_SR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_014_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_SR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_014_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_SR_mW_EB = new TH1D("ChargedIso_MR_900_R2_020_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_SR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_020_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_SR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_020_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_SR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_020_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_SR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_020_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_SR_mW_EB = new TH1D("ChargedIso_MR_900_R2_037_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_SR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_037_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_SR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_037_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_SR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_037_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_SR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_037_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_SR_mW_EB = new TH1D("ChargedIso_MR_900_R2_075_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_SR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_075_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_SR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_075_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_SR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_075_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_SR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_075_SR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_CR_mW_EB = new TH1D("ChargedIso_MR_900_R2_010_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_CR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_010_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_CR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_010_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_CR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_010_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_CR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_010_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_CR_mW_EB = new TH1D("ChargedIso_MR_900_R2_014_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_CR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_014_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_CR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_014_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_CR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_014_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_CR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_014_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_CR_mW_EB = new TH1D("ChargedIso_MR_900_R2_020_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_CR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_020_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_CR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_020_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_CR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_020_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_CR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_020_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_CR_mW_EB = new TH1D("ChargedIso_MR_900_R2_037_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_CR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_037_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_CR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_037_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_CR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_037_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_CR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_037_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_CR_mW_EB = new TH1D("ChargedIso_MR_900_R2_075_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_CR_mW_EB = new TH1D("ChargedIso_MR_1100_R2_075_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_CR_mW_EB = new TH1D("ChargedIso_MR_1400_R2_075_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_CR_mW_EB = new TH1D("ChargedIso_MR_1800_R2_075_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_CR_mW_EB = new TH1D("ChargedIso_MR_3000_R2_075_CR_mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_SR_mW_EE = new TH1D("ChargedIso_MR_900_R2_010_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_SR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_010_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_SR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_010_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_SR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_010_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_SR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_010_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_SR_mW_EE = new TH1D("ChargedIso_MR_900_R2_014_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_SR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_014_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_SR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_014_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_SR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_014_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_SR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_014_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_SR_mW_EE = new TH1D("ChargedIso_MR_900_R2_020_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_SR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_020_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_SR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_020_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_SR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_020_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_SR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_020_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_SR_mW_EE = new TH1D("ChargedIso_MR_900_R2_037_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_SR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_037_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_SR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_037_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_SR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_037_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_SR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_037_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_SR_mW_EE = new TH1D("ChargedIso_MR_900_R2_075_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_SR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_075_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_SR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_075_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_SR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_075_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_SR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_075_SR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_CR_mW_EE = new TH1D("ChargedIso_MR_900_R2_010_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_CR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_010_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_CR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_010_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_CR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_010_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_CR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_010_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_CR_mW_EE = new TH1D("ChargedIso_MR_900_R2_014_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_CR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_014_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_CR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_014_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_CR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_014_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_CR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_014_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_CR_mW_EE = new TH1D("ChargedIso_MR_900_R2_020_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_CR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_020_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_CR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_020_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_CR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_020_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_CR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_020_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_CR_mW_EE = new TH1D("ChargedIso_MR_900_R2_037_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_CR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_037_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_CR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_037_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_CR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_037_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_CR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_037_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_CR_mW_EE = new TH1D("ChargedIso_MR_900_R2_075_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_CR_mW_EE = new TH1D("ChargedIso_MR_1100_R2_075_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_CR_mW_EE = new TH1D("ChargedIso_MR_1400_R2_075_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_CR_mW_EE = new TH1D("ChargedIso_MR_1800_R2_075_CR_mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_CR_mW_EE = new TH1D("ChargedIso_MR_3000_R2_075_CR_mW_EE", ";Charged Isolation", 100, 0, 10);

  h_ChargedIso_MR_900_R2_010_SR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_010_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_SR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_010_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_SR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_010_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_SR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_010_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_SR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_010_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_SR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_014_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_SR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_014_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_SR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_014_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_SR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_014_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_SR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_014_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_SR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_020_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_SR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_020_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_SR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_020_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_SR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_020_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_SR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_020_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_SR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_037_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_SR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_037_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_SR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_037_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_SR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_037_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_SR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_037_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_SR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_075_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_SR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_075_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_SR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_075_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_SR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_075_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_SR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_075_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_CR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_010_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_CR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_010_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_CR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_010_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_CR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_010_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_CR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_010_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_CR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_014_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_CR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_014_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_CR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_014_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_CR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_014_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_CR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_014_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_CR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_020_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_CR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_020_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_CR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_020_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_CR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_020_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_CR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_020_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_CR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_037_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_CR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_037_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_CR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_037_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_CR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_037_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_CR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_037_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_CR_0mW_EB = new TH1D("ChargedIso_MR_900_R2_075_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_CR_0mW_EB = new TH1D("ChargedIso_MR_1100_R2_075_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_CR_0mW_EB = new TH1D("ChargedIso_MR_1400_R2_075_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_CR_0mW_EB = new TH1D("ChargedIso_MR_1800_R2_075_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_CR_0mW_EB = new TH1D("ChargedIso_MR_3000_R2_075_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_SR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_010_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_SR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_010_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_SR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_010_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_SR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_010_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_SR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_010_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_SR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_014_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_SR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_014_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_SR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_014_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_SR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_014_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_SR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_014_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_SR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_020_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_SR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_020_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_SR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_020_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_SR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_020_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_SR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_020_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_SR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_037_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_SR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_037_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_SR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_037_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_SR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_037_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_SR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_037_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_SR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_075_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_SR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_075_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_SR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_075_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_SR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_075_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_SR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_075_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_CR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_010_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_CR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_010_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_CR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_010_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_CR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_010_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_CR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_010_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_CR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_014_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_CR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_014_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_CR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_014_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_CR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_014_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_CR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_014_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_CR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_020_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_CR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_020_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_CR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_020_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_CR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_020_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_CR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_020_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_CR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_037_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_CR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_037_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_CR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_037_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_CR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_037_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_CR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_037_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_CR_1mW_EB = new TH1D("ChargedIso_MR_900_R2_075_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_CR_1mW_EB = new TH1D("ChargedIso_MR_1100_R2_075_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_CR_1mW_EB = new TH1D("ChargedIso_MR_1400_R2_075_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_CR_1mW_EB = new TH1D("ChargedIso_MR_1800_R2_075_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_CR_1mW_EB = new TH1D("ChargedIso_MR_3000_R2_075_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_SR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_010_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_SR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_010_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_SR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_010_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_SR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_010_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_SR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_010_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_SR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_014_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_SR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_014_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_SR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_014_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_SR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_014_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_SR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_014_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_SR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_020_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_SR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_020_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_SR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_020_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_SR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_020_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_SR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_020_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_SR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_037_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_SR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_037_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_SR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_037_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_SR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_037_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_SR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_037_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_SR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_075_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_SR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_075_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_SR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_075_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_SR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_075_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_SR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_075_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_CR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_010_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_CR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_010_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_CR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_010_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_CR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_010_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_CR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_010_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_CR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_014_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_CR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_014_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_CR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_014_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_CR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_014_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_CR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_014_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_CR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_020_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_CR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_020_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_CR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_020_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_CR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_020_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_CR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_020_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_CR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_037_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_CR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_037_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_CR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_037_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_CR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_037_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_CR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_037_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_CR_0mW_EE = new TH1D("ChargedIso_MR_900_R2_075_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_CR_0mW_EE = new TH1D("ChargedIso_MR_1100_R2_075_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_CR_0mW_EE = new TH1D("ChargedIso_MR_1400_R2_075_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_CR_0mW_EE = new TH1D("ChargedIso_MR_1800_R2_075_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_CR_0mW_EE = new TH1D("ChargedIso_MR_3000_R2_075_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_SR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_010_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_SR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_010_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_SR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_010_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_SR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_010_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_SR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_010_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_SR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_014_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_SR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_014_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_SR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_014_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_SR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_014_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_SR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_014_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_SR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_020_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_SR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_020_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_SR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_020_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_SR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_020_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_SR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_020_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_SR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_037_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_SR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_037_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_SR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_037_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_SR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_037_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_SR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_037_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_SR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_075_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_SR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_075_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_SR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_075_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_SR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_075_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_SR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_075_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_010_CR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_010_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_010_CR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_010_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_010_CR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_010_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_010_CR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_010_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_010_CR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_010_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_014_CR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_014_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_014_CR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_014_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_014_CR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_014_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_014_CR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_014_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_014_CR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_014_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_020_CR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_020_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_020_CR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_020_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_020_CR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_020_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_020_CR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_020_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_020_CR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_020_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_037_CR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_037_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_037_CR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_037_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_037_CR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_037_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_037_CR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_037_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_037_CR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_037_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_R2_075_CR_1mW_EE = new TH1D("ChargedIso_MR_900_R2_075_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_R2_075_CR_1mW_EE = new TH1D("ChargedIso_MR_1100_R2_075_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_R2_075_CR_1mW_EE = new TH1D("ChargedIso_MR_1400_R2_075_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_R2_075_CR_1mW_EE = new TH1D("ChargedIso_MR_1800_R2_075_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_R2_075_CR_1mW_EE = new TH1D("ChargedIso_MR_3000_R2_075_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);

  h_ChargedIso_MR_900_SR_0mW_EB = new TH1D("ChargedIso_MR_900_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_SR_0mW_EB = new TH1D("ChargedIso_MR_1100_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_SR_0mW_EB = new TH1D("ChargedIso_MR_1400_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_SR_0mW_EB = new TH1D("ChargedIso_MR_1800_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_SR_0mW_EB = new TH1D("ChargedIso_MR_3000_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_SR_0mW_EB = new TH1D("ChargedIso_R2_010_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_SR_0mW_EB = new TH1D("ChargedIso_R2_014_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_SR_0mW_EB = new TH1D("ChargedIso_R2_020_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_SR_0mW_EB = new TH1D("ChargedIso_R2_037_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_SR_0mW_EB = new TH1D("ChargedIso_R2_075_SR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_CR_0mW_EB = new TH1D("ChargedIso_MR_900_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_CR_0mW_EB = new TH1D("ChargedIso_MR_1100_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_CR_0mW_EB = new TH1D("ChargedIso_MR_1400_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_CR_0mW_EB = new TH1D("ChargedIso_MR_1800_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_CR_0mW_EB = new TH1D("ChargedIso_MR_3000_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_CR_0mW_EB = new TH1D("ChargedIso_R2_010_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_CR_0mW_EB = new TH1D("ChargedIso_R2_014_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_CR_0mW_EB = new TH1D("ChargedIso_R2_020_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_CR_0mW_EB = new TH1D("ChargedIso_R2_037_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_CR_0mW_EB = new TH1D("ChargedIso_R2_075_CR_0mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_SR_1mW_EB = new TH1D("ChargedIso_MR_900_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_SR_1mW_EB = new TH1D("ChargedIso_MR_1100_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_SR_1mW_EB = new TH1D("ChargedIso_MR_1400_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_SR_1mW_EB = new TH1D("ChargedIso_MR_1800_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_SR_1mW_EB = new TH1D("ChargedIso_MR_3000_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_SR_1mW_EB = new TH1D("ChargedIso_R2_010_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_SR_1mW_EB = new TH1D("ChargedIso_R2_014_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_SR_1mW_EB = new TH1D("ChargedIso_R2_020_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_SR_1mW_EB = new TH1D("ChargedIso_R2_037_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_SR_1mW_EB = new TH1D("ChargedIso_R2_075_SR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_CR_1mW_EB = new TH1D("ChargedIso_MR_900_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_CR_1mW_EB = new TH1D("ChargedIso_MR_1100_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_CR_1mW_EB = new TH1D("ChargedIso_MR_1400_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_CR_1mW_EB = new TH1D("ChargedIso_MR_1800_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_CR_1mW_EB = new TH1D("ChargedIso_MR_3000_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_CR_1mW_EB = new TH1D("ChargedIso_R2_010_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_CR_1mW_EB = new TH1D("ChargedIso_R2_014_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_CR_1mW_EB = new TH1D("ChargedIso_R2_020_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_CR_1mW_EB = new TH1D("ChargedIso_R2_037_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_CR_1mW_EB = new TH1D("ChargedIso_R2_075_CR_1mW_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_SR_0mW_EE = new TH1D("ChargedIso_MR_900_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_SR_0mW_EE = new TH1D("ChargedIso_MR_1100_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_SR_0mW_EE = new TH1D("ChargedIso_MR_1400_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_SR_0mW_EE = new TH1D("ChargedIso_MR_1800_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_SR_0mW_EE = new TH1D("ChargedIso_MR_3000_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_SR_0mW_EE = new TH1D("ChargedIso_R2_010_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_SR_0mW_EE = new TH1D("ChargedIso_R2_014_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_SR_0mW_EE = new TH1D("ChargedIso_R2_020_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_SR_0mW_EE = new TH1D("ChargedIso_R2_037_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_SR_0mW_EE = new TH1D("ChargedIso_R2_075_SR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_CR_0mW_EE = new TH1D("ChargedIso_MR_900_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_CR_0mW_EE = new TH1D("ChargedIso_MR_1100_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_CR_0mW_EE = new TH1D("ChargedIso_MR_1400_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_CR_0mW_EE = new TH1D("ChargedIso_MR_1800_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_CR_0mW_EE = new TH1D("ChargedIso_MR_3000_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_CR_0mW_EE = new TH1D("ChargedIso_R2_010_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_CR_0mW_EE = new TH1D("ChargedIso_R2_014_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_CR_0mW_EE = new TH1D("ChargedIso_R2_020_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_CR_0mW_EE = new TH1D("ChargedIso_R2_037_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_CR_0mW_EE = new TH1D("ChargedIso_R2_075_CR_0mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_SR_1mW_EE = new TH1D("ChargedIso_MR_900_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_SR_1mW_EE = new TH1D("ChargedIso_MR_1100_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_SR_1mW_EE = new TH1D("ChargedIso_MR_1400_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_SR_1mW_EE = new TH1D("ChargedIso_MR_1800_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_SR_1mW_EE = new TH1D("ChargedIso_MR_3000_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_SR_1mW_EE = new TH1D("ChargedIso_R2_010_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_SR_1mW_EE = new TH1D("ChargedIso_R2_014_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_SR_1mW_EE = new TH1D("ChargedIso_R2_020_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_SR_1mW_EE = new TH1D("ChargedIso_R2_037_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_SR_1mW_EE = new TH1D("ChargedIso_R2_075_SR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_CR_1mW_EE = new TH1D("ChargedIso_MR_900_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_CR_1mW_EE = new TH1D("ChargedIso_MR_1100_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_CR_1mW_EE = new TH1D("ChargedIso_MR_1400_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_CR_1mW_EE = new TH1D("ChargedIso_MR_1800_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_CR_1mW_EE = new TH1D("ChargedIso_MR_3000_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_CR_1mW_EE = new TH1D("ChargedIso_R2_010_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_CR_1mW_EE = new TH1D("ChargedIso_R2_014_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_CR_1mW_EE = new TH1D("ChargedIso_R2_020_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_CR_1mW_EE = new TH1D("ChargedIso_R2_037_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_CR_1mW_EE = new TH1D("ChargedIso_R2_075_CR_1mW_EE", ";Charged Isolation", 100, 0, 10);

  h_ChargedIso_MR_900_SR_1mTop_EB = new TH1D("ChargedIso_MR_900_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_SR_1mTop_EB = new TH1D("ChargedIso_MR_1100_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_SR_1mTop_EB = new TH1D("ChargedIso_MR_1400_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_SR_1mTop_EB = new TH1D("ChargedIso_MR_1800_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_SR_1mTop_EB = new TH1D("ChargedIso_MR_3000_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_SR_1mTop_EB = new TH1D("ChargedIso_R2_010_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_SR_1mTop_EB = new TH1D("ChargedIso_R2_014_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_SR_1mTop_EB = new TH1D("ChargedIso_R2_020_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_SR_1mTop_EB = new TH1D("ChargedIso_R2_037_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_SR_1mTop_EB = new TH1D("ChargedIso_R2_075_SR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_CR_1mTop_EB = new TH1D("ChargedIso_MR_900_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_CR_1mTop_EB = new TH1D("ChargedIso_MR_1100_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_CR_1mTop_EB = new TH1D("ChargedIso_MR_1400_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_CR_1mTop_EB = new TH1D("ChargedIso_MR_1800_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_CR_1mTop_EB = new TH1D("ChargedIso_MR_3000_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_CR_1mTop_EB = new TH1D("ChargedIso_R2_010_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_CR_1mTop_EB = new TH1D("ChargedIso_R2_014_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_CR_1mTop_EB = new TH1D("ChargedIso_R2_020_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_CR_1mTop_EB = new TH1D("ChargedIso_R2_037_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_CR_1mTop_EB = new TH1D("ChargedIso_R2_075_CR_1mTop_EB", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_SR_1mTop_EE = new TH1D("ChargedIso_MR_900_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_SR_1mTop_EE = new TH1D("ChargedIso_MR_1100_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_SR_1mTop_EE = new TH1D("ChargedIso_MR_1400_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_SR_1mTop_EE = new TH1D("ChargedIso_MR_1800_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_SR_1mTop_EE = new TH1D("ChargedIso_MR_3000_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_SR_1mTop_EE = new TH1D("ChargedIso_R2_010_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_SR_1mTop_EE = new TH1D("ChargedIso_R2_014_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_SR_1mTop_EE = new TH1D("ChargedIso_R2_020_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_SR_1mTop_EE = new TH1D("ChargedIso_R2_037_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_SR_1mTop_EE = new TH1D("ChargedIso_R2_075_SR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_900_CR_1mTop_EE = new TH1D("ChargedIso_MR_900_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1100_CR_1mTop_EE = new TH1D("ChargedIso_MR_1100_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1400_CR_1mTop_EE = new TH1D("ChargedIso_MR_1400_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_1800_CR_1mTop_EE = new TH1D("ChargedIso_MR_1800_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_MR_3000_CR_1mTop_EE = new TH1D("ChargedIso_MR_3000_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_010_CR_1mTop_EE = new TH1D("ChargedIso_R2_010_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_014_CR_1mTop_EE = new TH1D("ChargedIso_R2_014_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_020_CR_1mTop_EE = new TH1D("ChargedIso_R2_020_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_037_CR_1mTop_EE = new TH1D("ChargedIso_R2_037_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);
  h_ChargedIso_R2_075_CR_1mTop_EE = new TH1D("ChargedIso_R2_075_CR_1mTop_EE", ";Charged Isolation", 100, 0, 10);

  for (unsigned int i=0; i<=syst_nSyst; ++i) {
    std::stringstream histoname, title;
    title<<"Systematic variation #="<<i;

    histoname<<"ht_AK4_S_syst"<<i;
    vh_ht_AK4_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_S_syst"<<i;
    vh_jet1_pt_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_S_syst"<<i;
    vh_jet2_pt_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_S_syst"<<i;
    vh_jet3_pt_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_S_syst"<<i;
    vh_MR_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_S_syst"<<i;
    vh_MTR_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_S_syst"<<i;
    vh_R2_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_S_syst"<<i;
    vh_tau21_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_S_syst"<<i;
    vh_MET_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_S_syst"<<i;
    vh_R2_MR_S.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"AK8_jet1_pt_S_syst"<<i;
    vh_AK8_jet1_pt_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_S_syst"<<i;
    vh_AK8_jet1_eta_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_S_syst"<<i;
    vh_MET_phi_S.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_S_nj35_syst"<<i;
    vh_AK8_jet1_pt_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_S_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_S_nj35.nj35_syst"<<i;
    vh_MET_phi_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_S_nj6_syst"<<i;
    vh_AK8_jet1_pt_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_S_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_S_nj6.nj6_syst"<<i;
    vh_MET_phi_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_s_syst"<<i;
    vh_AK8_jet1_pt_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_s_syst"<<i;
    vh_AK8_jet1_eta_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_s_syst"<<i;
    vh_MET_phi_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_s_nj35_syst"<<i;
    vh_AK8_jet1_pt_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_s_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_s_nj35.nj35_syst"<<i;
    vh_MET_phi_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_s_nj6_syst"<<i;
    vh_AK8_jet1_pt_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_s_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_s_nj6.nj6_syst"<<i;
    vh_MET_phi_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_Q_syst"<<i;
    vh_AK8_jet1_pt_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_Q_syst"<<i;
    vh_AK8_jet1_eta_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_Q_syst"<<i;
    vh_MET_phi_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_Q_nj35_syst"<<i;
    vh_AK8_jet1_pt_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_Q_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_Q_nj35.nj35_syst"<<i;
    vh_MET_phi_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_Q_nj6_syst"<<i;
    vh_AK8_jet1_pt_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_Q_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_Q_nj6.nj6_syst"<<i;
    vh_MET_phi_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_q_syst"<<i;
    vh_AK8_jet1_pt_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_q_syst"<<i;
    vh_AK8_jet1_eta_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_q_syst"<<i;
    vh_MET_phi_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_q_nj35_syst"<<i;
    vh_AK8_jet1_pt_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_q_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_q_nj35.nj35_syst"<<i;
    vh_MET_phi_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_q_nj6_syst"<<i;
    vh_AK8_jet1_pt_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_q_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_q_nj6.nj6_syst"<<i;
    vh_MET_phi_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_T_syst"<<i;
    vh_AK8_jet1_pt_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_T_syst"<<i;
    vh_AK8_jet1_eta_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_T_syst"<<i;
    vh_MET_phi_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_T_nj35_syst"<<i;
    vh_AK8_jet1_pt_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_T_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_T_nj35.nj35_syst"<<i;
    vh_MET_phi_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_T_nj6_syst"<<i;
    vh_AK8_jet1_pt_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_T_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_T_nj6.nj6_syst"<<i;
    vh_MET_phi_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_W_syst"<<i;
    vh_AK8_jet1_pt_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_W_syst"<<i;
    vh_AK8_jet1_eta_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_W_syst"<<i;
    vh_MET_phi_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_W_nj35_syst"<<i;
    vh_AK8_jet1_pt_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_W_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_W_nj35.nj35_syst"<<i;
    vh_MET_phi_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_W_nj6_syst"<<i;
    vh_AK8_jet1_pt_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_W_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_W_nj6.nj6_syst"<<i;
    vh_MET_phi_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_Z_syst"<<i;
    vh_AK8_jet1_pt_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_Z_syst"<<i;
    vh_AK8_jet1_eta_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_Z_syst"<<i;
    vh_MET_phi_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_Z_nj35_syst"<<i;
    vh_AK8_jet1_pt_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_Z_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_Z_nj35.nj35_syst"<<i;
    vh_MET_phi_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_Z_nj6_syst"<<i;
    vh_AK8_jet1_pt_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_Z_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_Z_nj6.nj6_syst"<<i;
    vh_MET_phi_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_G_syst"<<i;
    vh_AK8_jet1_pt_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_G_syst"<<i;
    vh_AK8_jet1_eta_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_G_syst"<<i;
    vh_MET_phi_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_G_nj35_syst"<<i;
    vh_AK8_jet1_pt_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_G_nj35.nj35_syst"<<i;
    vh_AK8_jet1_eta_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_G_nj35.nj35_syst"<<i;
    vh_MET_phi_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));
    histoname.str("");
    histoname<<"AK8_jet1_pt_G_nj6_syst"<<i;
    vh_AK8_jet1_pt_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"AK8_jet1_eta_G_nj6.nj6_syst"<<i;
    vh_AK8_jet1_eta_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#eta_{jet1}").c_str(), 640,-3.2,3.2));
    histoname.str("");
    histoname<<"MET_phi_G_nj6.nj6_syst"<<i;
    vh_MET_phi_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";#phi_{MET}").c_str(),640,-3.2,3.2));

    histoname.str("");
    histoname<<"ht_AK4_s_syst"<<i;
    vh_ht_AK4_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_s_syst"<<i;
    vh_jet1_pt_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_s_syst"<<i;
    vh_jet2_pt_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_s_syst"<<i;
    vh_jet3_pt_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_s_syst"<<i;
    vh_MR_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_s_syst"<<i;
    vh_MTR_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_s_syst"<<i;
    vh_R2_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_s_syst"<<i;
    vh_tau21_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_s_syst"<<i;
    vh_MET_s.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_s_syst"<<i;
    vh_R2_MR_s.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_Q_syst"<<i;
    vh_ht_AK4_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_Q_syst"<<i;
    vh_jet1_pt_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_Q_syst"<<i;
    vh_jet2_pt_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_Q_syst"<<i;
    vh_jet3_pt_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_Q_syst"<<i;
    vh_MR_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_Q_syst"<<i;
    vh_MTR_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_Q_syst"<<i;
    vh_R2_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_Q_syst"<<i;
    vh_tau21_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_Q_syst"<<i;
    vh_MET_Q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_Q_syst"<<i;
    vh_R2_MR_Q.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_q_syst"<<i;
    vh_ht_AK4_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_q_syst"<<i;
    vh_jet1_pt_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_q_syst"<<i;
    vh_jet2_pt_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_q_syst"<<i;
    vh_jet3_pt_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_q_syst"<<i;
    vh_MR_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_q_syst"<<i;
    vh_MTR_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_q_syst"<<i;
    vh_R2_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_q_syst"<<i;
    vh_tau21_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_q_syst"<<i;
    vh_MET_q.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_q_syst"<<i;
    vh_R2_MR_q.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_T_syst"<<i;
    vh_ht_AK4_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_T_syst"<<i;
    vh_jet1_pt_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_T_syst"<<i;
    vh_jet2_pt_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_T_syst"<<i;
    vh_jet3_pt_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_T_syst"<<i;
    vh_MR_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_T_syst"<<i;
    vh_MTR_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_T_syst"<<i;
    vh_R2_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_T_syst"<<i;
    vh_tau21_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_T_syst"<<i;
    vh_MET_T.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_T_syst"<<i;
    vh_R2_MR_T.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_W_syst"<<i;
    vh_ht_AK4_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_W_syst"<<i;
    vh_jet1_pt_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_W_syst"<<i;
    vh_jet2_pt_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_W_syst"<<i;
    vh_jet3_pt_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_W_syst"<<i;
    vh_MR_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_W_syst"<<i;
    vh_MTR_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_W_syst"<<i;
    vh_R2_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_W_syst"<<i;
    vh_tau21_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_W_syst"<<i;
    vh_MET_W.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_W_syst"<<i;
    vh_R2_MR_W.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_Z_syst"<<i;
    vh_ht_AK4_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_Z_syst"<<i;
    vh_jet1_pt_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_Z_syst"<<i;
    vh_jet2_pt_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_Z_syst"<<i;
    vh_jet3_pt_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_Z_syst"<<i;
    vh_MR_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_Z_syst"<<i;
    vh_MTR_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_Z_syst"<<i;
    vh_R2_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_Z_syst"<<i;
    vh_tau21_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_Z_syst"<<i;
    vh_MET_Z.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_Z_syst"<<i;
    vh_R2_MR_Z.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
    histoname.str("");

    histoname<<"ht_AK4_G_syst"<<i;
    vh_ht_AK4_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_G_syst"<<i;
    vh_jet1_pt_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_G_syst"<<i;
    vh_jet2_pt_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_G_syst"<<i;
    vh_jet3_pt_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_G_syst"<<i;
    vh_MR_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_G_syst"<<i;
    vh_MTR_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_G_syst"<<i;
    vh_R2_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_G_syst"<<i;
    vh_tau21_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_G_syst"<<i;
    vh_MET_G.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_G_syst"<<i;
    vh_R2_MR_G.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));
    histoname.str("");

    histoname<<"ht_AK4_S_nj35_syst"<<i;
    vh_ht_AK4_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_S_nj35_syst"<<i;
    vh_jet1_pt_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_S_nj35_syst"<<i;
    vh_jet2_pt_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_S_nj35_syst"<<i;
    vh_jet3_pt_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_S_nj35_syst"<<i;
    vh_MR_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_S_nj35_syst"<<i;
    vh_MTR_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_S_nj35_syst"<<i;
    vh_R2_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_S_nj35_syst"<<i;
    vh_tau21_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_S_nj35_syst"<<i;
    vh_MET_S_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_S_nj35_syst"<<i;
    vh_R2_MR_S_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_s_nj35_syst"<<i;
    vh_ht_AK4_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_s_nj35_syst"<<i;
    vh_jet1_pt_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_s_nj35_syst"<<i;
    vh_jet2_pt_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_s_nj35_syst"<<i;
    vh_jet3_pt_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_s_nj35_syst"<<i;
    vh_MR_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_s_nj35_syst"<<i;
    vh_MTR_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_s_nj35_syst"<<i;
    vh_R2_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_s_nj35_syst"<<i;
    vh_tau21_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_s_nj35_syst"<<i;
    vh_MET_s_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_s_nj35_syst"<<i;
    vh_R2_MR_s_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_Q_nj35_syst"<<i;
    vh_ht_AK4_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_Q_nj35_syst"<<i;
    vh_jet1_pt_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_Q_nj35_syst"<<i;
    vh_jet2_pt_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_Q_nj35_syst"<<i;
    vh_jet3_pt_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_Q_nj35_syst"<<i;
    vh_MR_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_Q_nj35_syst"<<i;
    vh_MTR_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_Q_nj35_syst"<<i;
    vh_R2_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_Q_nj35_syst"<<i;
    vh_tau21_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_Q_nj35_syst"<<i;
    vh_MET_Q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_Q_nj35_syst"<<i;
    vh_R2_MR_Q_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_q_nj35_syst"<<i;
    vh_ht_AK4_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_q_nj35_syst"<<i;
    vh_jet1_pt_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_q_nj35_syst"<<i;
    vh_jet2_pt_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_q_nj35_syst"<<i;
    vh_jet3_pt_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_q_nj35_syst"<<i;
    vh_MR_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_q_nj35_syst"<<i;
    vh_MTR_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_q_nj35_syst"<<i;
    vh_R2_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_q_nj35_syst"<<i;
    vh_tau21_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_q_nj35_syst"<<i;
    vh_MET_q_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_q_nj35_syst"<<i;
    vh_R2_MR_q_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_T_nj35_syst"<<i;
    vh_ht_AK4_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_T_nj35_syst"<<i;
    vh_jet1_pt_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_T_nj35_syst"<<i;
    vh_jet2_pt_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_T_nj35_syst"<<i;
    vh_jet3_pt_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_T_nj35_syst"<<i;
    vh_MR_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_T_nj35_syst"<<i;
    vh_MTR_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_T_nj35_syst"<<i;
    vh_R2_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_T_nj35_syst"<<i;
    vh_tau21_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_T_nj35_syst"<<i;
    vh_MET_T_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_T_nj35_syst"<<i;
    vh_R2_MR_T_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_W_nj35_syst"<<i;
    vh_ht_AK4_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_W_nj35_syst"<<i;
    vh_jet1_pt_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_W_nj35_syst"<<i;
    vh_jet2_pt_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_W_nj35_syst"<<i;
    vh_jet3_pt_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_W_nj35_syst"<<i;
    vh_MR_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_W_nj35_syst"<<i;
    vh_MTR_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_W_nj35_syst"<<i;
    vh_R2_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_W_nj35_syst"<<i;
    vh_tau21_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_W_nj35_syst"<<i;
    vh_MET_W_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_W_nj35_syst"<<i;
    vh_R2_MR_W_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_Z_nj35_syst"<<i;
    vh_ht_AK4_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_Z_nj35_syst"<<i;
    vh_jet1_pt_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_Z_nj35_syst"<<i;
    vh_jet2_pt_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_Z_nj35_syst"<<i;
    vh_jet3_pt_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_Z_nj35_syst"<<i;
    vh_MR_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_Z_nj35_syst"<<i;
    vh_MTR_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_Z_nj35_syst"<<i;
    vh_R2_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_Z_nj35_syst"<<i;
    vh_tau21_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_Z_nj35_syst"<<i;
    vh_MET_Z_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_Z_nj35_syst"<<i;
    vh_R2_MR_Z_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_G_nj35_syst"<<i;
    vh_ht_AK4_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_G_nj35_syst"<<i;
    vh_jet1_pt_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_G_nj35_syst"<<i;
    vh_jet2_pt_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_G_nj35_syst"<<i;
    vh_jet3_pt_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_G_nj35_syst"<<i;
    vh_MR_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_G_nj35_syst"<<i;
    vh_MTR_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_G_nj35_syst"<<i;
    vh_R2_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_G_nj35_syst"<<i;
    vh_tau21_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_G_nj35_syst"<<i;
    vh_MET_G_nj35.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_G_nj35_syst"<<i;
    vh_R2_MR_G_nj35.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_S_nj6_syst"<<i;
    vh_ht_AK4_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_S_nj6_syst"<<i;
    vh_jet1_pt_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_S_nj6_syst"<<i;
    vh_jet2_pt_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_S_nj6_syst"<<i;
    vh_jet3_pt_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_S_nj6_syst"<<i;
    vh_MR_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_S_nj6_syst"<<i;
    vh_MTR_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_S_nj6_syst"<<i;
    vh_R2_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_S_nj6_syst"<<i;
    vh_tau21_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_S_nj6_syst"<<i;
    vh_MET_S_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_S_nj6_syst"<<i;
    vh_R2_MR_S_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_s_nj6_syst"<<i;
    vh_ht_AK4_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_s_nj6_syst"<<i;
    vh_jet1_pt_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_s_nj6_syst"<<i;
    vh_jet2_pt_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_s_nj6_syst"<<i;
    vh_jet3_pt_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_s_nj6_syst"<<i;
    vh_MR_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_s_nj6_syst"<<i;
    vh_MTR_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_s_nj6_syst"<<i;
    vh_R2_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_s_nj6_syst"<<i;
    vh_tau21_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_s_nj6_syst"<<i;
    vh_MET_s_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_s_nj6_syst"<<i;
    vh_R2_MR_s_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_Q_nj6_syst"<<i;
    vh_ht_AK4_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_Q_nj6_syst"<<i;
    vh_jet1_pt_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_Q_nj6_syst"<<i;
    vh_jet2_pt_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_Q_nj6_syst"<<i;
    vh_jet3_pt_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_Q_nj6_syst"<<i;
    vh_MR_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_Q_nj6_syst"<<i;
    vh_MTR_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_Q_nj6_syst"<<i;
    vh_R2_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_Q_nj6_syst"<<i;
    vh_tau21_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_Q_nj6_syst"<<i;
    vh_MET_Q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_Q_nj6_syst"<<i;
    vh_R2_MR_Q_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_q_nj6_syst"<<i;
    vh_ht_AK4_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_q_nj6_syst"<<i;
    vh_jet1_pt_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_q_nj6_syst"<<i;
    vh_jet2_pt_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_q_nj6_syst"<<i;
    vh_jet3_pt_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_q_nj6_syst"<<i;
    vh_MR_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_q_nj6_syst"<<i;
    vh_MTR_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_q_nj6_syst"<<i;
    vh_R2_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_q_nj6_syst"<<i;
    vh_tau21_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_q_nj6_syst"<<i;
    vh_MET_q_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_q_nj6_syst"<<i;
    vh_R2_MR_q_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_T_nj6_syst"<<i;
    vh_ht_AK4_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_T_nj6_syst"<<i;
    vh_jet1_pt_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_T_nj6_syst"<<i;
    vh_jet2_pt_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_T_nj6_syst"<<i;
    vh_jet3_pt_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_T_nj6_syst"<<i;
    vh_MR_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_T_nj6_syst"<<i;
    vh_MTR_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_T_nj6_syst"<<i;
    vh_R2_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_T_nj6_syst"<<i;
    vh_tau21_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_T_nj6_syst"<<i;
    vh_MET_T_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_T_nj6_syst"<<i;
    vh_R2_MR_T_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_W_nj6_syst"<<i;
    vh_ht_AK4_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_W_nj6_syst"<<i;
    vh_jet1_pt_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_W_nj6_syst"<<i;
    vh_jet2_pt_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_W_nj6_syst"<<i;
    vh_jet3_pt_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_W_nj6_syst"<<i;
    vh_MR_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_W_nj6_syst"<<i;
    vh_MTR_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_W_nj6_syst"<<i;
    vh_R2_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_W_nj6_syst"<<i;
    vh_tau21_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_W_nj6_syst"<<i;
    vh_MET_W_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_W_nj6_syst"<<i;
    vh_R2_MR_W_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_Z_nj6_syst"<<i;
    vh_ht_AK4_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_Z_nj6_syst"<<i;
    vh_jet1_pt_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_Z_nj6_syst"<<i;
    vh_jet2_pt_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_Z_nj6_syst"<<i;
    vh_jet3_pt_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_Z_nj6_syst"<<i;
    vh_MR_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_Z_nj6_syst"<<i;
    vh_MTR_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_Z_nj6_syst"<<i;
    vh_R2_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_Z_nj6_syst"<<i;
    vh_tau21_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_Z_nj6_syst"<<i;
    vh_MET_Z_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_Z_nj6_syst"<<i;
    vh_R2_MR_Z_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

    histoname.str("");
    histoname<<"ht_AK4_G_nj6_syst"<<i;
    vh_ht_AK4_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";H_{T}").c_str(), nbn_HT, bn_HT));
    histoname.str("");
    histoname<<"jet1_pt_G_nj6_syst"<<i;
    vh_jet1_pt_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet1}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet2_pt_G_nj6_syst"<<i;
    vh_jet2_pt_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet2}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"jet3_pt_G_nj6_syst"<<i;
    vh_jet3_pt_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";p_{T, jet3}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"MR_G_nj6_syst"<<i;
    vh_MR_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_MR, bn_MR));
    histoname.str("");
    histoname<<"MTR_G_nj6_syst"<<i;
    vh_MTR_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MTR_{AK4}").c_str(), 200, 0,2000));
    histoname.str("");
    histoname<<"R2_G_nj6_syst"<<i;
    vh_R2_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MR_{AK4}").c_str(), nbn_R2, bn_R2));
    histoname.str("");
    histoname<<"tau21_G_nj6_syst"<<i;
    vh_tau21_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";tau21").c_str(), 200, 0,1));
    histoname.str("");
    histoname<<"MET_G_nj6_syst"<<i;
    vh_MET_G_nj6.push_back(new TH1D(histoname.str().c_str(), (title.str()+";MET").c_str(), 400, 0,2000));
    histoname.str("");
    histoname<<"R2_MR_G_nj6_syst"<<i;
    vh_R2_MR_G_nj6.push_back(new TH2D(histoname.str().c_str(), (title.str()+";MR_{AK4};R2_{AK4}").c_str(),nbn_MR,bn_MR,nbn_R2, bn_R2));

  }


  // HT
  h_HT_pre = new TH1D("h_HT_pre", ";H_{T} [GeV]", nbnHT, bnHT);
  h_HT_pre_pass = new TH1D("h_HT_pre_pass", ";H_{T} [GeV]", nbnHT, bnHT);

  h_HT_j1pt_pre = new TH2D("h_HT_j1pt_pre", ";H_{T} [GeV];Leading AK8 jet p_{T} [GeV]", nbnHT, bnHT, nbnj1pt, bnj1pt);
  h_HT_j1pt_pre_pass = new TH2D("h_HT_j1pt_pre_pass", ";H_{T} [GeV];Leading AK8 jet p_{T} [GeV]", nbnHT, bnHT, nbnj1pt, bnj1pt);

  h_AK8Jet1Pt_W_fakerate = new TH1D("AK8Jet1Pt_W_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_W_fakerate = new TH1D("AK8Jet1Pt_no_W_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_mW_fakerate = new TH1D("AK8Jet1Pt_mW_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_mW_fakerate = new TH1D("AK8Jet1Pt_no_mW_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_m0bW_fakerate = new TH1D("AK8Jet1Pt_m0bW_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_m0bW_fakerate = new TH1D("AK8Jet1Pt_no_m0bW_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_aW_fakerate = new TH1D("AK8Jet1Pt_aW_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_aW_fakerate = new TH1D("AK8Jet1Pt_no_aW_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);

  h_AK8Jet1Pt_Top_fakerate = new TH1D("AK8Jet1Pt_Top_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_Top_fakerate = new TH1D("AK8Jet1Pt_no_Top_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_mTop_fakerate = new TH1D("AK8Jet1Pt_mTop_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_mTop_fakerate = new TH1D("AK8Jet1Pt_no_mTop_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_m0bTop_fakerate = new TH1D("AK8Jet1Pt_m0bTop_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_m0bTop_fakerate = new TH1D("AK8Jet1Pt_no_m0bTop_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_aTop_fakerate = new TH1D("AK8Jet1Pt_aTop_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_no_aTop_fakerate = new TH1D("AK8Jet1Pt_no_aTop_fakerate",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);

  h_AK8Jet1Pt_Eta_W_fakerate = new TH2D("AK8Jet1Pt_Eta_W_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_W_fakerate = new TH2D("AK8Jet1Pt_Eta_no_W_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_mW_fakerate = new TH2D("AK8Jet1Pt_Eta_mW_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_mW_fakerate = new TH2D("AK8Jet1Pt_Eta_no_mW_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_m0bW_fakerate = new TH2D("AK8Jet1Pt_Eta_m0bW_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_m0bW_fakerate = new TH2D("AK8Jet1Pt_Eta_no_m0bW_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_aW_fakerate = new TH2D("AK8Jet1Pt_Eta_aW_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_aW_fakerate = new TH2D("AK8Jet1Pt_Eta_no_aW_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);

  h_AK8Jet1Pt_Eta_Top_fakerate = new TH2D("AK8Jet1Pt_Eta_Top_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_Top_fakerate = new TH2D("AK8Jet1Pt_Eta_no_Top_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_mTop_fakerate = new TH2D("AK8Jet1Pt_Eta_mTop_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_mTop_fakerate = new TH2D("AK8Jet1Pt_Eta_no_mTop_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_m0bTop_fakerate = new TH2D("AK8Jet1Pt_Eta_m0bTop_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_m0bTop_fakerate = new TH2D("AK8Jet1Pt_Eta_no_m0bTop_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_aTop_fakerate = new TH2D("AK8Jet1Pt_Eta_aTop_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_Eta_no_aTop_fakerate = new TH2D("AK8Jet1Pt_Eta_no_aTop_fakerate",";p_{T, AK8 jet} [GeV];#eta_{AK8 jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);

  h_tau1_GenW_W = new TH1D("tau1_GenW_W",";#tau_{1}",50,0,1);
  h_tau1_GenW_no_W = new TH1D("tau1_GenW_no_W",";#tau_{1}",50,0,1);
  h_tau1_GenTop_Top = new TH1D("tau1_GenTop_Top",";#tau_{1}",50,0,1);
  h_tau1_GenTop_no_Top = new TH1D("tau1_GenTop_no_Top",";#tau_{1}",50,0,1);

  h_tau2_GenW_W = new TH1D("tau2_GenW_W",";#tau_{2}",50,0,1);
  h_tau2_GenW_no_W = new TH1D("tau2_GenW_no_W",";#tau_{2}",50,0,1);
  h_tau2_GenTop_Top = new TH1D("tau2_GenTop_Top",";#tau_{2}",50,0,1);
  h_tau2_GenTop_no_Top = new TH1D("tau2_GenTop_no_Top",";#tau_{2}",50,0,1);

  h_tau3_GenW_W = new TH1D("tau3_GenW_W",";#tau_{3}",50,0,1);
  h_tau3_GenW_no_W = new TH1D("tau3_GenW_no_W",";#tau_{3}",50,0,1);
  h_tau3_GenTop_Top = new TH1D("tau3_GenTop_Top",";#tau_{3}",50,0,1);
  h_tau3_GenTop_no_Top = new TH1D("tau3_GenTop_no_Top",";#tau_{3}",50,0,1);

  h_tau21_GenW_W = new TH1D("tau21_GenW_W",";#tau_{2}/#tau_{1}",50,0,1);
  h_tau21_GenW_no_W = new TH1D("tau21_GenW_no_W",";#tau_{2}/#tau_{1}",50,0,1);
  h_tau21_GenTop_Top = new TH1D("tau21_GenTop_Top",";#tau_{2}/#tau_{1}",50,0,1);
  h_tau21_GenTop_no_Top = new TH1D("tau21_GenTop_no_Top",";#tau_{2}/#tau_{1}",50,0,1);

  h_tau32_GenW_W = new TH1D("tau32_GenW_W",";#tau_{3}/#tau_{2}",50,0,1);
  h_tau32_GenW_no_W = new TH1D("tau32_GenW_no_W",";#tau_{3}/#tau_{2}",50,0,1);
  h_tau32_GenTop_Top = new TH1D("tau32_GenTop_Top",";#tau_{3}/#tau_{2}",50,0,1);
  h_tau32_GenTop_no_Top = new TH1D("tau32_GenTop_no_Top",";#tau_{3}/#tau_{2}",50,0,1);

  h_SubjetBTag_GenTop_Top    = new TH1D("SubjetBTag_GenTop_Top",   ";Subjet B tag",50,0,1);
  h_SubjetBTag_GenTop_no_Top = new TH1D("SubjetBTag_GenTop_no_Top",";Subjet B tag",50,0,1);

  h_AK8Jet1Pt_GenW_W = new TH1D("AK8Jet1Pt_GenW_W",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_GenW_no_W = new TH1D("AK8Jet1Pt_GenW_no_W",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_GenTop_Top = new TH1D("AK8Jet1Pt_GenTop_Top",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);
  h_AK8Jet1Pt_GenTop_no_Top = new TH1D("AK8Jet1Pt_GenTop_no_Top",";p_{T, AK8 jet} [GeV]",nbn_AK8J1pt,bn_AK8J1pt);

  h_AK8Jet1Pt_eta_GenW_W = new TH2D("AK8Jet1Pt_eta_GenW_W",";p_{T, AK8 jet} [GeV];#eta_{jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_eta_GenW_no_W = new TH2D("AK8Jet1Pt_eta_GenW_no_W",";p_{T, AK8 jet} [GeV];#eta_{jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_eta_GenTop_Top = new TH2D("AK8Jet1Pt_eta_GenTop_Top",";p_{T, AK8 jet} [GeV];#eta_{jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);
  h_AK8Jet1Pt_eta_GenTop_no_Top = new TH2D("AK8Jet1Pt_eta_GenTop_no_Top",";p_{T, AK8 jet} [GeV];#eta_{jet}",nbn_AK8J1pt,bn_AK8J1pt,nbn_eta,bn_eta);

}

//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_analysis_histos(DataStruct& data, const unsigned int& syst_index, const double& weight)
{
    double w = weight; // No scale factor applied
  if (syst_index == 0) {
    // syst_index should only be non-0 if settings.varySystematics is true
    // in case of studying systematics, one should fill a different histogram for each syst_index
    // this variable can be used to chose the correct vector element in case there is a vector of histograms
    // It makes sense, to cut on syst_index == 0, for all ordinary plots
    // syst_index == 0 always guarantees, there are no variations in any systematics
    
    // Check what common variables are available in AnalysisBase.h
    // There a good chance a lot of stuff is already calculated!
    // Especially common object selections or variables to cut on in Analysis

    //bool pass = data.hlt.AK8PFJet360_TrimMass30+data.hlt.PFHT800 == 0 ? false : true;
    //bool pass1 = (data.hlt.AK8PFJet450 == 1 || data.hlt.PFHT800 == 1 || data.hlt.PFHT900 == 1);
    //bool pass2 = (data.hlt.AK8PFJet450 == 1 || data.hlt.PFHT800 == 1);
    //bool pass3 = (data.hlt.AK8PFJet450 == 1 || data.hlt.PFHT900 == 1);

    /*
      Weight:
      They now include trigger efficiencies for MC by default
      w is the event weight without any scale factor applied
      Because scale factors are region dependend, then
      in order to apply them, one has to use the sf_weight[region] variable instead
      eg. sf_weight['S']
     */

    // Baseline cuts 
    // Additionally, let's apply the trigger selection
    // since the weight contains the trigger scaling in MC
    // no specific region, so don't apply scale factors
    // Especially for comparison plots with Ufuk/Fatma
    // Alternatively, could apply SF of the Signal regio

    //double w = sf_weight['S']; // Scale factors applied for the Signal region

    if (apply_all_cuts('P')) h_R2_MR->Fill(data.evt.MR, data.evt.R2, weight);

    if (apply_cut('S',"HLT")) {
      //h_HT_weight->Fill(AK4_Ht,1,w);
      h_njet   ->Fill(nJet,        w);
      h_nb     ->Fill(nMediumBTag, w);
      h_nw     ->Fill(nTightWTag,  w);
      h_ht_gen->Fill(data.evt.Gen_Ht,  w);  // in ntuple
      h_ht_AK4->Fill(AK4_Ht, w); // Calculated in AnalysisBase.h
      h_ht_AK8->Fill(AK8_Ht, w); // Calculated in AnalysisBase.h
    }

    w = sf_weight['Q'];
    if (apply_all_cuts('Q')) h_njetAK8_Q->Fill(nJetAK8,w);
    if (apply_all_cuts('Q')) h_njet_Q->Fill(nJet,w);
    if (apply_all_cuts_except('Q', "3Jet")) h_njet_Q_3Jet->Fill(nJet,w);
    w = sf_weight['S'];
    if (apply_all_cuts('S')) h_njetAK8_S->Fill(nJetAK8,w);
    if (apply_all_cuts('S')) h_njet_S->Fill(nJet,w);
    if (apply_all_cuts_except('S', "3Jet")) h_njet_S_3Jet->Fill(nJet,w);
    if (apply_all_cuts_except('S', "MR_R2")) h_njet_S_MRR2->Fill(nJet,w);
    w = sf_weight['T'];
    if (apply_all_cuts('T')) h_njetAK8_T->Fill(nJetAK8,w);
    if (apply_all_cuts_except('T', "3Jet")) h_njet_T_3Jet->Fill(nJet,w);
    if (apply_all_cuts('T')) h_njet_T->Fill(nJet,w);
    w = sf_weight['W'];
    if (apply_all_cuts('W')) h_njetAK8_W->Fill(nJetAK8,w);
    if (apply_all_cuts_except('W', "3Jet")) h_njet_W_3Jet->Fill(nJet,w);
    if (apply_all_cuts('W')) h_njet_W->Fill(nJet,w);

    // W boson tag fake rate
    w = sf_weight['Y'];
    if (apply_all_cuts_except('Y', "1W")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_W_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_no_m0bW_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_no_aW_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_Eta_no_W_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    h_AK8Jet1Pt_Eta_no_m0bW_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    h_AK8Jet1Pt_Eta_no_aW_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passTightWTag[i]) h_AK8Jet1Pt_W_fakerate->Fill(data.jetsAK8.Pt[i], w);
    if(passWMassTag[i])  h_AK8Jet1Pt_m0bW_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passTightWAntiTag[i])h_AK8Jet1Pt_aW_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passTightWTag[i]) h_AK8Jet1Pt_Eta_W_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w);
    if(passWMassTag[i])  h_AK8Jet1Pt_Eta_m0bW_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passTightWAntiTag[i])h_AK8Jet1Pt_Eta_aW_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
}}

    w = sf_weight['K'];
    if (apply_all_cuts_except('K', "1a0bTop")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_Top_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_no_m0bTop_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_no_aTop_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_Eta_no_Top_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    h_AK8Jet1Pt_Eta_no_m0bTop_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    h_AK8Jet1Pt_Eta_no_aTop_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passHadTopTag[i]) h_AK8Jet1Pt_Top_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passHadTop0BMassTag[i])h_AK8Jet1Pt_m0bTop_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passHadTop0BAntiTag[i])h_AK8Jet1Pt_aTop_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passHadTopTag[i]) h_AK8Jet1Pt_Eta_Top_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passHadTop0BMassTag[i])h_AK8Jet1Pt_Eta_m0bTop_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passHadTop0BAntiTag[i])h_AK8Jet1Pt_Eta_aTop_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    }}

    w = sf_weight['k'];
    if (apply_all_cuts_except('k', "1a0bTop")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_mW_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_Eta_no_mW_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passWMassTag[i])  h_AK8Jet1Pt_mW_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passWMassTag[i])  h_AK8Jet1Pt_Eta_mW_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    h_AK8Jet1Pt_no_mTop_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    h_AK8Jet1Pt_Eta_no_mTop_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    if(passHadTopMassTag[i])h_AK8Jet1Pt_mTop_fakerate->Fill(data.jetsAK8.Pt[i], w); 
    if(passHadTopMassTag[i])h_AK8Jet1Pt_Eta_mTop_fakerate->Fill(data.jetsAK8.Pt[i], std::abs(data.jetsAK8.Eta[i]), w); 
    }}

/*
    // W boson tag fake rate
    w = sf_weight['Y'];
    if (apply_all_cuts('Y')){
    for (size_t i=0; i<iTightWTag.size(); ++i) {
    h_AK8Jet1Pt_W_fakerate->Fill(data.jetsAK8.Pt[iTightWTag[i]], w);
    }}
    if (apply_all_cuts_except('Y', "1W")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_W_fakerate->Fill(data.jetsAK8.Pt[iJetAK8[i]], w); 
    }}
    
    // W boson mass tag fake rate
    w = sf_weight['U'];
    if (apply_all_cuts('U')) {
    for (size_t i=0; i<iWMassTag.size(); ++i) {
    h_AK8Jet1Pt_mW_fakerate->Fill(data.jetsAK8.Pt[iWMassTag[i]], w); 
    }}
    if (apply_all_cuts_except('U', "1mW")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_mW_fakerate->Fill(data.jetsAK8.Pt[iJetAK8[i]], w); 
    }}

    // W boson anti tag fake rate
    w = sf_weight['I'];
    if (apply_all_cuts('I')) {
    for (size_t i=0; i<iTightWAntiTag.size(); ++i) {
    h_AK8Jet1Pt_aW_fakerate->Fill(data.jetsAK8.Pt[iTightWAntiTag[i]], w); 
    }}
    if (apply_all_cuts_except('I', "1aW")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_aW_fakerate->Fill(data.jetsAK8.Pt[iJetAK8[i]], w); 
    }}

    // Top quark tag fake rate
    w = sf_weight['H'];
    if (apply_all_cuts('H')) {
    for (size_t i=0; i<iHadTopTag.size(); ++i) {
    h_AK8Jet1Pt_Top_fakerate->Fill(data.jetsAK8.Pt[iHadTopTag[i]], w); 
    }}
    if (apply_all_cuts_except('H', "1Top")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_Top_fakerate->Fill(data.jetsAK8.Pt[iJetAK8[i]], w); 
    }}

    // Top quark mass tag fake rate
    w = sf_weight['J'];
    if (apply_all_cuts('J')) {
    for (size_t i=0; i<iHadTop0BMassTag.size(); ++i) {
    h_AK8Jet1Pt_mTop_fakerate->Fill(data.jetsAK8.Pt[iHadTop0BMassTag[i]], w); 
    }}
    if (apply_all_cuts_except('J', "1m0bTop")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_mTop_fakerate->Fill(data.jetsAK8.Pt[iJetAK8[i]], w); 
    }}

    // Top quark anti tag fake rate
    w = sf_weight['K'];
    if (apply_all_cuts('K')) {
    for (size_t i=0; i<iHadTop0BAntiTag.size(); ++i) {
    h_AK8Jet1Pt_aTop_fakerate->Fill(data.jetsAK8.Pt[iHadTop0BAntiTag[i]], w); 
    }}
    if (apply_all_cuts_except('K', "1a0bTop")){
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    h_AK8Jet1Pt_no_aTop_fakerate->Fill(data.jetsAK8.Pt[iJetAK8[i]], w); 
    }}
*/
    // W enriched region
    w = sf_weight['W'];
    if (apply_all_cuts('W')) {
      h_ht_AK4_W->Fill(AK4_Ht, w);
      h_ht_AK8_W->Fill(AK8_Ht, w);
      h_jet1_pt_W->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_W->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_W->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_W->Fill(data.evt.MR, w);
      h_MTR_W->Fill(data.evt.MTR, w);
      h_R_W->Fill(data.evt.R, w);
      h_R2_W->Fill(data.evt.R2, w);
      h_tau21_W->Fill(tau21.at(0),w);
      h_MET_W->Fill(data.met.Pt.at(0),w);
      h_R2_MR_W->Fill(data.evt.MR, data.evt.R2, w);
      h_R2_MET_W->Fill(data.met.Pt.at(0), data.evt.R2, w);
      h_MR_MET_W->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_W->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_W->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_W->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_W->Fill(data.met.Phi.at(0),w);

      if(nJet>=3 && nJet<6){
        h_ht_AK4_W_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_W_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_W_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_W_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_W_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_W_nj35->Fill(data.evt.MR, w);
        h_MTR_W_nj35->Fill(data.evt.MTR, w);
        h_R_W_nj35->Fill(data.evt.R, w);
        h_R2_W_nj35->Fill(data.evt.R2, w);
        h_tau21_W_nj35->Fill(tau21.at(0),w);
        h_MET_W_nj35->Fill(data.met.Pt.at(0),w);
        h_R2_MR_W_nj35->Fill(data.evt.MR, data.evt.R2, w);
        h_AK8_jet1_pt_W_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
        h_AK8_jet1_eta_W_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
        h_MET_phi_W_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_W_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_W_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_W_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_W_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_W_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_W_nj6->Fill(data.evt.MR, w);
        h_MTR_W_nj6->Fill(data.evt.MTR, w);
        h_R_W_nj6->Fill(data.evt.R, w);
        h_R2_W_nj6->Fill(data.evt.R2, w);
        h_tau21_W_nj6->Fill(tau21.at(0),w);
        h_MET_W_nj6->Fill(data.met.Pt.at(0),w);
        h_R2_MR_W_nj6->Fill(data.evt.MR, data.evt.R2, w);
        h_AK8_jet1_pt_W_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
        h_AK8_jet1_eta_W_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
        h_MET_phi_W_nj6->Fill(data.met.Phi.at(0),w);
      }

      //h_HT_j1pt_W->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
    if (apply_all_cuts_except('W', "0b")) {
      h_ht_AK4_nonb_W->Fill(AK4_Ht, w);
      h_ht_AK8_nonb_W->Fill(AK8_Ht, w);
      h_jet1_pt_nonb_W->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nonb_W->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nonb_W->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nonb_W->Fill(data.evt.MR, w);
      h_MTR_nonb_W->Fill(data.evt.MTR, w);
      h_R_nonb_W->Fill(data.evt.R, w);
      h_R2_nonb_W->Fill(data.evt.R2, w);
      h_tau21_nonb_W->Fill(tau21.at(0),w);
      h_MET_nonb_W->Fill(data.met.Pt.at(0),w);

      if(nJet>=3 && nJet<6){
        h_ht_AK4_nonb_W_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_nonb_W_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_nonb_W_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_nonb_W_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_nonb_W_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_nonb_W_nj35->Fill(data.evt.MR, w);
        h_MTR_nonb_W_nj35->Fill(data.evt.MTR, w);
        h_R_nonb_W_nj35->Fill(data.evt.R, w);
        h_R2_nonb_W_nj35->Fill(data.evt.R2, w);
        h_tau21_nonb_W_nj35->Fill(tau21.at(0),w);
        h_MET_nonb_W_nj35->Fill(data.met.Pt.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_nonb_W_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_nonb_W_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_nonb_W_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_nonb_W_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_nonb_W_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_nonb_W_nj6->Fill(data.evt.MR, w);
        h_MTR_nonb_W_nj6->Fill(data.evt.MTR, w);
        h_R_nonb_W_nj6->Fill(data.evt.R, w);
        h_R2_nonb_W_nj6->Fill(data.evt.R2, w);
        h_tau21_nonb_W_nj6->Fill(tau21.at(0),w);
        h_MET_nonb_W_nj6->Fill(data.met.Pt.at(0),w);
      }
    }
    if (apply_all_cuts_except('W', "MT")) {
      h_ht_AK4_noMT_W->Fill(AK4_Ht, w);
      h_ht_AK8_noMT_W->Fill(AK8_Ht, w);
      h_jet1_pt_noMT_W->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_noMT_W->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_noMT_W->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_noMT_W->Fill(data.evt.MR, w);
      h_MTR_noMT_W->Fill(data.evt.MTR, w);
      h_R_noMT_W->Fill(data.evt.R, w);
      h_R2_noMT_W->Fill(data.evt.R2, w);
      h_tau21_noMT_W->Fill(tau21.at(0),w);
      h_MET_noMT_W->Fill(data.met.Pt.at(0),w);

      if(nJet>=3 && nJet<6){
        h_ht_AK4_noMT_W_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_noMT_W_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_noMT_W_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_noMT_W_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_noMT_W_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_noMT_W_nj35->Fill(data.evt.MR, w);
        h_MTR_noMT_W_nj35->Fill(data.evt.MTR, w);
        h_R_noMT_W_nj35->Fill(data.evt.R, w);
        h_R2_noMT_W_nj35->Fill(data.evt.R2, w);
        h_tau21_noMT_W_nj35->Fill(tau21.at(0),w);
        h_MET_noMT_W_nj35->Fill(data.met.Pt.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_noMT_W_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_noMT_W_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_noMT_W_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_noMT_W_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_noMT_W_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_noMT_W_nj6->Fill(data.evt.MR, w);
        h_MTR_noMT_W_nj6->Fill(data.evt.MTR, w);
        h_R_noMT_W_nj6->Fill(data.evt.R, w);
        h_R2_noMT_W_nj6->Fill(data.evt.R2, w);
        h_tau21_noMT_W_nj6->Fill(tau21.at(0),w);
        h_MET_noMT_W_nj6->Fill(data.met.Pt.at(0),w);
      }
    }

    // top enriched region
    w = sf_weight['T'];
    if (apply_all_cuts('T')) {
      h_ht_AK4_T->Fill(AK4_Ht, w);
      h_ht_AK8_T->Fill(AK8_Ht, w);
      h_jet1_pt_T->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_T->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_T->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_T->Fill(data.evt.MR, w);
      h_MTR_T->Fill(data.evt.MTR, w);
      h_R_T->Fill(data.evt.R, w);
      h_R2_T->Fill(data.evt.R2, w);
      h_tau21_T->Fill(tau21.at(0),w);
      h_MET_T->Fill(data.met.Pt.at(0),w);
      h_R2_MR_T->Fill(data.evt.MR, data.evt.R2, w);
      h_R2_MET_T->Fill(data.met.Pt.at(0), data.evt.R2, w);
      h_MR_MET_T->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_T->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_T->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_T->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_T->Fill(data.met.Phi.at(0),w);

      if(nJet>=3 && nJet<6){
      h_ht_AK4_T_nj35->Fill(AK4_Ht, w);
      h_ht_AK8_T_nj35->Fill(AK8_Ht, w);
      h_jet1_pt_T_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_T_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_T_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_T_nj35->Fill(data.evt.MR, w);
      h_MTR_T_nj35->Fill(data.evt.MTR, w);
      h_R_T_nj35->Fill(data.evt.R, w);
      h_R2_T_nj35->Fill(data.evt.R2, w);
      h_tau21_T_nj35->Fill(tau21.at(0),w);
      h_MET_T_nj35->Fill(data.met.Pt.at(0),w);
      h_R2_MR_T_nj35->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_T_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_T_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_T_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
      h_ht_AK4_T_nj6->Fill(AK4_Ht, w);
      h_ht_AK8_T_nj6->Fill(AK8_Ht, w);
      h_jet1_pt_T_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_T_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_T_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_T_nj6->Fill(data.evt.MR, w);
      h_MTR_T_nj6->Fill(data.evt.MTR, w);
      h_R_T_nj6->Fill(data.evt.R, w);
      h_R2_T_nj6->Fill(data.evt.R2, w);
      h_tau21_T_nj6->Fill(tau21.at(0),w);
      h_MET_T_nj6->Fill(data.met.Pt.at(0),w);
      h_R2_MR_T_nj6->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_T_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_T_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_T_nj6->Fill(data.met.Phi.at(0),w);
      }

      //h_HT_j1pt_T->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
    if (apply_all_cuts_except('T', "MT")) {
      h_ht_AK4_noMT_T->Fill(AK4_Ht, w);
      h_ht_AK8_noMT_T->Fill(AK8_Ht, w);
      h_jet1_pt_noMT_T->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_noMT_T->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_noMT_T->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_noMT_T->Fill(data.evt.MR, w);
      h_MTR_noMT_T->Fill(data.evt.MTR, w);
      h_R_noMT_T->Fill(data.evt.R, w);
      h_R2_noMT_T->Fill(data.evt.R2, w);
      h_tau21_noMT_T->Fill(tau21.at(0),w);
      h_MET_noMT_T->Fill(data.met.Pt.at(0),w);

      if(nJet>=3 && nJet<6){
      h_ht_AK4_noMT_T_nj35->Fill(AK4_Ht, w);
      h_ht_AK8_noMT_T_nj35->Fill(AK8_Ht, w);
      h_jet1_pt_noMT_T_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_noMT_T_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_noMT_T_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_noMT_T_nj35->Fill(data.evt.MR, w);
      h_MTR_noMT_T_nj35->Fill(data.evt.MTR, w);
      h_R_noMT_T_nj35->Fill(data.evt.R, w);
      h_R2_noMT_T_nj35->Fill(data.evt.R2, w);
      h_tau21_noMT_T_nj35->Fill(tau21.at(0),w);
      h_MET_noMT_T_nj35->Fill(data.met.Pt.at(0),w);
      }
      if(nJet>=6){
      h_ht_AK4_noMT_T_nj6->Fill(AK4_Ht, w);
      h_ht_AK8_noMT_T_nj6->Fill(AK8_Ht, w);
      h_jet1_pt_noMT_T_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_noMT_T_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_noMT_T_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_noMT_T_nj6->Fill(data.evt.MR, w);
      h_MTR_noMT_T_nj6->Fill(data.evt.MTR, w);
      h_R_noMT_T_nj6->Fill(data.evt.R, w);
      h_R2_noMT_T_nj6->Fill(data.evt.R2, w);
      h_tau21_noMT_T_nj6->Fill(tau21.at(0),w);
      h_MET_noMT_T_nj6->Fill(data.met.Pt.at(0),w);
      }
    }
    

    // Z enriched region
    w = sf_weight['Z'];
    if (apply_all_cuts('Z')) {
      h_ht_AK4_Z->Fill(AK4_Ht, w);
      h_ht_AK8_Z->Fill(AK8_Ht, w);
      h_jet1_pt_Z->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_Z->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_Z->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_Z->Fill(data.evt.MR, w);
      h_MTR_Z->Fill(data.evt.MTR, w);
      h_R_Z->Fill(data.evt.R, w);
      h_R2_Z->Fill(R2_ll, w);
      h_tau21_Z->Fill(tau21.at(0),w);
      h_MET_Z->Fill(data.met.Pt.at(0),w);
      h_R2_MR_Z->Fill(data.evt.MR, R2_ll, w);
      h_R2_MET_Z->Fill(data.met.Pt.at(0), R2_ll, w);
      h_MR_MET_Z->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_Z->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_Z->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_Z->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_Z->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        h_ht_AK4_Z_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_Z_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_Z_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_Z_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_Z_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_Z_nj35->Fill(data.evt.MR, w);
        h_R_Z_nj35->Fill(data.evt.R, w);
        h_MTR_Z_nj35->Fill(data.evt.MTR, w);
        h_R2_Z_nj35->Fill(R2_ll, w);
        h_tau21_Z_nj35->Fill(tau21.at(0),w);
        h_MET_Z_nj35->Fill(data.met.Pt.at(0),w);
        h_R2_MR_Z_nj35->Fill(data.evt.MR, R2_ll, w);
      h_AK8_jet1_pt_Z_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_Z_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_Z_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_Z_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_Z_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_Z_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_Z_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_Z_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_Z_nj6->Fill(data.evt.MR, w);
        h_R_Z_nj6->Fill(data.evt.R, w);
        h_MTR_Z_nj6->Fill(data.evt.MTR, w);
        h_R2_Z_nj6->Fill(R2_ll, w);
        h_tau21_Z_nj6->Fill(tau21.at(0),w);
        h_MET_Z_nj6->Fill(data.met.Pt.at(0),w);
        h_R2_MR_Z_nj6->Fill(data.evt.MR, R2_ll, w);
      h_AK8_jet1_pt_Z_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_Z_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_Z_nj6->Fill(data.met.Phi.at(0),w);
      }
    }

    // G enriched region
    w = sf_weight['G'];
    if (apply_all_cuts('G')) {
      h_ht_AK4_G->Fill(AK4_Ht, w);
      h_ht_AK8_G->Fill(AK8_Ht, w);
      h_jet1_pt_G->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_G->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_G->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_G->Fill(data.evt.MR, w);
      h_MTR_G->Fill(data.evt.MTR, w);
      h_R_G->Fill(data.evt.R, w);
      h_R2_G->Fill(R2_pho, w);
      h_tau21_G->Fill(tau21.at(0),w);
      h_MET_G->Fill(data.met.Pt.at(0),w);
      h_R2_MR_G->Fill(data.evt.MR, R2_pho, w);
      h_R2_MET_G->Fill(data.met.Pt.at(0), R2_pho, w);
      h_MR_MET_G->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_G->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_G->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_G->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_G->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        h_ht_AK4_G_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_G_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_G_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_G_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_G_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_G_nj35->Fill(data.evt.MR, w);
        h_R_G_nj35->Fill(data.evt.R, w);
        h_MTR_G_nj35->Fill(data.evt.MTR, w);
        h_R2_G_nj35->Fill(R2_pho, w);
        h_tau21_G_nj35->Fill(tau21.at(0),w);
        h_MET_G_nj35->Fill(data.met.Pt.at(0),w);
        h_R2_MR_G_nj35->Fill(data.evt.MR, R2_pho, w);
      h_AK8_jet1_pt_G_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_G_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_G_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_G_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_G_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_G_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_G_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_G_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_G_nj6->Fill(data.evt.MR, w);
        h_R_G_nj6->Fill(data.evt.R, w);
        h_MTR_G_nj6->Fill(data.evt.MTR, w);
        h_R2_G_nj6->Fill(R2_pho, w);
        h_tau21_G_nj6->Fill(tau21.at(0),w);
        h_MET_G_nj6->Fill(data.met.Pt.at(0),w);
        h_R2_MR_G_nj6->Fill(data.evt.MR, R2_pho, w);
      h_AK8_jet1_pt_G_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_G_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_G_nj6->Fill(data.met.Phi.at(0),w);
      }
    }

    // QCD enriched region
    w = sf_weight['q'];
    if (apply_all_cuts('q')) {
      h_ht_AK4_q->Fill(AK4_Ht, w);
      h_ht_AK8_q->Fill(AK8_Ht, w);
      h_jet1_pt_q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_q->Fill(data.evt.MR, w);
      h_MTR_q->Fill(data.evt.MTR, w);
      h_R_q->Fill(data.evt.R, w);
      h_R2_q->Fill(data.evt.R2, w);
      h_tau21_q->Fill(tau21.at(0),w);
      h_MET_q->Fill(data.met.Pt.at(0),w);
      h_R2_MR_q->Fill(data.evt.MR, data.evt.R2, w);
      h_R2_MET_q->Fill(data.met.Pt.at(0), data.evt.R2, w);
      h_MR_MET_q->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_q->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_q->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_q->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_q->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        h_ht_AK4_q_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_q_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_q_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_q_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_q_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_q_nj35->Fill(data.evt.MR, w);
        h_R_q_nj35->Fill(data.evt.R, w);
        h_MTR_q_nj35->Fill(data.evt.MTR, w);
        h_R2_q_nj35->Fill(data.evt.R2, w);
        h_tau21_q_nj35->Fill(tau21.at(0),w);
        h_MET_q_nj35->Fill(data.met.Pt.at(0),w);
        h_R2_MR_q_nj35->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_q_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_q_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_q_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_q_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_q_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_q_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_q_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_q_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_q_nj6->Fill(data.evt.MR, w);
        h_R_q_nj6->Fill(data.evt.R, w);
        h_MTR_q_nj6->Fill(data.evt.MTR, w);
        h_R2_q_nj6->Fill(data.evt.R2, w);
        h_tau21_q_nj6->Fill(tau21.at(0),w);
        h_MET_q_nj6->Fill(data.met.Pt.at(0),w);
        h_R2_MR_q_nj6->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_q_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_q_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_q_nj6->Fill(data.met.Phi.at(0),w);
      }
    }

    // QCD enriched region
    w = sf_weight['Q'];
    if (apply_all_cuts('Q')) {
      h_ht_AK4_Q->Fill(AK4_Ht, w);
      h_ht_AK8_Q->Fill(AK8_Ht, w);
      h_jet1_pt_Q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_Q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_Q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_Q->Fill(data.evt.MR, w);
      h_MTR_Q->Fill(data.evt.MTR, w);
      h_R_Q->Fill(data.evt.R, w);
      h_R2_Q->Fill(data.evt.R2, w);
      h_tau21_Q->Fill(tau21.at(0),w);
      h_MET_Q->Fill(data.met.Pt.at(0),w);
      h_R2_MR_Q->Fill(data.evt.MR, data.evt.R2, w);
      h_R2_MET_Q->Fill(data.met.Pt.at(0), data.evt.R2, w);
      h_MR_MET_Q->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_Q->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_Q->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_Q->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_Q->Fill(data.met.Phi.at(0),w);

      if(nJet>=3 && nJet<6){
      h_ht_AK4_Q_nj35->Fill(AK4_Ht, w);
      h_ht_AK8_Q_nj35->Fill(AK8_Ht, w);
      h_jet1_pt_Q_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_Q_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_Q_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_Q_nj35->Fill(data.evt.MR, w);
      h_MTR_Q_nj35->Fill(data.evt.MTR, w);
      h_R_Q_nj35->Fill(data.evt.R, w);
      h_R2_Q_nj35->Fill(data.evt.R2, w);
      h_tau21_Q_nj35->Fill(tau21.at(0),w);
      h_MET_Q_nj35->Fill(data.met.Pt.at(0),w);
      h_R2_MR_Q_nj35->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_Q_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_Q_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_Q_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
      h_ht_AK4_Q_nj6->Fill(AK4_Ht, w);
      h_ht_AK8_Q_nj6->Fill(AK8_Ht, w);
      h_jet1_pt_Q_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_Q_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_Q_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_Q_nj6->Fill(data.evt.MR, w);
      h_MTR_Q_nj6->Fill(data.evt.MTR, w);
      h_R_Q_nj6->Fill(data.evt.R, w);
      h_R2_Q_nj6->Fill(data.evt.R2, w);
      h_tau21_Q_nj6->Fill(tau21.at(0),w);
      h_MET_Q_nj6->Fill(data.met.Pt.at(0),w);
      h_R2_MR_Q_nj6->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_Q_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_Q_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_Q_nj6->Fill(data.met.Phi.at(0),w);
      }

      //h_HT_Q->Fill(AK4_Ht,w);
      //h_j1_pt_Q->Fill(data.jetsAK8.Pt[iJetAK8[0]],w);
      //h_HT_j1pt_Q->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
    if (apply_all_cuts_except('Q', "InvmDPhi0p3")) {
      h_ht_AK4_nodPhi_Q->Fill(AK4_Ht, w);
      h_ht_AK8_nodPhi_Q->Fill(AK8_Ht, w);
      h_jet1_pt_nodPhi_Q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nodPhi_Q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nodPhi_Q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nodPhi_Q->Fill(data.evt.MR, w);
      h_MTR_nodPhi_Q->Fill(data.evt.MTR, w);
      h_R_nodPhi_Q->Fill(data.evt.R, w);
      h_R2_nodPhi_Q->Fill(data.evt.R2, w);
      h_tau21_nodPhi_Q->Fill(tau21.at(0),w);
      h_MET_nodPhi_Q->Fill(data.met.Pt.at(0),w);
      h_dPhiRazor_nodPhi_Q->Fill(dPhiRazor,w);

      if(nJet>=3 && nJet<6){
      h_ht_AK4_nodPhi_Q_nj35->Fill(AK4_Ht, w);
      h_ht_AK8_nodPhi_Q_nj35->Fill(AK8_Ht, w);
      h_jet1_pt_nodPhi_Q_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nodPhi_Q_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nodPhi_Q_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nodPhi_Q_nj35->Fill(data.evt.MR, w);
      h_MTR_nodPhi_Q_nj35->Fill(data.evt.MTR, w);
      h_R_nodPhi_Q_nj35->Fill(data.evt.R, w);
      h_R2_nodPhi_Q_nj35->Fill(data.evt.R2, w);
      h_tau21_nodPhi_Q_nj35->Fill(tau21.at(0),w);
      h_MET_nodPhi_Q_nj35->Fill(data.met.Pt.at(0),w);
      h_dPhiRazor_nodPhi_Q_nj35->Fill(dPhiRazor,w);
      }
      if(nJet>=6){
      h_ht_AK4_nodPhi_Q_nj6->Fill(AK4_Ht, w);
      h_ht_AK8_nodPhi_Q_nj6->Fill(AK8_Ht, w);
      h_jet1_pt_nodPhi_Q_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nodPhi_Q_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nodPhi_Q_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nodPhi_Q_nj6->Fill(data.evt.MR, w);
      h_MTR_nodPhi_Q_nj6->Fill(data.evt.MTR, w);
      h_R_nodPhi_Q_nj6->Fill(data.evt.R, w);
      h_R2_nodPhi_Q_nj6->Fill(data.evt.R2, w);
      h_tau21_nodPhi_Q_nj6->Fill(tau21.at(0),w);
      h_MET_nodPhi_Q_nj6->Fill(data.met.Pt.at(0),w);
      h_dPhiRazor_nodPhi_Q_nj6->Fill(dPhiRazor,w);
      }
    }
    if (apply_all_cuts_except('Q', "0b")) {
      h_ht_AK4_nonb_Q->Fill(AK4_Ht, w);
      h_ht_AK8_nonb_Q->Fill(AK8_Ht, w);
      h_jet1_pt_nonb_Q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nonb_Q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nonb_Q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nonb_Q->Fill(data.evt.MR, w);
      h_MTR_nonb_Q->Fill(data.evt.MTR, w);
      h_R_nonb_Q->Fill(data.evt.R, w);
      h_R2_nonb_Q->Fill(data.evt.R2, w);
      h_tau21_nonb_Q->Fill(tau21.at(0),w);
      h_MET_nonb_Q->Fill(data.met.Pt.at(0),w);

      if(nJet>=3 && nJet<6){
      h_ht_AK4_nonb_Q_nj35->Fill(AK4_Ht, w);
      h_ht_AK8_nonb_Q_nj35->Fill(AK8_Ht, w);
      h_jet1_pt_nonb_Q_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nonb_Q_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nonb_Q_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nonb_Q_nj35->Fill(data.evt.MR, w);
      h_MTR_nonb_Q_nj35->Fill(data.evt.MTR, w);
      h_R_nonb_Q_nj35->Fill(data.evt.R, w);
      h_R2_nonb_Q_nj35->Fill(data.evt.R2, w);
      h_tau21_nonb_Q_nj35->Fill(tau21.at(0),w);
      h_MET_nonb_Q_nj35->Fill(data.met.Pt.at(0),w);
      }
      if(nJet>=6){
      h_ht_AK4_nonb_Q_nj6->Fill(AK4_Ht, w);
      h_ht_AK8_nonb_Q_nj6->Fill(AK8_Ht, w);
      h_jet1_pt_nonb_Q_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_nonb_Q_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_nonb_Q_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_nonb_Q_nj6->Fill(data.evt.MR, w);
      h_MTR_nonb_Q_nj6->Fill(data.evt.MTR, w);
      h_R_nonb_Q_nj6->Fill(data.evt.R, w);
      h_R2_nonb_Q_nj6->Fill(data.evt.R2, w);
      h_tau21_nonb_Q_nj6->Fill(tau21.at(0),w);
      h_MET_nonb_Q_nj6->Fill(data.met.Pt.at(0),w);
      }
    }
    if (apply_all_cuts_except('Q', "1aW")) {
      h_ht_AK4_notau21_Q->Fill(AK4_Ht, w);
      h_ht_AK8_notau21_Q->Fill(AK8_Ht, w);
      h_jet1_pt_notau21_Q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_notau21_Q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_notau21_Q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_notau21_Q->Fill(data.evt.MR, w);
      h_MTR_notau21_Q->Fill(data.evt.MTR, w);
      h_R_notau21_Q->Fill(data.evt.R, w);
      h_R2_notau21_Q->Fill(data.evt.R2, w);
      h_tau21_notau21_Q->Fill(tau21.at(0),w);
      h_MET_notau21_Q->Fill(data.met.Pt.at(0),w);

      if(nJet>=3 && nJet<6){
      h_ht_AK4_notau21_Q_nj35->Fill(AK4_Ht, w);
      h_ht_AK8_notau21_Q_nj35->Fill(AK8_Ht, w);
      h_jet1_pt_notau21_Q_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_notau21_Q_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_notau21_Q_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_notau21_Q_nj35->Fill(data.evt.MR, w);
      h_MTR_notau21_Q_nj35->Fill(data.evt.MTR, w);
      h_R_notau21_Q_nj35->Fill(data.evt.R, w);
      h_R2_notau21_Q_nj35->Fill(data.evt.R2, w);
      h_tau21_notau21_Q_nj35->Fill(tau21.at(0),w);
      h_MET_notau21_Q_nj35->Fill(data.met.Pt.at(0),w);
      }
      if(nJet>=6){
      h_ht_AK4_notau21_Q_nj6->Fill(AK4_Ht, w);
      h_ht_AK8_notau21_Q_nj6->Fill(AK8_Ht, w);
      h_jet1_pt_notau21_Q_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_notau21_Q_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_notau21_Q_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_notau21_Q_nj6->Fill(data.evt.MR, w);
      h_MTR_notau21_Q_nj6->Fill(data.evt.MTR, w);
      h_R_notau21_Q_nj6->Fill(data.evt.R, w);
      h_R2_notau21_Q_nj6->Fill(data.evt.R2, w);
      h_tau21_notau21_Q_nj6->Fill(tau21.at(0),w);
      h_MET_notau21_Q_nj6->Fill(data.met.Pt.at(0),w);
      }
    }


    // Signal' region
    w = sf_weight['s'];
    //if (apply_all_cuts_except('s', "HLT")) {
    if (apply_all_cuts('s')) {
      h_ht_AK4_s->Fill(AK4_Ht, w);
      h_ht_AK8_s->Fill(AK8_Ht, w);
      h_jet1_pt_s->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_s->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_s->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_s->Fill(data.evt.MR, w);
      h_R_s->Fill(data.evt.R, w);
      h_MTR_s->Fill(data.evt.MTR, w);
      h_R2_s->Fill(data.evt.R2, w);
      h_tau21_s->Fill(tau21.at(0),w);
      h_MET_s->Fill(data.met.Pt.at(0),w);
      h_R2_MR_s->Fill(data.evt.MR, data.evt.R2, w);
      h_R2_MET_s->Fill(data.met.Pt.at(0), data.evt.R2, w);
      h_MR_MET_s->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_s->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_s->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_s->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_s->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        h_ht_AK4_s_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_s_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_s_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_s_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_s_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_s_nj35->Fill(data.evt.MR, w);
        h_R_s_nj35->Fill(data.evt.R, w);
        h_MTR_s_nj35->Fill(data.evt.MTR, w);
        h_R2_s_nj35->Fill(data.evt.R2, w);
        h_tau21_s_nj35->Fill(tau21.at(0),w);
        h_MET_s_nj35->Fill(data.met.Pt.at(0),w);
        h_R2_MR_s_nj35->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_s_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_s_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_s_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_s_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_s_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_s_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_s_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_s_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_s_nj6->Fill(data.evt.MR, w);
        h_R_s_nj6->Fill(data.evt.R, w);
        h_MTR_s_nj6->Fill(data.evt.MTR, w);
        h_R2_s_nj6->Fill(data.evt.R2, w);
        h_tau21_s_nj6->Fill(tau21.at(0),w);
        h_MET_s_nj6->Fill(data.met.Pt.at(0),w);
        h_R2_MR_s_nj6->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_s_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_s_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_s_nj6->Fill(data.met.Phi.at(0),w);
      }
    }

    // Signal region
    w = sf_weight['S'];
    //if (apply_all_cuts_except('S', "HLT")) {
    if (apply_all_cuts('S')) {
      h_ht_AK4_S->Fill(AK4_Ht, w);
      h_ht_AK8_S->Fill(AK8_Ht, w);
      h_jet1_pt_S->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_S->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_S->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_S->Fill(data.evt.MR, w);
      h_R_S->Fill(data.evt.R, w);
      h_MTR_S->Fill(data.evt.MTR, w);
      h_R2_S->Fill(data.evt.R2, w);
      h_tau21_S->Fill(tau21.at(0),w);
      h_MET_S->Fill(data.met.Pt.at(0),w);
      h_R2_MR_S->Fill(data.evt.MR, data.evt.R2, w);
      h_R2_MET_S->Fill(data.met.Pt.at(0), data.evt.R2, w);
      h_MR_MET_S->Fill(data.met.Pt.at(0), data.evt.MR, w);
      h_AK8Jet1pT_MET_S->Fill(data.met.Pt.at(0), data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_pt_S->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_S->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_S->Fill(data.met.Phi.at(0),w);

      if(nJet>=3 && nJet<6){
        h_ht_AK4_S_nj35->Fill(AK4_Ht, w);
        h_ht_AK8_S_nj35->Fill(AK8_Ht, w);
        h_jet1_pt_S_nj35->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_S_nj35->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_S_nj35->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_S_nj35->Fill(data.evt.MR, w);
        h_R_S_nj35->Fill(data.evt.R, w);
        h_MTR_S_nj35->Fill(data.evt.MTR, w);
        h_R2_S_nj35->Fill(data.evt.R2, w);
        h_tau21_S_nj35->Fill(tau21.at(0),w);
        h_MET_S_nj35->Fill(data.met.Pt.at(0),w);
        h_R2_MR_S_nj35->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_S_nj35->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_S_nj35->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_S_nj35->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        h_ht_AK4_S_nj6->Fill(AK4_Ht, w);
        h_ht_AK8_S_nj6->Fill(AK8_Ht, w);
        h_jet1_pt_S_nj6->Fill(data.jetsAK4.Pt[iJet[0]], w);
        h_jet2_pt_S_nj6->Fill(data.jetsAK4.Pt[iJet[1]], w);
        h_jet3_pt_S_nj6->Fill(data.jetsAK4.Pt[iJet[2]], w);
        h_MR_S_nj6->Fill(data.evt.MR, w);
        h_R_S_nj6->Fill(data.evt.R, w);
        h_MTR_S_nj6->Fill(data.evt.MTR, w);
        h_R2_S_nj6->Fill(data.evt.R2, w);
        h_tau21_S_nj6->Fill(tau21.at(0),w);
        h_MET_S_nj6->Fill(data.met.Pt.at(0),w);
        h_R2_MR_S_nj6->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_S_nj6->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_S_nj6->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_S_nj6->Fill(data.met.Phi.at(0),w);
      }
      h_softDropMass->Fill(softDropMassW.at(0),w);
      h_GluinoLSPMass->Fill(data.evt.SUSY_Gluino_Mass,data.evt.SUSY_LSP_Mass,w);
      h_StopMass->Fill(data.evt.SUSY_Stop_Mass,w);
      h_GluinoMass->Fill(data.evt.SUSY_Gluino_Mass,w);
      h_LSPMass->Fill(data.evt.SUSY_LSP_Mass,w);

      //h_HT_S->Fill(AK4_Ht,w);
      //h_j1_pt_S->Fill(data.jetsAK8.Pt[iJetAK8[0]],w);
      //h_HT_j1pt_S->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }



    w = sf_weight['A'];
    if (apply_all_cuts('A')) {
      h_ht_AK4_topana_S->Fill(AK4_Ht, w);
      h_ht_AK8_topana_S->Fill(AK8_Ht, w);
      h_jet1_pt_topana_S->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_S->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_S->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_S->Fill(data.evt.MR, w);
      h_R_topana_S->Fill(data.evt.R, w);
      h_MTR_topana_S->Fill(data.evt.MTR, w);
      h_R2_topana_S->Fill(data.evt.R2, w);
      h_tau21_topana_S->Fill(tau21.at(0),w);
      h_MET_topana_S->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_S->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_topana_S->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_S->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_S->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['a'];
    if (apply_all_cuts('a')) {
      h_ht_AK4_topana_s->Fill(AK4_Ht, w);
      h_ht_AK8_topana_s->Fill(AK8_Ht, w);
      h_jet1_pt_topana_s->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_s->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_s->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_s->Fill(data.evt.MR, w);
      h_R_topana_s->Fill(data.evt.R, w);
      h_MTR_topana_s->Fill(data.evt.MTR, w);
      h_R2_topana_s->Fill(data.evt.R2, w);
      h_tau21_topana_s->Fill(tau21.at(0),w);
      h_MET_topana_s->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_s->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_topana_s->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_s->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_s->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['X'];
    if (apply_all_cuts('X')) {
      h_ht_AK4_topana_Q->Fill(AK4_Ht, w);
      h_ht_AK8_topana_Q->Fill(AK8_Ht, w);
      h_jet1_pt_topana_Q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_Q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_Q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_Q->Fill(data.evt.MR, w);
      h_R_topana_Q->Fill(data.evt.R, w);
      h_MTR_topana_Q->Fill(data.evt.MTR, w);
      h_R2_topana_Q->Fill(data.evt.R2, w);
      h_tau21_topana_Q->Fill(tau21.at(0),w);
      h_MET_topana_Q->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_Q->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_topana_Q->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_Q->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_Q->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['x'];
    if (apply_all_cuts('x')) {
      h_ht_AK4_topana_q->Fill(AK4_Ht, w);
      h_ht_AK8_topana_q->Fill(AK8_Ht, w);
      h_jet1_pt_topana_q->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_q->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_q->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_q->Fill(data.evt.MR, w);
      h_R_topana_q->Fill(data.evt.R, w);
      h_MTR_topana_q->Fill(data.evt.MTR, w);
      h_R2_topana_q->Fill(data.evt.R2, w);
      h_tau21_topana_q->Fill(tau21.at(0),w);
      h_MET_topana_q->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_q->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_topana_q->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_q->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_q->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['V'];
    if (apply_all_cuts('V')) {
      h_ht_AK4_topana_T->Fill(AK4_Ht, w);
      h_ht_AK8_topana_T->Fill(AK8_Ht, w);
      h_jet1_pt_topana_T->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_T->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_T->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_T->Fill(data.evt.MR, w);
      h_R_topana_T->Fill(data.evt.R, w);
      h_MTR_topana_T->Fill(data.evt.MTR, w);
      h_R2_topana_T->Fill(data.evt.R2, w);
      h_tau21_topana_T->Fill(tau21.at(0),w);
      h_MET_topana_T->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_T->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_topana_T->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_T->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_T->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['C'];
    if (apply_all_cuts('C')) {
      h_ht_AK4_topana_W->Fill(AK4_Ht, w);
      h_ht_AK8_topana_W->Fill(AK8_Ht, w);
      h_jet1_pt_topana_W->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_W->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_W->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_W->Fill(data.evt.MR, w);
      h_R_topana_W->Fill(data.evt.R, w);
      h_MTR_topana_W->Fill(data.evt.MTR, w);
      h_R2_topana_W->Fill(data.evt.R2, w);
      h_tau21_topana_W->Fill(tau21.at(0),w);
      h_MET_topana_W->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_W->Fill(data.evt.MR, data.evt.R2, w);
      h_AK8_jet1_pt_topana_W->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_W->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_W->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['F'];
    if (apply_all_cuts('F')) {
      h_ht_AK4_topana_Z->Fill(AK4_Ht, w);
      h_ht_AK8_topana_Z->Fill(AK8_Ht, w);
      h_jet1_pt_topana_Z->Fill(data.jetsAK4.Pt[iJet[0]], w);
      h_jet2_pt_topana_Z->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_jet3_pt_topana_Z->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_topana_Z->Fill(data.evt.MR, w);
      h_R_topana_Z->Fill(data.evt.R, w);
      h_MTR_topana_Z->Fill(data.evt.MTR, w);
      h_R2_topana_Z->Fill(R2_ll, w);
      h_tau21_topana_Z->Fill(tau21.at(0),w);
      h_MET_topana_Z->Fill(data.met.Pt.at(0),w);
      h_R2_MR_topana_Z->Fill(data.evt.MR, R2_ll, w);
      h_AK8_jet1_pt_topana_Z->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      h_AK8_jet1_eta_topana_Z->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      h_MET_phi_topana_Z->Fill(data.met.Phi.at(0),w);
      }

    w = sf_weight['z'];
		if (apply_all_cuts('z') && nPhotonSelect==1) {
			if(std::abs(data.pho.Eta[iPhotonSelect[0]])<=1.479) {
       if(nWMassTag == 0) {
         if(nDirectPromptPhoton > 0){
           h_MR_0mW_EB_DPP->Fill(data.evt.MR,w);
           h_R2_0mW_EB_DPP->Fill(R2_pho,w);
           h_R2_MR_0mW_EB_DPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nDirectPromptPhoton == 0 && nFragmentationPromptPhoton > 0){
           h_MR_0mW_EB_FPP->Fill(data.evt.MR,w);
           h_R2_0mW_EB_FPP->Fill(R2_pho,w);
           h_R2_MR_0mW_EB_FPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nFakePhoton > 0){
           h_MR_0mW_EB_FP->Fill(data.evt.MR,w);
           h_R2_0mW_EB_FP->Fill(R2_pho,w);
           h_R2_MR_0mW_EB_FP->Fill(data.evt.MR,R2_pho,w);
         }
       }
       else if(nWMassTag == 1) {
         if(nDirectPromptPhoton > 0){
           h_MR_1mW_EB_DPP->Fill(data.evt.MR,w);
           h_R2_1mW_EB_DPP->Fill(R2_pho,w);
           h_R2_MR_1mW_EB_DPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nDirectPromptPhoton == 0 && nFragmentationPromptPhoton > 0){
           h_MR_1mW_EB_FPP->Fill(data.evt.MR,w);
           h_R2_1mW_EB_FPP->Fill(R2_pho,w);
           h_R2_MR_1mW_EB_FPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nFakePhoton > 0){
           h_MR_1mW_EB_FP->Fill(data.evt.MR,w);
           h_R2_1mW_EB_FP->Fill(R2_pho,w);
           h_R2_MR_1mW_EB_FP->Fill(data.evt.MR,R2_pho,w);
         }
       }
      }
			if(std::abs(data.pho.Eta[iPhotonSelect[0]])>1.479) {
       if(nWMassTag == 0) {
         if(nDirectPromptPhoton > 0){
           h_MR_0mW_EE_DPP->Fill(data.evt.MR,w);
           h_R2_0mW_EE_DPP->Fill(R2_pho,w);
           h_R2_MR_0mW_EE_DPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nDirectPromptPhoton == 0 && nFragmentationPromptPhoton > 0){
           h_MR_0mW_EE_FPP->Fill(data.evt.MR,w);
           h_R2_0mW_EE_FPP->Fill(R2_pho,w);
           h_R2_MR_0mW_EE_FPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nFakePhoton > 0){
           h_MR_0mW_EE_FP->Fill(data.evt.MR,w);
           h_R2_0mW_EE_FP->Fill(R2_pho,w);
           h_R2_MR_0mW_EE_FP->Fill(data.evt.MR,R2_pho,w);
         }
       else if(nWMassTag == 1) {
         if(nDirectPromptPhoton > 0){
           h_MR_1mW_EE_DPP->Fill(data.evt.MR,w);
           h_R2_1mW_EE_DPP->Fill(R2_pho,w);
           h_R2_MR_1mW_EE_DPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nDirectPromptPhoton == 0 && nFragmentationPromptPhoton > 0){
           h_MR_1mW_EE_FPP->Fill(data.evt.MR,w);
           h_R2_1mW_EE_FPP->Fill(R2_pho,w);
           h_R2_MR_1mW_EE_FPP->Fill(data.evt.MR,R2_pho,w);
         }
         else if(nFakePhoton > 0){
           h_MR_1mW_EE_FP->Fill(data.evt.MR,w);
           h_R2_1mW_EE_FP->Fill(R2_pho,w);
           h_R2_MR_1mW_EE_FP->Fill(data.evt.MR,R2_pho,w);
         }
       }
     }
    }
  }

    if(apply_all_cuts('z') && nPhotonSelect==1) {
      if(std::abs(data.pho.Eta[iPhotonSelect[0]])<=1.479) {
       h_nmW_zinv_EB->Fill(nWMassTag,w);
       h_MR_zinv_EB->Fill(data.evt.MR,w);
       h_R2_zinv_EB->Fill(R2_pho,w);
       h_R2_MR_zinv_EB->Fill(data.evt.MR,R2_pho,w);
       if(nWMassTag == 0) {
       h_MR_zinv_0mW_EB->Fill(data.evt.MR,w);
       h_R2_zinv_0mW_EB->Fill(R2_pho,w);
       h_R2_MR_zinv_0mW_EB->Fill(data.evt.MR,R2_pho,w);
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_SR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
        }
       else if(nHadTopMassTag >= 1) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_SR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
				}
       else if(nWMassTag >= 1) {
       h_MR_zinv_1mW_EB->Fill(data.evt.MR,w);
       h_R2_zinv_1mW_EB->Fill(R2_pho,w);
       h_R2_MR_zinv_1mW_EB->Fill(data.evt.MR,R2_pho,w);
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_SR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
        }
        if(R2_pho >= 0.08 && R2_pho < 0.12){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.12 && R2_pho < 0.16){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.16 && R2_pho < 0.24){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.24 && R2_pho < 0.50){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.50 && R2_pho < 1.00){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_SR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
       }
      else if(std::abs(data.pho.Eta[iPhotonSelect[0]])>1.479) {
       h_nmW_zinv_EE->Fill(nWMassTag,w);
       h_MR_zinv_EE->Fill(data.evt.MR,w);
       h_R2_zinv_EE->Fill(R2_pho,w);
       h_R2_MR_zinv_EE->Fill(data.evt.MR,R2_pho,w);
       if(nWMassTag == 0) {
       h_MR_zinv_0mW_EE->Fill(data.evt.MR,w);
       h_R2_zinv_0mW_EE->Fill(R2_pho,w);
       h_R2_MR_zinv_0mW_EE->Fill(data.evt.MR,R2_pho,w);
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_SR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
       }
       else if(nHadTopMassTag >= 1) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_SR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
       }
       else if(nWMassTag >= 1) {
       h_MR_zinv_1mW_EE->Fill(data.evt.MR,w);
       h_R2_zinv_1mW_EE->Fill(R2_pho,w);
       h_R2_MR_zinv_1mW_EE->Fill(data.evt.MR,R2_pho,w);
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_SR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         }
        }
        if(R2_pho >= 0.08 && R2_pho < 0.12){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.12 && R2_pho < 0.16){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.16 && R2_pho < 0.24){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.24 && R2_pho < 0.50){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
        else if(R2_pho >= 0.50 && R2_pho < 1.00){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_SR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonSelect[0]],w);
        }
       }
      }

    if (apply_all_cuts('z') && nPhotonFake==1) {
      if(std::abs(data.pho.Eta[iPhotonFake[0]])<=1.479) {
       if(nWMassTag == 0) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_CR_0mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
        }
       else if(nHadTopMassTag >= 1) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_CR_1mTop_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
       else if(nWMassTag == 1) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_CR_1mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
        }
        if(R2_pho >= 0.08 && R2_pho < 0.12){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.12 && R2_pho < 0.16){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.16 && R2_pho < 0.24){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.24 && R2_pho < 0.50){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.50 && R2_pho < 1.00){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_CR_mW_EB->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
       }
      else if(std::abs(data.pho.Eta[iPhotonFake[0]])>1.479) {
       if(nWMassTag == 0) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_CR_0mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
        }
       else if(nHadTopMassTag >= 1) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_CR_1mTop_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
       else if(nWMassTag >= 1) {
         if(data.evt.MR >= 800 && data.evt.MR < 1000)  h_ChargedIso_MR_900_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12) h_ChargedIso_R2_010_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.12 && R2_pho < 0.16) h_ChargedIso_R2_014_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.16 && R2_pho < 0.24) h_ChargedIso_R2_020_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.24 && R2_pho < 0.50) h_ChargedIso_R2_037_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.50 && R2_pho < 1.00) h_ChargedIso_R2_075_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         if(R2_pho >= 0.08 && R2_pho < 0.12){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.12 && R2_pho < 0.16){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.16 && R2_pho < 0.24){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.24 && R2_pho < 0.50){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
         else if(R2_pho >= 0.50 && R2_pho < 1.00){
          if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
          else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_CR_1mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         }
        }
        if(R2_pho >= 0.08 && R2_pho < 0.12){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_010_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_010_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_010_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_010_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_010_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.12 && R2_pho < 0.16){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_014_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_014_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_014_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_014_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_014_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.16 && R2_pho < 0.24){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_020_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_020_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_020_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_020_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_020_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.24 && R2_pho < 0.50){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_037_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_037_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_037_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_037_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_037_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
        else if(R2_pho >= 0.50 && R2_pho < 1.00){
         if(data.evt.MR >= 800 && data.evt.MR < 1000) h_ChargedIso_MR_900_R2_075_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1000 && data.evt.MR < 1200) h_ChargedIso_MR_1100_R2_075_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1200 && data.evt.MR < 1600) h_ChargedIso_MR_1400_R2_075_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 1600 && data.evt.MR < 2000) h_ChargedIso_MR_1800_R2_075_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
         else if(data.evt.MR >= 2000 && data.evt.MR < 4000) h_ChargedIso_MR_3000_R2_075_CR_mW_EE->Fill(ChargedHadronIsoEACorr[iPhotonFake[0]],w);
        }
       }
      }


    // Trigger efficiencies
    // No weighting
 
    // I guess these histots are for data only (trigger efficiencies)
    // One could just use the normal weight  for that ( = 1)
    // Since MC already has trigger efficiency weight applied already
    // one could use simply 1 as the weight there also
    // N-1 weights are not currently supported

    //w = isData ? 1 : sf_weight['S'];
    w = 1;


bool iscomFFsim = TString(sample).Contains("TTJets_madgraphMLM");
  if(iscomFFsim){ 
    for (size_t i=0; i<data.jetsAK8.size; ++i) {
    if (apply_cut('g', "genW" )){
      h_tau1_GenW_no_W->Fill(data.jetsAK8.tau1Puppi[i],w);
      h_tau2_GenW_no_W->Fill(data.jetsAK8.tau2Puppi[i],w);
      h_tau3_GenW_no_W->Fill(data.jetsAK8.tau3Puppi[i],w);
      h_tau21_GenW_no_W->Fill(tau21[i],w);
      h_tau32_GenW_no_W->Fill(tau32[i],w);
      h_AK8Jet1Pt_GenW_no_W->Fill(data.jetsAK8.Pt[i],w);
      h_AK8Jet1Pt_eta_GenW_no_W->Fill(data.jetsAK8.Pt[i],std::abs(data.jetsAK8.Eta[i]),w);
    }
    if (apply_cut('B', "genTop" )){
      h_tau1_GenTop_no_Top->Fill(data.jetsAK8.tau1Puppi[i],w);
      h_tau2_GenTop_no_Top->Fill(data.jetsAK8.tau2Puppi[i],w);
      h_tau3_GenTop_no_Top->Fill(data.jetsAK8.tau3Puppi[i],w);
      h_tau21_GenTop_no_Top->Fill(tau21[i],w);
      h_tau32_GenTop_no_Top->Fill(tau32[i],w);
      h_SubjetBTag_GenTop_no_Top->Fill(data.jetsAK8.maxSubjetCSVv2[i],w);
      h_AK8Jet1Pt_GenTop_no_Top->Fill(data.jetsAK8.Pt[i],w);
      h_AK8Jet1Pt_eta_GenTop_no_Top->Fill(data.jetsAK8.Pt[i],std::abs(data.jetsAK8.Eta[i]),w);
    }
    }
    if (apply_all_cuts_except('g', "WTag")){
    for (size_t i=0; i<iGenMassW.size(); ++i) {
      h_tau1_GenW_W->Fill(data.jetsAK8.tau1Puppi[iGenMassW[i]],w);
      h_tau2_GenW_W->Fill(data.jetsAK8.tau2Puppi[iGenMassW[i]],w);
      h_tau3_GenW_W->Fill(data.jetsAK8.tau3Puppi[iGenMassW[i]],w);
      h_tau21_GenW_W->Fill(tau21[iGenMassW[i]],w);
      h_tau32_GenW_W->Fill(tau32[iGenMassW[i]],w);
    }
    }
    if (apply_all_cuts('g')){
    for (size_t i=0; i<iGenHadW.size(); ++i) {
      h_AK8Jet1Pt_GenW_W->Fill(data.jetsAK8.Pt[iGenHadW[i]],w);
      h_AK8Jet1Pt_eta_GenW_W->Fill(data.jetsAK8.Pt[iGenHadW[i]],std::abs(data.jetsAK8.Eta[iGenHadW[i]]),w);
    }
    }
    if (apply_all_cuts_except('B', "TopTag")){
    for (size_t i=0; i<iGenMassTop.size(); ++i) {
      h_tau1_GenTop_Top->Fill(data.jetsAK8.tau1Puppi[iGenMassTop[i]],w);
      h_tau2_GenTop_Top->Fill(data.jetsAK8.tau2Puppi[iGenMassTop[i]],w);
      h_tau3_GenTop_Top->Fill(data.jetsAK8.tau3Puppi[iGenMassTop[i]],w);
      h_tau21_GenTop_Top->Fill(tau21[iGenMassTop[i]],w);
      h_tau32_GenTop_Top->Fill(tau32[iGenMassTop[i]],w);
      h_SubjetBTag_GenTop_Top->Fill(data.jetsAK8.maxSubjetCSVv2[iGenMassTop[i]],w);
    }
    }
    if (apply_all_cuts('B')){
    for (size_t i=0; i<iGenTop.size(); ++i) {
      h_AK8Jet1Pt_GenTop_Top->Fill(data.jetsAK8.Pt[iGenTop[i]],w);
      h_AK8Jet1Pt_eta_GenTop_Top->Fill(data.jetsAK8.Pt[iGenTop[i]],std::abs(data.jetsAK8.Eta[iGenTop[i]]),w);
    }
    }
  }
    
  w=1;

  bool istriggerData = TString(sample).Contains("SingleElectron");
  if(istriggerData && (data.hlt.Ele23_WPLoose_Gsf==1||data.hlt.Ele27_WPTight_Gsf==1)&&nEleTight==1&&nMuVeto==0){ 
    h_HT_pre->Fill(AK4_Ht,w);
    h_HT_j1pt_pre->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
  
    if(apply_cut('S',"HLT")){
      h_HT_pre_pass->Fill(AK4_Ht,w);
      h_HT_j1pt_pre_pass->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
  }

  istriggerData = TString(sample).Contains("SingleMuon");
  if(istriggerData && (data.hlt.IsoMu24==1||data.hlt.IsoTkMu24==1)&&nMuTight==1&&nEleVeto==0){ 
    h_HT_pre->Fill(AK4_Ht,w);
    h_HT_j1pt_pre->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
  
    if(apply_cut('S',"HLT")){
      h_HT_pre_pass->Fill(AK4_Ht,w);
      h_HT_j1pt_pre_pass->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
  }

  istriggerData = TString(sample).Contains("MET");
  if(istriggerData && data.hlt.PFMET120_PFMHT120_IDTight==1&&data.met.Pt[0]>200&&nLepVeto==0&&data.evt.NIsoTrk==0){ 
    h_HT_pre->Fill(AK4_Ht,w);
    h_HT_j1pt_pre->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
  
    if(apply_cut('S',"HLT")){
      h_HT_pre_pass->Fill(AK4_Ht,w);
      h_HT_j1pt_pre_pass->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
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
    h_njet_pre[syst_index]->Fill(nJet,        w);
    h_nb_pre[syst_index]->Fill(nMediumBTag, w);
    h_nw_pre[syst_index]->Fill(nTightWTag,  w);
    h_j1_pt_pre[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]],w);
    h_j2_pt_pre[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
    h_j3_pt_pre[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
    h_MR_pre[syst_index]->Fill(data.evt.MR, w);
    h_R2_pre[syst_index]->Fill(data.evt.R2, w);
    h_tau21_pre[syst_index]->Fill(tau21.at(0),w);
    h_MET_pre[syst_index]->Fill(data.met.Pt.at(0),w);
    h_softDropMass_pre[syst_index]->Fill(softDropMassW.at(0),w);
    h_HT_pre[syst_index]->Fill(AK4_Ht,w);
    h_HT_j1pt_pre[syst_index]->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);

    if(pass){
      h_njet_pre_pass[syst_index]->Fill(nJet,        w);
      h_nb_pre_pass[syst_index]->Fill(nMediumBTag, w);
      h_nw_pre_pass[syst_index]->Fill(nTightWTag,  w);
      h_j1 pt_pre_pass[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]],w);
      h_j2_pt_pre_pass[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      h_j3_pt_pre_pass[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      h_MR_pre_pass[syst_index]->Fill(data.evt.MR, w);
      h_R2_pre_pass[syst_index]->Fill(data.evt.R2, w);
      h_tau21_pre_pass[syst_index]->Fill(tau21.at(0),w);
      h_MET_pre_pass[syst_index]->Fill(data.met.Pt.at(0),w);
      h_softDropMass_pre_pass[syst_index]->Fill(softDropMassW.at(0),w);
      h_HT_pre_pass[syst_index]->Fill(AK4_Ht,w);
      h_HT_j1pt_pre_pass[syst_index]->Fill(AK4_Ht,data.jetsAK8.Pt[iJetAK8[0]],w);
    }
*/

  // Vary systematics and save each variation into a different historgam
  // Switch on settings.varySystematics to be effective
  //if (apply_all_cuts('S')) vh_jet1_pt[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);

    w = sf_weight['S'];
    if (apply_all_cuts('S')) {
      vh_ht_AK4_S[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_S[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_S[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_S[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_S[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_S[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_S[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_S[syst_index]->Fill(tau21.at(0),w);
      vh_MET_S[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_S[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_S[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_S[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_S[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_S_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_S_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_S_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_S_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_S_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_S_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_S_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_S_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_S_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_S_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_S_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_S_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_S_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_S_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_S_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_S_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_S_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_S_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_S_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_S_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_S_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_S_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_S_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_S_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_S_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_S_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['s'];
    if (apply_all_cuts('s')) {
      vh_ht_AK4_s[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_s[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_s[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_s[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_s[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_s[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_s[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_s[syst_index]->Fill(tau21.at(0),w);
      vh_MET_s[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_s[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_s[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_s[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_s[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_s_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_s_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_s_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_s_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_s_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_s_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_s_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_s_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_s_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_s_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_s_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_s_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_s_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_s_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_s_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_s_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_s_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_s_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_s_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_s_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_s_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_s_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_s_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_s_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_s_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_s_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['Q'];
    if (apply_all_cuts('Q')) {
      vh_ht_AK4_Q[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_Q[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_Q[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_Q[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_Q[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_Q[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_Q[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_Q[syst_index]->Fill(tau21.at(0),w);
      vh_MET_Q[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_Q[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_Q[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_Q[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_Q[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_Q_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_Q_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_Q_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_Q_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_Q_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_Q_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_Q_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_Q_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_Q_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_Q_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_Q_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_Q_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_Q_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_Q_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_Q_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_Q_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_Q_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_Q_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_Q_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_Q_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_Q_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_Q_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_Q_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_Q_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_Q_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_Q_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['q'];
    if (apply_all_cuts('q')) {
      vh_ht_AK4_q[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_q[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_q[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_q[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_q[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_q[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_q[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_q[syst_index]->Fill(tau21.at(0),w);
      vh_MET_q[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_q[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_q[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_q[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_q[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_q_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_q_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_q_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_q_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_q_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_q_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_q_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_q_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_q_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_q_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_q_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_q_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_q_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_q_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_q_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_q_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_q_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_q_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_q_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_q_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_q_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_q_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_q_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_q_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_q_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_q_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['T'];
    if (apply_all_cuts('T')) {
      vh_ht_AK4_T[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_T[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_T[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_T[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_T[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_T[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_T[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_T[syst_index]->Fill(tau21.at(0),w);
      vh_MET_T[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_T[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_T[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_T[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_T[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_T_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_T_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_T_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_T_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_T_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_T_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_T_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_T_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_T_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_T_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_T_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_T_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_T_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_T_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_T_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_T_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_T_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_T_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_T_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_T_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_T_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_T_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_T_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_T_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_T_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_T_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['W'];
    if (apply_all_cuts('W')) {
      vh_ht_AK4_W[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_W[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_W[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_W[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_W[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_W[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_W[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_W[syst_index]->Fill(tau21.at(0),w);
      vh_MET_W[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_W[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_W[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_W[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_W[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_W_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_W_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_W_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_W_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_W_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_W_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_W_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_W_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_W_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_W_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_W_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_W_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_W_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_W_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_W_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_W_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_W_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_W_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_W_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_W_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_W_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_W_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_W_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_W_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_W_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_W_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['Z'];
    if (apply_all_cuts('Z')) {
      vh_ht_AK4_Z[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_Z[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_Z[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_Z[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_Z[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_Z[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_Z[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_Z[syst_index]->Fill(tau21.at(0),w);
      vh_MET_Z[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_Z[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_Z[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_Z[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_Z[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_Z_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_Z_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_Z_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_Z_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_Z_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_Z_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_Z_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_Z_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_Z_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_Z_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_Z_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_Z_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_Z_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_Z_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_Z_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_Z_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_Z_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_Z_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_Z_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_Z_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_Z_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_Z_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_Z_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_Z_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_Z_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_Z_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }

    w = sf_weight['G'];
    if (apply_all_cuts('G')) {
      vh_ht_AK4_G[syst_index]->Fill(AK4_Ht, w);
      vh_jet1_pt_G[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
      vh_jet2_pt_G[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
      vh_jet3_pt_G[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
      vh_MR_G[syst_index]->Fill(data.evt.MR, w);
      vh_MTR_G[syst_index]->Fill(data.evt.MTR, w);
      vh_R2_G[syst_index]->Fill(data.evt.R2, w);
      vh_tau21_G[syst_index]->Fill(tau21.at(0),w);
      vh_MET_G[syst_index]->Fill(data.met.Pt.at(0),w);
      vh_R2_MR_G[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_G[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_G[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_G[syst_index]->Fill(data.met.Phi.at(0),w);
      if(nJet>=3 && nJet<6){
        vh_ht_AK4_G_nj35[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_G_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_G_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_G_nj35[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_G_nj35[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_G_nj35[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_G_nj35[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_G_nj35[syst_index]->Fill(tau21.at(0),w);
        vh_MET_G_nj35[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_G_nj35[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_G_nj35[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_G_nj35[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_G_nj35[syst_index]->Fill(data.met.Phi.at(0),w);
      }
      if(nJet>=6){
        vh_ht_AK4_G_nj6[syst_index]->Fill(AK4_Ht, w);
        vh_jet1_pt_G_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[0]], w);
        vh_jet2_pt_G_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[1]], w);
        vh_jet3_pt_G_nj6[syst_index]->Fill(data.jetsAK4.Pt[iJet[2]], w);
        vh_MR_G_nj6[syst_index]->Fill(data.evt.MR, w);
        vh_MTR_G_nj6[syst_index]->Fill(data.evt.MTR, w);
        vh_R2_G_nj6[syst_index]->Fill(data.evt.R2, w);
        vh_tau21_G_nj6[syst_index]->Fill(tau21.at(0),w);
        vh_MET_G_nj6[syst_index]->Fill(data.met.Pt.at(0),w);
        vh_R2_MR_G_nj6[syst_index]->Fill(data.evt.MR, data.evt.R2, w);
      vh_AK8_jet1_pt_G_nj6[syst_index]->Fill(data.jetsAK8.Pt[iJetAK8[0]], w);
      vh_AK8_jet1_eta_G_nj6[syst_index]->Fill(data.jetsAK8.Eta[iJetAK8[0]], w);
      vh_MET_phi_G_nj6[syst_index]->Fill(data.met.Phi.at(0),w);
      }
    }
}

// Methods used by SmartHistos (Plotter)
// Can leave them empty
void
Analysis::define_histo_options(const double& w, const DataStruct& d, const unsigned int& syst_nSyst,
			       const unsigned int& syst_index, bool runOnSkim=false)
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
