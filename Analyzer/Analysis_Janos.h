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

// Cut variables
static size_t cut_index;
std::map<char, unsigned int> cutbits;

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
  }
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

  // Signal skim
  //return signal_selection(data);
}

//_______________________________________________________
//          Define Analysis event selection cuts
//     Can define all sorts of Signal/Control regions

enum SCuts { S_3Jet, S_MR_R2, S_1W,    S_1b, S_0Ele, S_0Mu, S_mDPhi};
enum WCuts { W_3Jet, W_MR_R2, W_1Wpre, W_0b, W_1Lep,        W_mDPhi,    W_MT};
enum TCuts { T_3Jet, T_MR_R2, T_1W,    T_1b, T_1Lep,        T_mDPhi,    T_MT};
enum QCuts { Q_3Jet, Q_MR_R2, Q_1W,    Q_0b, Q_0Ele, Q_0Mu, Q_mDPhiInv};

void
Analysis::define_selections(const DataStruct& d)
{
  analysis_cuts.clear();

  // Define here cuts that are common in all Signal/Control regions
  // MET Filters, etc. are already applied in AnalysisBase.h, See baseline_cuts
  baseline_cuts.push_back({ .name="Skim_1JetAK8",    .func = []    { return nJetAK8>=1;                      }}); // Similar to pt>200, one AK8 jet has pt>170
  baseline_cuts.push_back({ .name="Skim_2Jet",       .func = []    { return nJet>=2;                         }});
  
  // S: Signal region
  analysis_cuts['S'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                         }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['S'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08; }});
  analysis_cuts['S'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                   }});
  analysis_cuts['S'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                  }});
  analysis_cuts['S'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['S'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
  //analysis_cuts['S'].push_back({ .name="0TauTrk",    .func = []    { return;  }});
  //analysis_cuts['S'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['S'].push_back({ .name="mDPhi>=0p4", .func = []    { return minDeltaPhi>=0.4;                }}); // Decreased it to the AK4 cone size (from 0.5)
  
  // W: W enriched control sample
  analysis_cuts['W'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                         }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['W'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08; }});
  analysis_cuts['W'].push_back({ .name="1Wpre",      .func = []    { return nWPreTag>=1;                     }});
  analysis_cuts['W'].push_back({ .name="0b",         .func = []    { return nLooseBTag==1;                   }});
  analysis_cuts['W'].push_back({ .name="1Lep",       .func = []    { return nLepTight==1;                    }});
  //analysis_cuts['W'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['W'].push_back({ .name="mDPhi>=0p4",  .func = []   { return minDeltaPhi>=0.4;                }}); // Decreased it to the AK4 cone size (from 0.5)
  analysis_cuts['W'].push_back({ .name="30<=MT<100",  .func = []   { return MT>=30 && MT<100;                }});
  
  // T: Top enriched control sample
  analysis_cuts['T'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                         }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['T'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08; }});
  analysis_cuts['T'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                   }});
  analysis_cuts['T'].push_back({ .name="1b",         .func = []    { return nMediumBTag>=1;                  }});
  analysis_cuts['T'].push_back({ .name="1Lep",       .func = []    { return nLepTight==1;                    }});
  //analysis_cuts['T'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['T'].push_back({ .name="mDPhi>=0p4", .func = []    { return minDeltaPhi>=0.4;                }}); // Decreased it to the AK4 cone size (from 0.5)
  analysis_cuts['T'].push_back({ .name="MT<100",     .func = []    { return MT<100;                          }});
  
  // Q: QCD enriched control sample
  analysis_cuts['Q'].push_back({ .name="3Jet",       .func = []    { return nJet>=3;                         }}); // Separate cut, so one can exclude (N-1)
  analysis_cuts['Q'].push_back({ .name="MR_R2",      .func = [&d]  { return d.evt.MR>=800 && d.evt.R2>=0.08; }});
  analysis_cuts['Q'].push_back({ .name="1W",         .func = []    { return nTightWTag>=1;                   }});
  analysis_cuts['Q'].push_back({ .name="0b",         .func = []    { return nLooseBTag==0;                   }});
  analysis_cuts['Q'].push_back({ .name="0Ele",       .func = []    { return nEleVeto==0;                     }});
  analysis_cuts['Q'].push_back({ .name="0Mu",        .func = []    { return nMuVeto==0;                      }});
  //analysis_cuts['Q'].push_back({ .name="0TauTrk",    .func = []    { return;  }});
  //analysis_cuts['Q'].push_back({ .name="mDPhiHat",   .func = []    { return;  }});
  analysis_cuts['Q'].push_back({ .name="mDPhi<0.25", .func = []    { return minDeltaPhi<0.25;                }}); // Decreased it to 0.25 (from 0.3)
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!

bool
Analysis::signal_selection(const DataStruct& data) {
  return apply_all_cuts('S');
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
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst, const unsigned int& syst_index, std::string dirname, bool runOnSkim=false)
{
  sh.SetHistoWeights({ [&weight] { return weight; } });

  // Keep this to be able to use analysis cuts
  define_selections(d);

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
			   "TTJets_HT-600to800_ext1", "TTJets_HT-800to1200_ext1", "TTJets_HT-1200to2500_ext1", "TTJets_HT-2500toInf" 
			 } });
  
  bkg_ttbars.push_back({ .postfix="TTJets_madgraph_FullSim",  .legend="t#bar{t} (madgraphMLM, FullSim)", .color="901",/*Pink*/   .dirs={ "TTJets_madgraphMLM-pythia8" } });
  bkg_ttbars.push_back({ .postfix="TTJets_madgraph_FastSim",  .legend="t#bar{t} (madgraphMLM, FastSim)", .color="903",/*DPink*/  .dirs={ "TTJets_madgraphMLM_FastSim" } });
  bkg_ttbars.push_back({ .postfix="TTJets_amcatnlo",          .legend="t#bar{t} (aMC@NLO FxFx)",         .color="617",/*Magent*/ .dirs={ "TTJets_amcatnloFXFX-pythia8_reHLT" } });

  bkg_ttbars.push_back({ .postfix="TT_amcatnlo",              .legend="t#bar{t} (aMC@NLO)",              .color="619",/*DMagen*/ .dirs={ "TT_amcatnlo-pythia8_ext1" } });
  bkg_ttbars.push_back({ .postfix="TT_powheg_pythia8",        .legend="t#bar{t} (powheg, pythia8)",      .color="633",/*Red*/    .dirs={ 
			   /*"TT_powheg-pythia8_evtgen", "TT_powheg-pythia8_ext3_reHLT",*/ "TT_powheg-pythia8_ext4" // Use only last one, it has 182M events
			 } });
  bkg_ttbars.push_back({ .postfix="TT_powheg_herwigpp",       .legend="t#bar{t} (powheg, herwigpp)",     .color="803",/*DOran*/  .dirs={ "TT_powheg-herwigpp" } });
  std::vector<Sample> bkg_ttbar_hlt;
  bkg_ttbar_hlt.push_back({ .postfix="TT_powheg_pythia8",        .legend="t#bar{t} (powheg, pythia8)",   .color="633",/*Red*/    .dirs={ "TT_powheg-pythia8_ext3_reHLT" } });

  std::vector<Sample> bkg_nonttbars;
  //bkg_nonttbars.push_back({ .postfix="QCD",         .legend="QCD",                 .color=  "4",/*Blue*/    
  //      		      .dirs={ 
  //      		      "QCD_HT100to200",  "QCD_HT200to300",   "QCD_HT300to500",   "QCD_HT500to700",
  //      		      "QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"
  //			    } });
  //bkg_nonttbars.push_back({ .postfix="VToQQ",  .legend="W/Z/#gamma#rightarrowqq",  .color="862",/*Azure*/ .dirs={ 
  //      		      "ZJetsToQQ_HT600toInf_reHLT", "WJetsToQQ_HT180", "DYJetsToQQ_HT180"
  //			    } });
  //bkg_nonttbars.push_back({ .postfix="TTX",         .legend="t#bar{t} + W/Z/#gamma/t#bar{t}", .color="803",/*Brown*/   .dirs={ 
  //      		      "TTWJetsToLNu", "TTWJetsToQQ",
  //      		      "TTZToLLNuNu", "TTZToQQ",
  //      		      "TTGJets",
  //      		      "TTTT"
  //      		    } });
  //bkg_nonttbars.push_back({ .postfix="Diboson",     .legend="Diboson",            .color="804",/*DOrange*/ .dirs={ 
  //      		      "WWTo2L2Nu", "WWToLNuQQ", "WWToLNuQQ_ext1", "WWTo4Q",
  //      		      "WZTo3LNu", "WZTo2L2Q", "WZTo1L1Nu2Q", "WZTo1L3Nu",
  //      		      "ZZTo4L", "ZZTo2L2Nu", "ZZTo2L2Q", "ZZTo2Q2Nu", "ZZTo4Q"
  //      		    } });
  //bkg_nonttbars.push_back({ .postfix="Triboson",     .legend="Triboson",          .color="805",/*DOrange*/ .dirs={ "WWW", "WWZ", "WZZ", "ZZZ" } });

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

  bkg_nonttbars.push_back({ .postfix="Multijet",   .legend="multijet (QCD, W/Z/#gamma#rightarrowqq)", .color="619",/*DMagen*/ .dirs={ 
			      "QCD_HT100to200",  "QCD_HT200to300",   "QCD_HT300to500",   "QCD_HT500to700",
			      "QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"
			      "ZJetsToQQ_HT600toInf_reHLT", "WJetsToQQ_HT180", "DYJetsToQQ_HT180" // V Multi
			      "WWTo4Q", "ZZTo4Q" // VV Multi
			    } });
  bkg_nonttbars.push_back({ .postfix="WToLNu",     .legend="W#rightarrowl#nu+jets",                   .color="418",/*Green*/  .dirs={ 
			      "WJetsToLNu_HT-100To200_ext1", "WJetsToLNu_HT-200To400", "WJetsToLNu_HT-400To600", "WJetsToLNu_HT-600To800"
			      "WJetsToLNu_HT-800To1200", "WJetsToLNu_HT-1200To2500_ext1", "WJetsToLNu_HT-2500ToInf",
			    } });
  bkg_nonttbars.push_back({ .postfix="ZToNuNu",    .legend="Z#rightarrow#nu#nu+jets",                 .color="401",/*Yellow*/ .dirs={ 
			      "ZJetsToNuNu_HT-100To200", "ZJetsToNuNu_HT-200To400", "ZJetsToNuNu_HT-400To600", "ZJetsToNuNu_HT-600To800", 
			      "ZJetsToNuNu_HT-800To1200", "ZJetsToNuNu_HT-1200To2500", "ZJetsToNuNu_HT-2500ToInf"
			    } });
  bkg_nonttbars.push_back({ .postfix="Multiboson", .legend="VV + VVV",                                .color="601",/*Blue*/   .dirs={ 
			      "WWTo2L2Nu", "WWToLNuQQ", "WWToLNuQQ_ext1",
        		      "WZTo3LNu", "WZTo2L2Q", "WZTo1L1Nu2Q", "WZTo1L3Nu",
        		      "ZZTo4L", "ZZTo2L2Nu", "ZZTo2L2Q", "ZZTo2Q2Nu",
			      "WWW", "WWZ", "WZZ", "ZZZ"
			    } });
  bkg_nonttbars.push_back({ .postfix="TTX",        .legend="t#bar{t} + W/Z/#gamma/t#bar{t}",          .color="843",/*DTeal*/  .dirs={ 
			      "TTWJetsToLNu", "TTWJetsToQQ",
			      "TTZToLLNuNu", "TTZToQQ",
			      "TTGJets",
			      "TTTT"
			    } });
  bkg_nonttbars.push_back({ .postfix="Top",        .legend="single top",                              .color="433",/*Cyan*/   .dirs={ 
			      "ST_s-channel_4f_leptonDecays",
			      "ST_t-channel_top_4f_inclusiveDecays", "ST_t-channel_antitop_4f_inclusiveDecays",
			      "ST_tW_top_5f_inclusiveDecays",        "ST_tW_antitop_5f_inclusiveDecays"
			    } });
  bkg_nonttbars.push_back({ .postfix="DYToLL",     .legend="Z/#gamma#rightarrowll",                   .color="803",/*Brown*/  .dirs={ 
			      "DYJetsToLL_M-5to50_HT-100to200", "DYJetsToLL_M-5to50_HT-200to400", "DYJetsToLL_M-5to50_HT-400to600", "DYJetsToLL_M-5to50_HT-600toInf",
			      "DYJetsToLL_M-50_HT-100to200",    "DYJetsToLL_M-50_HT-200to400",    "DYJetsToLL_M-50_HT-400to600",    "DYJetsToLL_M-50_HT-600toInf",
			    } });

  std::vector<Sample> bkg_all, bkg_selected;
  bkg_all.insert(bkg_all.end(), bkg_ttbars.begin(), bkg_ttbars.end());
  bkg_all.insert(bkg_all.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());
  bkg_selected.push_back(bkg_ttbars[5]);   // powheg
  bkg_selected.insert(bkg_selected.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());

  std::vector<Sample> data_all, data_selected, single_ele, single_mu, met;
  data_all.push_back({ .postfix="Data",      .legend="Data",             .color="1", .dirs={
			 "JetHT_Run2016B_PRv2", "JetHT_Run2016C_PRv2", "JetHT_Run2016D_PRv2", 
			 "JetHT_Run2016E_PRv2", "JetHT_Run2016F_PRv1", "JetHT_Run2016G_PRv1"
		       } });
  data_all.push_back({ .postfix="SingleEle", .legend="Data (SingleEle)", .color="1", .dirs={
			 "SingleElectron_Run2016B_PRv2", "SingleElectron_Run2016C_PRv2", "SingleElectron_Run2016D_PRv2", 
			 "SingleElectron_Run2016E_PRv2", "SingleElectron_Run2016F_PRv1", "SingleElectron_Run2016G_PRv1"
		       } });
  data_all.push_back({ .postfix="SingleMu",  .legend="Data (SingleMu)",  .color="1", .dirs={
			 "SingleMuon_Run2016B_PRv2", "SingleMuon_Run2016C_PRv2", "SingleMuon_Run2016D_PRv2", 
			 "SingleMuon_Run2016E_PRv2", "SingleMuon_Run2016F_PRv1", "SingleMuon_Run2016G_PRv1"
		       } });
  data_all.push_back({ .postfix="MET",       .legend="Data (MET)",       .color="1", .dirs={
			 "MET_Run2016B_PRv2", "MET_Run2016C_PRv2", "MET_Run2016D_PRv2", 
			 "MET_Run2016E_PRv2", "MET_Run2016F_PRv1", "MET_Run2016G_PRv1"
		       } });
  data_selected.push_back(data_all[0]);
  single_ele.push_back(data_all[1]);
  single_mu.push_back(data_all[2]);
  met.push_back(data_all[3]);

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
  for (int i=3; i<5; ++i) signal_stop   .push_back(signal_all[i]);
  
  //"T5ttttDeg (M_{#tilde{g}}=1TeV)","1",/*Black*/
  //"T1tttt (M_{#tilde{g}}=1.5TeV)", "862",/*Azure*/
  //"T1tttt (M_{#tilde{g}}=1.5TeV, M_{#tilde{#chi}^{0}}=100GeV)", "841",/*Teal*/     
  //"T1tttt (M_{#tilde{g}}=1.2TeV, M_{#tilde{#chi}^{0}}=800GeV)", "843",/*DarkTeal*/ 
  //"T5ttttDeg (M_{#tilde{g}}=1TeV, 2,3-body)",                   "12", /*DarkGrey*/ 
  //"T2ttDeg (M_{#tilde{t}}=350GeV)",                             "434",/*Cyan*/     

  // Sample postfixes
  static const PostfixOptions all_samples_opt=get_pf_opts_({data_all, bkg_all, signal_all}, dirname);
  sh.AddNewPostfix("AllSamples", [] { return all_samples_opt.index; }, all_samples_opt.postfixes, all_samples_opt.legends, all_samples_opt.colors);

  static const PostfixOptions plot_samples_opt=get_pf_opts_({data_selected, signal_fastsim, bkg_selected}, dirname);
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
  static const PostfixOptions background_opt = get_pf_opts_({background}, dirname);
  sh.AddNewPostfix("Background",  [] { return background_opt.index; }, background_opt.postfixes, background_opt.legends, background_opt.colors);

  static const PostfixOptions gluino_signalscans_opt = get_pf_opts_({signal_gluino}, dirname);
  sh.AddNewPostfix("GluinoSignalScans",  [] { return gluino_signalscans_opt.index; }, gluino_signalscans_opt.postfixes, gluino_signalscans_opt.legends, gluino_signalscans_opt.colors);

  static const PostfixOptions stop_signalscans_opt = get_pf_opts_({signal_stop}, dirname);
  sh.AddNewPostfix("StopSignalScans",  [] { return stop_signalscans_opt.index; }, stop_signalscans_opt.postfixes, stop_signalscans_opt.legends, stop_signalscans_opt.colors);

  static const PostfixOptions background_signal_opt = get_pf_opts_({background, signal_selected}, dirname);
  sh.AddNewPostfix("Background,Signal", [&d] { 
		     if (background_signal_opt.index==1) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1;
		     }
		     return background_signal_opt.index;
		   }, background_signal_opt.postfixes, background_signal_opt.legends, "633,601");

  static const PostfixOptions signals_opt = get_pf_opts_({signal_all}, dirname);
  sh.AddNewPostfix("Signals",  [&d] { 
		     // Select gluino/stop mass to give ~1k events with 40 fb^-1
		     if (signals_opt.index<3) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1; // T1ttt/T5ttcc/T5tttt
		     } else if (signals_opt.index==3) {
		       if (d.evt.SUSY_Stop_Mass  != 850 || d.evt.SUSY_LSP_Mass != 100) return (size_t)-1; // T2tt - Same as FullSim point
		     }
		     return signals_opt.index; 
		   }, signals_opt.postfixes, signals_opt.legends, signals_opt.colors);

  static const PostfixOptions signals_background_opt = get_pf_opts_({signal_all, background}, dirname);
  sh.AddNewPostfix("Signals,Background",  [&d] { 
		     // Select gluino/stop mass to give ~1k events with 40 fb^-1
		     if (signals_background_opt.index<3) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1; // T1ttt/T5ttcc/T5tttt
		     } else if (signals_background_opt.index==3) {
		       if (d.evt.SUSY_Stop_Mass  != 850 || d.evt.SUSY_LSP_Mass != 100) return (size_t)-1; // T2tt - Same as FullSim point
		     }
		     return signals_background_opt.index; 
		   }, signals_background_opt.postfixes, signals_background_opt.legends, signals_background_opt.colors);

  static const PostfixOptions signals_ttbar_opt = get_pf_opts_({signal_all, bkg_ttbars}, dirname);
  sh.AddNewPostfix("Signals,TT",  [&d] { 
		     // Select gluino/stop mass to give ~1k events with 40 fb^-1
		     if (signals_ttbar_opt.index<3) {
		       if (d.evt.SUSY_Gluino_Mass!=1400 || d.evt.SUSY_LSP_Mass != 300) return (size_t)-1; // T1ttt/T5ttcc/T5tttt
		     } else if (signals_ttbar_opt.index==3) {
		       if (d.evt.SUSY_Stop_Mass  != 850 || d.evt.SUSY_LSP_Mass != 100) return (size_t)-1; // T2tt - Same as FullSim point
		     }
		     return signals_ttbar_opt.index; 
		   }, signals_ttbar_opt.postfixes, signals_ttbar_opt.legends, signals_ttbar_opt.colors);

  static const PostfixOptions mgluinopoints_opt = get_pf_opts_({signal_gluino}, dirname);
  sh.AddNewPostfix("MGluinoPoints",  [&d] { 
		     if (mgluinopoints_opt.index==(size_t)-1) return (size_t)-1;
		     else {
		       if (d.evt.SUSY_LSP_Mass != 300) return (size_t)-1;
		       if      (d.evt.SUSY_Gluino_Mass== 800) return (size_t)0;
		       else if (d.evt.SUSY_Gluino_Mass==1000) return (size_t)1;
		       else if (d.evt.SUSY_Gluino_Mass==1200) return (size_t)2;
		       else if (d.evt.SUSY_Gluino_Mass==1400) return (size_t)3;
		       else if (d.evt.SUSY_Gluino_Mass==1600) return (size_t)4;
		       else return (size_t)-1;
		     }
		   }, "Mlsp300_Mglu[800to1600++200]", "M_{#tilde{#chi}^{0}}=300GeV, M_{#tilde{g}}=[0.8to1.6++0.2]TeV", col5_green_to_red);

  static const PostfixOptions mstoppoints_opt = get_pf_opts_({signal_stop}, dirname);
  sh.AddNewPostfix("MStopPoints",  [&d] { 
		     if (mstoppoints_opt.index==(size_t)-1) return (size_t)-1;
		     else {
		       if (d.evt.SUSY_LSP_Mass != 100) return (size_t)-1;
		       if      (d.evt.SUSY_Stop_Mass== 700) return (size_t)0;
		       else if (d.evt.SUSY_Stop_Mass== 800) return (size_t)1;
		       else if (d.evt.SUSY_Stop_Mass== 900) return (size_t)2;
		       else if (d.evt.SUSY_Stop_Mass==1000) return (size_t)3;
		       else if (d.evt.SUSY_Stop_Mass==1100) return (size_t)4;
		       else return (size_t)-1;
		     }
		   }, "Mlsp100_Mstop[800to1600++200]", "M_{#tilde{#chi}^{0}}=100GeV, M_{#tilde{t}}=[0.7to1.1++0.1]TeV", col5_green_to_red);

  static const PostfixOptions data_mc_opt = get_pf_opts_({data_selected, background}, dirname);
  sh.AddNewPostfix("Data,MC",  [] { return data_mc_opt.index; }, data_mc_opt.postfixes, data_mc_opt.legends, "1,633");

  static const PostfixOptions single_lep_opt = get_pf_opts_({single_ele, single_mu}, dirname);
  sh.AddNewPostfix("SingleEle,SingleMu", [] { return single_lep_opt.index; }, single_lep_opt.postfixes, single_lep_opt.legends, "1,633");

  static const PostfixOptions triggers_opt = get_pf_opts_({data_selected, single_ele, single_mu, met, bkg_ttbar_hlt}, dirname);
  sh.AddNewPostfix("Triggers", [&d]()
                    {
                      if (triggers_opt.index==0) {
        		//bool Pass_any_PFHT = //(d.hlt.AK8PFJet360_TrimMass30==1) ||
        		//  (d.hlt.PFHT400==1) || (d.hlt.PFHT475==1) || (d.hlt.PFHT600==1) || 
        		//  (d.hlt.PFHT650==1) || (d.hlt.PFHT800==1) ||
        		//  (d.hlt.PFHT550_4JetPt50==1) ||(d.hlt.PFHT650_4JetPt50==1) || (d.hlt.PFHT750_4JetPt50==1);
        		//// JetHT Dataset: Pass any low threshold HT Trigger
                        //if (Pass_any_PFHT) return (size_t)0;
                        return (size_t)0;
                      } else if (triggers_opt.index==1) { // SingleElectron
                        if (d.hlt.Ele27_WPTight_Gsf==1&&nEleTight==1) return (size_t)1;
                      } else if (triggers_opt.index==2) { // SingleMuon
                        if (d.hlt.IsoMu24==1&&nMuTight==1) return (size_t)2;
                      } else if (triggers_opt.index==3) { // MET
        		if (d.hlt.MET200==1&&d.met.Pt[0]>300) return (size_t)3;
                      } else if (triggers_opt.index==4) { // Background MC
                        return (size_t)4;
                      }
                      return (size_t)-1; 
                    }, "JetHT;Ele27WPTight;IsoMu24;MET200;MC", "JetHT (All events);SingleEle (Ele27, Tight);SingleMu (IsoMu24, Tight);MET (MET>300);Simulation", "1,417,601,618,633");
  //sh.AddNewPostfix("Triggers", [&d] { return triggers_opt.index; },
  //      	   "JetHT;SingleEle;SingleMu;MET;MC", "Data: JetHT;Data: SingleEle;Data: SingleMu;Data: MET;MC: t#bar{t}", "1,417,601,618,633");
  
  // Systematics Postfixes
  sh.AddNewPostfix("Syst", [&syst_index] { return syst_index; }, std::string(";Syst[0to")+std::to_string(syst_nSyst)+"]", std::string(";systematics [0to")+std::to_string(syst_nSyst)+"]", "1-999");
  if (syst_nSyst>998) utils::error("Error: Too large number of systematics, define more colors!");



  // Cut Postfixes
  sh.AddNewPostfix("PassBaselineCuts", [] { return 0; }, "Pass0Cuts", "No cuts", "1");
  all_cuts.push_back("PassBaselineCuts");
  for (const auto& region : analysis_cuts) {
    std::string cutflow_str = "";
    sh.AddNewPostfix("PassAllCuts"+std::string(1,region.first), [this,region] { return apply_all_cuts(region.first) ? 0 : (size_t)-1; },
		     "PassAllCuts"+std::string(1,region.first), "Pass all cuts in "+std::string(1,region.first), "1");
    for (size_t i=0, n=region.second.size(); i<n; ++i) {
      // Cuts in order 1-N: "PassNCuts[search region]"
      sh.AddNewPostfix("Pass"+std::to_string(i+1)+"Cuts"+std::string(1,region.first), [this,i,region] { return apply_ncut(region.first, i) ? 0 : (size_t)-1; },
		       "Pass"+std::to_string(i+1)+"Cuts"+std::string(1,region.first), "Cuts up to "+region.second[i].name+" in "+std::string(1,region.first), "1");
      all_cuts.push_back("Pass"+std::to_string(i+1)+"Cuts"+std::string(1,region.first));
      cutflow_str += region.second[i].name+std::string(1,region.first)+";";
      // N-1 Cuts: "PassAllCuts[search region]Excl[cut]"
      sh.AddNewPostfix("PassAllCuts"+std::string(1,region.first)+"Excl"+region.second[i].name, [this,i,region] { 
			 unsigned int mask = (1<<region.second.size())-1 - (1<<i); 
			 return ((cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
		       "PassAllCuts"+std::string(1,region.first)+"Excl"+region.second[i].name, region.second[i].name+" (N-1) in "+std::string(1,region.first), "1");
      // N-2 Cuts: "PassAllCuts[search region]Excl[cut1][cut2]"
      for (size_t j=i+1, n=region.second.size(); j<n; ++j)
	sh.AddNewPostfix("PassAllCuts"+std::string(1,region.first)+"Excl"+region.second[i].name+region.second[j].name, [this,i,j,region] { 
			   unsigned int mask = (1<<region.second.size())-1 - (1<<i) - (1<<j); 
			   return ((cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
			 "PassAllCuts"+std::string(1,region.first)+"Excl"+region.second[i].name+region.second[j].name, region.second[i].name+", "+region.second[j].name+" (N-2) in "+std::string(1,region.first), "1");
    }
    // Stackable Cut Histos: "CutFlow"
    sh.AddNewPostfix("CutFlow"+std::string(1,region.first), [this,region] { for (size_t i=0, n=region.second.size(); i<n; ++i) if (!region.second[i].func()) return i; return region.second.size(); }, 
		     cutflow_str+"PassAll"+std::string(1,region.first), cutflow_str+std::string(1,region.first)+" region", col10+col10);
  }
  sh.AddNewPostfix("PassTriggerPreSelection", [this] { return apply_cuts('W', {W_MR_R2, W_1Wpre})?0:(size_t)-1; }, "PassTriggerPreSelection", "Pass trigger pre-selection", "1");
  sh.AddNewPostfix("PassTriggerPReSelNoMass", [this] { return apply_cut('W', W_MR_R2)?0:(size_t)-1; },      "PassTriggerPreSelNoMass", "Loose Razor cut only", "1");
  sh.AddNewPostfix("OtherUnisoLep",           [this] { return min(nLepVetoNoIso - nLepTight,(unsigned int)1); }, "NoOtherUnisoLep;OtherUnisoLep", "0 other unisol. lepton;#geq1 other unisol. lepton", "418,633");
  sh.AddNewPostfix("IsLepIsolated",           [this] { return isLepIsolated;                                  }, "LeptonNotIsolated;LeptonIsolated", "Lepton not isolated;Lepton is isolated", "633,418");
  sh.AddNewPostfix("AllLepIsolated",          [this] { return allVetoLepIsolated;                             }, "NotAllLepIsoLated;AllLepIsolated", "Not all lepton isolated;All lepton isolated", "633,418");
  
  // Individual Cuts implemented as Postfixes
  // Triggers
  sh.AddNewPostfix("PassHLT",      [this,&dirname,&d] { 
		     if (dirname.find("SingleElectron")!=std::string::npos) return (size_t)-1;
		     else if (dirname.find("SingleMuon")!=std::string::npos) return (size_t)-1;
		     else if (dirname.find("MET")!=std::string::npos) return (size_t)-1;
		     //else if (dirname.find("JetHT")!=std::string::npos) return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1 ? 0 : (size_t)-1;
		     else if (dirname.find("JetHT")!=std::string::npos) return (size_t)0;
		     else return (size_t)0;
		 //}, "PassHLT",      "HLT cut", "1");
		   }, "PassHLT",      "", "1");
  sh.AddNewPostfix("PFHT475",         [&d] { if (d.hlt.PFHT475==-9999) return (size_t)-1; else return (size_t)d.hlt.PFHT475; }, "NoPassHLT_PFHT475;PassHLT_PFHT475", "Do not pass HLT_PFHT475;Pass HLT_PFHT475", "633;418");



  // AK4 Jet Postfixes
  sh.AddNewPostfix("Jets",    [&d] {  size_t i=itJet[d.jetsAK4Puppi.it];        return (i<4)?i:(size_t)-1; }, "Jet[1to5]",  "1st Jet;2nd Jet;3rd Jet;[4to5]th Jet", col5_red_to_green);
  sh.AddNewPostfix("BTags",   [&d] {  size_t i=itMediumBTag[d.jetsAK4Puppi.it]; return (i<4)?i:(size_t)-1; }, "BTag[1to5]", "1st b;2nd b;3rd b;[4to5]th b",         col5_red_to_green);



  // AK8 Jet Postfixes
  sh.AddNewPostfix("JetsAK8",  [&d] {  size_t i=itJetAK8[d.jetsAK8Puppi.it];    return (i<4)?i:(size_t)-1; }, "Jet[1to4]",     "1st Jet;2nd Jet;3rd Jet;4th Jet",                     col4_red_to_cyan);
  sh.AddNewPostfix("WPreTags", [&d] {  size_t i=itWPreTag[d.jetsAK8Puppi.it];   return (i<4)?i:(size_t)-1; }, "WPreTag[1to4]", "1st W-pretag;2nd W-pretag;3rd W-pretag;4th W-pretag", col4_red_to_cyan);
  sh.AddNewPostfix("Ws",       [&d] {  size_t i=itTightWTag[d.jetsAK8Puppi.it]; return (i<4)?i:(size_t)-1; }, "W[1to4]",       "1st W;2nd W;3rd W;4th W",                             col4_red_to_cyan);
  sh.AddNewPostfix("Jet1AK8Pt>450",  [&d] {  return d.jetsAK8Puppi.Pt[iJetAK8[0]]>450 ? 0 : (size_t)-1; },          "Jet1AK8_Pt450",  "1st jet p_{T} (AK8) > 450", "1");
  sh.AddNewPostfix("Jet1AK8Pt>500",  [&d] {  return d.jetsAK8Puppi.Pt[iJetAK8[0]]>500 ? 0 : (size_t)-1; },          "Jet1AK8_Pt500",  "1st jet p_{T} (AK8) > 500", "1");
  sh.AddNewPostfix("Jet1AK8Mass>65", [&d] {  return d.jetsAK8Puppi.softDropMass[iJetAK8[0]]>65 ? 0 : (size_t)-1; }, "Jet1AK8_Mass65", "1st jet M_{SD} (AK8) > 65", "1");


  // Event
  sh.AddNewPostfix("RBins",          [&d] { return (size_t)((d.evt.R>=0.1)+(d.evt.R>=0.2)+(d.evt.R>=0.4)); }, "R0to0p1;R0p1to0p2;R0p2to0p4;R0p4", "0.0<R<0.1;0.1<R<0.2;0.2<R<0.4;R>=0.4", "1,4,418,633");



  // --------------------------------------------------------------------
  //                         Fill Parameters
  // --------------------------------------------------------------------



  // Bins
  std::vector<double> E   = {0, 100, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 5000, 10000};
  std::vector<double> Pt  = {0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000, 5000};
  std::vector<double> M   = {0, 10, 20, 30, 40, 50, 65, 75, 85, 95, 105, 120, 135, 150, 165, 180, 195, 210, 230, 260, 300, 500, 1000};
  std::vector<double> MW  = {65, 75, 85, 95, 105};
  std::vector<double> PtF = {0, 300, 400, 600, 1000, 2000, 5000};
  std::vector<double> CSV = {0, 0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 };
  std::vector<double> R   = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.7, 0.8, 1.0, 1.2, 2.0 };
  std::vector<double> MR  = {0, 600, 800, 1000, 1200, 1600, 2000, 4000, 10000};
  std::vector<double> R2  = {0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.5, 1.0, 5.0};
  std::vector<double> HT  = {0, 100, 200, 300, 400, 500, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500, 2000, 4000, 10000};
  sh.AddNewFillParam("Bin", { .nbin=1,   .bins={0,1}, .fill=[&d] { return 0; }, .axis_title="Bin"}); // For averages/counts

  // AK4 Jets
  /* Variables for ID
  sh.AddNewFillParam("JetPhotonE",           { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.PhotonEnergy[d.jetsAK4Puppi.it];              }, .axis_title="Jet Photon Energy (GeV)"});
  sh.AddNewFillParam("JetElectronE",         { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.ElectronEnergy[d.jetsAK4Puppi.it];            }, .axis_title="Jet Electron Energy (GeV)"});
  sh.AddNewFillParam("JetMuonE",             { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.MuonEnergy[d.jetsAK4Puppi.it];                }, .axis_title="Jet Muon Energy (GeV)"});
  sh.AddNewFillParam("JetChargedMuE",        { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.ChargeMuEnergy[d.jetsAK4Puppi.it];            }, .axis_title="Jet Charged Mu Energy (GeV)"});
  sh.AddNewFillParam("JetChargedEmE",        { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.chargedEmEnergy[d.jetsAK4Puppi.it];           }, .axis_title="Jet Charged Em Energy (GeV)"});
  sh.AddNewFillParam("JetChargedHadronE",    { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.chargedHadronEnergy[d.jetsAK4Puppi.it];       }, .axis_title="Jet Charged Hadron Energy (GeV)"});
  sh.AddNewFillParam("JetNeutralHadronE",    { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.neutralHadronEnergy[d.jetsAK4Puppi.it];       }, .axis_title="Jet Neutral Hadron Energy (GeV)"});
  sh.AddNewFillParam("JetNeutralEmE",        { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.neutralEmEnergy[d.jetsAK4Puppi.it];           }, .axis_title="Jet Neutral Em Energy (GeV)"});
  sh.AddNewFillParam("JetHFHadronE",         { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.HFHadronEnergy[d.jetsAK4Puppi.it];            }, .axis_title="Jet HF Hadron Energy (GeV)"});
  sh.AddNewFillParam("JetHFEME",             { .nbin= E.size()-1, .bins=E,         .fill=[&d] { return d.jetsAK4Puppi.HFEMEnergy[d.jetsAK4Puppi.it];                }, .axis_title="Jet HF EM Energy (GeV)"});
  sh.AddNewFillParam("JetNumberOfDaughters", { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4Puppi.numberOfDaughters[d.jetsAK4Puppi.it];         }, .axis_title="Jet Number Of Daughters"});
  sh.AddNewFillParam("JetPhotonMult",        { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4Puppi.photonMultiplicity[d.jetsAK4Puppi.it];        }, .axis_title="Jet Photon Multiplicity"});
  sh.AddNewFillParam("JetElectronMult",      { .nbin=  10, .bins={     0,     10}, .fill=[&d] { return d.jetsAK4Puppi.electronMultiplicity[d.jetsAK4Puppi.it];      }, .axis_title="Jet Electron Multiplicity"});
  sh.AddNewFillParam("JetMuonMult",          { .nbin=  10, .bins={     0,     10}, .fill=[&d] { return d.jetsAK4Puppi.muonMultiplicity[d.jetsAK4Puppi.it];          }, .axis_title="Jet Muon Multiplicity"});
  sh.AddNewFillParam("JetNeutralMult",       { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4Puppi.neutralMultiplicity[d.jetsAK4Puppi.it];       }, .axis_title="Jet Neutral Multiplicity"});
  sh.AddNewFillParam("JetChargedMult",       { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4Puppi.chargedMultiplicity[d.jetsAK4Puppi.it];       }, .axis_title="Jet Charged Multiplicity"});
  sh.AddNewFillParam("JetChargedHadronMult", { .nbin= 120, .bins={     0,    120}, .fill=[&d] { return d.jetsAK4Puppi.ChargedHadronMultiplicity[d.jetsAK4Puppi.it]; }, .axis_title="Jet Charged Hadron Multiplicity"});
  sh.AddNewFillParam("JetNeutralHadronMult", { .nbin=  40, .bins={     0,     40}, .fill=[&d] { return d.jetsAK4Puppi.neutralHadronMultiplicity[d.jetsAK4Puppi.it]; }, .axis_title="Jet Neutral Hadron Multiplicity"});
  sh.AddNewFillParam("JetHFHadronMult",      { .nbin=  80, .bins={     0,     80}, .fill=[&d] { return d.jetsAK4Puppi.HFHadronMultiplicity[d.jetsAK4Puppi.it];      }, .axis_title="Jet HF Hadron Multiplicity"});
  sh.AddNewFillParam("JetHFEMMult",          { .nbin=  50, .bins={     0,     50}, .fill=[&d] { return d.jetsAK4Puppi.HFEMMultiplicity[d.jetsAK4Puppi.it];          }, .axis_title="Jet HF EM Multiplicity"});
  */
  sh.AddNewFillParam("JetPtBins",            { .nbin=Pt.size()-1,   .bins=Pt,      .fill=[&d] { return d.jetsAK4Puppi.Pt[d.jetsAK4Puppi.it];           }, .axis_title="Jet p_{T} (GeV)", .def_range={0,1000} });
  sh.AddNewFillParam("JetPtFewBins",         { .nbin=PtF.size()-1,  .bins=PtF,     .fill=[&d] { return d.jetsAK4Puppi.Pt[d.jetsAK4Puppi.it];           }, .axis_title="Jet p_{T} (GeV)", .def_range={0,1000} });
  sh.AddNewFillParam("JetPtOneBin",          { .nbin=   1, .bins={   400,   5000}, .fill=[&d] { return d.jetsAK4Puppi.Pt[d.jetsAK4Puppi.it];           }, .axis_title="Jet p_{T} (GeV)"});
  sh.AddNewFillParam("JetPt",                { .nbin= 500, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK4Puppi.Pt[d.jetsAK4Puppi.it];           }, .axis_title="Jet p_{T} (GeV)", .def_range={0,1000} });
  sh.AddNewFillParam("JetEta",               { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK4Puppi.Eta[d.jetsAK4Puppi.it];          }, .axis_title="Jet #eta",        .def_range={-3,3}});
  sh.AddNewFillParam("JetPhi",               { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK4Puppi.Phi[d.jetsAK4Puppi.it];          }, .axis_title="Jet #phi"});
  sh.AddNewFillParam("JetCSV",               { .nbin= 101, .bins={     0,   1.01}, .fill=[&d] { return d.jetsAK4Puppi.CSVv2[d.jetsAK4Puppi.it];        }, .axis_title="Jet CSV"});
  // BJets
  sh.AddNewFillParam("BJetPtBins",           { .nbin=Pt.size()-1,   .bins=Pt,      .fill=[&d] { return d.jetsAK4Puppi.Pt[d.jetsAK4Puppi.it];           }, .axis_title="B-jet p_{T} (GeV)", .def_range={0,1000} });
  sh.AddNewFillParam("BJetPt",               { .nbin= 500, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK4Puppi.Pt[d.jetsAK4Puppi.it];           }, .axis_title="B-jet p_{T} (GeV)", .def_range={0,1000} });
  sh.AddNewFillParam("BJetEta",              { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK4Puppi.Eta[d.jetsAK4Puppi.it];          }, .axis_title="B-jet #eta",        .def_range={-3,3}});
  sh.AddNewFillParam("BJetPhi",              { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK4Puppi.Phi[d.jetsAK4Puppi.it];          }, .axis_title="B-jet #phi"});



  // AK8 Jets
  sh.AddNewFillParam("JetAK8PtOneBin",       { .nbin=   1, .bins={   400,   5000}, .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="AK8 jet p_{T} (GeV)"});
  sh.AddNewFillParam("JetAK8PtFewBins",      { .nbin=PtF.size()-1, .bins=PtF,      .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("JetAK8PtBins",         { .nbin=Pt.size()-1 , .bins=Pt,       .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("JetAK8Pt",             { .nbin= 500, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("JetAK8Eta",            { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8Puppi.Eta[d.jetsAK8Puppi.it];          }, .axis_title="AK8 jet #eta",        .def_range={-3,3}});
  sh.AddNewFillParam("JetAK8Phi",            { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8Puppi.Phi[d.jetsAK8Puppi.it];          }, .axis_title="AK8 jet #phi"});
  sh.AddNewFillParam("JetAK8SoftDropMass",   { .nbin= 200, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8Puppi.softDropMass[d.jetsAK8Puppi.it]; }, .axis_title="AK8 jet soft-drop mass (GeV)", .def_range={0,150}});
  /*
  sh.AddNewFillParam("JetAK8PrunedMass",     { .nbin= 400, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8Puppi.prunedMass[d.jetsAK8Puppi.it];   }, .axis_title="AK8 jet pruned mass (GeV)",    .def_range={0,150}});
  sh.AddNewFillParam("JetAK8FilteredMass",   { .nbin= 400, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8Puppi.filteredMass[d.jetsAK8Puppi.it]; }, .axis_title="AK8 jet filtered mass (GeV)",  .def_range={0,150}});
  sh.AddNewFillParam("JetAK8TrimmedMass",    { .nbin= 400, .bins={     0,   2000}, .fill=[&d] { return d.jetsAK8Puppi.trimmedMass[d.jetsAK8Puppi.it];  }, .axis_title="AK8 jet trimmed mass (GeV)",   .def_range={0,150}});
  */
  sh.AddNewFillParam("JetAK8Tau1",           { .nbin= 100, .bins={     0,      1}, .fill=[&d] { return d.jetsAK8Puppi.tau1[d.jetsAK8Puppi.it];         }, .axis_title="AK8 jet #tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau2",           { .nbin= 100, .bins={     0,      1}, .fill=[&d] { return d.jetsAK8Puppi.tau2[d.jetsAK8Puppi.it];         }, .axis_title="AK8 jet #tau_{2}"});
  sh.AddNewFillParam("JetAK8Tau3",           { .nbin= 100, .bins={     0,      1}, .fill=[&d] { return d.jetsAK8Puppi.tau3[d.jetsAK8Puppi.it];         }, .axis_title="AK8 jet #tau_{3}"});
  sh.AddNewFillParam("JetAK8Tau21",          { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return tau21[d.jetsAK8Puppi.it];                       }, .axis_title="AK8 jet #tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau31",          { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return tau31[d.jetsAK8Puppi.it];                       }, .axis_title="AK8 jet #tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetAK8Tau32",          { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return tau32[d.jetsAK8Puppi.it];                       }, .axis_title="AK8 jet #tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("MaxAK8SubjetCSV",      { .nbin= 101, .bins={     0,   1.01}, .fill=[&d] { return maxSubjetCSV[d.jetsAK8Puppi.it];                }, .axis_title="Max. AK8 subjet CSV"});
  sh.AddNewFillParam("MaxAK8SubJetCSVBins",  { .nbin=CSV.size()-1, .bins=CSV,      .fill=[&d] { return maxSubjetCSV[d.jetsAK8Puppi.it];                }, .axis_title="Max. AK8 subjet CSV"});
  // WPreTags
  sh.AddNewFillParam("WPreTagPtBins",        { .nbin=Pt.size()-1 , .bins=Pt,       .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="W-pretagged AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("WPreTagPt",            { .nbin= 500, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="W-pretagged AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("WPreTagEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8Puppi.Eta[d.jetsAK8Puppi.it];          }, .axis_title="W-pretagged AK8 jet #eta",        .def_range={-3,3}});
  sh.AddNewFillParam("WPreTagPhi",           { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8Puppi.Phi[d.jetsAK8Puppi.it];          }, .axis_title="W-pretagged AK8 jet #phi"});
  sh.AddNewFillParam("WPreTagTau21",         { .nbin=  50, .bins={     0,      1}, .fill=[&d] { return tau21[d.jetsAK8Puppi.it];                       }, .axis_title="W-pretagged AK8 jet #tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("WPreTagMass",          { .nbin=MW.size()-1, .bins=MW,        .fill=[&d] { return d.jetsAK8Puppi.softDropMass[d.jetsAK8Puppi.it]; }, .axis_title="W-pretagged AK8 jet M_{Soft-Drop} (GeV)"});
  // WPreTags
  sh.AddNewFillParam("WTagPtBins",           { .nbin=Pt.size()-1 , .bins=Pt,       .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="W-tagged AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("WTagPt",               { .nbin= 500, .bins={     0,  10000}, .fill=[&d] { return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];           }, .axis_title="W-tagged AK8 jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("WTagEta",              { .nbin=  40, .bins={    -4,      4}, .fill=[&d] { return d.jetsAK8Puppi.Eta[d.jetsAK8Puppi.it];          }, .axis_title="W-tagged AK8 jet #eta",        .def_range={-3,3}});
  sh.AddNewFillParam("WTagPhi",              { .nbin=  16, .bins={-3.142,  3.142}, .fill=[&d] { return d.jetsAK8Puppi.Phi[d.jetsAK8Puppi.it];          }, .axis_title="W-tagged AK8 jet #phi"});

  // Leptons
  sh.AddNewFillParam("ElePt",                { .nbin= 400, .bins={    0,    2000}, .fill=[&d] { return d.ele.Pt[d.ele.it];                             }, .axis_title="Electron p_{T} (GeV)"});
  sh.AddNewFillParam("MuPt",                 { .nbin= 400, .bins={    0,    2000}, .fill=[&d] { return d.mu.Pt[d.mu.it];                               }, .axis_title="Muon p_{T} (GeV)"});
  /*
  sh.AddNewFillParam("EleDRJet",             { .nbin=  60, .bins={    0,       6}, .fill=[&d] { return d.evt.EleDRJet[d.ele.it];                       }, .axis_title="#DeltaR (e, jet)"});
  sh.AddNewFillParam("EleRelPtJet",          { .nbin=  50, .bins={    0,     500}, .fill=[&d] { return d.evt.EleRelPtJet[d.ele.it];                    }, .axis_title="p_{T}^{rel} (e, jet) (GeV)"});
  sh.AddNewFillParam("MuDRJet",              { .nbin=  60, .bins={    0,       6}, .fill=[&d] { return d.evt.MuDRJet[d.mu.it];                         }, .axis_title="#DeltaR (#mu, jet)"});
  sh.AddNewFillParam("MuRelPtJet",           { .nbin=  50, .bins={    0,     500}, .fill=[&d] { return d.evt.MuRelPtJet[d.mu.it];                      }, .axis_title="p_{T}^{rel} (#mu, jet) (GeV)"});
  */



  // Event
  // Cuts
  // Object counts
  sh.AddNewFillParam("NVtx",                 { .nbin= 100, .bins={    0,     100}, .fill=[&d] { return d.evt.NGoodVtx;          }, .axis_title="N_{Vertices}",         .def_range={0,50}});
  sh.AddNewFillParam("NJet",                 { .nbin=  50, .bins={    0,      50}, .fill=[&d] { return nJet;                    }, .axis_title="N_{Jet}",              .def_range={0,10}});
  sh.AddNewFillParam("NJetAK8",              { .nbin=  10, .bins={    0,      10}, .fill=[&d] { return nJetAK8;                 }, .axis_title="N_{AK8 jet}",          .def_range={0,5}});
  sh.AddNewFillParam("NBTag",                { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nMediumBTag;             }, .axis_title="N_{b}",                .def_range={0,5}});
  sh.AddNewFillParam("NLooseBTag",           { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nLooseBTag;              }, .axis_title="N_{b, loose tag}",     .def_range={0,5}});
  sh.AddNewFillParam("NTightBTag",           { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nTightBTag;              }, .axis_title="N_{b, tight tag}",     .def_range={0,5}});
  sh.AddNewFillParam("NWPreTag",             { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nWPreTag;                }, .axis_title="N_{W, pre-tag}",       .def_range={0,5}});
  sh.AddNewFillParam("NWTag",                { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nTightWTag;              }, .axis_title="N_{W}",                .def_range={0,5}});
  sh.AddNewFillParam("NLooseWTag",           { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nLooseWTag;              }, .axis_title="N_{W, loose tag}",     .def_range={0,5}});
  sh.AddNewFillParam("NHadTopTag",           { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nHadTopTag;              }, .axis_title="N_{top (had.)}",       .def_range={0,5}});
  sh.AddNewFillParam("NLepVeto",             { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nLepVeto;                }, .axis_title="N_{lepton, Veto}",     .def_range={0,5}});
  sh.AddNewFillParam("NEleVeto",             { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nEleVeto;                }, .axis_title="N_{ele, Veto}",        .def_range={0,5}});
  sh.AddNewFillParam("NMuVeto",              { .nbin=  20, .bins={    0,      20}, .fill=[&d] { return nMuVeto;                 }, .axis_title="N_{muon, Veto}",       .def_range={0,5}});
  sh.AddNewFillParam("NLep",                 { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nLepTight;               }, .axis_title="N_{lepton}",           .def_range={0,5}});
  sh.AddNewFillParam("NEle",                 { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nEleTight;               }, .axis_title="N_{ele}",              .def_range={0,5}});
  sh.AddNewFillParam("NMu",                  { .nbin=   5, .bins={    0,       5}, .fill=[&d] { return nMuTight;                }, .axis_title="N_{muon}",             .def_range={0,5}});
  // Razor
  sh.AddNewFillParam("R",                    { .nbin=  40, .bins={    0,     2.0}, .fill=[&d] { return d.evt.R;                 }, .axis_title="R",                    .def_range={0,1}});
  sh.AddNewFillParam("RFine",                { .nbin= 200, .bins={    0,     2.0}, .fill=[&d] { return d.evt.R;                 }, .axis_title="R",                    .def_range={0,1}});
  sh.AddNewFillParam("RBins",                { .nbin=R.size()-1, .bins=R,          .fill=[&d] { return d.evt.R;                 }, .axis_title="R",                    .def_range={0,1}});
  sh.AddNewFillParam("MTR",                  { .nbin= 200, .bins={    0,    4000}, .fill=[&d] { return d.evt.MTR;               }, .axis_title="M_{T}^{R} (GeV)",      .def_range={0,1000}});
  //sh.AddNewFillParam("MR",                   { .nbin= 100, .bins={    0,   10000}, .fill=[&d] { return d.evt.MR;                }, .axis_title="M_{R} (GeV)",          .def_range={0,2000}});
  //sh.AddNewFillParam("R2",                   { .nbin=  80, .bins={    0,     4.0}, .fill=[&d] { return d.evt.R2;                }, .axis_title="R^{2}",                .def_range={0,1}});
  sh.AddNewFillParam("MR",                   { .nbin=MR.size()-1, .bins=MR,        .fill=[&d] { return d.evt.MR;                }, .axis_title="M_{R} (GeV)",          .def_range={0,4000}});
  sh.AddNewFillParam("R2",                   { .nbin=R2.size()-1, .bins=R2,        .fill=[&d] { return d.evt.R2;                }, .axis_title="R^{2}",                .def_range={0,1}});
  // HT
  sh.AddNewFillParam("HT",                   { .nbin=HT.size()-1, .bins=HT,        .fill=[&d] { return AK4Puppi_Ht;             }, .axis_title="H_{T} (GeV)",          .def_range={0, 2000}});
  sh.AddNewFillParam("GenHT",                { .nbin=HT.size()-1, .bins=HT,        .fill=[&d] { return d.evt.Gen_Ht;            }, .axis_title="H_{T}^{Gen} (GeV)",    .def_range={0, 2000}});
  sh.AddNewFillParam("AK8HT",                { .nbin=HT.size()-1, .bins=HT,        .fill=[&d] { return AK8Puppi_Ht;             }, .axis_title="H_{T}^{AK8} (GeV)",    .def_range={0, 2000}});
  // MET
  sh.AddNewFillParam("MET",                  { .nbin= 200, .bins={    0,    4000}, .fill=[&d] { return d.met.Pt[0];             }, .axis_title="#slash{p}_{T}) (GeV)", .def_range={0,1000}});
  // DPhi
  sh.AddNewFillParam("MinDeltaPhi",          { .nbin=  64, .bins={    0,     3.2}, .fill=[]   { return minDeltaPhi;             }, .axis_title="Min(|#Delta#phi(Jet_{1-3}, MET)|"});
  // SUSY
  sh.AddNewFillParam("MGluino",              { .nbin= 121, .bins={-12.5, 3012.5 }, .fill=[&d] { return d.evt.SUSY_Gluino_Mass;  }, .axis_title="M_{#tilde{g}} (GeV)",        .def_range={550,2350}});
  sh.AddNewFillParam("MStop",                { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[&d] { return d.evt.SUSY_Stop_Mass;    }, .axis_title="M_{#tilde{s}} (GeV)",        .def_range={  0,1650}});
  sh.AddNewFillParam("MLSP",                 { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[&d] { return d.evt.SUSY_LSP_Mass;     }, .axis_title="M_{#tilde{#chi}^{0}} (GeV)", .def_range={  0,1650}});
  // AK8 JetN
  sh.AddNewFillParam("Jet1AK8Mass",         { .nbin=M.size()-1, .bins=M,           .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8Puppi.softDropMass[iJetAK8[0]]; }, .axis_title="Leading AK8 jet M_{Soft-Drop} (GeV)",    .def_range={0, 300}});
  sh.AddNewFillParam("Jet2AK8Mass",         { .nbin=M.size()-1, .bins=M,           .fill=[&d] { return (nJetAK8<2) ? -9999. : d.jetsAK8Puppi.softDropMass[iJetAK8[1]]; }, .axis_title="Subleading AK8 jet M_{Soft-Drop} (GeV)", .def_range={0, 300}});
  sh.AddNewFillParam("Jet1AK8Pt",           { .nbin=Pt.size()-1, .bins=Pt,         .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8Puppi.Pt[iJetAK8[0]];           }, .axis_title="Leading AK8 jet p_{T} (GeV)",    .def_range={0, 2000}});
  sh.AddNewFillParam("Jet2AK8Pt",           { .nbin=Pt.size()-1, .bins=Pt,         .fill=[&d] { return (nJetAK8<2) ? -9999. : d.jetsAK8Puppi.Pt[iJetAK8[1]];           }, .axis_title="Subleading AK8 jet p_{T} (GeV)", .def_range={0, 2000}});
  sh.AddNewFillParam("Jet1AK8Eta",          { .nbin=   80, .bins={   -4,       4}, .fill=[&d] { return (nJetAK8<1) ? -9999. : d.jetsAK8Puppi.Eta[iJetAK8[0]];          }, .axis_title="Leading AK8 jet #eta",    .def_range={-3, 3}});
  sh.AddNewFillParam("Jet2AK8Eta",          { .nbin=   80, .bins={   -4,       4}, .fill=[&d] { return (nJetAK8<2) ? -9999. : d.jetsAK8Puppi.Eta[iJetAK8[1]];          }, .axis_title="Subleading AK8 jet #eta", .def_range={-3, 3}});
  sh.AddNewFillParam("Jet1AK8Tau32",        { .nbin=   50, .bins={    0,       1}, .fill=[&d] { return (nJetAK8<1) ? -9999. : tau32[iJetAK8[0]];                       }, .axis_title="Leading AK8 jet #tau_{32}"});
  sh.AddNewFillParam("Jet2AK8Tau32",        { .nbin=   50, .bins={    0,       1}, .fill=[&d] { return (nJetAK8<2) ? -9999. : tau32[iJetAK8[1]];                       }, .axis_title="Subleading AK8 jet #tau_{32}"});
  sh.AddNewFillParam("Jet1AK8BTagCSV",      { .nbin=  101, .bins={    0,    1.01}, .fill=[&d] { return (nJetAK8<1) ? -9999. : maxSubjetCSV[iJetAK8[0]];                }, .axis_title="Leading AK8 jet - Max. Subjet CSV",    .def_range={0,1}});
  sh.AddNewFillParam("Jet2AK8BTagCSV",      { .nbin=  101, .bins={    0,    1.01}, .fill=[&d] { return (nJetAK8<2) ? -9999. : maxSubjetCSV[iJetAK8[1]];                }, .axis_title="Subleading AK8 jet - Max. Subjet CSV", .def_range={0,1}});
  sh.AddNewFillParam("WPreTag1Mass",        { .nbin=MW.size()-1, .bins=MW,         .fill=[&d] { return (nWPreTag<1)? -9999. : d.jetsAK8Puppi.softDropMass[iWPreTag[0]];}, .axis_title="W-pretagged AK8 jet M_{Soft-Drop} (GeV)"});
  // SPECIAL
  // Special Y/Z axis parameters:
  sh.AddSpecial({ .name="Counts",                         .name_plus_1d="Syst",                          .axis="Counts (Incl Syst Unc)",              .axis_plus_1d="Systematics variation index"});
  /*
  sh.AddSpecial({ .name="MergedTopFraction",              .name_plus_1d="IsGenTopMerged",                .axis="Fraction of Merged Tops",             .axis_plus_1d="Gen. b and W Merged in R<0.8 (bool)"});
  sh.AddSpecial({ .name="JetFindingEfficiency",           .name_plus_1d="HasJet",                        .axis="Jet finding Efficiency",              .axis_plus_1d="Found AK8 (Puppi) jet (bool)"});
  sh.AddSpecial({ .name="TopFindingEfficiency",           .name_plus_1d="HasHadTopTaggedJet",            .axis="Top finding Efficiency",              .axis_plus_1d="Has hadronic top-tagged jet (bool)"});
  sh.AddSpecial({ .name="TopTagEfficiency",               .name_plus_1d="JetIsHadTopTagged",             .axis="Top-tagging Efficiency",              .axis_plus_1d="Jet is hadronic top-tagged (bool)"});
  sh.AddSpecial({ .name="MisTagRate",                     .name_plus_1d="JetHasNoGenTop",                .axis="Mis-tag Rate",                        .axis_plus_1d="Jet is mis-matched (bool)"});
  */
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360",             .name_plus_1d="HLT_AK8PFJet360_TrimMass30",    .axis="#epsilon_{HLT_AK8PFJet360_TrimMass30}",                   .axis_plus_1d="HLT_AK8PFJet360_TrimMass30"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet450",             .name_plus_1d="HLT_Ak8PFJet450",               .axis="#epsilon_{HLT_AK8PFJet450}",                              .axis_plus_1d="HLT_AK8PFJet450"});
  sh.AddSpecial({ .name="HLTEff_AK8PFHT700_TrimMass50",   .name_plus_1d="HLT_AK8PFHT700_TrimMass50",     .axis="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}",          .axis_plus_1d="HLT_AK8PFHT700_TrimR0p1PT0p03Mass50"});
  sh.AddSpecial({ .name="HLTEff_PFHT750_4JetPt50",        .name_plus_1d="HLT_PFHT750_4JetPt50",          .axis="#epsilon_{HLT_PFHT750_4JetPt50}",                         .axis_plus_1d="HLT_PFHT750_4JetPt50"});
  sh.AddSpecial({ .name="HLTEff_PFHT800",                 .name_plus_1d="HLT_PFHT800",                   .axis="#epsilon_{HLT_PFHT800}",                                  .axis_plus_1d="HLT_PFHT800"});
  sh.AddSpecial({ .name="HLTEff_AK8DiPFJet250_200",       .name_plus_1d="HLT_AK8DiPFJet250_200",         .axis="#epsilon_{HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20}", .axis_plus_1d="HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20"});
  sh.AddSpecial({ .name="HLTEff_Rsq0p25",                 .name_plus_1d="HLT_Rsq0p25",                   .axis="#epsilon_{HLT_Rsq0p25}",                                  .axis_plus_1d="HLT_Rsq0p25"});
  sh.AddSpecial({ .name="HLTEff_RsqMR270_Rsq0p09_MR200",  .name_plus_1d="HLT_RsqMR270_Rsq0p09_MR200",    .axis="#epsilon_{HLT_RsqMR270_Rsq0p09_MR200}",                   .axis_plus_1d="HLT_RsqMR270_Rsq0p09_MR200"});
  sh.AddSpecial({ .name="HLTEff_AK8PFHT700orPFHT800",     .name_plus_1d="HLT_AK8PFHT700_or_PFHT800",     .axis="#epsilon_{HLT_AK8PFHT700 OR HLT_PFHT800}",                .axis_plus_1d="HLT_AK8PFHT700 OR HLT_PFHT800"});
  sh.AddSpecial({ .name="HLTEff_PFJet450orPFHT800",       .name_plus_1d="HLT_PFJet450_or_PFHT800",       .axis="#epsilon_{HLT_PFJet450 OR HLT_PFHT800}",                  .axis_plus_1d="HLT_PFJet450 OR HLT_PFHT800"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet450orPFHT800",    .name_plus_1d="HLT_AK8PFJet450_or_PFHT800",    .axis="#epsilon_{HLT_AK8PFJet450 OR HLT_PFHT800}",               .axis_plus_1d="HLT_AK8PFJet450 OR HLT_PFHT800"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360orPFHT800",    .name_plus_1d="HLT_AK8PFJet360_or_PFHT800",    .axis="#epsilon_{HLT_AK8PFJet360 OR HLT_PFHT800}",               .axis_plus_1d="HLT_AK8PFJet360 OR HLT_PFHT800"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360orAK8PFHT700", .name_plus_1d="HLT_AK8PFJet360_or_AK8PFHT700", .axis="#epsilon_{HLT_AK8PFJet360 OR HLT_AK8PFHT700}",            .axis_plus_1d="HLT_AK8PFJet360 OR HLT_AK8PFHT700"});

  sh.AddNewFillParam("Counts",                         { .nbin= 1+syst_nSyst, .bins={-0.5, syst_nSyst+0.5}, .fill=[&syst_index] { return syst_index; }, .axis_title="Counts (Incl Syst Unc)"});
  /*
  sh.AddNewFillParam("MergedTopFraction",              { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return top_Children_Within_Cone[d.jetsAK8Puppi.it]; }, .axis_title="Fraction of Merged Tops" });
  sh.AddNewFillParam("JetFindingEfficiency",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return has_Matched_Jet[d.gen.it];                   }, .axis_title="Jet finding Efficiency" });
  sh.AddNewFillParam("TopFindingEfficiency",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return has_Matched_Tagged_Jet[d.gen.it];            }, .axis_title="Top finding Efficiency" });
  sh.AddNewFillParam("TopTagEfficiency",               { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return passHadTopTag[d.jetsAK8Puppi.it];            }, .axis_title="Top-tagging Efficiency" });
  */
  sh.AddNewFillParam("HLTEff_AK8PFJet360",             { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet360_TrimMass30;                   }, .axis_title="#epsilon_{HLT_AK8PFJet360_TrimMass30}",                   .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet450",             { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet450;                              }, .axis_title="#epsilon_{HLT_AK8PFJet450}",                              .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFHT700_TrimMass50",   { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50;          }, .axis_title="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}",          .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_PFHT750_4JetPt50",        { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.PFHT750_4JetPt50;                         }, .axis_title="#epsilon_{HLT_PFHT750_4JetPt50}",                         .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_PFHT800",                 { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.PFHT800;                                  }, .axis_title="#epsilon_{HLT_PFHT800}",                                  .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8DiPFJet250_200",       { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8DiPFJet250_200_TrimMass30_BTagCSV_p20; }, .axis_title="#epsilon_{HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20}", .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_Rsq0p25",                 { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.Rsq0p25;                                  }, .axis_title="#epsilon_{HLT_Rsq0p25}",                                  .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_RsqMR270_Rsq0p09_MR200",  { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.RsqMR270_Rsq0p09_MR200;                   }, .axis_title="#epsilon_{HLT_RsqMR270_Rsq0p09_MR200}",                   .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFHT700orPFHT800",     { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1 || d.hlt.PFHT800==1;                }, .axis_title="#epsilon_{HLT_AK8PFHT700 OR HLT_PFHT800}",  .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_PFJet450orPFHT800",       { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.PFJet450==1 || d.hlt.PFHT800==1;                                       }, .axis_title="#epsilon_{HLT_PFJet450 OR HLT_PFHT800}",    .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet450orPFHT800",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet450==1 || d.hlt.PFHT800==1;                                    }, .axis_title="#epsilon_{HLT_AK8PFJet450 OR HLT_PFHT800}", .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet360orPFHT800",    { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet360_TrimMass30==1 || d.hlt.PFHT800==1;                         }, .axis_title="#epsilon_{HLT_AK8PFJet360 OR HLT_PFHT800}",     .def_range={0,1} });
  sh.AddNewFillParam("HLTEff_AK8PFJet360orAK8PFHT700", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[&d] { return d.hlt.AK8PFJet360_TrimMass30==1 || d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; }, .axis_title="#epsilon_{HLT_AK8PFJet360 OR HLT_AK8PFHT700}",  .def_range={0,1} });

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
      + W-pretags (mass tag only)  - num, dau1 vs dau2 pt, dau1 vs dau2 mass, massdrop(/pt), yasym, mdr(/pt) vs yasym
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
  sh.AddHistoType("AK4");
  sh.AddHistoType("b");
  sh.AddHistoType("b loose");
  sh.AddHistoType("AK8");
  sh.AddHistoType("Wpre");
  sh.AddHistoType("W");
  sh.AddHistoType("ele");
  sh.AddHistoType("ele veto");
  sh.AddHistoType("mu");
  sh.AddHistoType("mu veto");
  sh.AddHistoType("evt");
  sh.AddHistoType("syst");

  // Histo options
  std::string d = "HISTE1";
  std::string o_stk = "LogSumw2Stack5AddRatio";
  std::vector<double> r_stk = {0,0, 1.01e-3,1e10, 0.55,0.9};
  std::string Stack = "StackPlot";

  // -------------------------------------------------------------------------
  //                              Selected AK4 Jets
  
  // Event
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"PassBaselineCuts",    "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"Pass2CutsS",          "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"Pass3CutsS",          "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"Pass4CutsS",          "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"Pass5CutsS",          "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"Pass6CutsS",          "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"PassAllCutsSExcl3Jet","PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJet",               .pfs={Stack,"PassAllCutsS",        "PassHLT"},   .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  
  for (const auto& cut : {"Pass1CutsS","PassAllCutsS","PassAllCutsT","PassAllCutsW","PassAllCutsQ"}) {
    // All jets
    sh.AddHistos("AK4",  { .fill="JetPtBins",          .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetPt",              .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetEta",             .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetPhi",             .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetCSV",             .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    // pt order
    sh.AddHistos("AK4",  { .fill="JetPtBins",          .pfs={Stack,cut,"PassHLT","Jets"},     .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetPt",              .pfs={Stack,cut,"PassHLT","Jets"},     .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetEta",             .pfs={Stack,cut,"PassHLT","Jets"},     .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetPhi",             .pfs={Stack,cut,"PassHLT","Jets"},     .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK4",  { .fill="JetCSV",             .pfs={Stack,cut,"PassHLT","Jets"},     .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  }
  
  // -------------------------------------------------------------------------
  //                                 Selected bs
  
  // Event  
  for (const auto& cut : {"PassBaselineCuts","Pass1CutsS","Pass2CutsS","Pass3CutsS","PassAllCutsSExcl1b","PassAllCutsTExcl1b","PassAllCutsWExcl0b","PassAllCutsQExcl0b"}) {
    sh.AddHistos("evt",  { .fill="NBTag",              .pfs={Stack,cut,"PassHLT"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NBTag",              .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="NBTag",              .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
  }
  
  for (const auto& cut : {"Pass1CutsS","PassAllCutsS","PassAllCutsT"}) {
    // b jets (cut applied)
    sh.AddHistos("b",    { .fill="BJetPtBins",         .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("b",    { .fill="BJetPt",             .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("b",    { .fill="BJetEta",            .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("b",    { .fill="BJetPhi",            .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    // pt order
    sh.AddHistos("b",    { .fill="BJetPtBins",         .pfs={Stack,cut,"PassHLT","BTags"},    .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("b",    { .fill="BJetPt",             .pfs={Stack,cut,"PassHLT","BTags"},    .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("b",    { .fill="BJetEta",            .pfs={Stack,cut,"PassHLT","BTags"},    .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("b",    { .fill="BJetPhi",            .pfs={Stack,cut,"PassHLT","BTags"},    .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  }
  
  // -------------------------------------------------------------------------
  //                              Selected AK8 Jets
  
  // Event
  sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,"PassBaselineCuts","PassHLT"},      .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,"Pass1CutsS",      "PassHLT"},      .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,"Pass2CutsS",      "PassHLT"},      .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,"Pass3CutsS",      "PassHLT"},      .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,"Pass6CutsS",      "PassHLT"},      .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  sh.AddHistos("evt",  { .fill="NJetAK8",            .pfs={Stack,"PassAllCutsS",    "PassHLT"},      .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  
  // All jets
  for (const auto& cut : {"Pass1CutsS","PassAllCutsSExcl1W","PassAllCutsTExcl1W"}) {
    sh.AddHistos("AK8",  { .fill="JetAK8PtBins",       .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Pt",           .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Eta",          .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Phi",          .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8SoftDropMass", .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Tau21",        .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Tau32",        .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="MaxAK8SubjetCSV",    .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    // pt order
    sh.AddHistos("AK8",  { .fill="JetAK8PtBins",       .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Pt",           .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Eta",          .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Phi",          .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8SoftDropMass", .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Tau21",        .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="JetAK8Tau32",        .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("AK8",  { .fill="MaxAK8SubjetCSV",    .pfs={Stack,cut,"PassHLT","JetsAK8"},  .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  }
  
  // -------------------------------------------------------------------------
  //                              Selected W Pretags
  
  // Event
  for (const auto& cut : {"Pass1CutsS","Pass2CutsS","PassAllCutsSExcl1W","PassAllCutsTExcl1W","PassAllCutsQExcl1W","PassAllCutsWExcl1Wpre"}) {
    sh.AddHistos("evt",  { .fill="NWPreTag",           .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  }
  
  // W pre-tags (cut must be applied)
  for (const auto& cut : {"PassAllCutsS","PassAllCutsT","PassAllCutsQ","PassAllCutsW"}) {
    sh.AddHistos("Wpre", { .fill="WPreTagPtBins",      .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagPt",          .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagEta",         .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagPhi",         .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagTau21",       .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    // pt order
    sh.AddHistos("Wpre", { .fill="WPreTagPtBins",      .pfs={Stack,cut,"PassHLT","WPreTags"}, .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagPt",          .pfs={Stack,cut,"PassHLT","WPreTags"}, .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagEta",         .pfs={Stack,cut,"PassHLT","WPreTags"}, .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagPhi",         .pfs={Stack,cut,"PassHLT","WPreTags"}, .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("Wpre", { .fill="WPreTagTau21",       .pfs={Stack,cut,"PassHLT","WPreTags"}, .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  } 
  
  // -------------------------------------------------------------------------
  //                                 Selected Ws
  
  // Event
  for (const auto& cut : {"Pass1CutsS","Pass2CutsS","PassAllCutsSExcl1W","PassAllCutsTExcl1W","PassAllCutsQExcl1W","PassAllCutsWExcl1Wpre"}) {
    sh.AddHistos("evt",  { .fill="NWTag",              .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NWTag",              .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="NWTag",              .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
  }
  
  // W tags (cut applied)
  for (const auto& cut : {"PassAllCutsS","PassAllCutsT","PassAllCutsQ"}) {
    sh.AddHistos("W",    { .fill="WTagPtBins",         .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("W",    { .fill="WTagPt",             .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("W",    { .fill="WTagEta",            .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("W",    { .fill="WTagPhi",            .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    // pt order
    sh.AddHistos("W",    { .fill="WTagPtBins",         .pfs={Stack,cut,"PassHLT","Ws"},       .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("W",    { .fill="WTagPt",             .pfs={Stack,cut,"PassHLT","Ws"},       .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("W",    { .fill="WTagEta",            .pfs={Stack,cut,"PassHLT","Ws"},       .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("W",    { .fill="WTagPhi",            .pfs={Stack,cut,"PassHLT","Ws"},       .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  }
  
  
  // -------------------------------------------------------------------------
  //                            Selected tops (loose)
  
  // Event
  for (const auto& cut : {"Pass1CutsS","Pass2CutsS","PassAllCutsSExcl1W1b","PassAllCutsTExcl1W1b","PassAllCutsQExcl1W0b","PassAllCutsWExcl1Wpre0b"}) {
    sh.AddHistos("evt",  { .fill="NHadTopTag",         .pfs={Stack,cut,"PassHLT"},            .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NHadTopTag",         .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="NHadTopTag",         .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
  }
  
  
  // -------------------------------------------------------------------------
  //                                  Leptons
  
  // Event
  for (const auto& cut : {"PassAllCutsTExcl1Lep","PassAllCutsWExcl1Lep"}) {
    sh.AddHistos("evt",  { .fill="NEle",               .pfs={Stack,cut,"PassHLT"},           .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NMu",                .pfs={Stack,cut,"PassHLT"},           .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="NLep",               .pfs={Stack,cut,"PassHLT"},           .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  }
  for (const auto& cut : {"PassAllCutsSExcl0Mu","PassAllCutsQExcl0Mu"})
    sh.AddHistos("evt",  { .fill="NMuVeto",            .pfs={Stack,cut,"PassHLT"},           .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  for (const auto& cut : {"PassAllCutsSExcl0Ele","PassAllCutsQExcl0Ele"})
    sh.AddHistos("evt",  { .fill="NEleVeto",           .pfs={Stack,cut,"PassHLT"},           .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  for (const auto& cut : {"PassAllCutsSExcl0Ele0Mu","PassAllCutsQExcl0Ele0Mu"})
    sh.AddHistos("evt",  { .fill="NLepVeto",           .pfs={Stack,cut,"PassHLT"},           .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
  
  
  // -------------------------------------------------------------------------
  //                                  Gen Info
  
  
  // -------------------------------------------------------------------------
  //                                   Trigger
  

  for (const auto& iso : {"OtherUnisoLep","IsLepIsolated","AllLepIsolated"}) for (const auto& cut : {"PassTriggerPReSelNoMass","PassTriggerPreSelection","PassAllCutsSExcl1W","PassAllCutsSExclMR_R2","PassAllCutsS","PassAllCutsT","PassAllCutsQ","PassAllCutsW"}) {
    // (AK8)HT triggers
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_R2_vs_MR",              .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_WPreTag1Mass_vs_AK8HT", .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Mass_vs_AK8HT",  .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_WPreTag1Mass_vs_HT",    .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Mass_vs_HT",     .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_AK8HT",                 .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_HT",                    .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Bin",                   .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Pt",             .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_Jet1AK8Mass",           .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_TrimMass50_vs_WPreTag1Mass",          .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_R2_vs_MR",                            .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_WPreTag1Mass_vs_HT",                  .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,1, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_Jet1AK8Mass_vs_HT",                   .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,1, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_HT",                                  .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_Bin",                                 .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_Jet1AK8Pt",                           .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_Jet1AK8Mass",                         .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_WPreTag1Mass",                        .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    
    //sh.AddHistos("evt", { .fill="HLTEff_PFHT750_4JetPt50_vs_R2_vs_MR",                   .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_PFHT750_4JetPt50_vs_HT",                         .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_PFHT750_4JetPt50_vs_Bin",                        .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    
    // AK8/B2G triggers
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_R2_vs_MR",                        .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt_vs_Jet1AK8Mass",        .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt_vs_WPreTag1Mass",       .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_WPreTag1Mass",                    .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt",                       .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt",                       .pfs={"Triggers",cut,iso,"Jet1AK8Pt>450"},  .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Pt",                       .pfs={"Triggers",cut,iso,"Jet1AK8Pt>500"},  .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Mass",                     .pfs={"Triggers",cut,iso},                  .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Jet1AK8Mass",                     .pfs={"Triggers",cut,iso,"Jet1AK8Mass>65"}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_HT",                              .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,iso,"Jet1AK8Pt>450"},  .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,iso,"Jet1AK8Pt>500"},  .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,iso,"Jet1AK8Mass>65"}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,iso,"Jet1AK8Mass>65","Jet1AK8Pt>450"}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet360_vs_Bin",                             .pfs={"Triggers",cut,iso,"Jet1AK8Mass>65","Jet1AK8Pt>500"}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_R2_vs_MR",                        .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_WPreTag1Mass",                    .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_HT",                              .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFJet450_vs_Bin",                             .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    
    //sh.AddHistos("evt", { .fill="HLTEff_AK8DiPFJet250_200_vs_R2_vs_MR",                  .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8DiPFJet250_200_vs_WPreTag1Mass",              .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });
    //sh.AddHistos("evt", { .fill="HLTEff_AK8DiPFJet250_200_vs_Bin",                       .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,0, 0.45,0.45} });

    // Suggestion: AK8PFHT700 OR PFHT800
    //for (auto trigger_comb : std::vector<std::string>({"HLTEff_AK8PFHT700orPFHT800", "HLTEff_PFJet450orPFHT800", "HLTEff_AK8PFJet450orPFHT800", "HLTEff_AK8PFJet360orPFHT800","HLTEff_AK8PFJet360orAK8PFHT700"})) {
    for (auto trigger_comb : std::vector<std::string>({"HLTEff_AK8PFJet450orPFHT800", "HLTEff_AK8PFJet360orPFHT800"})) {
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_R2_vs_MR",                .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      //sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt_vs_AK8HT",         .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      //sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass_vs_AK8HT",       .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      //sh.AddHistos("evt", { .fill=trigger_comb+"_vs_WPreTag1Mass_vs_AK8HT",      .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt_vs_HT",         .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass_vs_HT",       .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_WPreTag1Mass_vs_HT",      .pfs={"Triggers",cut,iso}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,0, 0,0, 0,1} });
      //sh.AddHistos("evt", { .fill=trigger_comb+"_vs_AK8HT",                   .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_HT",                      .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Pt",               .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Jet1AK8Mass",             .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_WPreTag1Mass",            .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
      sh.AddHistos("evt", { .fill=trigger_comb+"_vs_Bin",                     .pfs={"Triggers",cut,iso}, .cuts={}, .draw="PE1",  .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    }
  }
  
  
  // -------------------------------------------------------------------------
  //                           Signal Event Variables
  
  for (const auto& cut : {"Pass1CutsS","Pass2CutsS","Pass3CutsS","Pass4CutsS","Pass5CutsS","Pass6CutsS","PassAllCutsSExclMR_R2","PassAllCutsS"}) {
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={Stack,cut,"PassHLT"},                .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={Stack,cut,"PassHLT"},                .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={Stack,cut,"PassHLT"},                .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={Stack,cut,"PassHLT"},                .cuts={},.draw=d,.opt=o_stk,.ranges=r_stk});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={cut,"PassHLT","Signals,Background"}, .cuts={},.draw="COLZ",.opt="Sumw2",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT_vs_MR",           .pfs={cut,"PassHLT","Signals,Background"}, .cuts={},.draw="COLZ",.opt="Sumw2",.ranges={}});
    // MGlunio/MStop plots
    sh.AddHistos("evt",  { .fill="MLSP_vs_MGluino",    .pfs={cut,"PassHLT","GluinoSignalScans"},  .cuts={},.draw="COLZ",.opt="Sumw2",.ranges={}});
    sh.AddHistos("evt",  { .fill="MLSP_vs_MStop",      .pfs={cut,"PassHLT","StopSignalScans"},    .cuts={},.draw="COLZ",.opt="Sumw2",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="HT",                 .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="MET",                .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="MR",                 .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={"MGluinoPoints",cut,"PassHLT","GluinoSignalScans"}, .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2",                 .pfs={"MStopPoints",  cut,"PassHLT","StopSignalScans"},   .cuts={},.draw=d,.opt="Sumw2Norm",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={cut,"PassHLT","GluinoSignalScans","MGluinoPoints"}, .cuts={},.draw="COLZ",.opt="Sumw2",.ranges={}});
    sh.AddHistos("evt",  { .fill="R2_vs_MR",           .pfs={cut,"PassHLT","StopSignalScans",  "MStopPoints"},   .cuts={},.draw="COLZ",.opt="Sumw2",.ranges={}});
  }
  
  
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
      sh.AddHistos("evt",        { .fill="NJet",            .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("evt",        { .fill="NHadTopTag",      .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"PassHLT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPhi",          .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPhi",          .pfs={Stack,cut,"PassHLT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetEta",          .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetEta",          .pfs={Stack,cut,"PassHLT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetTau32",        .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetTau32",        .pfs={Stack,cut,"PassHLT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetSoftDropMass", .pfs={Stack,cut,"PassHLT"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetSoftDropMass", .pfs={Stack,cut,"PassHLT","PtOrder"},     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      // Subjet BTags
      sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,"NSubjetBTag",cut,"PassHLT"},                       .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,"NSubjetBTag",cut,"PassHLT","PtOrder"},       .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"PassHLT","JetPassSubjetBTag"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetPt",           .pfs={Stack,cut,"PassHLT","JetPassSubjetBTag","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetBTagCSV",      .pfs={Stack,cut,"PassHLT"},                                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("AK8",      { .fill="JetBTagCSV",      .pfs={Stack,cut,"PassHLT","PtOrder"},                     .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      
      // Cutflow
      sh.AddHistos("cutflow", { .fill="Cutflow", .pfs={Stack,cut,"PassHLT"}, .cuts={"PassAnaSelection"}, .draw=d, .opt=o_stk, .ranges={0,0, 0,0} });
      sh.AddHistos("cutflow", { .fill="Cutflow", .pfs={"Background", cut,"PassHLT"}, .cuts={"PassAnaSelection"}, .draw=d, .opt="Sumw2Log",     .ranges={0,0, 0,0} });
    }
    
    // N-1 plots
    sh.AddHistos("evt", { .fill="Jet1Eta",     .pfs={Stack,"AllCutsExclJet1Eta",   "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet2Eta",     .pfs={Stack,"AllCutsExclJet2Eta",   "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet1Pt",      .pfs={Stack,"AllCutsExclJet1Pt",    "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet2Pt",      .pfs={Stack,"AllCutsExclJet2Pt",    "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet1Mass",    .pfs={Stack,"AllCutsExclJet1Mass",  "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet2Mass",    .pfs={Stack,"AllCutsExclJet2Mass",  "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet1BTagCSV", .pfs={Stack,"AllCutsExclJet1BTag",  "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet2BTagCSV", .pfs={Stack,"AllCutsExclJet2BTag",  "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet1Tau32",   .pfs={Stack,"AllCutsExclJet1Tau32", "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="Jet2Tau32",   .pfs={Stack,"AllCutsExclJet2Tau32", "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="DPhi",        .pfs={Stack,"AllCutsExclDeltaPhi",  "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt", { .fill="R",           .pfs={Stack,"AllCutsExclR",         "PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    
    // --------------------------------------------------------------------------
    //                                 Trigger

    sh.AddHistos("evt", { .fill="GenHt",      .pfs={"AllSamples"},                           .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("evt", { .fill="AK4PuppiHt", .pfs={"AllSamples"},                           .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("evt", { .fill="AK8PuppiHt", .pfs={"AllSamples"},                           .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("evt", { .fill="GenHt",      .pfs={"AllSamples","PassHLT"},              .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("evt", { .fill="AK4PuppiHt", .pfs={"AllSamples","PassHLT"},              .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("evt", { .fill="AK8PuppiHt", .pfs={"AllSamples","PassHLT"},              .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("evt", { .fill="GenHt",      .pfs={"AllSamples","Pass12Cuts","PassHLT"},  .cuts={},  .draw=d, .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });

    // Event composition - under different cuts
    for (const auto& cut : all_cuts) {
      sh.AddHistos("evt", { .fill="SumPt",      .pfs={Stack,cut},              .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("evt", { .fill="SumPt",      .pfs={Stack,cut,"PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("evt", { .fill="AK8PuppiHt", .pfs={Stack,cut},              .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("evt", { .fill="AK8PuppiHt", .pfs={Stack,cut,"PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("evt", { .fill="AK4PuppiHt", .pfs={Stack,cut},              .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
      sh.AddHistos("evt", { .fill="AK4PuppiHt", .pfs={Stack,cut,"PassHLT"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });

      sh.AddHistos("evt", { .fill="MetPt",       .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
      sh.AddHistos("evt", { .fill="MinDeltaPhi", .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
      sh.AddHistos("evt", { .fill="MR",      .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
      sh.AddHistos("evt", { .fill="MTR",     .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
      sh.AddHistos("evt", { .fill="R",       .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("evt", { .fill="DPhi",    .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("evt", { .fill="DEta",    .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("evt", { .fill="DR",      .pfs={Stack,cut,"PassHLT"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    }
    sh.AddHistos("evt",   { .fill="SumPt",      .pfs={Stack,"Pass12Cuts","PassHLT","PassDPhiCut"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt",   { .fill="AK8PuppiHt", .pfs={Stack,"Pass12Cuts","PassHLT","PassDPhiCut"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("evt",   { .fill="AK4PuppiHt", .pfs={Stack,"Pass12Cuts","PassHLT","PassDPhiCut"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });

    // N-1 plot for all cuts
    

    // --------------------------------------
    //               Efficiencies

    // 1D plots - sumpt - cutflow
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass11Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass11Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK8PuppiHt", .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass11Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_AK4PuppiHt", .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });

    // Check other triggers
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass5Cuts"},  .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass7Cuts"},  .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass9Cuts"},  .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass11Cuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass12Cuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass12Cuts","PassDPhiCut"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass5Cuts"},  .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass7Cuts"},  .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass9Cuts"},  .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass11Cuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass12Cuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass12Cuts","PassDPhiCut"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });

    // 1 Bin
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0.95,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0.95,1} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPtOneBin",          .pfs={"Triggers","Pass12Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("evt", { .fill="HLTEff_PFHT800_vs_SumPtOneBin",          .pfs={"Triggers","Pass12Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });

    // 2D plots
    // sumpt vs mass1
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass1Cuts"},                .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass3Cuts"},                .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts"},                .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass7Cuts"},                .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass9Cuts"},                .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass11Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("evt", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass12Cuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });

  
    // Main variables (Shape, area, Data-MC agreement)
    // Hadronic top selection
    //   // No Cut
    //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
    //   sh.AddHistos("AK8", { .fill="JetTau2",         .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetTau3",         .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} }); // Note
    //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note
    //   // Apply 1 Cut
    //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack,"JetMassCut"},  .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack,"JetTau32Cut"}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack,"JetPtCut"},    .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack,"JetMassCut"},  .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack,"JetTau32Cut"}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack,"JetPtCut"},    .cuts={},  .draw=d, .opt=o_stk, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
    //   // Apply 2 Cuts (N-1)
    //   sh.AddHistos("AK8", { .fill="JetPt",           .pfs={Stack,"JetMassCut","JetTau32Cut"}, .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //   sh.AddHistos("AK8", { .fill="JetTau32",        .pfs={Stack,"JetPtCut","JetMassCut"},    .cuts={},  .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });  // Note
    //   sh.AddHistos("AK8", { .fill="JetSoftDropMass", .pfs={Stack,"JetTau32Cut","JetPtCut"},   .cuts={},  .draw=d, .opt=o_stk, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note

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
    sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} }); // No Cut
    sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau1",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau2",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau3",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau21",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau31",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    //sh.AddHistos("AK8", { .fill="JetNSubJets",                           .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} }); // 1 Cut
    sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("AK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetMassCut","JetTau32Cut"}, .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} }); // 2 Cuts
    sh.AddHistos("AK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw=d, .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("AK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={},  .draw=d, .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",           .pfs={"JetGenTruth","Signals,Background"},             .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 2D No Cut
    sh.AddHistos("AK8", { .fill="JetSoftDropMass_vs_JetTau32",           .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 1 Cut

    // Fraction of merged sub-jets
    sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT"},                                   .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT","JetMatchedGenTopType"},            .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"Signals,TT"},                                   .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("AK8", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"Signals,TT","JetMatchedGenTopType"},            .cuts={"JetHasMatchedGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
  
    // Top Tag/Finding Efficiency
    sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT"},              .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });
    sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT","TopSizeCut"}, .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });
    sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"Signals,TT"},              .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });
    sh.AddHistos("AK8", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"Signals,TT","TopSizeCut"}, .cuts={"JetHasMatchedGenTop"}, .draw=d, .opt="", .ranges={0,0, 0,1} });

    // --------------------------------------------------------------------------
    //                              Gen particles

    // Jet Finding Efficiency
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });  // Note
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw=d,     .opt="",        .ranges={0,0, 0,1} });

    // --------------------------------------------------------------------------
    //                                   HLT

    // // Trigger Efficiencies
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPt",       .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={800,2000, 0.8,1.05, 0.15,0.95} });
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPtOneBin", .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtBins",        .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,2000,   0.8,1.05, 0.15,0.95} });
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtOneBin",      .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });

    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPt",       .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={800,2000, 0.8,1.05, 0.15,0.95} });
    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPtOneBin", .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });
    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtBins",        .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,2000,   0.8,1.05, 0.15,0.95} });
    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtOneBin",      .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });

    // --------------------------------------------------------------------------
    //                                MET/Razor

    // MET, Razor Variables
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts","RBands"},           .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={Stack,"Tau32Cuts","RBands","DPhiBands"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    //sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    //sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MR",               .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MR",               .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MR",               .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MTR",              .pfs={Stack},                              .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MTR",              .pfs={Stack,"Tau32Cuts"},                    .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MTR",              .pfs={Stack,"Tau32Cuts","DPhiBands"},          .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1e-3,1e6} });
  
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
  
    // Signal selection (Apply loose pretag selection)
    sh.AddHistos("baseline events",   { .fill="R",                 .pfs={Stack},                     .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={Stack},                     .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="R",                 .pfs={Stack,"Tau32Cuts"},           .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={Stack,"Tau32Cuts"},           .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="R",                 .pfs={Stack,"Tau32Cuts","DPhiBands"}, .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={Stack,"Tau32Cuts","RBands"},  .cuts={}, .draw=d, .opt=o_stk, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  
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
    sh.AddHistos("baseline events",   { .fill="MLSP_vs_MGluino", .pfs={"DPhiBands","Tau32Cuts","RBands","Signals"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={} });
  } // End all plots

  // --------------------------------------------------------------------------
  //                             Systematics

  // Data-MC Comparison - background composition
  std::vector<std::string> cuts = { "Pass5Cuts", "Pass11Cuts", "Pass12Cuts", "Pass13Cuts", "Pass15Cuts", "Pass16Cuts" };
  for (auto cut : cuts) { 
    sh.AddHistos("all events syst",   { .fill="Counts_vs_NJet",            .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("all events syst",   { .fill="Counts_vs_NHadTopTag",      .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPt",           .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPt",           .pfs={Stack,cut,"PassHLT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPhi",          .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPhi",          .pfs={Stack,cut,"PassHLT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetEta",          .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetEta",          .pfs={Stack,cut,"PassHLT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetTau32",        .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetTau32",        .pfs={Stack,cut,"PassHLT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetSoftDropMass", .pfs={Stack,cut,"PassHLT"},                 .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetSoftDropMass", .pfs={Stack,cut,"PassHLT","PtOrder"}, .cuts={},  .draw=d, .opt=o_stk, .ranges=r_stk });
  }

  // Background estimation
  //sh.AddHistos("baseline events",      { .fill="HtAll_vs_DPhiFine_vs_RFine",      .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  sh.AddHistos("baseline events",      { .fill="AK8PuppiHt_vs_DPhiFine_vs_RFine", .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
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
    while(d.jetsAK4Puppi.Loop()) if (passLooseJet   [d.jetsAK4Puppi.it]) sh.Fill("AK4");
    while(d.jetsAK8Puppi.Loop()) if (passTightJetAK8[d.jetsAK8Puppi.it]) sh.Fill("AK8");
    while(d.jetsAK4Puppi.Loop()) if (passMediumBTag [d.jetsAK4Puppi.it]) sh.Fill("b");
    while(d.jetsAK4Puppi.Loop()) if (passLooseBTag  [d.jetsAK4Puppi.it]) sh.Fill("b loose");
    while(d.jetsAK8Puppi.Loop()) if (passWPreTag    [d.jetsAK8Puppi.it]) sh.Fill("Wpre");
    while(d.jetsAK8Puppi.Loop()) if (passTightWTag  [d.jetsAK8Puppi.it]) sh.Fill("W");
    while(d.ele.Loop())          if (passEleTight   [d.ele.it])          sh.Fill("ele");
    while(d.ele.Loop())          if (passEleVeto    [d.ele.it])          sh.Fill("ele veto");
    while(d.mu.Loop())           if (passMuTight    [d.mu.it])           sh.Fill("mu");
    while(d.mu.Loop())           if (passMuVeto     [d.mu.it])           sh.Fill("mu veto");
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
