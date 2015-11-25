#define TEST 0
#if TEST == 1
#define NTHSTAT 10
#else
#define NTHSTAT 1
#endif

#define VTX_REWEIGHT_BKG_EST 0

#define R_TYPE 1 // 1: AK8, 2: AK4, 3: Top-Pair

// ----------- DEFAULT SETTINGS -----------

#define NEW_TOP_DEF 1

/* 
   Top Tagging working points:
   https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging#13_TeV_working_points
   e(B) = 3% SD WP2 	3% 	  	SD(beta=0,z=0.1, R=0.8) [110,210] 	tau32 < 0.75, max SD subjet b-discriminant > 0.39 

   --> Use WP2
   --> Use SoftDrop with PUPPI

*/
#define SD_MASS_CUT_LOW 110
#define SD_MASS_CUT_HIGH 210
#define TAU32_CUT_NEW 0.75
#define PT_CUT_NEW 400
#define R_CUT_NEW 0.4
#define R_CUT_LOW 0.2
#define DPHI_CUT_NEW 2.7

#define PRUNED_MASS_CUT_LOW 140
#define TAU32_CUT_OLD 0.75
#define PT_CUT_OLD 400
#define R_CUT_OLD 0.4
#define DPHI_CUT_OLD 2.8

// NNLO k-factors, QCD Scale
#define QCD_SCALE_FACTOR 1 // Do not rescale QCD

// Integrated Luminositied in uits of fb^-1
#define Data_IntLumi_invfb  0.01634

// Cross section weighting for MC
#define MC_Expected_IntLumi_invfb 4.0

// Vertex reweighting
#define MAX_NVTX 100
#define VTX_REWEIGHT_MC 1

#define HIGHEST_PT_JETS_ONLY 1

// ------------------------------------------


#include <cstdlib>
#include <unistd.h>
#include <vector>
#include "../interface/Samples.h"
#include "../interface/SmartHistos.h"
#include "../plugins/B2GTreeReader.cc"
#include "../plugins/B2GTreeLooper.cc"

int main(int argc, char* argv[]) {
  bool debug = 0;
  // Get arguments from shell
  std::vector<std::string> filelist_fromshell;
  std::string inputfile="";
  std::string outputfile="plots.root";
  // -o <output file> option:
  // Specify the output root filename
  
  // Rest of the arguments are treated as files added
  // Note:
  // If using postfixes with the v.pf_file_add variable
  // each added file will increase this variable so when using *
  // add ""-s so instead of the shell TChain will parse the argument
  bool is_i = false;
  bool is_o = false;
  Samples samples;
  
  for(int i=1; i<argc; i++) {
    std::string arg = argv[i];
    if (arg[0]=='-'&&arg.size()==2) { 
      is_i = (arg[1]=='i'); 
      is_o = (arg[1]=='o'); 
    } else if (is_i) {
      inputfile=arg; is_i=0; 
    } else if (is_o) {
      outputfile=arg; is_o=0; 
    } else {
      filelist_fromshell.push_back(arg);
    }
  }
  bool Run = inputfile.size()==0;
  bool test = TEST;
  bool all = 0;
  const size_t ndata = TEST ? 1 : all ? 6 : 2;
  const size_t nsignal = TEST ? 0 : all ? 6 : 2;
  const size_t nttbar = TEST ? 1 : all ? 7 : 7;
  if (filelist_fromshell.size()) {
    samples.AddSample("test", "test", "1", { { .dir="", .scale_factor=1 } });
  } else {
    if (test) {
      //samples.AddSample("test", "T5ttttDeg (#tilde{g} #rightarrow t + #tilde{t}_{4-body decay} )", { { .dir="../../../B2GTTreeNtupleExtra_susy.root", .xsec_pb=0.0460525 } });
      //samples.AddSample("TTJets_NLO",         "t#bar{t}+jets (NLO)", "633", { { .dir="/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Jun16_edm_Jun03/TTJets_amcatnloFXFX-pythia8/*.root", .xsec_pb=831.76 } }); // Red
      //samples.AddSample("TTJets_LO_madgraph", "t#bar{t}+jets (LO)",  "601", { { .dir="/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Jun16_edm_Jun03/TTJets_madgraphMLM-pythia8/*.root", .xsec_pb=831.76 } }); // Blue
      // Phys14 comparison
      // std::string sample_dir = "/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Jul24_edm_Jun03/";
      // samples.AddSample("TTJets_RunII", "t#bar{t}+jets", "633", { { .dir=sample_dir+"TTJets_madgraphMLM-pythia8/*.root", .xsec_pb=831.76 } }); // Red
      // samples.AddSample("WJets_RunII",  "W+Jets", "417", // Green
      //   		{ { .dir=sample_dir+"WJetsToLNu_HT-100To200/*.root", .xsec_pb=2234.9 },
      //   		  { .dir=sample_dir+"WJetsToLNu_HT-200To400/*.root", .xsec_pb=580.1 }, 
      //   		  { .dir=sample_dir+"WJetsToLNu_HT-400To600/*.root", .xsec_pb=68.4 },
      //       	          { .dir=sample_dir+"WJetsToLNu_HT-600ToInf/*.root", .xsec_pb=23.14 } });
      // samples.AddSample("TTJets", "t#bar{t}+jets - PHYS14", "901", { { .dir="", .xsec_pb=1 } }); // Pink
      // samples.AddSample("WJets", "W+Jets - PHYS14", "841", { { .dir="", .xsec_pb=1 } }); // Teal
      // test
      //samples.AddSample("Data",         "Data", "1", { { .dir="/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/CMSSW_7_4_7_patch2/src/B2GTTreeNtupleExtra_Data.root", .xsec_pb=1 } });
      //samples.AddSample("JetHT",          "Data",                        "1",/*Black*/   { { .dir="/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Oct03_edm_Oct02/JetHT_2015D.root", .xsec_pb=0} });
      //samples.AddSample("TTJetsMadgraph", "t#bar{t}+Jets (Madpgraph)", "632",/*Red*/     { { .dir="/data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple/Oct03_edm_Oct02/TTJets_madgraph.root", .xsec_pb=0} });
    } else {
      // Colors
      // 400 kYellow    800 kOrange
      // 416 kGreen	820 kSpring
      // 432 kCyan	840 kTeal
      // 600 kBlue	860 kAzure
      // 616 kMagenta	880 kViolet
      // 632 kRed    	900 kPink
      std::string sample_dir = "Nov18_ntuple/";
      // First data is plotted on most plots
      if (ndata>=1) samples.AddSample("JetHT_25ns",      "Data",                              "1",/*Black*/   { { .dir=sample_dir+"JetHT_25ns_*.root", .scale_factor=1} });
      if (ndata>=2) samples.AddSample("MET_25ns",        "MET Data",                          "1",/*Black*/   { { .dir=sample_dir+"JetHT_25ns_*.root", .scale_factor=1} });
      if (ndata>=3) samples.AddSample("SingleEle_25ns",  "SingleElectron Data",             "803",/*DBrown*/  { { .dir=sample_dir+"SingleElectron_25ns_*.root", .scale_factor=1} });
      if (ndata>=4) samples.AddSample("SingleMuon_25ns", "SingleMuon Data",                 "434",/*DCyan*/   { { .dir=sample_dir+"SingleMuon_25ns_*.root", .scale_factor=1} });
      if (nsignal>=1) samples.AddSample("T1tttt",         "T1tttt",                         "862",/*Azure*/   { { .dir=sample_dir+"T1tttt_mGluino_*.root", .scale_factor=1} });
      if (nsignal>=2) samples.AddSample("T1tttt_FullSim", "T1tttt FullSim"                  "902",/*Pink*/    { { .dir=sample_dir+"T1tttt_FullSim_mGluino-1500_mLSP-100.root", .scale_factor=1} });
      //if (nsignal>=2) samples.AddSample("T5ttttDeg_mGo1000_4bodydec",     "T5ttttDeg (M_{#tilde{g}}=1TeV)",                       "1",/*Black*/    { { .dir=sample_dir+"T5ttttDeg_mGo1000_4bodydec.root", .scale_factor=1} });
      if (nttbar>=1) samples.AddSample("TTJetsHT600",     "t#bar{t}+Jets (HT600)",          "634",/*DRed*/   { { .dir=sample_dir+"TTJets_HT-*.root",             .scale_factor=1} });
      if (nttbar>=2) samples.AddSample("TTHerwig",        "t#bar{t} Herwig",                "801",/*Orange*/ { { .dir=sample_dir+"TT_herwigpp.root",             .scale_factor=1} });
      if (nttbar>=3) samples.AddSample("TTPowhegPythia8", "t#bar{t} Powheg-Pythia8",        "901",/*Pink*/   { { .dir=sample_dir+"TT_powheg_pythia8*.root",      .scale_factor=1} }); //!!ext
      if (nttbar>=4) samples.AddSample("TTJetsNLO",       "t#bar{t}+Jets MC@NLO",           "617",/*Magent*/ { { .dir=sample_dir+"TTJets_amcatnloFXFX.root",     .scale_factor=1} });
      if (nttbar>=5) samples.AddSample("TTNLO",           "t#bar{t} MC@NLO",                "619",/*DMagen*/ { { .dir=sample_dir+"TT_amcatnlo.root",             .scale_factor=1} });
      if (nttbar>=6) samples.AddSample("TTJetsMadgraph",  "t#bar{t}+Jets (Madpgraph)",      "632",/*Red*/    { { .dir=sample_dir+"TTJets_madgraph.root",         .scale_factor=1} });
      if (nttbar>=7) samples.AddSample("TTJetsMadgraphFS","t#bar{t}+Jets (Madpgraph - FS)", "903",/*DPink*/ { { .dir=sample_dir+"TTJets_madgraph_FastSim.root", .scale_factor=1} });
      samples.AddSample("WJets",          "W+Jets",                                    "418",/*Green*/  { { .dir=sample_dir+"WJetsToLNu_HT*.root",          .scale_factor=1} });
      samples.AddSample("ZJetsToNuNu",    "Z+Jets (Z#rightarrow#nu#nu)",               "401",/*Yellow*/ { { .dir=sample_dir+"ZJetsToNuNu_HT*.root",         .scale_factor=1} });
      samples.AddSample("DYJetsToLL",     "DY+Jets (Z/#gamma#rightarrowll)",           "403",/*DYellow*/{ { .dir=sample_dir+"DYJetsToLL_M-50_HT*.root",     .scale_factor=1} });
      samples.AddSample("QCD",            "QCD",                                         "4",/*Blue*/   { { .dir=sample_dir+"QCD_HT*.root",                 .scale_factor=QCD_SCALE_FACTOR} });
      samples.AddSample("TTV",            "t#bar{t}+V",                                "803",/*Brown*/  { { .dir=sample_dir+"TTZToQQ.root",                 .scale_factor=1},
                                                                                                          { .dir=sample_dir+"TTZToLLNuNu.root",             .scale_factor=1},
                                                                                                          { .dir=sample_dir+"TTWJetsToQQ.root",             .scale_factor=1},
                                                                                                          { .dir=sample_dir+"TTWJetsToLNu.root",            .scale_factor=1} });
      samples.AddSample("Diboson",        "Diboson",                                   "804",/*DOrange*/{ { .dir=sample_dir+"WWTo4Q_4f.root",               .scale_factor=1},
                                                                                                          { .dir=sample_dir+"WWTo1L1Nu2Q.root",             .scale_factor=1},
                                                                                                          { .dir=sample_dir+"WWTo2L2Nu.root",               .scale_factor=1},
                                                                                                          { .dir=sample_dir+"WZTo1L1Nu2Q.root",             .scale_factor=1},
                                                                                                          { .dir=sample_dir+"WZTo1L3Nu.root",               .scale_factor=1},
                                                                                                          { .dir=sample_dir+"WZTo2L2Q.root",                .scale_factor=1},
                                                                                                          { .dir=sample_dir+"WZTo3LNu.root",                .scale_factor=1},
                                                                                                          { .dir=sample_dir+"ZZTo4Q.root",                  .scale_factor=1},
                                                                                                          { .dir=sample_dir+"ZZTo2Q2Nu.root",               .scale_factor=1},
                                                                                                          { .dir=sample_dir+"ZZTo2L2Nu.root",               .scale_factor=1},
                                                                                                          { .dir=sample_dir+"ZZTo2L2Q.root",                .scale_factor=1},
													  { .dir=sample_dir+"ZZTo4L.root",                  .scale_factor=1} });
      samples.AddSample("SingleTop",      "Single top",                                "403",/*DYellow*/{ { .dir=sample_dir+"ST_s-channel_4f_leptonDecays.root",     .scale_factor=1},
                                                                                                          { .dir=sample_dir+"ST_t-channel_4f_leptonDecays*.root",    .scale_factor=1}, //!!ext
                                                                                                          { .dir=sample_dir+"ST_tW_top_5f_inclusiveDecays.root",     .scale_factor=1},
													  { .dir=sample_dir+"ST_tW_antitop_5f_inclusiveDecays.root", .scale_factor=1} }); 
      //samples.AddSample("TZQ",            "top+Z+Q",                                   "882",/*DViolet*/{ { .dir=sample_dir+"tZq_ll.root",                  .scale_factor=3.763389},
      //                                                                                                    { .dir=sample_dir+"tZq_nunu.root",                .scale_factor=3.820400} });
      // Zero event in: R>0.4
      //samples.AddSample("ZJetsToQQ",      "Z+Jets (Z#rightarrowq#bar{q})",       "403",/*DYellow*/ { { .dir=sample_dir+"ZJetsToQQ_HT600toInf.root", .scale_factor=1} });
      //samples.AddSample("GJets",          "#gamma+Jets",                         "617",/*Magenta*/ { { .dir=sample_dir+"GJets_HT-40ToInf.root", .scale_factor=1} });
      
      /* 
	 NLO Corr factors - Have to renormalize NLO events to the true number of events due to negative weights
	 Corr factor is: sum ( abs(evt_weight) ) / sum ( evt_weight ) - For unfiltered events
	 
         TT_herwigpp_TuneEE5C                    Fraction of negative weights: 21.92%   Corr factor: 1.780401
         TTJets_amcatnloFXFX                     Fraction of negative weights: 33.42%   Corr factor: 3.015442
         TT_amcatnlo                             Fraction of negative weights: 25.43%   Corr factor: 2.035141
	 
         TTZToQQ                                 Fraction of negative weights: 26.57%   Corr factor: 2.133763
         TTZToLLNuNu                             Fraction of negative weights: 26.76%   Corr factor: 2.151468
         TTWJetsToQQ                             Fraction of negative weights: 24.20%   Corr factor: 1.937964
         TTWJetsToLNu                            Fraction of negative weights: 24.33%   Corr factor: 1.947693
	 
         WZTo1L1Nu2Q                             Fraction of negative weights: 20.97%   Corr factor: 1.722645
         WZTo1L3Nu                               Fraction of negative weights: 22.33%   Corr factor: 1.807001
         WZTo2L2Q                                Fraction of negative weights: 20.06%   Corr factor: 1.669946
         ZZTo4Q                                  Fraction of negative weights: 19.06%   Corr factor: 1.616203
         ZZTo2Q2Nu                               Fraction of negative weights: 19.30%   Corr factor: 1.628800
         ZZTo2L2Q                                Fraction of negative weights: 18.43%   Corr factor: 1.583833
	 
         tZq_ll                                  Fraction of negative weights: 36.71%   Corr factor: 3.763389
         tZq_nunu                                Fraction of negative weights: 36.91%   Corr factor: 3.820400
	 
         ST_s-channel_4f_leptonDecays            Fraction of negative weights: 18.84%   Corr factor: 1.604867
      */
    }
  }
  std::vector<size_t > dir_to_index = samples.GetDirToIndex();
  std::vector<double> sample_scale_factors = samples.GetKFactors();
  if (debug) std::cout<<"Samples ok\n";
  
  // Initialize TreeReader
  B2GTreeReader reader;
  
  // Class to Loop on files and read the Trees
  B2GTreeLooper looper(NTHSTAT,1);
  
  // Data variable
  Data d;
  
  // Histogram storage class
  SmartHistos sh;
  sh.AddHistoType("mu");
  sh.AddHistoType("ele");
  sh.AddHistoType("jetsAK4");
  sh.AddHistoType("jetsAK8");
  sh.AddHistoType("met");
  sh.AddHistoType("evt");
  sh.AddHistoType("gen");
  
  // --------------------------------------------------------------------
  //                             CUTS
  // --------------------------------------------------------------------
  
  // Sample
  
  // Ele
  sh.AddNewCut("GoodEle",       [&d](){ return d.ele.Pt[d.ele.it] > 35 && fabs(d.ele.Eta[d.ele.it]) < 2.5 && d.ele.isTight[d.ele.it] > 0; });
  sh.AddNewCut("EleJetCombMass>90", [&d](){ return d.evt.EleJetCombMass[d.ele.it] > 90; });
  // Mu
  sh.AddNewCut("GoodMu",        [&d](){ return d.mu.Pt[d.mu.it] > 45 && fabs(d.mu.Eta[d.mu.it]) < 2.1 && d.mu.IsTightMuon[d.mu.it] > 0; });
  sh.AddNewCut("MuJetCombMass>90",  [&d](){ return d.evt.MuJetCombMass[d.mu.it] > 90; });
  
  // Jets
  //sh.AddNewCut("AK4Highest2Jet",   [&d](){ return d.jetsAK4.size>=2 && d.jetsAK4.it<2; });
  //sh.AddNewCut("AK4Highest3Jet",   [&d](){ return d.jetsAK4.size>=3 && d.jetsAK4.it<3; });
  //sh.AddNewCut("Highest2Jet",   [&d](){ return d.jetsAK8.size>=2 && d.jetsAK8.it<2; });
  //sh.AddNewCut("Highest3Jet",   [&d](){ return d.jetsAK8.size>=3 && d.jetsAK8.it<3; });
  //sh.AddNewCut("HadTop",           [&d](){ return d.jetsAK8.Pt[d.jetsAK8.it] > PT_CUT_OLD && d.jetsAK8.prunedMass[d.jetsAK8.it] > PRUNED_MASS_CUT_LOW && (d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it]) < TAU32_CUT_OLD; });
  //sh.AddNewCut("HadTopNoPtCut",    [&d](){ return d.jetsAK8.prunedMass[d.jetsAK8.it] > PRUNED_MASS_CUT_LOW && (d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it]) < TAU32_CUT_OLD; });
  //sh.AddNewCut("HadTopNoTauCut",   [&d](){ return d.jetsAK8.Pt[d.jetsAK8.it] > PT_CUT_OLD && d.jetsAK8.prunedMass[d.jetsAK8.it] > PRUNED_MASS_CUT_LOW; });
  //sh.AddNewCut("HadTopNoMassCut",  [&d](){ return d.jetsAK8.Pt[d.jetsAK8.it] > PT_CUT_OLD && (d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it]) < TAU32_CUT_OLD; });
  //sh.AddNewCut("HadTopNoMassNoTauCut", [&d](){ return d.jetsAK8.Pt[d.jetsAK8.it] > PT_CUT_OLD; });
  
  // Gen Jets
  sh.AddNewCut("JetHasMatchedGenTop",  [&d](){ return d.jetsAK8.HasNearGenTop[d.jetsAK8.it]==1; });
  sh.AddNewCut("GenLepTop",            [&d](){ return d.evt.JetMatchedGenTopType[d.jetsAK8.it]==1; });
  sh.AddNewCut("JetIsHadTopTagged",    [&d](){ return d.evt.JetIsHadTopTagged[d.jetsAK8.it]; });
  
  // Gen Particles
  sh.AddNewCut("IsGenTop",         [&d](){ return d.evt.IsGenTop[d.gen.it]; });
  
  // Events
  //sh.AddNewCut("NTop==2",          [&d](){ return d.evt.NTopHad+d.evt.NTopLep==2; });
  //sh.AddNewCut("NTopLep==1",       [&d](){ return d.evt.NTopLep==1; });
  //sh.AddNewCut("NLepTight==0",     [&d](){ return (d.evt.nmu+d.evt.neletight)==0; });
  //sh.AddNewCut("NLepVeto==0",      [&d](){ return (d.evt.nmuveto+d.evt.neleveto)==0; });
  sh.AddNewCut("NLepTight==1",     [&d](){ return (d.evt.neletight+d.evt.nmu)==1; });
  sh.AddNewCut("NEleTight==1",     [&d](){ return d.evt.neletight==1; });
  sh.AddNewCut("NMuTight==1",      [&d](){ return d.evt.nmu==1; });
  
  //sh.AddNewCut("NGenLepFromTop==0",[&d](){ return d.evt.NGenLepFromTop==0; });

  sh.AddNewCut("NTopHadPreTag>=2",   [&d](){ return d.evt.NTopHadPreTag>=2; });
  
  sh.AddNewCut("HadTrigger",       [&d](){ return d.evt.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50; });
  sh.AddNewCut("LowHTTrig",        [&d](){ return d.evt.HLT_PFHT600; });
  sh.AddNewCut("NoDataSignalRegion", [&looper,&dir_to_index,&d,&ndata](){ return dir_to_index[looper.it_sample]<ndata ? !(d.evt.NTopHad==2 && fabs(d.evt.TTHadDPhi)<DPHI_CUT_OLD && d.evt.R>R_CUT_OLD ): 1; });
  sh.AddNewCut("AllFilters",       [&d](){
		 return d.evt.NGoodVtx>0 // d.evt.Flag_goodVertices - Doesn't work currently
		   && d.evt.Flag_HBHEIsoNoiseFilterResult
		   && d.evt.Flag_HBHENoiseFilterResultRun2Loose
		   && d.evt.Flag_CSCTightHaloFilter // Will need to rerun 2015 filter for Dec Jamboree
		   && d.evt.Flag_eeBadScFilter
		   // Exclude bad HCal runs
		   && d.evt.RunNumber!=256729 && d.evt.RunNumber!=256734
		   ;});
  sh.AddNewCut("ExclBadHcal", [&d](){ return d.evt.RunNumber!=256729 && d.evt.RunNumber!=256734; });
  sh.AddNewCut("BaselineCuts",      [&d,&looper,&dir_to_index,&ndata](){
		 return
		   // Trigger
		   d.evt.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50
		   // Event filters
		   && d.evt.NGoodVtx>0 // d.evt.Flag_goodVertices - Doesn't work currently
		   && d.evt.Flag_HBHEIsoNoiseFilterResult
		   && d.evt.Flag_HBHENoiseFilterResultRun2Loose
		   && d.evt.Flag_CSCTightHaloFilter // Will need to rerun 2015 filter for Dec Jamboree
		   && d.evt.Flag_eeBadScFilter
		   // At least 2 jets pre-tagged (pass pt and mass cut)
		   && d.evt.NTopHadPreTag>=2
		   // Blinding
#if NEW_TOP_DEF == 1
		   // Exclude signal region in Data (NTophad==2, |DPhi_tt|<2.7, R>0.4)
		   && (dir_to_index[looper.it_sample]<ndata ? !(d.evt.NTopHad==2 && fabs(d.evt.TTHadDPhi)<DPHI_CUT_NEW && d.evt.R>R_CUT_NEW) : 1)
#else
		   // Exclude signal region in Data (NTophad==2, |DPhi_tt|<2.8, R>0.4)
		   && (dir_to_index[looper.it_sample]<ndata ? !(d.evt.NTopHad==2 && fabs(d.evt.TTHadDPhi)<DPHI_CUT_OLD && d.evt.R>R_CUT_OLD) : 1)
#endif
		   // Exclude bad HCal runs
		   && d.evt.RunNumber!=256729 && d.evt.RunNumber!=256734
		   ;});
  
  if (debug) std::cout<<"Cut definitions ok\n";
  
  // --------------------------------------------------------------------
  //                         Postfixes
  //                  (Alternatives to cuts)
  // --------------------------------------------------------------------
  
  // Colors
  // 400 kYellow  800 kOrange
  // 416 kGreen	  820 kSpring
  // 432 kCyan	  840 kTeal
  // 600 kBlue	  860 kAzure
  // 616 kMagenta 880 kViolet
  // 632 kRed     900 kPink
  
  std::string col3_red_to_blue = "633,618,601,"; // red, purple, blue
  std::string col4_cyan_to_red = "434,601,618,633,"; // Cyan, blue, purple, red
  std::string col5_green_to_red = "418,434,601,618,633,"; // green, cyan, blue, purple, red
  std::string col5_red_to_green = "633,618,601,434,418,"; // red, , purple, blue, cyan, green
  std::string col6_rainbow_dark = "601,434,418,402,633,618,"; // blue, cyan, green, yellow, red, purple
  std::string col8 = "1,601,434,418,402,807,633,618,"; // above plus black and orange
  std::string col12 = "1,4,6,2,800,402,417,433,9,618,633,924,"; // black, blue, magenta, red, orange, darker yellow, darker green, darker cyan, blue-purple, dark purple, dark red
  std::string col12_rainbow = "402,416,433,600,617,632,802,813,833,863,883,892,"; // Go around first bright and then dark colors
  
  // Samples
  std::string Samples_PFs = samples.GetPFNames();
  std::string Samples_Lat = samples.GetLatexNames();
  std::string Samples_Cols = samples.GetColors();
  
  sh.AddNewPostfix("Directories",                 [&looper](){ return looper.it_sample; }, "[0to50]", "Sample [0to50]", "1-51");
  sh.AddNewPostfix("AllSamples",                  [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
  sh.AddNewPostfix("StackSamples",                [&looper,&dir_to_index,&ndata,&nsignal,&nttbar](){
		     return
		       /* Only Include Selected Data sample   */ dir_to_index[looper.it_sample]<ndata && dir_to_index[looper.it_sample]!=0 ? -1 :
		       /* Only Include 1st two Signal samples */ (dir_to_index[looper.it_sample]>ndata+1)&&(dir_to_index[looper.it_sample]<ndata+nsignal) ? -1 :
		       /* Only Include 1st TTbar sample       */ (dir_to_index[looper.it_sample]>ndata+nsignal)&&(dir_to_index[looper.it_sample]<ndata+nsignal+nttbar) ? -1 :
		       dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
  if (test) {
    sh.AddNewPostfix("Signals,TT,NonTT",   [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
    sh.AddNewPostfix("Signals,TT",   [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
    sh.AddNewPostfix("Signals,Background", [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
    sh.AddNewPostfix("Background", [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
    sh.AddNewPostfix("Data,MC", [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
    sh.AddNewPostfix("Background,Signal", [&looper,&dir_to_index](){ return dir_to_index[looper.it_sample]; }, Samples_PFs, Samples_Lat, Samples_Cols);
  } else {
    sh.AddNewPostfix("Signals,TT,NonTT",
		     [&looper,&dir_to_index,&ndata,&nsignal,&nttbar]() {
		       return
			 /* Exclude Data */ dir_to_index[looper.it_sample]<ndata ? -1 :
			 dir_to_index[looper.it_sample]<ndata+nsignal+nttbar ? dir_to_index[looper.it_sample] : nsignal+nttbar; 
		     },
		     std::string(Samples_PFs).replace(Samples_PFs.find("WJets"),5,"NonTT"),
		     std::string(Samples_Lat).replace(Samples_Lat.find("W+Jets"),6,"Non-ttbar"),
		     Samples_Cols);
    sh.AddNewPostfix("Signals,TT",
		     [&looper,&dir_to_index,&ndata,&nsignal,&nttbar]() {
		       return
			 /* Exclude Data */ dir_to_index[looper.it_sample]<ndata ? -1 :
			 dir_to_index[looper.it_sample]<ndata+nsignal+nttbar ? dir_to_index[looper.it_sample] : -1; 
		     }, 
		     std::string(Samples_PFs).replace(Samples_PFs.find("WJets"),5,"NonTT"),
		     std::string(Samples_Lat).replace(Samples_Lat.find("W+Jets"),6,"Non-ttbar"),
		     Samples_Cols);
    sh.AddNewPostfix("Signals,Background", 
		     [&looper,&dir_to_index,&ndata,&nsignal,&nttbar]() {
		       return
			 /* Exclude Data                  */ dir_to_index[looper.it_sample]<ndata ? -1 :
			 /* Only Include 1st TTbar sample */ (dir_to_index[looper.it_sample]>ndata+nsignal)&&(dir_to_index[looper.it_sample]<ndata+nsignal+nttbar) ? -1 :
			 dir_to_index[looper.it_sample]<ndata+nsignal ? dir_to_index[looper.it_sample] : ndata+nsignal;
		     },
		     std::string(Samples_PFs).replace(Samples_PFs.find("TTJetsHT600"),11,"Background"),
		     std::string(Samples_Lat).replace(Samples_Lat.find("t#bar{t}+Jets (HT600)"),21,"Background"),
		     Samples_Cols);
    sh.AddNewPostfix("Background",
		     [&looper,&dir_to_index,&ndata,&nsignal,&nttbar]() {
		       return
			 /* Exclude Data and Signal       */ dir_to_index[looper.it_sample]<ndata+nsignal ? -1 :
			 /* Only Include 1st TTbar sample */ (dir_to_index[looper.it_sample]<ndata+nsignal+nttbar)&&(dir_to_index[looper.it_sample]>ndata+nsignal) ? -1 : 0;
		     }, "Background", "Background", "1");
    sh.AddNewPostfix("Data,MC",
		     [&looper,&dir_to_index,&ndata,&nsignal,&nttbar]() {
		       return
			 /* Only Include Selected Data sample */ dir_to_index[looper.it_sample]<ndata && dir_to_index[looper.it_sample]!=0 ? -1 :
			 /* Exclude Signal                    */ (dir_to_index[looper.it_sample]<ndata+nsignal)&&(dir_to_index[looper.it_sample]>=ndata) ? -1 :
			 /* Only Include 1st TTbar sample     */ (dir_to_index[looper.it_sample]<ndata+nsignal+nttbar)&&(dir_to_index[looper.it_sample]>ndata+nsignal) ? -1 :
			 dir_to_index[looper.it_sample]==0 ? 0 : 1;
		     },
		     "Data;MC", "Data;MC", "1,633");
    sh.AddNewPostfix("Background,Signal",
		     [&looper,&dir_to_index,&ndata,&nsignal,&nttbar]() {
		       return
			 /* Exclude Data                   */ dir_to_index[looper.it_sample]<ndata ? -1 :
			 /* Only Include 1st Signal sample */ (dir_to_index[looper.it_sample]<ndata+nsignal)&&(dir_to_index[looper.it_sample]>ndata) ? -1 :
			 /* Only Include 1st TTbar  sample */ (dir_to_index[looper.it_sample]<ndata+nsignal+nttbar)&&(dir_to_index[looper.it_sample]>ndata+nsignal) ? -1 :
			 dir_to_index[looper.it_sample]==ndata ? 1 : 0;
		     },
		     "Background;T5ttttDeg", "MC - Background;MC - T5ttttDeg Signal", "633,601");
    sh.AddNewPostfix("ABCD,Signals",           [&looper,&dir_to_index,&ndata,&nsignal,&nttbar,&d](){
		       return
			 /* Exclude Data                   */ dir_to_index[looper.it_sample]<ndata ? -1 :
			 /* Only Include 1st TTbar  sample */ (dir_to_index[looper.it_sample]<ndata+nsignal+nttbar)&&(dir_to_index[looper.it_sample]>ndata+nsignal) ? -1 :
			 /* Signals                        */ dir_to_index[looper.it_sample]<ndata+nsignal ? dir_to_index[looper.it_sample]-ndata+4 :
			 /* Background - ABCD              */
			 d.evt.R<R_CUT_LOW || d.evt.NTopHad>2 ? -1 : // Exclude low R region
			 (d.evt.R>R_CUT_NEW) + (d.evt.NTopHad==2)*2;
		     },
		     "A;B;C;D;T5ttttDeg_mGo1000_4bodydec;T1tttt_mGo1500_mChi100", 
		     "A: 0.2<R<0.4, N_{top-tag} = 0,1;B: R > 0.4, N_{top-tag} = 0,1;C: 0.2<R<0.4, N_{top-tag} = 2;D (SR): R > 0.4, N_{top-tag} = 2;"
		     "T5ttttDeg (M_{#tilde{g}}=1TeV);T1tttt (M_{#tilde{g}}=1.5TeV)",
		     "858,602,628,634,419,402");
  }
  
  // Jets
  //sh.AddNewPostfix("AK4JetsPtOrdered", [&d](){ return d.jetsAK4.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("JetsPtOrdered",    [&d](){ return d.jetsAK8.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("NSubJet",          [&d](){ return (size_t)d.jetsAK8.nSubJets[d.jetsAK8.it]; }, "NSubJet[0to4]", "N_{subjet}=[0to4]", "1-5");
#if NEW_TOP_DEF == 1
  sh.AddNewPostfix("JetPtCut",      [&d](){ return (size_t)(d.jetsAK8.Pt[d.jetsAK8.it] > PT_CUT_NEW); }, "PtBelow400;PtAbove400", "p_{T} < 400;p_{T} > 400", "2,3");
  sh.AddNewPostfix("JetMassCut",    [&d](){ return (size_t)(d.jetsAK8.softDropMass[d.jetsAK8.it] > SD_MASS_CUT_LOW && d.jetsAK8.softDropMass[d.jetsAK8.it] < SD_MASS_CUT_HIGH); }, "OutMassWindow;InMassWindow", "M_{Soft-drop} < 110 or M_{Soft-drop} > 210;110 < M_{Soft-drop} < 210", "2,3");
  sh.AddNewPostfix("JetTau32Cut",   [&d](){ return (size_t)(d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it] < TAU32_CUT_NEW); }, "Tau32Above0p75;Tau32Below0p75", "#tau_{3}/#tau_{2} > 0.75;#tau_{3}/#tau_{2} < 0.75", "2,3");
#else
  sh.AddNewPostfix("JetPtCut",         [&d](){ return (size_t)(d.jetsAK8.Pt[d.jetsAK8.it] > PT_CUT_OLD); }, "PtBelow400;PtAbove400", "p_{T} < 400;p_{T} > 400", "2,3");
  sh.AddNewPostfix("JetMassCut",       [&d](){ return (size_t)(d.jetsAK8.prunedMass[d.jetsAK8.it] > PRUNED_MASS_CUT_LOW); }, "MassBelow140;MassAbove140", "M_{pruned} < 140;M_{pruned} > 140", "2,3");
  sh.AddNewPostfix("JetTau32Cut",      [&d](){ return (size_t)(d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it] < TAU32_CUT_OLD); }, "Tau32Above0p75;Tau32Below0p75", "#tau_{3}/#tau_{2} > 0.75;#tau_{3}/#tau_{2} < 0.75", "2,3");
#endif
  
  sh.AddNewPostfix("JetGenTruth",          [&d](){ return (size_t)(d.jetsAK8.HasNearGenTop[d.jetsAK8.it]==0 ? 0 : d.jetsAK8.NearGenTopIsHadronic[d.jetsAK8.it]!=1 ? 1 : d.jetsAK8.DRNearGenWFromTop[d.jetsAK8.it]<0.6&&d.jetsAK8.DRNearGenBFromTop[d.jetsAK8.it]<0.6 ? 3 : 2); }, "NotTop;SemiLepTop;NonMergedHadTop;MergedHadTop", "Non-top jet;Semi-leptonic top;Non-Merged hadronic top;Merged hadronic top", "1,601,633,418");
  sh.AddNewPostfix("JetMatchedGenTopType", [&d](){ return (size_t)(d.jetsAK8.NearGenTopIsHadronic[d.jetsAK8.it]!=-9999 ? d.jetsAK8.NearGenTopIsHadronic[d.jetsAK8.it] : -1); }, "MatchedGenTopLep;MatchedGenTopHad", "Semi-leptonic top;Hadronic top", "4,2");
  sh.AddNewPostfix("TopSizeCut",             [&d](){ return (size_t)(d.jetsAK8.HasNearGenTop[d.jetsAK8.it]==1 ? d.jetsAK8.DRNearGenWFromTop[d.jetsAK8.it]<0.6&&d.jetsAK8.DRNearGenBFromTop[d.jetsAK8.it]<0.6 : -1); }, "TopSizeAbove0p6;TopSizeBelow0p6", "Non-merged top;Merged top", "2,4");
  sh.AddNewPostfix("JetPassToptag",    [&d](){ return (size_t)(d.jetsAK8.Pt[d.jetsAK8.it]>PT_CUT_NEW && d.jetsAK8.softDropMass[d.jetsAK8.it] > SD_MASS_CUT_LOW && d.jetsAK8.softDropMass[d.jetsAK8.it] < SD_MASS_CUT_HIGH && d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it] < TAU32_CUT_NEW); }, "FailTopTag;PassTopTag", "Non top-tagged AK8 jet;Top-tagged AK8 jet", "2,3");
  
  // Event
  sh.AddNewPostfix("RBins",          [&d](){ return (size_t)((d.evt.R>0.1)+(d.evt.R>0.2)+(d.evt.R>0.4)); }, "R0to0p1;R0p1to0p2;R0p2to0p4;R0p4", "0.0<R<0.1;0.1<R<0.2;0.2<R<0.4;R>0.4", "1,4,418,633");
  //sh.AddNewPostfix("RBins0p1",       [&d](){ return (size_t)(d.evt.R<0.4 ? d.evt.R/0.1 : 4); }, "R0to0p1;R0p1to0p2;R0p2to0p3;R0p3to0p4;R0p4", "0.0<R<0.1;0.1<R<0.2;0.2<R<0.3;0.3<R<0.4;R>0.4", "1,4,418,401,807,633,618");
  sh.AddNewPostfix("TTHadMRBins",    [&d](){ return (size_t)(d.evt.TTHadMR > 5000 ? -1 : d.evt.TTHadMR/500); }, "MR[250to4750++500]", "MR_{tt} [250to4750++500]#pm250", "1,4,418,401,807,633,618,1,4,418,401,807,633,618");
#if NEW_TOP_DEF == 1
  sh.AddNewPostfix("RBands",         [&d](){ return d.evt.R < R_CUT_LOW ? -1 : d.evt.R > R_CUT_NEW; }, "RBelow0p4;RAbove0p4", "0.2 < R < 0.4;R > 0.4", "4,2");
  sh.AddNewPostfix("DPhiBands",      [&d](){ return fabs(d.evt.TTHadDPhi) > DPHI_CUT_NEW; }, "DPhiBelow2p7;DPhiAbove2p7", "#Delta#phi_{t#bar{t}} < 2.7;#Delta#phi_{t#bar{t}} > 2.7", "2,4");
  sh.AddNewPostfix("NTopBands",      [&d](){ return d.evt.NTopHad>2 ? -1 : (d.evt.NTopHad>1); }, "0To1HadTop;2HadTop", "N_{top-tag} = 0,1;N_{top-tag} = 2", "4,2");
  sh.AddNewPostfix("NTopHad",     [&d](){ return d.evt.NTopHad; },   "0HadTopTag;1HadTopTag;2HadTopTag;3HadTopTag;4HadTopTag", "N_{top-tag,hadronic}=0;N_{top-tag,hadronic}=1;N_{top-tag,hadronic}=2;N_{top-tag,hadronic}=3;N_{top-tag,hadronic}=4", col5_green_to_red);
  sh.AddNewPostfix("NTopHadPreTag",     [&d](){ return d.evt.NTopHadPreTag; },  "0HadTopPreTag;1HadTopPreTag;2HadTopPreTag;3HadTopPreTag;4HadTopPreTag", "N_{top-like,hadronic}=0;N_{top-like,hadronic}=1;N_{top-like,hadronic}=2;N_{top-like,hadronic}=3;N_{top-like,hadronic}=4", col5_green_to_red);
  sh.AddNewPostfix("ABCD",           [&d](){ return 
		       d.evt.R<R_CUT_LOW || d.evt.NTopHad>2 ? -1 : // Exclude low R region
		       (d.evt.R>R_CUT_NEW) + (d.evt.NTopHad==2)*2;
		   }, "A;B;C;D", "A: 0.2<R<0.4, N_{top-tag} = 0,1;B: R > 0.4, N_{top-tag} = 0,1;C: 0.2<R<0.4, N_{top-tag} = 2;D (SR): R > 0.4, N_{top-tag} = 2", "858,602,628,634");
# else
  sh.AddNewPostfix("RBands",         [&d](){ return d.evt.R < R_CUT_LOW ? -1 : d.evt.R > R_CUT_OLD; }, "RBelow0p4;RAbove0p4", "0.2 < R < 0.4;R > 0.4", "4,2");
  sh.AddNewPostfix("DPhiBands",      [&d](){ return fabs(d.evt.TTHadDPhi) > DPHI_CUT_OLD; },  "DPhiBelow2p8;DPhiAbove2p8", "#Delta#phi_{t#bar{t}} < 2.8;#Delta#phi_{t#bar{t}} > 2.8", "2,4");
  sh.AddNewPostfix("NTopBands",      [&d](){ return d.evt.NTopHad>2 ? -1 : (d.evt.NTopHad>1); }, "0To1HadTop;2HadTop", "N_{top-tag} = 0,1;N_{top-tag} = 2", "4,2");
  sh.AddNewPostfix("NTopHad",     [&d](){ return d.evt.NTopHad; },      "0HadTopTag;1HadTopTag;2HadTopTag;3HadTopTag;4HadTopTag", "N_{top-tag,hadronic}=0;N_{top-tag,hadronic}=1;N_{top-tag,hadronic}=2;N_{top-tag,hadronic}=3;N_{top-tag,hadronic}=4", col5_green_to_red);
  sh.AddNewPostfix("NTopHadPreTag",     [&d](){ return d.evt.NTopHadPreTag; },     "0HadTopPreTag;1HadTopPreTag;2HadTopPreTag;3HadTopPreTag;4HadTopPreTag", "N_{top-like,hadronic}=0;N_{top-like,hadronic}=1;N_{top-like,hadronic}=2;N_{top-like,hadronic}=3;N_{top-like,hadronic}=4", col5_green_to_red);
  sh.AddNewPostfix("ABCD",           [&d](){ return 
		       d.evt.R<R_CUT_LOW || d.evt.NTopHad>2 ? -1 : // Exclude low R region
		       (d.evt.R>R_CUT_OLD) + (d.evt.NTopHad==2)*2;
		   }, "A;B;C;D", "A: 0.2<R<0.4, <2 Top tag;B: R > 0.4, <2 Top tag;C: 0.2<R<0.4, 2 Top tag;D (Signal): R > 0.4, 2 Top tag", "602,413,618,629");
#endif
  sh.AddNewPostfix("CutHtAll",        [&d](){ return d.evt.HtAll > 1200; },  "HtAllBelow1200;HtAllAbove1200", "H_{T,all} < 1200;H_{T,all} > 1200", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1450",       [&d](){ return d.evt.HtAll > 1450; },  "HtAllBelow1450;HtAllAbove1450", "H_{T,all} < 1450;H_{T,all} > 1450", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1500",       [&d](){ return d.evt.HtAll > 1500; },  "HtAllBelow1500;HtAllAbove1500", "H_{T,all} < 1500;H_{T,all} > 1500", "4,2");
  sh.AddNewPostfix("NGenLepFromTop",   [&d](){ return (size_t)d.evt.NGenLepFromTop; }, "FullHad;[1to4]LepTop", "[0to4]l (e/#mu, from top)", "1-5");
  
  // Trigger
  sh.AddNewPostfix("PFHT475",         [&d](){ return (size_t)d.evt.HLT_PFHT475; }, "NoPassHLT_PFHT475;PassHLT_PFHT475", "Do not pass HLT_PFHT475;Pass HLT_PFHT475", "633;418");
  //sh.AddNewPostfix("LowPFHT",         [&d](){ return size_t(d.evt.HLT_PFHT350||d.evt.HLT_PFHT400||d.evt.HLT_PFHT475); }, "NoPassLowPFHT;PassLowPFHT", "Do not pass HLT_PFHT(350,400,475);Pass HLT_PFHT(350,400,475)", "633;418");
  
  // Gen Particles
  sh.AddNewPostfix("GenTopType",      [&d](){ return (size_t)(d.evt.GenTopType[d.gen.it]!=-9999 ? d.evt.GenTopType[d.gen.it] : -1); }, "GenTopHad;GenTopLep", "Hadronic top;Semi-leptonic top", "2,4");
  
  if (debug) std::cout<<"Postfix definitions ok\n";
  
  // --------------------------------------------------------------------
  //                         Fill Parameters
  // --------------------------------------------------------------------
  
  // Define histo parameters and filling variable
  // X/Y/Z - axis parameters:

  // Samples
  sh.AddNewFillParam("Sample",            { .nbin= 50,  .bins={   0,     50}, .fill=[&looper](){ return looper.it_sample;     }, .axis_title="iSample"});
  
  // Muons
  sh.AddNewFillParam("MuEnergy",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.mu.E[d.mu.it];              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuPt",              { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.mu.Pt[d.mu.it];             }, .axis_title="Muon p_{T} (GeV/c)"});
  sh.AddNewFillParam("MuDRJet",           { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.MuDRJet[d.mu.it];    }, .axis_title="#DeltaR (#mu, jet)"});
  sh.AddNewFillParam("MuRelPtJet",        { .nbin=  50, .bins={   0,    500}, .fill=[&d](){ return d.evt.MuRelPtJet[d.mu.it]; }, .axis_title="p_{T}^{rel} (#mu, jet) (GeV/c)"});
  sh.AddNewFillParam("MuJetCombMass",     { .nbin=  80, .bins={   0,   2000}, .fill=[&d](){ return d.evt.MuJetCombMass[d.mu.it]; }, .axis_title="Mass_{#mu+jet comb.} (GeV/c^{2})"});
  
  // Electrons
  sh.AddNewFillParam("EleEnergy",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.ele.E[d.ele.it];              }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.ele.Pt[d.ele.it];             }, .axis_title="Electron p_{T} (GeV/c)"});
  sh.AddNewFillParam("EleDRJet",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.EleDRJet[d.ele.it];   }, .axis_title="#DeltaR (e, jet)"});
  sh.AddNewFillParam("EleRelPtJet",       { .nbin=  50, .bins={   0,    500}, .fill=[&d](){ return d.evt.EleRelPtJet[d.ele.it];}, .axis_title="p_{T}^{rel} (e, jet) (GeV/c)"});
  sh.AddNewFillParam("EleJetCombMass",    { .nbin=  80, .bins={   0,   2000}, .fill=[&d](){ return d.evt.EleJetCombMass[d.ele.it]; }, .axis_title="Mass_{e+jet comb.} (GeV/c^{2})"});
  
  // MET
  sh.AddNewFillParam("MetPt",             { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.met.Pt;                   }, .axis_title="MET p_{T} (GeV/c)"});
  
  // AK4 Jets
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.E[d.jetsAK8.it];          }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.Pt[d.jetsAK8.it];         }, .axis_title="AK4-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("AK4JetMass",        { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.Mass[d.jetsAK8.it];       }, .axis_title="AK4-jet Mass (GeV/c^{2})"});
  
  // Jets (AK8)
  sh.AddNewFillParam("JetEnergy",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.E[d.jetsAK8.it];                    }, .axis_title="AK8-jet Energy (GeV)"});
  sh.AddNewFillParam("JetPt",              { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.Pt[d.jetsAK8.it];                   }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtBins",          { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return d.jetsAK8.Pt[d.jetsAK8.it]; }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtFewBins",       { .nbin=   5, .bins={0, 300, 400, 600, 1000, 2000}, .fill=[&d](){ return d.jetsAK8.Pt[d.jetsAK8.it]; }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetPtOneBin",        { .nbin=   1, .bins={ 400,    5000}, .fill=[&d](){ return d.jetsAK8.Pt[d.jetsAK8.it]; }, .axis_title="AK8-jet p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetEta",             { .nbin=  40, .bins={  -4,      4},                  .fill=[&d](){ return d.jetsAK8.Eta[d.jetsAK8.it]; }, .axis_title="AK8-jet #eta"});
  sh.AddNewFillParam("JetPhi",             { .nbin=  16, .bins={-3.1416, 3.1416}, .fill=[&d](){ return d.jetsAK8.Phi[d.jetsAK8.it]; }, .axis_title="AK8-jet #phi"});
  sh.AddNewFillParam("JetNeutralHadronMultiplicity", { .nbin=  20, .bins={0,  20}, .fill=[&d](){ return d.jetsAK8.neutralHadronMultiplicity[d.jetsAK8.it]; }, .axis_title="AK8-jet Neutral Hadron Multiplicity"});
  sh.AddNewFillParam("JetChargedHadronMultiplicity", { .nbin=  50, .bins={0, 100}, .fill=[&d](){ return d.jetsAK8.ChargedHadronMultiplicity[d.jetsAK8.it]; }, .axis_title="AK8-jet Charged Hadron Multiplicity"});
  
  sh.AddNewFillParam("JetMass",            { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.Mass[d.jetsAK8.it];                 }, .axis_title="AK8-jet Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetPrunedMass",      { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.prunedMass[d.jetsAK8.it];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetPrunedMassCoarse",{ .nbin= 200, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.prunedMass[d.jetsAK8.it];           }, .axis_title="AK8-jet Pruned Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetFilteredMass",    { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.filteredMass[d.jetsAK8.it];         }, .axis_title="AK8-jet Filtered Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetTrimmedMass",     { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.trimmedMass[d.jetsAK8.it];          }, .axis_title="AK8-jet Trimmed Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetSoftDropMass",    { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.softDropMass[d.jetsAK8.it];         }, .axis_title="AK8-jet Soft-drop Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetTopMass",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.topMass[d.jetsAK8.it];              }, .axis_title="AK8-jet Top Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetMinMass",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.minmass[d.jetsAK8.it];              }, .axis_title="AK8-jet Min. Subjet-pair Mass (GeV/c^{2})"});
  sh.AddNewFillParam("JetNSubJets",        { .nbin=  11, .bins={-0.5,   10.5}, .fill=[&d](){ return d.jetsAK8.nSubJets[d.jetsAK8.it];             }, .axis_title="AK8-jet N_{subjet}"});
  sh.AddNewFillParam("JetTau1",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau1[d.jetsAK8.it];                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("JetTau2",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau2[d.jetsAK8.it];                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("JetTau3",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau3[d.jetsAK8.it];                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("JetTau21",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau2[d.jetsAK8.it]/d.jetsAK8.tau1[d.jetsAK8.it];   }, .axis_title="#tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("JetTau31",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau1[d.jetsAK8.it];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetTau32",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8.tau3[d.jetsAK8.it]/d.jetsAK8.tau2[d.jetsAK8.it];   }, .axis_title="#tau_{3}/#tau_{2}"});
  sh.AddNewFillParam("JetDRLep",           { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.DRJetLep[d.jetsAK8.it];         }, .axis_title="#DeltaR (lepton, jet)"});
  sh.AddNewFillParam("JetRelPtLep",        { .nbin= 100, .bins={   0,    500}, .fill=[&d](){ return d.evt.RelPtJetLep[d.jetsAK8.it];      }, .axis_title="p_{T}^{rel} (lepton, jet) [GeV/c]"});
  
  sh.AddNewFillParam("JetMatchedGenTopPt",        { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.PtNearGenTop[d.jetsAK8.it];         }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMatchedGenTopPtCoarse",  { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8.PtNearGenTop[d.jetsAK8.it];         }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMatchedGenTopPtBins",    { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return d.jetsAK8.PtNearGenTop[d.jetsAK8.it]; }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("JetMatchedGenTopJetDR",     { .nbin=  50, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenTop[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. top, jet)"});
  sh.AddNewFillParam("JetMatchedGenTopJetDRFine", { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenTop[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. top, jet)"});
  sh.AddNewFillParam("GenBJetDR",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenBFromTop[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. b, jet)"});
  sh.AddNewFillParam("GenBJetDRFine",      { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenBFromTop[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. b, jet)"});
  sh.AddNewFillParam("GenWJetDR",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenWFromTop[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, jet)"});
  sh.AddNewFillParam("GenWJetDRFine",      { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenWFromTop[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, jet)"});
  sh.AddNewFillParam("GenLepJetDR",        { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenLepFromSLTop[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. lep, jet)"});
  sh.AddNewFillParam("GenLepJetDRFine",    { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.jetsAK8.DRNearGenLepFromSLTop[d.jetsAK8.it];      }, .axis_title="#DeltaR (Gen. lep, jet)"});
  //sh.AddNewFillParam("GenWGenBDR",         { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenWGenBDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, Gen. b)"});
  //sh.AddNewFillParam("GenWGenBDRFine",     { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenWGenBDR[d.jetsAK8.it];	  }, .axis_title="#DeltaR (Gen. W, Gen. b)"});
  //sh.AddNewFillParam("GenLepGenBDR",       { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenLepGenBDR[d.jetsAK8.it];     }, .axis_title="#DeltaR (Gen. lep, Gen. b)"});
  //sh.AddNewFillParam("GenLepGenBDRFine",   { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return d.evt.GenLepGenBDR[d.jetsAK8.it];     }, .axis_title="#DeltaR (Gen. lep, Gen. b)"});
  
  sh.AddNewFillParam("MaxSubJetCSV",       { .nbin=  9, .bins={ 0, 0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 }, .fill=[&d](){ return d.evt.maxSubjetCSV[d.jetsAK8.it]; }, .axis_title="Maximum Subjet CSV"});
  
  // Special Y/Z axis parameters:
  // Define how to calculate them in SmartHistos!
  // SmartHistos::init_() for extra histo name and axis title
  // SmartHistos::calc_spec_[1,2]d() for the method
  sh.AddNewFillParam("MergedTopFraction",  { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.JetMatchedGenTopIsMerged[d.jetsAK8.it]; }, .axis_title="Fraction of Merged Tops" });   
  sh.AddNewFillParam("TopTagEfficiency",   { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.JetIsHadTopTagged[d.jetsAK8.it]; }, .axis_title="Top-tagging Efficiency" });
  sh.AddNewFillParam("MisTagRate",         { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return (size_t)(!d.evt.JetHasMatchedGenTop[d.jetsAK8.it]); }, .axis_title="Mis-tag Rate" });
  
  // Event variables
  sh.AddNewFillParam("NJetSelected",        { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopHad+d.evt.NTopLep;     }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NJetHadronic",        { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopHad;                   }, .axis_title="N_{Hadronic AK8-jet}"});
  sh.AddNewFillParam("NJetLeptonic",        { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{Leptonic AK8-jet}"});
  sh.AddNewFillParam("NLep",                { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NLep;                      }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NLepTight",           { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neletight+d.evt.nmu;       }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NMu",                 { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nmu;                       }, .axis_title="N_{muon}"});
  sh.AddNewFillParam("NEle",                { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nele;                      }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NEleTight",           { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neletight;                 }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NLepVeto",            { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neleveto+d.evt.nmuveto;    }, .axis_title="N_{lepton}"});
  sh.AddNewFillParam("NMuVeto",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nmuveto;                   }, .axis_title="N_{muon}"});
  sh.AddNewFillParam("NEleVeto",            { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neleveto;                  }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NTopHadPreTag",         { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NTopHadPreTag;               }, .axis_title="N_{jet} (p_{T}>400, M>100)"});
  sh.AddNewFillParam("NTopLep",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{top, leptonic}"});
  sh.AddNewFillParam("NTopHad",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NTopHad;                   }, .axis_title="N_{top, hadronic}"});
  sh.AddNewFillParam("TTHadMR",             { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMRCoarse",       { .nbin=  20, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMTR",            { .nbin=  50, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMTR;                  }, .axis_title="M_{T,t#bar{t}}^{R} (GeV/c)"});
  sh.AddNewFillParam("DPhiFine",            { .nbin=  64, .bins={   0,     3.2}, .fill=[&d](){ return fabs(d.evt.TTHadDPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DPhiBins",            { .nbin=   9, .bins={ 0, 0.5, 1.0, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.15 }, .fill=[&d](){ return fabs(d.evt.TTHadDPhi); }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DPhi",                { .nbin=  16, .bins={   0,     3.2}, .fill=[&d](){ return fabs(d.evt.TTHadDPhi);           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("TTHadSumPt",          { .nbin=   7, .bins={ 600, 800, 850, 900, 1000, 1200, 2000, 10000}, .fill=[&d](){ return d.evt.TTHadSumPt; }, .axis_title="p_{T, 1st AK8 jet} + p_{T, 2nd AK8 jet}"});
  sh.AddNewFillParam("TTHadSumPtOneBin",    { .nbin=   1, .bins={ 800,   10000}, .fill=[&d](){ return d.evt.TTHadSumPt;                }, .axis_title="p_{T, 1st AK8 jet} + p_{T, 2nd AK8 jet}"});
  sh.AddNewFillParam("TTHadDEta",           { .nbin=  50, .bins={   0,       5}, .fill=[&d](){ return d.evt.TTHadDEta;                 }, .axis_title="#Delta#eta_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadDR",             { .nbin=  60, .bins={   0,       6}, .fill=[&d](){ return d.evt.TTHadDR;                   }, .axis_title="#DeltaR_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadPz",             { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadPz;                   }, .axis_title="P_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadDPz",            { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadDPz;                  }, .axis_title="#DeltaP_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadHz",             { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadHz;                   }, .axis_title="H_{Z,t#bar{t}} (GeV/c)"});
  sh.AddNewFillParam("TTHadMass",           { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMass;                 }, .axis_title="M_{t#bar{t}} (GeV/c^{2})"});
  sh.AddNewFillParam("TTHadR",              { .nbin=  24, .bins={   0,     1.2}, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadRFine",          { .nbin= 120, .bins={   0,     1.2}, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadR2",             { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.TTHadR2;                   }, .axis_title="R_{t#bar{t}}^{2}"});
  sh.AddNewFillParam("R",                   { .nbin=  24, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("RFine",               { .nbin= 120, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R"});
  sh.AddNewFillParam("RBins",               { .nbin=  14, .bins={ 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2 }, .fill=[&d](){ return d.evt.R; }, .axis_title="R"});
  sh.AddNewFillParam("R2",                  { .nbin=  32, .bins={   0,     1.6}, .fill=[&d](){ return d.evt.R2;                        }, .axis_title="R^{2}"});
  sh.AddNewFillParam("MR",                  { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.MR;                        }, .axis_title="M_{R} (GeV/c)"});
  sh.AddNewFillParam("MTR",                 { .nbin= 100, .bins={   0,    2000}, .fill=[&d](){ return d.evt.MTR;                       }, .axis_title="M_{T}^{R} (GeV/c)"});
  sh.AddNewFillParam("HtTopFraction",       { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.HtTopFr;                   }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("HtExFraction",        { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.HtExFr;                    }, .axis_title="H_{T,extra}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV/c)"});
  sh.AddNewFillParam("Ht",                  { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.Ht;                        }, .axis_title="H_{T} (GeV/c)"});
  sh.AddNewFillParam("HtBins",              { .nbin=  19, .bins={ 0, 200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1500, 2000, 3000, 10000}, .fill=[&d](){ return d.evt.Ht; }, .axis_title="H_{T,AK8} (GeV/c)"});
  sh.AddNewFillParam("HtOneBin",            { .nbin=   1, .bins={800, 10000}, .fill=[&d](){ return d.evt.Ht;                        }, .axis_title="H_{T,AK8} (GeV/c)"});
  sh.AddNewFillParam("HtAllCoarse",         { .nbin=  20, .bins={   0,    6000}, .fill=[&d](){ return d.evt.HtAll;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtAll",               { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HtAll;                     }, .axis_title="H_{T}+H_{T,leptonic}+#slash{p}_{T} (GeV/c)"});
  sh.AddNewFillParam("HtTop",               { .nbin=  25, .bins={   0,    5000}, .fill=[&d](){ return d.evt.HtTop;                     }, .axis_title="H_{T,tops} (GeV/c)"});
  sh.AddNewFillParam("HtEx",                { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HtEx;                      }, .axis_title="H_{T,extra} (GeV/c)"});
  sh.AddNewFillParam("NVertices",           { .nbin= MAX_NVTX+1, .bins={-0.5,   MAX_NVTX+0.5}, .fill=[&d](){ return d.evt.npv;                       }, .axis_title="N_{Vertices}"});
  sh.AddNewFillParam("NVerticesReweighted", { .nbin= MAX_NVTX+1, .bins={-0.5,   MAX_NVTX+0.5}, .fill=[&d](){ return d.evt.npv;                       }, .axis_title="Reweighted N_{Vertices}"});
  
  sh.AddNewFillParam("FlaggoodVertices",                                { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_goodVertices;                       }, .axis_title="goodVertices Filter"});
  sh.AddNewFillParam("FlaggoodVerticesFix",                             { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.NGoodVtx>0;                              }, .axis_title="N_{PV,ndof>4,|z|<24,|#rho|<2} > 0"});
  sh.AddNewFillParam("FlageeBadScFilter",                               { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_eeBadScFilter;                      }, .axis_title="eeBadScFilter Filter"});
  sh.AddNewFillParam("FlagCSCTightHaloFilter",                          { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_CSCTightHaloFilter;                 }, .axis_title="CSCTightHaloFilter Filter"});
  sh.AddNewFillParam("FlagHBHEIsoNoiseFilterResult",                    { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_HBHEIsoNoiseFilterResult;           }, .axis_title="HBHEIsoNoiseFilterResult Filter"});
  sh.AddNewFillParam("FlagHBHENoiseFilterResult",                       { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_HBHENoiseFilterResult;              }, .axis_title="HBHENoiseFilterResult Filter"});
  sh.AddNewFillParam("FlagHBHENoiseFilterResultRun1",                   { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_HBHENoiseFilterResultRun1;          }, .axis_title="HBHENoiseFilterResultRun1 Filter"});
  sh.AddNewFillParam("FlagHBHENoiseFilterResultRun2Loose",              { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_HBHENoiseFilterResultRun2Loose;     }, .axis_title="HBHENoiseFilterResultRun2Loose Filter"});
  sh.AddNewFillParam("FlagHBHENoiseFilterResultRun2Tight",              { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_HBHENoiseFilterResultRun2Tight;     }, .axis_title="HBHENoiseFilterResultRun2Tight Filter"});
  sh.AddNewFillParam("FlagEcalDeadCellTriggerPrimitiveFilter",          { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.evt.Flag_EcalDeadCellTriggerPrimitiveFilter; }, .axis_title="EcalDeadCellTriggerPrimitiveFilter Filter"});
  
  // Trigger Efficiencies
  sh.AddNewFillParam("HLTEfficiencyAK8PFJet360TrimMass30",          { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.HLT_AK8PFJet360_TrimMass30;          }, .axis_title="#epsilon_{HLT_AK8PFJet360_TrimMass30}" });
  sh.AddNewFillParam("HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50", { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50; }, .axis_title="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}" });
  sh.AddNewFillParam("HLTEfficiencyPFHT750_4Jet",                   { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.HLT_PFHT750_4Jet;                    }, .axis_title="#epsilon_{HLT_PFHT750_4Jet}" });
  sh.AddNewFillParam("HLTEfficiencyPFHT350",                        { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.HLT_PFHT350;                         }, .axis_title="#epsilon_{HLT_PFHT350}" });
  sh.AddNewFillParam("HLTEfficiencyPFHT800",                        { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.HLT_PFHT800;                         }, .axis_title="#epsilon_{HLT_PFHT800}" });
  sh.AddNewFillParam("HLTEfficiencyPFHT900",                        { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.HLT_PFHT900;                         }, .axis_title="#epsilon_{HLT_PFHT900}" });
  
  // Gen Particles
  sh.AddNewFillParam("GenTopPtBins",      { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return d.evt.IsGenTop[d.gen.it] ? d.gen.Pt[d.gen.it] : -9999; }, .axis_title="Gen. top p_{T} (GeV/c)"});
  sh.AddNewFillParam("GenTopPtFewBins",   { .nbin=  4, .bins={0, 300, 400, 800, 2000}, .fill=[&d](){ return d.evt.IsGenTop[d.gen.it] ? d.gen.Pt[d.gen.it] : -9999; }, .axis_title="Gen. top p_{T} (GeV/c)"});
  
  // Special Y/Z axis parameters:
  sh.AddNewFillParam("JetFindingEfficiency", { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.GenTopHasMatchedJet[d.gen.it]; }, .axis_title="Jet finding Efficiency" });
  sh.AddNewFillParam("TopFindingEfficiency", { .nbin=  2, .bins={ -0.5, .high= 1.5}, .fill=[&d](){ return d.evt.GenTopHasMatchedTopTagJet[d.gen.it]; }, .axis_title="Top finding Efficiency" });
  
  if (debug) std::cout<<"Fill parameters ok\n";
  
  // --------------------------------------------------------------------------
  //                           Histogram Definitions
  
  // Weights and reweighting
  // Set Histogram weight (empty = 1)
  sh.SetHistoWeights({[&ndata,&looper,&dir_to_index,&d,&sample_scale_factors](){ return 
			  dir_to_index[looper.it_sample]<ndata ? 1 : 
			  /* Cross section weights: */ ( Data_IntLumi_invfb * d.evt.weight )
			  /* K-factors    */           * sample_scale_factors[looper.it_sample]
			  ;}});
  
  // Filters
  sh.AddHistos("evt",   { .fill="FlaggoodVertices",                       .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlaggoodVerticesFix",                    .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlageeBadScFilter",                      .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagCSCTightHaloFilter",                 .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagHBHEIsoNoiseFilterResult",           .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagHBHENoiseFilterResult",              .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagHBHENoiseFilterResultRun1",          .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagHBHENoiseFilterResultRun2Loose",     .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagHBHENoiseFilterResultRun2Tight",     .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="FlagEcalDeadCellTriggerPrimitiveFilter", .pfs={"AllSamples"}, .cuts={"HadTrigger","ExclBadHcal"}, .draw="HIST", .opt="Sumw2Norm", .ranges={0,0, 1e-3,1e5} });
  
  // Vertex/Pileup reweighting
  sh.AddHistos("evt",   { .fill="NVertices",                 .pfs={"AllSamples"}, .cuts={"HadTrigger","AllFilters"}, .draw="HIST", .opt="Sumw2LogStack4AddRatio", .ranges={0,50, 1.01e-3,1e5, 0.62,0.9} });
  
  if (debug) std::cout<<"Vertex plots ok\n";
  
#if VTX_REWEIGHT_MC == 1
  TFile *f_vtx = TFile::Open("ROOT_output/RunII_25nsDataNoHCal_AllMCNLO_VertexReweightingInput_fullstat.root");
  std::vector<std::string> samp = samples.GetListOfNames();
  std::vector<TH1D*> v_h;
  for (size_t i=0; i<samp.size(); ++i) {
    TH1D* h = (TH1D*)f_vtx->Get((std::string("NVertices/")+samp[i]).c_str());
    if (h->Integral()>0) h->Scale(1/h->Integral());
    v_h.push_back(h);
  }
  std::vector<std::vector<double> > vtx_weights;
  for (size_t imc=ndata; imc<samp.size(); ++imc) {
    vtx_weights.push_back(std::vector<double>());
    for (int nvtx=0; nvtx<=MAX_NVTX; ++nvtx)
      vtx_weights[imc-ndata].push_back(v_h[imc]->GetBinContent(nvtx+1)>0 ? v_h[0]->GetBinContent(nvtx+1) / v_h[imc]->GetBinContent(nvtx+1) : 0);
  }
  f_vtx->Close();
  
  // Set weights for all Histos, except Background estimate, where MC is scaled to MC_Expected_IntLumi_invfb
  sh.SetHistoWeights({[&ndata,&looper,&dir_to_index,&d,&sample_scale_factors,&vtx_weights](){ return
			  dir_to_index[looper.it_sample]<ndata ? 1 :
        		  /* Cross section weights: */ ( Data_IntLumi_invfb * d.evt.weight )
			  /* K-factors    */           * sample_scale_factors[looper.it_sample]
        		  /* Vertex weights:        */ * (d.evt.npv<=MAX_NVTX ? vtx_weights[dir_to_index[looper.it_sample]-ndata][d.evt.npv] : 0)
			  ;}});
  sh.AddHistos("evt",   { .fill="NVerticesReweighted", .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"}, .draw="HIST", .opt="Sumw2LogStack4AddRatio", .ranges={0,50, 1.01e-3,1e5, 0.62,0.9} });
#endif
  
  if (debug) std::cout<<"Reweighting ok\n";
  
  // --------------------------------------------------------------------------
  //                            Leptonic Top Selection
  
  //sh.AddHistos("evt",   { .fill="NTopLep",      .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NLep",         .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NEle",         .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NMu",          .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NLepTight",    .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NEleTight",    .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NLepVeto",     .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NEleVeto",     .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("evt",   { .fill="NMuVeto",      .pfs={"AllSamples"},           .cuts={}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,1.1} });
  //sh.AddHistos("jetsAK8", { .fill="JetDRLep",                .pfs={"AllSamples"},                 .cuts={"NLepTight==1"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetRelPtLep",             .pfs={"AllSamples"},                 .cuts={"NLepTight==1"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetRelPtLep_vs_JetDRLep", .pfs={"AllSamples"},                 .cuts={"NLepTight==1"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetRelPtLep_vs_JetDRLep", .pfs={"AllSamples","NSubJet"},       .cuts={"NLepTight==1"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("jetsAK8", { .fill="JetRelPtLep_vs_JetDRLep", .pfs={"AllSamples","JetsPtOrdered"}, .cuts={"NLepTight==1"}, .draw="COLZ", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("ele",   { .fill="EleJetCombMass",            .pfs={"AllSamples"},                 .cuts={"NEleTight==1","GoodEle"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("ele",   { .fill="EleDRJet",                  .pfs={"AllSamples"},                 .cuts={"NEleTight==1","GoodEle","EleJetCombMass>90"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("ele",   { .fill="EleRelPtJet",               .pfs={"AllSamples"},                 .cuts={"NEleTight==1","GoodEle","EleJetCombMass>90"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("mu",    { .fill="MuJetCombMass",             .pfs={"AllSamples"},                 .cuts={"NMuTight==1","GoodMu"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("mu",    { .fill="MuDRJet",                   .pfs={"AllSamples"},                 .cuts={"NMuTight==1","GoodMu","MuJetCombMass>90"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  //sh.AddHistos("mu",    { .fill="MuRelPtJet",                .pfs={"AllSamples"},                 .cuts={"NMuTight==1","GoodMu","MuJetCombMass>90"}, .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0, 0,0} });
  
  // --------------------------------------------------------------------------
  //                            Hadronic Top Selection
  //                       AK8 Jets, N-1 Cuts, Efficiencies
  
  // Since 03/10
  //sh.AddHistos("evt",   { .fill="NTopHad",   .pfs={"AllSamples"},      .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 1e-3,1e6} });
  //sh.AddHistos("evt",   { .fill="Sample",    .pfs={"NTopBands"},         .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 1e-3,1e6} });
  //sh.AddHistos("evt",   { .fill="Ht",        .pfs={"Directories"},     .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 1e-3,1e6} });
  //sh.AddHistos("evt",   { .fill="HtAll",     .pfs={"Directories"},     .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Log", .ranges={0,0, 1e-3,1e6} });
  
  // Jet cut variables Distributions
  // 1D
  // No Cut
  // sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau1",                               .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau2",                               .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau3",                               .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau21",                              .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau31",                              .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                       .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetNSubJets",                           .pfs={"Signals,TT,NonTT"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // // Apply 1 Cut
  // sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"Signals,TT,NonTT","JetMassCut"},  .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"Signals,TT,NonTT","JetMassCut"},  .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,NonTT","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,NonTT","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // // Apply 2 Cuts (N-1)
  // sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"Signals,TT,NonTT","JetMassCut","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"Signals,TT,NonTT","JetPtCut","JetMassCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // // Apply 3 Cuts
  // sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut","JetMassCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  // sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"Signals,TT,NonTT","JetTau32Cut","JetPtCut","JetMassCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  
  // Main variables (Shape, area, Data-MC agreement)
  // Hadronic top selection
  // No Cut
  sh.AddHistos("jetsAK8", { .fill="JetPt",           .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={200,2000, 1.01e-3,1e5, 0.55,0.9} }); // Note
  sh.AddHistos("jetsAK8", { .fill="JetTau2",         .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTau3",         .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",        .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} }); // Note
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass", .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",   .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",  .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass", .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  // Apply 1 Cut
  sh.AddHistos("jetsAK8", { .fill="JetPt",           .pfs={"StackSamples","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={200,2000, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetPt",           .pfs={"StackSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={200,2000, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",        .pfs={"StackSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",        .pfs={"StackSamples","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass", .pfs={"StackSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass", .pfs={"StackSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",   .pfs={"StackSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",   .pfs={"StackSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",  .pfs={"StackSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",  .pfs={"StackSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass", .pfs={"StackSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass", .pfs={"StackSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  // Apply 2 Cuts (N-1)
  sh.AddHistos("jetsAK8", { .fill="JetPt",           .pfs={"StackSamples","JetMassCut","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={200,2000, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",        .pfs={"StackSamples","JetPtCut","JetMassCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });  // Note
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass", .pfs={"StackSamples","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",   .pfs={"StackSamples","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",  .pfs={"StackSamples","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass", .pfs={"StackSamples","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
  
  // 2D
  // No Cut
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMassCoarse_vs_JetTau32",     .pfs={"AllSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMassCoarse_vs_JetPt",        .pfs={"AllSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetTau32",       .pfs={"AllSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetPt",          .pfs={"AllSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetPt_vs_JetTau3",                      .pfs={"AllSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });
  // Apply 1 Cut
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMassCoarse_vs_JetTau32",     .pfs={"AllSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMassCoarse_vs_JetPt",        .pfs={"AllSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetTau32",       .pfs={"AllSamples","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetPt",          .pfs={"AllSamples","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
  sh.AddHistos("jetsAK8", { .fill="JetPt_vs_JetTau32",                     .pfs={"AllSamples","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });  
  
  // Same plots, but use Gen Particle Truth
  sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // No Cut
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau1",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau2",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau3",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau21",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetNSubJets",                           .pfs={"JetGenTruth","Signals,Background"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // 1 Cut
  sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetMassCut","JetTau32Cut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // 2 Cuts
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                         .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",                        .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetTopMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut","JetMassCut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} }); // 3 Cuts
  sh.AddHistos("jetsAK8", { .fill="JetMinMass",                            .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut","JetMassCut"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMassCoarse_vs_JetTau32",     .pfs={"JetGenTruth","Signals,Background"},             .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 2D No Cut
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMassCoarse_vs_JetTau32",     .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 1 Cut
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetTau32",       .pfs={"JetGenTruth","Signals,Background"},             .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 2D No Cut
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMassCoarse_vs_JetTau32",       .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={"HadTrigger","AllFilters"},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 1 Cut
  
  //+NOT NEEDED+ // Jets with Matched GenTop (or constituents) are found in cone
  //+NOT NEEDED+ // --> Switch to Top to Jet matching (calc Efficiency)
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopPtCoarse",                     .pfs={"Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt",                                        .pfs={"Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt_vs_JetMatchedGenTopPtCoarse",            .pfs={"Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="Norm", .ranges={0,1000, 0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopPtCoarse",                     .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt",                                        .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt_vs_JetMatchedGenTopPtCoarse",            .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="Norm", .ranges={0,1000, 0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopPtCoarse",                     .pfs={"TopSizeCut","JetMatchedGenTopType","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt",                                        .pfs={"TopSizeCut","JetMatchedGenTopType","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt_vs_JetMatchedGenTopPtCoarse",            .pfs={"TopSizeCut","JetMatchedGenTopType","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="Norm", .ranges={0,1000, 0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"TopSizeCut","JetMatchedGenTopType","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopPtCoarse",                     .pfs={"JetMatchedGenTopType","TopSizeCut","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt",                                        .pfs={"JetMatchedGenTopType","TopSizeCut","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetPt_vs_JetMatchedGenTopPtCoarse",            .pfs={"JetMatchedGenTopType","TopSizeCut","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="Norm", .ranges={0,1000, 0,1000, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgJetPt_vs_JetMatchedGenTopPtCoarse",         .pfs={"JetMatchedGenTopType","TopSizeCut","Signals,Background"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,1000, 0,1000} });
  //+NOT NEEDED+ 
  //+NOT NEEDED+ // Distances - Samples
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopJetDR",                                  .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenBJetDR",                                              .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenWJetDR",                                              .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenWGenBDR",                                             .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenLepJetDR",                                            .pfs={"Signals,Background","TopSizeCut"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenLepGenBDR",                                           .pfs={"Signals,Background","TopSizeCut"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopJetDR_vs_JetMatchedGenTopPtBins",        .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenBJetDR_vs_JetMatchedGenTopPtBins",                    .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenWJetDR_vs_JetMatchedGenTopPtBins",                    .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenWGenBDR_vs_JetMatchedGenTopPtBins",                   .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, } });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenLepJetDR_vs_JetMatchedGenTopPtBins",                  .pfs={"Signals,Background","TopSizeCut"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenLepGenBDR_vs_JetMatchedGenTopPtBins",                 .pfs={"Signals,Background","TopSizeCut"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgJetMatchedGenTopJetDRFine_vs_JetMatchedGenTopPtBins", .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgGenBJetDRFine_vs_JetMatchedGenTopPtBins",             .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgGenWJetDRFine_vs_JetMatchedGenTopPtBins",             .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="AvgGenWGenBDRFine_vs_JetMatchedGenTopPtBins",            .pfs={"Signals,Background","TopSizeCut","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgGenLepJetDRFine_vs_JetMatchedGenTopPtBins",           .pfs={"Signals,Background","TopSizeCut"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="AvgGenLepGenBDRFine_vs_JetMatchedGenTopPtBins",          .pfs={"Signals,Background","TopSizeCut"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ // Distances - Merged/Type
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopJetDR",                                  .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenBJetDR",                                              .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenWJetDR",                                              .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenWGenBDR",                                             .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenLepJetDR",                                            .pfs={"TopSizeCut","Signals,Background"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenLepGenBDR",                                           .pfs={"TopSizeCut","Signals,Background"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="Norm", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="JetMatchedGenTopJetDR_vs_JetMatchedGenTopPtBins",        .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenBJetDR_vs_JetMatchedGenTopPtBins",                    .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenWJetDR_vs_JetMatchedGenTopPtBins",                    .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenWGenBDR_vs_JetMatchedGenTopPtBins",                   .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, } });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="GenLepJetDR_vs_JetMatchedGenTopPtBins",                  .pfs={"TopSizeCut","Signals,Background"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="GenLepGenBDR_vs_JetMatchedGenTopPtBins",                 .pfs={"TopSizeCut","Signals,Background"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="COLZ", .opt="", .ranges={0,0, 0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgJetMatchedGenTopJetDRFine_vs_JetMatchedGenTopPtBins", .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgGenBJetDRFine_vs_JetMatchedGenTopPtBins",             .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgGenWJetDRFine_vs_JetMatchedGenTopPtBins",             .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="AvgGenWGenBDRFine_vs_JetMatchedGenTopPtBins",            .pfs={"TopSizeCut","Signals,Background","JetMatchedGenTopType"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="AvgGenLepJetDRFine_vs_JetMatchedGenTopPtBins",           .pfs={"TopSizeCut","Signals,Background"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  //+NOT NEEDED+ //sh.AddHistos("jetsAK8", { .fill="AvgGenLepGenBDRFine_vs_JetMatchedGenTopPtBins",          .pfs={"TopSizeCut","Signals,Background"},                        .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop","GenLepTop"}, .draw="HISTE1",     .opt="", .ranges={0,0, 0,0} });
  
  // Jet Finding Efficiency
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT"},              .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });  // Note
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT"},              .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT","GenTopType"}, .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT","GenTopType"}, .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT"},              .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT"},              .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT","GenTopType"}, .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT","GenTopType"}, .cuts={"HadTrigger","AllFilters","IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
  // Fraction of merged sub-jets
  sh.AddHistos("jetsAK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT"},                                   .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT","JetMatchedGenTopType"},            .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"Signals,TT"},                                   .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"Signals,TT","JetMatchedGenTopType"},            .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
  // Top Tag/Finding Efficiency
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT"},              .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT","TopSizeCut"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"Signals,TT"},              .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
  sh.AddHistos("jetsAK8", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"Signals,TT","TopSizeCut"}, .cuts={"HadTrigger","AllFilters","JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
  
  //+NOT NEEDED+ // Check it only in TTJets samples (For all Bkg the rate is high)
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtBins",       .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters","JetIsHadTopTagged"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtOneBin",     .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters","JetIsHadTopTagged"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtFewBins",    .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtFewBins",    .pfs={"JetMassCut","Signals,TT"},                      .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtFewBins",    .pfs={"JetTau32Cut","Signals,TT"},                     .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtFewBins",    .pfs={"JetTau32Cut","JetMassCut","Signals,TT"},          .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPtFewBins",    .pfs={"JetMassCut","JetTau32Cut","Signals,TT"},          .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetTau32",        .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetTau32",        .pfs={"JetMassCut","Signals,TT"},                      .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetTau32",        .pfs={"JetPtCut","Signals,TT"},                        .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetTau32",        .pfs={"JetPtCut","JetMassCut","Signals,TT"},             .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetTau32",        .pfs={"JetMassCut","JetPtCut","Signals,TT"},             .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetSoftDropMass", .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetSoftDropMass", .pfs={"JetTau32Cut","Signals,TT"},                     .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetSoftDropMass", .pfs={"JetPtCut","Signals,TT"},                        .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetSoftDropMass", .pfs={"JetPtCut","JetTau32Cut","Signals,TT"},            .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetSoftDropMass", .pfs={"JetTau32Cut","JetPtCut","Signals,TT"},            .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPrunedMass",   .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPrunedMass",   .pfs={"JetTau32Cut","Signals,TT"},                     .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPrunedMass",   .pfs={"JetPtCut","Signals,TT"},                        .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPrunedMass",   .pfs={"JetPtCut","JetTau32Cut","Signals,TT"},            .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetPrunedMass",   .pfs={"JetTau32Cut","JetPtCut","Signals,TT"},            .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"Signals,TT"},                                 .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetMassCut","Signals,TT"},                      .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetTau32Cut","Signals,TT"},                     .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetPtCut","Signals,TT"},                        .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetMassCut","JetTau32Cut","Signals,TT"},          .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetMassCut","JetPtCut","Signals,TT"},             .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetTau32Cut","JetPtCut","Signals,TT"},            .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_MaxSubJetCSV",    .pfs={"JetMassCut","JetTau32Cut","JetPtCut","Signals,TT"}, .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetMinMass",      .pfs={"JetMassCut","JetTau32Cut","JetPtCut","Signals,TT"}, .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetTopMass",      .pfs={"JetMassCut","JetTau32Cut","JetPtCut","Signals,TT"}, .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  //+NOT NEEDED+ sh.AddHistos("jetsAK8", { .fill="MisTagRate_vs_JetNSubJets",     .pfs={"JetMassCut","JetTau32Cut","JetPtCut","Signals,TT"}, .cuts={"HadTrigger","AllFilters"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
  if (debug) std::cout<<"Jet histos ok\n";
  
  // --------------------------------------------------------------------------
  //                             Trigger Selection
  
  // Trigger Efficiencies
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_TTHadSumPt",       .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={800,2000, 0.8,1.05, 0.4,0.55} });
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_TTHadSumPtOneBin", .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={0,0,      0.8,1.05, 0.4,0.55} });
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_HtBins",           .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={0,2000,   0.8,1.05, 0.4,0.55} });
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_HtOneBin",         .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={0,0,      0.8,1.05, 0.4,0.55} });
  
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_TTHadSumPt",       .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={800,2000, 0.8,1.05, 0.4,0.55} });
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_TTHadSumPtOneBin", .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={0,0,      0.8,1.05, 0.4,0.55} });
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_HtBins",           .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={0,2000,   0.8,1.05, 0.4,0.55} });
  sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFHT700TrimR0p1PT0p03Mass50_vs_HtOneBin",         .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="Sumw2AddRatio", .ranges={0,0,      0.8,1.05, 0.4,0.55} });
  
  //sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFJet360TrimMass30_vs_HtBins",            .pfs={"Data,Background,Signals"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="", .ranges={0,1500, 0,1} });
  //sh.AddHistos("evt",   { .fill="HLTEfficiencyPFHT800_vs_HtBins",                          .pfs={"Data,Background,Signals"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="", .ranges={0,1500, 0,1} });
  //sh.AddHistos("evt",   { .fill="HLTEfficiencyPFHT900_vs_HtBins",                          .pfs={"Data,Background,Signals"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="", .ranges={0,1500, 0,1} });
  //sh.AddHistos("evt",   { .fill="HLTEfficiencyAK8PFJet360TrimMass30_vs_HtOneBin",          .pfs={"Data,Background,Signals"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="", .ranges={0,1500, 0,1} });
  //sh.AddHistos("evt",   { .fill="HLTEfficiencyPFHT800_vs_HtOneBin",                        .pfs={"Data,Background,Signals"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="", .ranges={0,1500, 0,1} });
  //sh.AddHistos("evt",   { .fill="HLTEfficiencyPFHT900_vs_HtOneBin",                        .pfs={"Data,Background,Signals"}, .cuts={"AllFilters","LowHTTrig","NTopHadPreTag>=2"}, .draw="PE1", .opt="", .ranges={0,1500, 0,1} });
  
  if (debug) std::cout<<"HLT histos ok\n";
  
  // --------------------------------------------------------------------------
  //                             Data-MC Comparison
  //                                Stack plots
  
  // Jets - Apply all top-tag cuts
  sh.AddHistos("jetsAK8", { .fill="JetPt",                        .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={200,2000, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTau32",                     .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetSoftDropMass",              .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,300,  1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetPrunedMass",                .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,300,  1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetTrimmedMass",               .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,300,  1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetFilteredMass",              .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,300,  1.01e-2,1e6, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetEta",                       .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetPhi",                       .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetEta",                       .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetPhi",                       .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetNeutralHadronMultiplicity", .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetChargedHadronMultiplicity", .pfs={"StackSamples"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetNeutralHadronMultiplicity", .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("jetsAK8", { .fill="JetChargedHadronMultiplicity", .pfs={"StackSamples","JetPassToptag"}, .cuts={"HadTrigger","AllFilters"},  .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={  0,0,    1.01e-3,1e5, 0.55,0.9} });
  
  // MET, Razor Variables
  sh.AddHistos("met", { .fill="MetPt",            .pfs={"StackSamples"},                              .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
  sh.AddHistos("met", { .fill="MetPt",            .pfs={"StackSamples","NTopBands"},                    .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("met", { .fill="MetPt",            .pfs={"StackSamples","NTopBands","DPhiBands"},          .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("met", { .fill="MetPt",            .pfs={"StackSamples","NTopBands","RBands"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("met", { .fill="MetPt",            .pfs={"StackSamples","NTopBands","RBands","DPhiBands"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
  sh.AddHistos("evt", { .fill="TTHadMR",          .pfs={"StackSamples"},                              .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="TTHadMR",          .pfs={"StackSamples","NTopBands"},                    .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="TTHadMR",          .pfs={"StackSamples","NTopBands","DPhiBands"},          .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MR",               .pfs={"StackSamples"},                              .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MR",               .pfs={"StackSamples","NTopBands"},                    .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MR",               .pfs={"StackSamples","NTopBands","DPhiBands"},          .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MTR",              .pfs={"StackSamples"},                              .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MTR",              .pfs={"StackSamples","NTopBands"},                    .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  sh.AddHistos("evt", { .fill="MTR",              .pfs={"StackSamples","NTopBands","DPhiBands"},          .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogStack4Sumw2AddRatio", .ranges={0,0, 1e-3,1e6} });
  
  // 2D Correlation plots
  sh.AddHistos("evt", { .fill="R_vs_DPhi",        .pfs={"AllSamples"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="R_vs_DPhi",        .pfs={"AllSamples","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadR_vs_DPhi",   .pfs={"AllSamples"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadR_vs_DPhi",   .pfs={"AllSamples","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MR_vs_DPhi",       .pfs={"AllSamples"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MR_vs_DPhi",       .pfs={"AllSamples","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadMR_vs_DPhi",  .pfs={"AllSamples"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadMR_vs_DPhi",  .pfs={"AllSamples","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_DPhi",      .pfs={"AllSamples"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_DPhi",      .pfs={"AllSamples","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"AllSamples"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"AllSamples","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"AllSamples","DPhiBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"AllSamples","NTopBands","DPhiBands"},          .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="R_vs_DPhi",        .pfs={"Background"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="R_vs_DPhi",        .pfs={"Background","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadR_vs_DPhi",   .pfs={"Background"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadR_vs_DPhi",   .pfs={"Background","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MR_vs_DPhi",       .pfs={"Background"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MR_vs_DPhi",       .pfs={"Background","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadMR_vs_DPhi",  .pfs={"Background"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="TTHadMR_vs_DPhi",  .pfs={"Background","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_DPhi",      .pfs={"Background"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_DPhi",      .pfs={"Background","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"Background"},                                  .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"Background","NTopBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"Background","DPhiBands"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="MTR_vs_MR",        .pfs={"Background","NTopBands","DPhiBands"},          .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgR_vs_DPhi",        .pfs={"NTopBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgTTHadR_vs_DPhi",   .pfs={"NTopBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMR_vs_DPhi",       .pfs={"NTopBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgTTHadMR_vs_DPhi",  .pfs={"NTopBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_DPhi",      .pfs={"NTopBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"NTopBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"NTopBands","DPhiBands","AllSamples"},          .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","AllSamples"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","NTopBands","AllSamples"},          .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgR_vs_DPhi",        .pfs={"NTopBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgTTHadR_vs_DPhi",   .pfs={"NTopBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMR_vs_DPhi",       .pfs={"NTopBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgTTHadMR_vs_DPhi",  .pfs={"NTopBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_DPhi",      .pfs={"NTopBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"NTopBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"NTopBands","DPhiBands","Background"},          .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Background"},                      .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","NTopBands","Background"},          .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  
  // Signal selection (Apply loose pretag selection)
  sh.AddHistos("evt",   { .fill="R",                 .pfs={"StackSamples"},                     .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt",   { .fill="DPhi",              .pfs={"StackSamples"},                     .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt",   { .fill="R",                 .pfs={"StackSamples","NTopBands"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt",   { .fill="DPhi",              .pfs={"StackSamples","NTopBands"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt",   { .fill="R",                 .pfs={"StackSamples","NTopBands","DPhiBands"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  sh.AddHistos("evt",   { .fill="DPhi",              .pfs={"StackSamples","NTopBands","RBands"},  .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="LogSumw2Stack4AddRatio", .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  
  if (debug) std::cout<<"Data-MC histos ok\n";
  
  // --------------------------------------------------------------------------
  //                              Signal Selection
  //                            Background Estimation
  //                          Razor variables, Deltaphi
  
  // Set weights for all Histos, except Background estimate, where MC is scaled to MC_Expected_IntLumi_invfb
#if VTX_REWEIGHT_BKG_EST == 1
  sh.SetHistoWeights({[&ndata,&looper,&dir_to_index,&d,&sample_scale_factors,&vtx_weights](){ return
#else
  sh.SetHistoWeights({[&ndata,&looper,&dir_to_index,&d,&sample_scale_factors](){ return
#endif
			  dir_to_index[looper.it_sample]<ndata ? MC_Expected_IntLumi_invfb/Data_IntLumi_invfb :
        		  /* Cross section weights: */ ( MC_Expected_IntLumi_invfb * d.evt.weight )
			  /* K-factors    */           * sample_scale_factors[looper.it_sample]
#if VTX_REWEIGHT_BKG_EST == 1
        		  /* Vertex weights:        */ * ( d.evt.npv<=MAX_NVTX ? vtx_weights[dir_to_index[looper.it_sample]-ndata][d.evt.npv] : 0 )
#endif
			  ;}});
  
  // 3D Plots to get best signal cuts (Maximize Smin) --> input for B2GAnalyzer
  sh.AddHistos("evt",   { .fill="MR_vs_DPhiFine_vs_RFine",           .pfs={"AllSamples","NTopBands"}, .cuts={"BaselineCuts"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="MR_vs_DPhiFine_vs_TTHadRFine",      .pfs={"AllSamples","NTopBands"}, .cuts={"BaselineCuts"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="TTHadMR_vs_DPhiFine_vs_RFine",      .pfs={"AllSamples","NTopBands"}, .cuts={"BaselineCuts"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("evt",   { .fill="HtAll_vs_DPhiFine_vs_RFine",        .pfs={"AllSamples","NTopBands"}, .cuts={"BaselineCuts"}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  
  // R plots - Distributions for All samples, All cut combinations, Ratios
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"DPhiBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"NTopBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"NTopBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  // ABCD regions - DPhi<2.7, define regions by: R and Ntop
  sh.AddHistos("evt",   { .fill="R",         .pfs={"ABCD,Signals","DPhiBands"},      .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="R",         .pfs={"ABCD","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="R",         .pfs={"ABCD","DPhiBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="R",         .pfs={"NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  sh.AddHistos("evt",   { .fill="R",         .pfs={"NTopBands","DPhiBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
  sh.AddHistos("evt",   { .fill="RFine",     .pfs={"NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RFine",     .pfs={"NTopBands","DPhiBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"NTopBands","DPhiBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"DPhiBands","NTopBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="RBins",     .pfs={"DPhiBands","NTopBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  
  // DPhi plots - Distributions for All samples, All cut combinations, Ratios
  sh.AddHistos("evt",   { .fill="DPhi",      .pfs={"ABCD,Signals","DPhiBands"},      .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBands", "AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBands", "Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBins",  "AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBins",  "Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"NTopBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"NTopBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  // Alternative ABCD regions - Ntop==2, define regions by: R and DPhi
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBands",   "NTopBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBands",   "NTopBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBins",    "NTopBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"RBins",    "NTopBands","Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"NTopBands","RBins",    "AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("evt",   { .fill="DPhiBins",  .pfs={"NTopBands","RBins",    "Background"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  
  //sh.AddHistos("evt",   { .fill="R2_vs_TTHadMRCoarse", .pfs={"AllSamples"},           .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2_vs_TTHadMRCoarse", .pfs={"DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2_vs_TTHadMRCoarse", .pfs={"NTopBands","AllSamples"},           .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2_vs_TTHadMRCoarse", .pfs={"NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"StackSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2LogStack4", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"NTopBands","AllSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"StackSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2LogStack4", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"NTopBands","AllSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"TTHadMRBins","AllSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"TTHadMRBins","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"TTHadMRBins","NTopBands","AllSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="R2",      .pfs={"TTHadMRBins","NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"RBins","NTopBands","AllSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"RBins","NTopBands","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"RBins","AllSamples"},           .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  //sh.AddHistos("evt",   { .fill="TTHadMRCoarse", .pfs={"RBins","DPhiBands","AllSamples"}, .cuts={"BaselineCuts"}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  
  if (debug) std::cout<<"Signal/Background histos ok\n";
  
  std::cout<<"-----------------------------------------------------------------\n";
  //std::cout<<"Creating the following plots:\n"; sh.PrintNames();
  //std::cout<<"-----------------------------------------------------------------\n";
  
  // Merge files
  //sh.Load("ROOT_output/RunII_Data_MC_StackPlots_fullstat.root");
  
  TFile *file;
  if (Run) {
    if (filelist_fromshell.size()) {
      std::cout<<"Adding "<<filelist_fromshell.size()<<" files from the shell arguments.\n";
      for (size_t i=0; i<filelist_fromshell.size(); ++i) looper.AddFile(filelist_fromshell[i], !i);
    } else {
      std::vector<std::string> dirs = samples.GetListOfDirectories();
      for ( std::string dir : dirs ) if (dir!="") looper.AddFile(dir);
    }
    if (debug) cout<<"Start ok\n";
    while (looper.LoopOnSamples()) {
      if (debug) cout<<"Sample ok\n";
      while (looper.LoopOnFiles()) {
	if (debug) cout<<"File ok\n";
	TFile *curr_file = looper.CurrentFile();
	reader.Load_Tree(*curr_file,looper.TreeName());
	if (debug) cout<<"TreeReader ok\n";
	while(looper.LoopOnEntries()) {
	  reader.GetEntry(looper.CurrentEntry());
	  d = reader.data;
	  if (debug) cout<<"GetEntry ok\n";
	  d.CalculateAllVariables();
	  if (debug) cout<<"CalcVar ok\n";
	  // loop on objects and fill their histos
	  while(d.mu.Loop())      sh.Fill("mu");  if (debug) cout<<"Fill Muons ok\n";
	  while(d.ele.Loop())     sh.Fill("ele"); if (debug) cout<<"Fill Electrons ok\n";
	  while(d.jetsAK4.Loop()) sh.Fill("jetsAK4"); if (debug) cout<<"Fill AK4Jets ok\n";
	  while(d.jetsAK8.Loop()) sh.Fill("jetsAK8"); if (debug) cout<<"Fill AK8Jets ok\n";
	  while(d.gen.Loop()) sh.Fill("gen"); if (debug) cout<<"Fill Gen ok\n";
	  sh.Fill("met"); if (debug) cout<<"Fill MET ok\n";
	  sh.Fill("evt"); if (debug) cout<<"Fill Evt ok\n";
	}
	curr_file->Close();
      }
    }
    std::cout<<std::endl;
  } else {
    std::cout<<"Loading Histos from file: "<<inputfile<<std::endl;
    sh.Load(inputfile.c_str());
  }
  std::cout<<"Finished ..."<<std::endl;
  std::cout<<"Writing Histograms to File: "<<outputfile<<std::endl;
  file = new TFile(outputfile.c_str(),"recreate");
  sh.DrawPlots();
  sh.Write();
  file->Close();
  
  return 0;
}
