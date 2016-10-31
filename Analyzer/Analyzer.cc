//-----------------------------------------------------------------------------
// File:        Analyzer.cc
// Created:     24-Nov-2015
// Author:      Janos Karancsi
//-----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "settings_Viktor.h" // Define all Analysis specific settings here

using namespace std;

int main(int argc, char** argv) {
  // List names in filenames from which the code can decide if it is data or signal
  // For the rest it's assumed it's background MC
  // if .txt file is given as input then from the directory name we can already tell
  std::vector<std::string> vname_data = { "JetHT", "MET", "SingleEle", "SingleMu", 
					  "2015B", "2015C", "2015D", 
					  "2016B", "2016C", "2016D", "2016E", 
					  "2016F", "2016G", "2016H", "2016I", 
					  "2016J", "2016K", "2016L", "2016M" };
  std::vector<std::string> vname_signal = { "T1tttt", "T5ttcc", "T5tttt", "T2tt" }; // SMS

  // ------------------------------
  // -- Parse command line stuff --
  // ------------------------------

  // Get file list and histogram filename from command line
  utils::commandLine cmdline;
  utils::decodeCommandLine(argc, argv, cmdline, vname_data, vname_signal);

  itreestream stream(cmdline.fileNames, settings.treeName);
  if ( !stream.good() ) utils::error("unable to open ntuple file(s)");

  if ( cmdline.isData ) cout << "Running on Data." << endl;
  else if ( cmdline.isBkg ) cout << "Running on Background MC." << endl;
  else cout << "Running on Signal MC." << endl;

  // Get number of events to be read
  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

  if (cmdline.quickTest) {
    cout << "quickTest (cmdline): true"<< endl;
    cout << "--> Doing a quick test on 1/100 statitics"<< endl;
    nevents /= 100;
  }

  // Select variables to be read
  DataStruct data;
  settings.selectVariables(stream, data);

  /*
	 Notes:
	
	 1. Use
	   ofile = outputFile(cmdline.outputfile, stream)
	
	 to skim events to output file in addition to writing out histograms.
	
	 2. Use
	   ofile.addEvent(event-weight)
	
	 to specify that the current event is to be added to the output file.
	 If omitted, the event-weight is defaulted to 1.
	
	 3. Use
	    ofile.count(cut-name, event-weight)
	
	 to keep track, in the count histogram, of the number of events
	 passing a given cut. If omitted, the event-weight is taken to be 1.
	 If you want the counts in the count histogram to appear in a given
	 order, specify the order, before entering the event loop, as in
	 the example below
	 
	    ofile.count("NoCuts", 0)
		ofile.count("GoodEvent", 0)
		ofile.count("Vertex", 0)
		ofile.count("MET", 0)
  */

  // Constuct the Analysis (specified in settings.h)
  Analysis ana;

  // ---------------------------------------------------------------------------
  // -- output file                                                           --
  // ---------------------------------------------------------------------------

  utils::outputFile* ofile;
  if ( settings.saveSkimmedNtuple ) {
    cout << "saveSkimmedNtuple (settings): true" << endl;
    ofile = new utils::outputFile(cmdline.outputFileName, stream);
  } else {
    ofile = new utils::outputFile(cmdline.outputFileName);
    cout << "saveSkimmedNtuple (settings): false" << endl;
  }
  TDirectory* out_dir = gDirectory;

  if (cmdline.noPlots) {
    cout << "noPlots (cmdline): true"<< endl;
    cout << "--> Will not save analysis histos"<< endl;
  }

  // ---------------------------------------------------------------------------
  // -- Read systematics file                                                 --
  // ---------------------------------------------------------------------------

  // Initialize all systemtics variables (0 = default/mean, no variation)
  struct Systematics {
    unsigned int index = 0;
    unsigned int nSyst = 0;
    std::vector<double> nSigmaLumi        = std::vector<double>(1,0);
    std::vector<double> nSigmaPU          = std::vector<double>(1,0);
    std::vector<double> nSigmaTrigger     = std::vector<double>(1,0);
    std::vector<double> nSigmaJEC         = std::vector<double>(1,0);
    std::vector<double> nSigmaWTagSF      = std::vector<double>(1,0);
    std::vector<double> nSigmaBTagSF      = std::vector<double>(1,0);
    //std::vector<double> nSigmaHadTopTagSF = std::vector<double>(1,0);
    std::vector<double> nSigmaHT          = std::vector<double>(1,0);
    std::vector<double> nSigmaAlphaS      = std::vector<double>(1,0);
    std::vector<double> nSigmaScale       = std::vector<double>(1,0);
    std::vector<unsigned int> numScale    = std::vector<unsigned int>(1,0);
    std::vector<unsigned int> numPdf      = std::vector<unsigned int>(1,0);
  } syst;

  if (settings.varySystematics) {
    cout << "varySystematics (settings): true" << endl;
    cout << "systematicsFileName (settings): " << settings.systematicsFileName << endl;
    std::ifstream systFile(settings.systematicsFileName.c_str());
    if ( !systFile.good() ) utils::error("unable to open systematics file: " + settings.systematicsFileName);

    // Read all nSigmas, nums
    double dbl = 0;
    unsigned int uint = 0;
    std::string line;
    std::cout<<"Systematics read from file:"<<std::endl;
    while ( std::getline(systFile, line) ) {
      ++syst.nSyst;
      std::stringstream nth_line;
      nth_line<<line;
      nth_line>>dbl; syst.nSigmaLumi.push_back(dbl);
      nth_line>>dbl; syst.nSigmaPU.push_back(dbl);
      nth_line>>dbl; syst.nSigmaTrigger.push_back(dbl);
      nth_line>>dbl; syst.nSigmaJEC.push_back(dbl);
      nth_line>>dbl; syst.nSigmaWTagSF.push_back(dbl);
      nth_line>>dbl; syst.nSigmaBTagSF.push_back(dbl);
      //nth_line>>dbl; syst.nSigmaHadTopTagSF.push_back(dbl);
      nth_line>>dbl; syst.nSigmaHT.push_back(dbl);
      nth_line>>dbl; syst.nSigmaAlphaS.push_back(dbl);
      nth_line>>dbl; syst.nSigmaScale.push_back(dbl);
      nth_line>>uint; syst.numScale.push_back(uint);
      nth_line>>uint; syst.numPdf.push_back(uint);
      std::cout<<" line "<<syst.nSyst<<": "<<line<<std::endl;
    }
    std::cout<<std::endl;

    /*
    if ( cmdline.numSyst <= 0 ) utils::error("varySystematics true, but command line argument numSyst=<positive integer> was not given");
    cout << "numSyst (cmdline): " << cmdline.numSyst << endl;

    // read nth line
    // error if end of file reached
    std::string line;
    for (int i=1; i<=cmdline.numSyst; ++i) if ( ! std::getline(systFile, line) ) 
      utils::error("systematics file contains less lines than numSyst");
    std::stringstream nth_line;
    nth_line<<line;

    // Assign the systematic sigmas:
    nth_line>>syst.nSigmaPU;
    nth_line>>syst.nSigmaJEC;
    nth_line>>syst.nSigmaAlphaS;
    nth_line>>syst.nSigmaScale;
    nth_line>>syst.numScale;
    nth_line>>syst.numPdf;
    cout << " nSigmaPU     = " << syst.nSigmaPU << endl;
    cout << " nSigmaJEC    = " << syst.nSigmaJEC << endl;
    cout << " nSigmaAlphaS = " << syst.nSigmaAlphaS << endl;
    cout << " nSigmaScale  = " << syst.nSigmaScale << endl;
    cout << " numScale     = " << syst.numScale << endl;
    cout << " numPdf       = " << syst.numPdf << endl;
    */
    
  } else {
    cout << "varySystematics (settings): false" << endl;
  }


  // ---------------------------------------------------------------------------
  // -- Read and apply JSON file (Data)                                       --
  // ---------------------------------------------------------------------------

  std::map<int, std::map<int, bool> > json_run_ls;
  if (settings.useJSON) {
    cout << "useJSON (settings): true" << endl;
    cout << "jsonFileName (settings): " << settings.jsonFileName << endl;
    std::ifstream jsonFile(settings.jsonFileName.c_str());
    if ( !jsonFile.good() ) utils::error("unable to open systematics file: " + settings.jsonFileName);

    std::string line;
    int run, ls_low, ls_high;
    while ( std::getline(jsonFile, line, ' ') ) {
      if (TString(line).Contains("\"")) {
	std::stringstream run_str;
	run_str<<line.substr(line.find("\"")+1, 6);
	run_str>>run;
      } else {
	if (TString(line).Contains("[")) {
	  std::stringstream ls_low_str;
	  while (line.find("[")!=std::string::npos) line.erase(line.find("["),1);
	  ls_low_str<<line;
	  ls_low_str>>ls_low;
	} else if (TString(line).Contains("]")) {
	  std::stringstream ls_high_str;
	  ls_high_str<<line;
	  ls_high_str>>ls_high;
	  for (int ls=ls_low; ls<=ls_high; ++ls)
	    json_run_ls[run].insert(std::pair<int, bool>(ls, 1));
	}
      }
    }
  } else {
    cout << "useJSON (settings): false" << endl;
  }


  // ---------------------------------------------------------------------------
  // -- Declare histograms                                                    --
  // ---------------------------------------------------------------------------

  // Histogram weight
  double w;
  TH1::SetDefaultSumw2();

  // This section is moved and defined in the Analysis class
  // this class specifies what histograms/cuts/methods etc. you want to implement
  // This is useful if someone wants to do quick study/define other search region etc.
  // But also, common methods in all anaylsis are defined in common/AnalysisBase.*

  if (!cmdline.noPlots)
    ana.define_histo_options(w, data, syst.nSyst, syst.index, cmdline.dirname, settings.runOnSkim);

  ana.init_common_histos();
  if (!cmdline.noPlots)
    ana.init_analysis_histos(syst.nSyst, syst.index);

  // --------------------------------------------------------------
  // -- Calculate the normalization factor for the event weights --
  // -- The original MC weight will be divided by this quantity  --
  // --------------------------------------------------------------

  cout << endl;
  double weightnorm = 1;
  int signal_index = -1;
  if ( cmdline.isBkg ) {
    cout << "intLumi (settings): " << settings.intLumi << endl; // given in settings.h

    cout << "useXSecFileForBkg (settings): " << ( settings.useXSecFileForBkg ? "true" : "false" ) << endl; // given in settings.h

    double xsec = 0;
    if (settings.useXSecFileForBkg) {
      cout << "xSecFileName (settings): " << settings.xSecFileName << endl; // given in settings.h
      xsec = ana.get_xsec_from_txt_file(settings.xSecFileName, cmdline.dirname); // xSecFileName given in settings.h
      cout << "xsec (txt file): " << xsec << endl;      
    } else {
      xsec = ana.get_xsec_from_ntuple(cmdline.fileNames, settings.treeName); // treename given in settings.h
      cout << "xsec (ntuple): " << xsec << endl;
    }
    if ( xsec==0 ) return 1;

    double totweight = ana.get_totweight_from_ntuple(cmdline.fileNames, settings.totWeightHistoName); // weight histo name given in settings.h
    cout << "totweight (ntuple): " << totweight << endl;

    weightnorm = (settings.intLumi*xsec)/totweight;
    cout << "weightnorm (calc): " << weightnorm << endl;
  } else if ( cmdline.isSignal ) {
    cout << "intLumi (settings): " << settings.intLumi << endl; // given in settings.h

    cout << "Normalization variables:" << endl;
    ana.calc_weightnorm_histo_from_ntuple(cmdline.fileNames, settings.intLumi, vname_signal,
					  settings.totWeightHistoNamesSignal); // histo names given in settings.h

    // Find the index of the current signal
    if (cmdline.fileNames.size()>0) for (size_t i=0, n=vname_signal.size(); i<n; ++i) 
      if (cmdline.fileNames[0].find(vname_signal[i])!=std::string::npos&&signal_index==-1) 
	signal_index = (i>=3); // 0: T1tttt, T5ttcc, T5tttt; 1: T2tt
  }

  // ---------------------------------------
  // --- ScaleFectors/Reweighting        ---
  // ---------------------------------------
  
  // Pile-up reweighting
  if ( !cmdline.isData && settings.doPileupReweighting ) {
    cout << "doPileupReweighting (settings): true" << endl;
    ana.init_pileup_reweightin(settings.pileupDir, settings.mcPileupHistoName, cmdline.fileNames);
  } else cout << "doPileupReweighting (settings): false" << endl;

  // Scale factors
  cout << "applyWTagSF (settings): " << ( settings.applyWTagSF ? "true" : "false" ) << endl;
  cout << "applyBTagSF (settings): " << ( settings.applyBTagSF ? "true" : "false" ) << endl;
  //cout << "applyHadTopTagSF (settings): " << ( settings.applyHadTopTagSF ? "true" : "false" ) << endl;

  // Scale QCD to match data in a QCD dominated region
  cout << "scaleQCD (settings): " << ( settings.scaleQCD ? "true" : "false" ) << endl;

  // Jet Pt Reweighting
  cout << "doHTReweighting (settings): " << ( settings.doHTReweighting ? "true" : "false" ) << endl;

  // ------------------------------------------------------------------------------
  // -- Define the order of cuts (and corresponding bins in the counts histogram --
  // ------------------------------------------------------------------------------

  // Define cuts that are common in all analyses
  // Given in common/AnalysisBase.h
  ana.define_preselections(data, cmdline.isData, cmdline.isSignal);

  // Define cuts that specific to this analysis
  // Given in [Name]_Analysis.h specified in setting.h
  ana.define_selections(data);

  // Define bin order for counts histogram
  ofile->count("NoCuts", 0);
  cout << endl;
  cout << "Number of events counted after applying" << endl;
  cout << "- Baseline cuts (common for all analysis):" << endl;
  for (auto cut : ana.baseline_cuts) {
    ofile->count(cut.name, 0);
    cout << "  "<<cut.name << endl;
  }
  cout << endl;  
  cout << "- Analysis specific cuts:\n";
  for (auto cut : ana.analysis_cuts) {
    ofile->count(cut.name, 0);
    cout << "  "<<cut.name << endl;
  }

  ofile->count("Signal", 0); // Dont worry, we blind data ;)

  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  cout << endl;
  cout << "Start looping on events ..." << endl;
  for(int entry=0; entry < nevents; ++entry) {

    // Read event into memory
    stream.read(entry);

    if ( entry%100000==0 ) cout << entry << " events analyzed." << endl;

    if ( cmdline.isData ) {
      w = 1;
      syst.index = 0;

      // Only analyze events that are in the JSON file
      if (settings.useJSON ? json_run_ls[data.evt.RunNumber][data.evt.LumiBlock] : 1) {

	// Calculate variables that do not exist in the ntuple
	ana.calculate_common_variables(data, syst.index);
	ana.calculate_variables(data, syst.index);

	// Save counts (after each cuts)
	ofile->count("NoCuts", w);
	bool pass_all_cuts = true;
	for (auto cut : ana.baseline_cuts) if (pass_all_cuts)
	  if ( ( pass_all_cuts = cut.func() ) ) ofile->count(cut.name, w);

	// If option (saveSkimmedNtuple) is specified save all 
	// skimmed events selected by the analysis to the output file
	// tree is copied and current weight is saved as "eventWeight"
	if ( settings.saveSkimmedNtuple ) if (ana.pass_skimming(data)) ofile->addEvent(w);

	// _______________________________________________________
	//                  BLINDING DATA

	// Before doing anything serious (eg. filling any histogram)
	// Define Signal region and blind it!!!!!!!!!!!!!!!!!!!!!!!!
	bool DATA_BLINDED = ! ( cmdline.isData && ana.signal_selection(data) );

	/*
	  Some more warning to make sure :)
	  _   _   _   ____    _        _____   _   _   _____    _____   _   _    _____   _   _   _ 
	  | | | | | | |  _ \  | |      |_   _| | \ | | |  __ \  |_   _| | \ | |  / ____| | | | | | |
	  | | | | | | | |_) | | |        | |   |  \| | | |  | |   | |   |  \| | | |  __  | | | | | |
	  | | | | | | |  _ <  | |        | |   | . ` | | |  | |   | |   | . ` | | | |_ | | | | | | |
	  |_| |_| |_| | |_) | | |____   _| |_  | |\  | | |__| |  _| |_  | |\  | | |__| | |_| |_| |_|
	  (_) (_) (_) |____/  |______| |_____| |_| \_| |_____/  |_____| |_| \_|  \_____| (_) (_) (_)
	*/

	//________________________________________________________
	//

	if (DATA_BLINDED) {

	  // Apply analysis cuts and fill histograms
	  // These are all defined in [Name]_Analysis.cc (included from settings.h)
	  // You specify there also which cut is applied for each histo
	  // But all common baseline cuts are alreay applied above
	  if (!cmdline.noPlots  && pass_all_cuts) ana.fill_analysis_histos(data, syst.index, w);

	  // Save counts for the analysis cuts
	  for (auto cut : ana.analysis_cuts) if (pass_all_cuts)
	    if ( ( pass_all_cuts = cut.func() ) ) ofile->count(cut.name, w);

	  // Count Signal events
	  if ( pass_all_cuts && ana.signal_selection(data) ) ofile->count("Signal", w);

	} // end Blinding

      } // end JSON file cut

    } // End DATA
    else {
      // Background and Signal MCs

      // Loop and vary systematics
      for (syst.index = 0; syst.index <= (settings.varySystematics ? syst.nSyst : 0); ++syst.index) {
	w = 1;

	// Event weights
	// Lumi normalization
	// Signals are binned so we get the total weight separately for each bin
	if (cmdline.isSignal) {
	  int bin = vh_weightnorm_signal[signal_index]->FindBin(signal_index ? data.evt.SUSY_Stop_Mass : data.evt.SUSY_Gluino_Mass, data.evt.SUSY_LSP_Mass);
	  weightnorm = vh_weightnorm_signal[signal_index]->GetBinContent(bin);
	}
	// Normalize to chosen luminosity, also consider symmeteric up/down variation in lumi uncertainty
	w *= ana.get_syst_weight(data.evt.Gen_Weight*weightnorm, settings.lumiUncertainty, syst.nSigmaLumi[syst.index]);

	// Pileup reweighting
	if (syst.index == 0) h_nvtx->Fill(data.evt.NGoodVtx, w);
	if ( settings.doPileupReweighting ) {
	  w *= ana.get_pileup_weight(data.pu.NtrueInt, syst.nSigmaPU[syst.index]);
	  if (syst.index == 0) h_nvtx_rw->Fill(data.evt.NGoodVtx, w);
	}

	// Trigger efficiency scale factor
	w *= ana.get_syst_weight(settings.triggerEffScaleFactor, settings.triggerEffUncertainty, syst.nSigmaTrigger[syst.index]);

	// Jet Energy Scale uncertainty
	// Rescale jet 4-momenta
	ana.rescale_jets(data, syst.index, syst.nSigmaJEC[syst.index]);

	// Scale factors
	if (settings.applyWTagSF)
	  w *= ana.get_w_tagging_sf(data, syst.nSigmaWTagSF[syst.index]);

	if (settings.applyBTagSF)
	  w *= ana.get_b_tagging_sf(data, syst.nSigmaBTagSF[syst.index]);

	//if (settings.applyHadTopTagSF)
	//  w *= ana.get_top_tagging_sf(data, syst.nSigmaHadTopTagSF[syst.index]);

	// Scale QCD to match data in QCD dominated region
	if (TString(cmdline.dirname).Contains("QCD")) {
	  // Scale factor
	  // value obtained with ROOT macro: scripts/CalcQCDNormFactor.C
	  if (settings.scaleQCD)
	    w *= settings.useJSON ? 0.776458 : 0.785087; // Golden/Silver JSON

	  // HT reweighting
	  if (settings.doHTReweighting)
	    w *= ana.get_ht_weight(data, syst.nSigmaHT[syst.index]);
	}

	// Theory weights
	// LHE weight variations
	// More info about them here:
	// https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L165-L237

	// Alpha_s variations
	// A set of two weights
	// Only stored for NLO, otherwise vector size==0
	// If vector was not filled (LO samples), not doing any weighting
	if ( data.syst_alphas.Weights.size() == 2 )
	  w *= ana.get_alphas_weight(data.syst_alphas.Weights, syst.nSigmaAlphaS[syst.index], data.evt.LHA_PDF_ID);

	// Scale variations
	// A set of six weights, unphysical combinations excluded
	// If numScale=0 is specified, not doing any weighting
	if ( syst.numScale[syst.index] >= 1 && syst.numScale[syst.index] <= 3 )
	  w *= ana.get_scale_weight(data.syst_scale.Weights, syst.nSigmaScale[syst.index], syst.numScale[syst.index]);

	// PDF weights
	// A set of 100 weights for the nominal PDF
	// If numPdf=0 is specified, not doing any weighting
	if ( syst.numPdf[syst.index] >= 1 && syst.numPdf[syst.index] <= data.syst_pdf.Weights.size() )
	  w *= data.syst_pdf.Weights[syst.numPdf[syst.index]-1];
	else if ( syst.numPdf[syst.index] > data.syst_pdf.Weights.size() )
	  utils::error("numPdf (syst) specified is larger than the number of PDF weights in the ntuple");

	// Analysis specific weights (comes last, as things may depend on jet energy)
	w *= ana.get_analysis_weight(data);

	// Calculate variables that do not exist in the ntuple
	ana.calculate_common_variables(data, syst.index);
	ana.calculate_variables(data, syst.index);

	// Save counts (after each cuts)
	bool pass_all_cuts = true;
	if (syst.index == 0) ofile->count("NoCuts", w);

	// Cuts that are likely to be implemented in all analyses
	// eg. MET filters, baseline event selections etc.
	for (auto cut : ana.baseline_cuts) if (pass_all_cuts) {
	  pass_all_cuts = cut.func();
	  if (pass_all_cuts && syst.index==0) ofile->count(cut.name, w);
	}

	// If option (saveSkimmedNtuple) is specified save all 
	// skimmed events selected by the analysis to the output file
	// tree is copied and current weight is saved as "eventWeight"
	if ( settings.saveSkimmedNtuple && syst.index==0 ) if (ana.pass_skimming(data)) ofile->addEvent(w);

	// Apply analysis cuts and fill histograms
	// These are all defined in [Name]_Analysis.cc (included from settings.h)
	// You specify there also which cut is applied for each histo
	// But all common baseline cuts will be already applied above
	if (!cmdline.noPlots && pass_all_cuts) ana.fill_analysis_histos(data, syst.index, w);

	// Save counts after each analysis cut
	for (auto cut : ana.analysis_cuts) if (pass_all_cuts) {
	  pass_all_cuts = cut.func();
	  if (pass_all_cuts && syst.index==0) ofile->count(cut.name, w);
	}

	// Count remaining signal events
	if ( pass_all_cuts && syst.index==0 && ana.signal_selection(data) ) ofile->count("Signal", w);

      } // end systematics loop
    } // end Background/Signal MC

  } // end event loop

  stream.close();
  out_dir->cd();
  if (!cmdline.noPlots)
    ana.save_analysis_histos();
  ofile->close();
  return 0;
}
