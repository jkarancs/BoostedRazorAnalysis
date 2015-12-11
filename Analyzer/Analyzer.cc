//-----------------------------------------------------------------------------
// File:        Analyzer.cc
// Created:     24-Nov-2015
// Author:      Janos Karancsi
//-----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "common/utils.h" // Helper functions
#include "settings.h"     // Define all Analysis specific settings here

using namespace std;

int main(int argc, char** argv) {
  // ------------------------------
  // -- Parse command line stuff --
  // ------------------------------

  // Get file list and histogram filename from command line
  utils::commandLine cmdline;
  utils::decodeCommandLine(argc, argv, cmdline);

  itreestream stream(cmdline.filenames, settings.treeName);      
  if ( !stream.good() ) utils::error("unable to open ntuple file(s)");                         

  if ( cmdline.isdata ) cout << "Running on Data." << endl;
  else cout << "Running on MC." << endl;

  // Get number of events to be read
  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

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

  // --------------------------------------------------------------
  // -- Calculate the normalization factor for the event weights --
  // -- The original MC weight will be divided by this quantity  --
  // --------------------------------------------------------------
  
  cout << endl;
  double weightnorm = 1 ;
  if ( !cmdline.isdata ) {
    cout << "intLumi (settings): " << settings.intLumi << endl; // given in settings.h

    double xsec = ana.get_xsec_from_ntuple(cmdline.filenames, settings.treeName); // treename given in settings.h
    if ( xsec==0 ) return 1;
    cout << "xsec (ntuple): " << xsec << endl;

    double totweight = ana.get_totweight_from_ntuple(cmdline.filenames, settings.totWeightHistoName); // weight histo name given in settings.h
    cout << "totweight (ntuple): " << totweight << endl;

    weightnorm = (xsec*settings.intLumi)/totweight;
    cout << "weightnorm (calc): " << weightnorm << endl;
  }
  
  // ---------------------------------------
  // --- Pileup Reweighting              ---
  // ---------------------------------------

  // TString pileupname = "../data/pileup/pileup_weights.root"; // default
  // TFile* fpileup = TFile::Open(pileupname);
  // if ( !fpileup ) {
  //   cout << "Could not find pileup weights root file... Where did you put it??" << endl;
  //   return 1;
  // }
  // TH1D* h_pileup = (TH1D*)fpileup->Get("pileup_weight");

  if ( !cmdline.isdata && settings.doPileupReweighting ) {
    cout << "doPileupReweighting (settings): true" << endl;
  } else cout << "doPileupReweighting (settings): false" << endl;

  utils::outputFile* ofile;
  if ( settings.saveSkimmedNtuple ) {
    cout << "saveSkimmedNtuple (settings): true" << endl;
    ofile = new utils::outputFile(cmdline.outputfilename, stream);
  } else {
    ofile = new utils::outputFile(cmdline.outputfilename);
    cout << "saveSkimmedNtuple (settings): false" << endl;
  }

  // ---------------------------------------------------------------------------
  // -- Declare histograms                                                    --
  // ---------------------------------------------------------------------------

  TH1::SetDefaultSumw2();

  // This section is moved and defined in the Analysis class
  // this class specifies what histograms/cuts/methods etc. you want to implement
  // This is useful if someone wants to do quick study/define other search region etc.
  // But also, common methods in all anaylsis are defined in common/AnalysisBase.*

  ana.declare_histograms();

  // ------------------------------------------------------------------------------
  // -- Define the order of cuts (and corresponding bins in the counts histogram --
  // ------------------------------------------------------------------------------

  // Define cuts that are common in all analyses
  // Given in common/AnalysisBase.h
  ana.define_preselections(data);

  // Define cuts that specific to this analysis
  // Given in [Name]_Analysis.h specified in setting.h
  ana.define_selections(data);

  // Define bin order for counts histogram
  ofile->count("NoCuts",   0);
  for (auto cut : ana.baseline_cuts) ofile->count(cut.name, 0);
  for (auto cut : ana.analysis_cuts) ofile->count(cut.name, 0);
  cout << endl;
  cout << "Number of events counted after applying" << endl;
  cout << "- Baseline cuts (common for all analysis):" << endl;
  for (auto cut : ana.baseline_cuts) cout << "  "<<cut.name << endl;

  cout << endl;  
  cout << "- Analysis specific cuts:\n";
  for (auto cut : ana.analysis_cuts) cout << "  "<<cut.name << endl;

  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  cout << endl;
  cout << "Start looping on events ..." << endl;
  for(int entry=0; entry < nevents; ++entry) {
    // Read event into memory
    stream.read(entry);
    //data.CalculateAllVariables();

    if ( entry%100000==0 ) cout << entry << " events analyzed." << endl;

    // Correctly normalized (to luminosity) event weight
    double w = 1.;
    if ( data.evt.Gen_Weight!=0 ) w = data.evt.Gen_Weight*weightnorm;

    ofile->count("NoCuts", w);

    // Apply preselections and save counts
    bool pass_all = true;
    for (auto cut : ana.baseline_cuts) if (pass_all)
      if ( ( pass_all = cut.func() ) ) ofile->count(cut.name, w);
    
    // Apply analysis cuts and fill histograms
    // These are all defined in [Name]_Analysis.cc (included from settings.h)
    // You specify there also which cut is applied for each histo
    // But all common baseline cuts are alreay applied above
    if (pass_all) ana.fill_histograms(data);

    // Save counts for the analysis cuts
    for (auto cut : ana.analysis_cuts) if (pass_all)
      if ( ( pass_all = cut.func() ) ) ofile->count(cut.name, w);
    
    // If option (saveSkimmedNtuple) is specified
    // save all events selected by the analysis to the output file
    // tree is copied and current weight is saved as "eventWeight"
    if ( settings.saveSkimmedNtuple ) ofile->addEvent(w);

  } // end event loop

  //fpileup->Close();
  stream.close();
  ofile->close();
  return 0;
}
