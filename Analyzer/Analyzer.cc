//-----------------------------------------------------------------------------
// File:        Analyzer.cc
// Created:     24-Nov-2015
// Author:      Janos Karancsi
//-----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <vector>
//#include "../interface/Samples.h"
//#include "../interface/SmartHistos.h"
//#include "../plugins/B2GTreeReader.cc"
//#include "../plugins/B2GTreeLooper.cc"

#include "Analyzer.h"

using namespace std;

int main(int argc, char** argv) {
  // ------------------------------                                                         
  // -- Parse command line stuff --                                                         
  // ------------------------------                                                         
  
  // Get file list and histogram filename from command line
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);
  
  // Get names of ntuple files to be processed and open chain of ntuples
  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "B2GTTreeMaker/B2GTree");          
  if ( !stream.good() ) error("unable to open ntuple file(s)");                             
  
  // Get number of events to be read
  int nevents = stream.size();      
  cout << "Number of events: " << nevents << endl;
  
  // Select variables to be read
  selectVariables(stream);
  
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
  
  double xsect = cmdline.xsect;
  double totweight = cmdline.totweight; 
  double lumi = cmdline.lumi;
  
  string sample = "";
  if ( argc > 6 ) sample = string(argv[6]);
  string Pileup = "";
  if ( argc > 7 ) Pileup = string(argv[7]);
  
  //bool doPileupReweighting = false;
  if (sample != "Data" && Pileup == "Pileup_True"){
    //doPileupReweighting = true;
    cout << "Will do pileup reweighting" << endl;
  }
  
  // ---------------------------------------
  // --- Get the correct pileup histogram --
  // ---------------------------------------
  
  // TString pileupname = "../data/pileup/pileup_weights.root"; // default
  // TFile* fpileup = TFile::Open(pileupname);
  // if (!fpileup){
  //   cout << "Could not find pileup weights root file... Where did you put it??" << endl;
  //   return 1;
  // }
  // TH1D* h_pileup = (TH1D*)fpileup->Get("pileup_weight");
  
  // --------------------------------------------------------------
  // -- Calculate the normalization factor for the event weights --
  // -- The original MC weight will be divided by this quantity  --
  // --------------------------------------------------------------
  
  double weightnorm = 1.;
  if (xsect != -1 && totweight != -1 && lumi != -1) { // i.e. is not data
    weightnorm = (xsect*lumi)/totweight;
  }
  
  cout << "lumi: " << lumi << endl;
  cout << "xsect: " << xsect << endl;
  cout << "totweight: " << totweight << endl;
  cout << "weightnorm: " << weightnorm << endl;
  
  double w = 1; // Event weight
  
  outputFile ofile(cmdline.outputfilename);
  
  // ------------------------------------------------------
  // -- Define the order of bins in the counts histogram --
  // ------------------------------------------------------
  
  ofile.count("NoCuts",   0);
  ofile.count("Cleaning", 0);
  
  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------
  
  for(int entry=0; entry < nevents; ++entry) {
    // Read event into memory
    stream.read(entry);
    data.CalculateAllVariables();
    
    ofile.count("NoCuts", w);
    
    // -------------------------------------------------------------------------
    // -- Get rid of the noise before start filling ANY histogram             --
    // -- by applying the filters:                                            --
    // -------------------------------------------------------------------------
    
    if (!(data.evt.NGoodVtx>0)) continue;                          // data.evt.Flag_goodVertices - Doesn't work currently
    if (!(data.evt.Flag_HBHEIsoNoiseFilterResult)) continue;
    if (!(data.evt.Flag_HBHENoiseFilterResultRun2Loose)) continue;
    if (!(data.evt.Flag_CSCTightHaloFilter)) continue;             // Will need to rerun 2015 filter for Dec Jamboree
    if (!(data.evt.Flag_eeBadScFilter)) continue;
    
    ofile.count("Cleaning", w);
    
  } // end event loop
  
  //fpileup->Close();
  stream.close();
  ofile.close();
  return 0;
}
