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

//#define __TEST__

#include "common/Analyzer_cmd.h"

using namespace std;

#ifndef __TEST__
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
  DataStruct data;
  selectVariables(stream, data);

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

  // ---------------------------------------
  // --- Get Lumi weighting information  ---
  // ---------------------------------------

  string fsample = "";
  if ( argc > 6 ) fsample = string(argv[6]);

  double lumi = cmdline.lumi;
  cout << "lumi: " << lumi << endl;
  double xsect = cmdline.xsect;
  if (xsect>0) cout << "xsect: " << xsect << endl;
  else if (fsample.find("Bkg")==0) {
    // Get cross section value from the ntuple
    float evt_XSec=0, prev_XSec=0;
    for (auto filename : filenames) {
      TTree* tree = (TTree*)TFile::Open(filename.c_str())->Get("B2GTTreeMaker/B2GTree");
      tree->GetBranch("evt_XSec")->SetAddress(&evt_XSec);
      tree->GetEntry(0);
      if (prev_XSec!=0&&prev_XSec!=evt_XSec) {
	cout << "!! Error !! Analyzer - Files added with different cross-sections. Add them separately!" << endl;
	return 0;
      }
      prev_XSec = evt_XSec;
    }
    xsect = evt_XSec;
    cout << "xsect (taken from ntuple): " << xsect << endl;
  }
  double totweight = cmdline.totweight;
  if (totweight>0) cout << "totweight: " << totweight << endl;
  else if (fsample.find("Bkg")==0) {
    // Get total number of event histo
    // (counting negative and positive weights)
    TH1D *h = new TH1D("NEventNoFilter","",2,-1,1); // 1st bin: negative weighted events, 2nd bin: positive
    float evt_Gen_Weight;
    for (auto filename : filenames) {
      TFile *f = TFile::Open(filename.c_str());
      h->Add((TH1D*)f->Get("EventCounter/NEventNoFilter"));
      // Also Get GenInfoProduct weight (+- constant)
      TTree* tree = (TTree*)f->Get("B2GTTreeMaker/B2GTree");
      tree->GetBranch("evt_Gen_Weight")->SetAddress(&evt_Gen_Weight);
      tree->GetEntry(0);
    }
    // total weight = ( N_pos - N_min ) * | gen_weight |
    totweight = ( h->GetBinContent(2) - h->GetBinContent(1) ) * fabs(evt_Gen_Weight);
    delete h;
    cout << "totweight (taken from ntuple): " << totweight << endl;
  }

  // --------------------------------------------------------------
  // -- Calculate the normalization factor for the event weights --
  // -- The original MC weight will be divided by this quantity  --
  // --------------------------------------------------------------
  double weightnorm = 1.;
  if (xsect > 0 && totweight > 0 && lumi > 0) { // i.e. is not data
    weightnorm = (xsect*lumi)/totweight;
  }
  cout << "weightnorm: " << weightnorm << endl;

  // ---------------------------------------
  // --- Pileup Reweighting              ---
  // ---------------------------------------

  // TString pileupname = "../data/pileup/pileup_weights.root"; // default
  // TFile* fpileup = TFile::Open(pileupname);
  // if (!fpileup){
  //   cout << "Could not find pileup weights root file... Where did you put it??" << endl;
  //   return 1;
  // }
  // TH1D* h_pileup = (TH1D*)fpileup->Get("pileup_weight");

  string Pileup = "";
  if ( argc > 7 ) Pileup = string(argv[7]);

  //bool doPileupReweighting = false;
  if (fsample.find("Data")!=0 && Pileup == "Pileup_True"){
    //doPileupReweighting = true;
    cout << "Will do pileup reweighting" << endl;
  }

  outputFile ofile(cmdline.outputfilename, stream);

  // ---------------------------------------------------------------------------
  // -- Declare histograms                                                    --
  // ---------------------------------------------------------------------------

  TH1::SetDefaultSumw2();

  // Save the total weight for cross check purposes
  TH1D* h_totweight = new TH1D("h_totweight", "h_totweight", 1, 1, 2);

  // ------------------------------------------------------
  // -- Define the order of bins in the counts histogram --
  // ------------------------------------------------------

  ofile.count("NoCuts",   0);
  ofile.count("Ntuple", 0);
  ofile.count("Cleaning", 0);
  ofile.count("Jet1Pt", 0);
  ofile.count("Jet2Pt", 0);
  ofile.count("Jet1MassSD", 0);
  ofile.count("Jet2MassSD", 0);
  ofile.count("Trigger", 0);

  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry < nevents; ++entry) {
    // Read event into memory
    stream.read(entry);
    //data.CalculateAllVariables();

    if (entry%100000==0) cout<<entry<<endl;

    // Count events and get the total weight contibuted by the event
    h_totweight->Fill(1, data.evt.Gen_Weight);

    // Event weight
    double w = 1.;
    if (data.evt.Gen_Weight != 0) {
      w = data.evt.Gen_Weight*weightnorm;
    }

    ofile.count("NoCuts", w);

    if (!(data.jetsAK8.size>=2)) continue;
    if (!(data.jetsAK8.Pt[0]>350)) continue;
    if (!(data.jetsAK8.Pt[1]>350)) continue;
    ofile.count("Ntuple", w);

    // -------------------------------------------------------------------------
    // -- Get rid of the noise before start filling ANY histogram             --
    // -- by applying the filters:                                            --
    // -------------------------------------------------------------------------

    data.evt.NGoodVtx = 0;
    for (int iVtx=0; iVtx<data.evt.vtx_size; ++iVtx) {
      if (data.evt.vtx_ndof[iVtx]<4) continue;
      if (fabs(data.evt.vtx_z[iVtx])>24) continue;
      if (fabs(data.evt.vtx_rho[iVtx])>2) continue;
      ++data.evt.NGoodVtx;
    }

    if (!(data.evt.NGoodVtx>0)) continue;                          // data.evt.Flag_goodVertices - Doesn't work currently
    if (!(data.evt.Flag_HBHEIsoNoiseFilterResult)) continue;
    if (!(data.evt.Flag_HBHENoiseFilterResultRun2Loose)) continue;
    if (!(data.evt.Flag_CSCTightHaloFilter)) continue;             // Will need to rerun 2015 filter for Dec Jamboree
    if (!(data.evt.Flag_eeBadScFilter)) continue;
    ofile.count("Cleaning", w);

    if (!(data.jetsAK8.Pt[0]>400)) continue;
    ofile.count("Jet1Pt", w);

    if (!(data.jetsAK8.Pt[1]>400)) continue;
    ofile.count("Jet2Pt", w);

    if (!(data.jetsAK8.softDropMass[0]>110)) continue;
    if (!(data.jetsAK8.softDropMass[0]<210)) continue;
    ofile.count("Jet1MassSD", w);

    if (!(data.jetsAK8.softDropMass[1]>110)) continue;
    if (!(data.jetsAK8.softDropMass[1]<210)) continue;
    ofile.count("Jet2MassSD", w);

    if (!(data.evt.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50==1))
    ofile.count("Trigger", w);

    // Skim events into output file
    ofile.addEvent(w);

  } // end event loop

  //fpileup->Close();
  stream.close();
  ofile.close();
  return 0;
}
#endif
