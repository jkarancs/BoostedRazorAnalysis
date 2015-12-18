#ifndef ANALYZER_UTILS_H
#define ANALYZER_UTILS_H
//-----------------------------------------------------------------------------
// File:        Analyzer.h
// Description: Analyzer header for ntuples created by B2GTTreeMaker
// Created:     Tue Nov 24 2015
// Author:      Sezen Sekmen, Janos Karancsi
//-----------------------------------------------------------------------------
// -- System

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "treestream.h"
#include "treestream.cc"

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"

//-----------------------------------------------------------------------------
// -- Utilities
//-----------------------------------------------------------------------------
namespace utils {
  
  void
  error(std::string message)
  {
    std::cout << "** error ** " << message << std::endl;
    exit(0);
  }
  
  std::string 
  strip(std::string line)
  {
    int l = line.size();
    if ( l == 0 ) return std::string("");
    int n = 0;
    while (((line[n] == 0)    ||
  	  (line[n] == ' ' ) ||
  	  (line[n] == '\n') ||
  	  (line[n] == '\t')) && n < l) n++;
  
    int m = l-1;
    while (((line[m] == 0)    ||
  	  (line[m] == ' ')  ||
  	  (line[m] == '\n') ||
  	  (line[m] == '\t')) && m > 0) m--;
    return line.substr(n,m-n+1);
  }
  
  std::string
  nameonly(std::string filename)
  {
    int i = filename.rfind("/");
    int j = filename.rfind(".");
    if ( j < 0 ) j = filename.size();
    return filename.substr(i+1,j-i-1);
  }
  //-----------------------------------------------------------------------------
  struct outputFile
  {
    outputFile(std::string filename)
      : filename_(filename),
        file_(new TFile(filename_.c_str(), "recreate")),
        tree_(0),
        b_weight_(0),
        entry_(0),
        SAVECOUNT_(50000)
    {
      file_->cd();
      hist_ = new TH1D("counts", "", 1,0,1);
      hist_->SetBit(TH1::kCanRebin);
      hist_->SetStats(0);
    }
    
    outputFile(std::string filename, itreestream& stream, int savecount=50000) 
      : filename_(filename),
        file_(new TFile(filename.c_str(), "recreate")),
        tree_(stream.tree()->CloneTree(0)),
        b_weight_(tree_->Branch("eventWeight", &weight_, "eventWeight/D")),
        entry_(0),
        SAVECOUNT_(savecount)
    {
      std::cout << "events will be skimmed to file "
  	      << filename_ << std::endl;
      file_->cd();
      hist_ = new TH1D("counts", "", 1,0,1);
      hist_->SetBit(TH1::kCanRebin);
      hist_->SetStats(0);
    }
    
    void addEvent(double weight=1)
    {
      if ( tree_ == 0 ) return;
      
      weight_ = weight;
      file_   = tree_->GetCurrentFile();
      file_->cd();
      tree_->Fill();
      
      entry_++;
      if ( entry_ % SAVECOUNT_ == 0 )
        tree_->AutoSave("SaveSelf");
    }
    
    void count(std::string cond, double w=1.0)
    {
      hist_->Fill(cond.c_str(), w);
    }
    
    void close()
    {
      std::cout << "==> histograms saved to file " << filename_ << std::endl;
      if ( tree_ != 0 )
        {
  	std::cout << "==> events skimmed to file " << filename_ << std::endl;
  	file_ = tree_->GetCurrentFile();
        }
      file_->cd();
      //file_->Write("", TObject::kWriteDelete);
      file_->Write();
      file_->ls();
      file_->Close();
    }
    
    std::string filename_;  
    TFile* file_;
    TTree* tree_;
    TH1D*  hist_;
    TBranch* b_weight_;
    double     weight_;
    int    entry_;
    int    SAVECOUNT_;
  };
  
  struct commandLine
  {
    std::string progName;
    std::string outputFileName;         // first non option argument
    std::vector<std::string> fileNames; // second and rest non optional arguments
    bool isData;                        // determined automatically from input file names
    std::string systematicsFileName;    // needs to be given only if doSystetmatics settings true
    int numSyst;                        // needs to be given only if doSystetmatics settings true
  };
  
  // Read ntuple fileNames from file list
  std::vector<std::string>
  getFilenames(std::string filelist)
  {
    std::vector<std::string> v;
    if (filelist.find(".root")!=std::string::npos) {
      v.push_back(filelist);
    } else {
      std::ifstream stream(filelist.c_str());
      if ( !stream.good() ) error("unable to open file: " + filelist);
      
      // Get list of ntuple files to be processed
      
      std::string filename;
      while ( stream >> filename )
        if ( strip(filename) != "" ) v.push_back(filename);
    }
    return v;
  }
  
  void
  decodeCommandLine(int argc, char** argv, commandLine& cl)
  {
    cl.progName = std::string(argv[0]);

    if ( argc < 3 ) error("<output filename> and <input file list> was not given!");

    // decide whether input is data
    cl.isData = false;
    int n_data_arg = 0;
    int n_input = 0;

    // specify the nth line of the systematic file you want to consider
    // if settings.doSystematics is true this number needs to be > 0
    cl.systematicsFileName = "";
    cl.numSyst = 0;

    for (int iarg=1; iarg<argc; ++iarg) {
      std::string arg = argv[iarg];
      // look for optional arguments (argument has "=" in it)
      size_t f = arg.find("=");
      if (f!=std::string::npos) {
	std::string option=arg.substr(0, f);
	std::stringstream value;
	value<<arg.substr(f+1, arg.size()-f-1);
	// reading option
	if (option=="numSyst") value>>cl.numSyst;
	if (option=="systematicsFileName") value>>cl.systematicsFileName;
      } else {
	if (cl.outputFileName=="") {
	  // 1st non-optional (i.e xxx=yyy) command line argument is output ifle
	  cl.outputFileName = arg;
	  std::cout<<"output file is: "<<cl.outputFileName<<std::endl;
	} else {
	  // 2nd and rest non-optional argument is input file(s)
	  if (arg.find(".txt")!=std::string::npos) {
	    // if txt file, read it's contents
	    std::vector<std::string> list = getFilenames(arg);
	    cl.fileNames.insert(cl.fileNames.end(), list.begin(), list.end());
	    if (arg.find("filelists/data")!=std::string::npos) n_data_arg++;
	  } else {
	    // Otherwise add all root files to input file list
	    if (std::string(arg).find(".root")!=std::string::npos) cl.fileNames.push_back(arg);
	    else error(std::string("argument ")+arg+" is not a root file or a list (.txt file)!");
	    if (std::string(arg).find("ns_2015")!=std::string::npos) n_data_arg++; // This part is in all data file names
	  }
	  ++n_input;
	}
      }
    }

    // check number of input file arguments containing data-like strings
    if (n_data_arg != 0 && n_data_arg != n_input) error(" Data mixed with MC!");
    else cl.isData = (n_data_arg == n_input);
  }
}

#endif
