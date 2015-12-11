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
    std::string progname;
    std::string outputfilename;
    std::vector<std::string> filenames;
    bool isdata;
  };
  
  // Read ntuple filenames from file list
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
    cl.progname = std::string(argv[0]);

    if ( argc < 3 ) error("<output filename> and <input file list> was not given!");

    // 1st command line argument
    cl.outputfilename = std::string(argv[1]);
    std::cout<<"output file is: "<<cl.outputfilename<<std::endl;

    cl.isdata = false;

    int n_data_arg = 0;
    // if 2nd argument is txt file
    if (std::string(argv[2]).find(".txt")!=std::string::npos&&argc==3) {
      cl.filenames = getFilenames(std::string(argv[2]));
      if (std::string(argv[2]).find("filelists/data")!=std::string::npos) n_data_arg++;
    }
    // otherwise add all of the rest arguments to the list
    else {
      cl.filenames = std::vector<std::string>();
      for (int i=2; i<argc; ++i) {
        if (std::string(argv[i]).find(".root")!=std::string::npos) {
  	cl.filenames.push_back(std::string(argv[i]));
  	if (std::string(argv[i]).find("ns_2015")!=std::string::npos) n_data_arg++; // This part is in all data file names
        }
        else error(std::string("argument ")+std::string(argv[i])+" is not a root file!");
      }
    }
    if (n_data_arg != 0 && n_data_arg != argc-2) error(" Data mixed with MC!");
    else cl.isdata = (n_data_arg == argc-2);
  }
}

#endif
