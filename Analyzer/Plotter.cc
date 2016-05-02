//-----------------------------------------------------------------------------
// File:        Analyzer.cc
// Created:     18-Jan-2016
// Author:      Janos Karancsi
//-----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <vector>

#include "common/utils.h" // Helper functions
#include "common/DataStruct.h"
#include "FullHad_Analysis.h"

int main(int argc, char** argv) {
  std::vector<std::string> vname_data = {};
  std::vector<std::string> vname_signal = {};

  // Get file list and histogram filename from command line
  utils::commandLine cmdline;
  utils::decodeCommandLine(argc, argv, cmdline, vname_data, vname_signal);

  utils::outputFile* ofile;
  ofile = new utils::outputFile(cmdline.outputFileName);
  TDirectory* out_dir = gDirectory;

  DataStruct data;
  double w = 1;
  Analysis ana;

  ana.define_histo_options(w, data, cmdline.dirname);
  ana.init_common_histos();
  ana.init_analysis_histos();

  for (auto in_file : cmdline.fileNames) ana.load_analysis_histos(in_file);

  out_dir->cd();
  ana.save_analysis_histos(1);
  ofile->close();
  
  return 0;
}
