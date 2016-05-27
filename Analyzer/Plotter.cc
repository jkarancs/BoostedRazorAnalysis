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
#include "settings_Janos.h" // Define all Analysis specific settings here

int main(int argc, char** argv) {
  std::vector<std::string> vname_data = {};
  std::vector<std::string> vname_signal = {};

  // Get file list and histogram filename from command line
  utils::commandLine cmdline;
  utils::decodeCommandLine(argc, argv, cmdline, vname_data, vname_signal);

  // Read systematics file (only need number of lines)
  struct Systematics {
    unsigned int index = 0;
    unsigned int nSyst = 0;
  } syst;
  if (settings.varySystematics) {
    std::ifstream systFile(settings.systematicsFileName.c_str());
    if ( !systFile.good() ) utils::error("unable to open systematics file: " + settings.systematicsFileName);
    std::string line;
    while ( std::getline(systFile, line) ) ++syst.nSyst;
  }

  utils::outputFile* ofile;
  ofile = new utils::outputFile(cmdline.outputFileName);
  TDirectory* out_dir = gDirectory;

  DataStruct data;
  double w = 1;
  Analysis ana;

  ana.define_histo_options(w, data, syst.nSyst, syst.index, cmdline.dirname, settings.runOnSkim);
  ana.init_common_histos();
  ana.init_analysis_histos(syst.nSyst, syst.index);

  for (auto in_file : cmdline.fileNames) ana.load_analysis_histos(in_file);

  out_dir->cd();
  ana.save_analysis_histos(1);
  ofile->close();
  
  return 0;
}
