import ROOT
#File = ROOT.TFile.Open("/data/jkarancs/CMSSW/ntuple/B2GTTreeNtupleExtra_MC_25ns_80X_QCD.root")
#File = ROOT.TFile.Open("/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/Validation/CMSSW_8_0_20/src/B2GTTreeNtupleExtra_MC_25ns_80X_QCD.root")
File = ROOT.TFile.Open("/data/jkarancs/CMSSW/SusyAnalysis/Ntuples/Validation/CMSSW_8_0_24_patch1/src/B2GTTreeNtupleExtra_MC_80X.root")
tree = File.Get("B2GTTreeMaker/B2GTree")

def printvars( varname, vartype, prefix, isvector, keep_prefix ):
    if varname.startswith(prefix):
        short = varname if keep_prefix else varname[len(prefix):]
        if "Keys" in varname:
            print "    std::vector<std::vector<int> > "+short+";"
        elif isvector:
            if "i" in vartype:
                print "    std::vector<unsigned int> "+short+";"
            elif "I" in vartype:
                print "    std::vector<int> "+short+";"
            elif "l" in vartype:
                print "    std::vector<long> "+short+";"
            elif "F" in vartype:
                print "    std::vector<float> "+short+";"
            elif "D" in vartype:
                print "    std::vector<double> "+short+";"
        else:
            if "i" in vartype:
                print "    unsigned int "+short+";"
            elif "I" in vartype:
                print "    int "+short+";"
            elif "l" in vartype:
                print "    long "+short+";"
            elif "F" in vartype:
                print "    float "+short+";"
            elif "D" in vartype:
                print "    double "+short+";"
    return

def printinits( varname, vartype, prefix, isvector, keep_prefix ):
    if varname.startswith(prefix):
        short = varname if keep_prefix else varname[len(prefix):]
        if isvector or "Keys" in varname:
            print "      init_vec("+short+");"
        else:
            if "i" in vartype:
                print "      "+short+"=9999;"
            elif "l" in vartype:
                print "      "+short+"=9999;"
            elif "I" in vartype:
                print "      "+short+"=NOVAL_I;"
            elif "F" in vartype:
                print "      "+short+"=NOVAL_F;"
            elif "D" in vartype:
                print "      "+short+"=NOVAL_F;"
    return

def printclass( tree, classname, name, prefix_list, add_loop ):
    print "  class "+classname+" {"
    print "  public:"
    print "    "+classname+"() { init(); };"
    print "    "
    # print variables
    for prefix in prefix_list:
        for branch in tree.GetListOfBranches(): printvars( branch.GetName(), str(branch.GetTitle())[-1:], prefix, ("[" in branch.GetTitle()), (prefix_list[0] != prefix) )
    print "    "
    if add_loop:
        print "    unsigned int it;"
        print "    "
    print "    void init() {"
    if add_loop: print "      it = -1;"
    # initialize variables
    for prefix in prefix_list:
        for branch in tree.GetListOfBranches(): printinits( branch.GetName(), str(branch.GetTitle())[-1:], prefix, ("[" in branch.GetTitle()), (prefix_list[0] != prefix) )
    print "    }"
    if add_loop:
        print "    "
        print "    bool Loop() {"
        print "      ++it;"
        print "      if (it<size) {"
        print "        return 1;"
        print "      } else {"
        print "        it=-1;"
        print "        return 0;"
        print "      }"
        print "    }"
    print "    "
    print "  } "+name+";"
    print "  "
    return


# Now print the contents for common/DataStruct.h
print """#ifndef DataStruct_h
#define DataStruct_h

#define NOVAL_I -9999
#define NOVAL_F -9999.0

#include <vector>
#include <iostream>

inline void init_vec(std::vector<int>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=-9999; }
inline void init_vec(std::vector<float>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=-9999; }
inline void init_vec(std::vector<double>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=-9999; }
inline void init_vec(std::vector<std::vector<int> >& vec) { vec.resize(1); vec[0].resize(500); for (int i=0; i<500; ++i) vec[0][i]=-9999; }

class DataStruct {
public:
  DataStruct() {};
  ~DataStruct() {};
"""

printclass( tree, "EventData",           "evt",             ["evt_", "SUSY_"],   0 )
printclass( tree, "METData",             "met",             ["met_"],            0 )
printclass( tree, "PuppiMETData",        "puppimet",        ["puppimet_"],       0 )
printclass( tree, "PileupData",          "pu",              ["pu_"],             0 )
printclass( tree, "VertexData",          "vtx",             ["vtx_"],            0 )
printclass( tree, "SystScaleData",       "syst_scale",      ["scale_"],          1 )
printclass( tree, "SystPDFData",         "syst_pdf",        ["pdf_"],            1 )
printclass( tree, "SystAlphaSData",      "syst_alphas",     ["alphas_"],         1 )
printclass( tree, "SystMETUncData",      "syst_met",        ["metsyst_"],        1 )
printclass( tree, "SystPuppiMETUncData", "syst_puppimet",   ["puppimetsyst_"],   1 )
printclass( tree, "FilterData",          "filter",          ["Flag_"],           0 )
printclass( tree, "HLTData",             "hlt",             ["HLT_"],            0 )
printclass( tree, "GenVars",             "gen",             ["gen_"],            1 )
printclass( tree, "ElectronVars",        "ele",             ["el_"],             1 )
printclass( tree, "MuonVars",            "mu",              ["mu_"],             1 )
printclass( tree, "AK4JetVars",          "jetsAK4",         ["jetAK4CHS_"],      1 )
printclass( tree, "AK8JetVars",          "jetsAK8",         ["jetAK8CHS_"],      1 )
printclass( tree, "AK8SubjetVars",       "subjetsAK8",      ["subjetAK8CHS_"],   1 )
printclass( tree, "AK8GenJetVars",       "genjetsAK8",      ["genjetAK8SD_"],    1 )

print """
};

#endif
"""
