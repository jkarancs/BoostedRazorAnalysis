#include <iostream>
#include <functional>
#include <map>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

#include "DataStruct.h"

class AnalysisBase
{
public:
  AnalysisBase();
  ~AnalysisBase();

  double get_xsec_from_ntuple(const std::vector<std::string>&, const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const std::string&);

  void define_preselections(const DataStruct&);

  typedef struct Cut { std::string name; std::function<bool()> func; } Cut;
  std::vector<Cut> baseline_cuts;
};


//_______________________________________________________
//                       Constructor
AnalysisBase::AnalysisBase() { }


//_______________________________________________________
//                       Destructor
AnalysisBase::~AnalysisBase() { }


//_______________________________________________________
//           Read cross-section from ntuple
double
AnalysisBase::get_xsec_from_ntuple(const std::vector<std::string>& filenames, const std::string& treename)
{
  float evt_XSec=0, prev_XSec=0;
  for (auto filename : filenames) {
    TTree* tree = (TTree*)TFile::Open(filename.c_str())->Get(treename.c_str());
    tree->GetBranch("evt_XSec")->SetAddress(&evt_XSec);
    tree->GetEntry(0);
    if (prev_XSec!=0&&prev_XSec!=evt_XSec) {
      cout << "!! Error !! Analysis - Files added with different cross-sections. Please, add them separately!" << endl;
      return 0;
    }
    prev_XSec = evt_XSec;
  }
  
  return evt_XSec;
}

  
//_______________________________________________________
//          Read total weight from ntuple histos
double
AnalysisBase::get_totweight_from_ntuple(const std::vector<std::string>& filenames, const std::string& histoname)
{
  // Merging totweight histos
  TH1D *h = new TH1D("totweight","",1,0,1); // 1st bin: negative weighted events, 2nd bin: positive
  for (auto filename : filenames)
    h->Add((TH1D*)TFile::Open(filename.c_str())->Get(histoname.c_str()));
  double totweight = h->GetBinContent(1);
  delete h;
  return totweight;
}

//_______________________________________________________
//                 Define baseline cuts
void
AnalysisBase::define_preselections(const DataStruct& data)
{ 
  // Apply the same cuts as it is in the ntuple - Only for check
  // cut is an std::function, which we can define easily with a lambda function
  baseline_cuts.push_back({ .name="ntuple_filter", .func = [&data]() { 
			      // Define cut function here:
			      if ( !(data.jetsAK8.size>=2) ) return 0;
			      if ( !(data.jetsAK8.Pt[0]>350) ) return 0;
			      if ( !(data.jetsAK8.Pt[1]>350) ) return 0;
			      return 1;
			    } });
  
  // Recommended event filters by MET group
  // In some cases we should use txt files
  //
  // Select at least one good vertex (z<24, rho<2, ndof>=4)
  // NGoodVtx defined in:
  // https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L272-L275
  baseline_cuts.push_back({ .name="met_filter_NVtx",            .func = [&data](){ return data.evt.NGoodVtx>0; } });
  // baseline_cuts.push_back({ .name="NGoodVtx",         .func = [&data](){ return data.evt.Flag_goodVertices; } }); // Did not work correctly in MiniAOD

  // Other filters (From MiniAOD)
  baseline_cuts.push_back({ .name="met_filter_HBHE_Iso",        .func = [&data](){ return data.evt.Flag_HBHEIsoNoiseFilterResult; } });
  baseline_cuts.push_back({ .name="met_filter_HBHE_Run2_Loose", .func = [&data](){ return data.evt.Flag_HBHENoiseFilterResultRun2Loose; } });
  baseline_cuts.push_back({ .name="met_filter_EE_Bad_Sc",       .func = [&data](){ return data.evt.Flag_eeBadScFilter; } });
  baseline_cuts.push_back({ .name="met_filter_CSC_Halo_Tight",  .func = [&data](){ return data.evt.Flag_CSCTightHaloFilter; } });  
}
