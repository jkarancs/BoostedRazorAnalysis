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

  typedef struct Cut { std::string name; std::function<bool()> func; } Cut;
  std::vector<Cut> baseline_cuts;

  // Functions used by the Analyzer
  void define_preselections(const DataStruct&);

  void declare_common_histos();

  double get_xsec_from_ntuple(const std::vector<std::string>&, const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const std::string&);

  void init_pileup_reweightin(const std::string&, const std::string&, const std::vector<std::string>&);

  double get_pileup_weight(const int&, const bool&, const double&);

  double get_alphas_weight(const std::vector<float>&, const double&, const int&);

  double get_scale_weight(const std::vector<float>&, const double&, const unsigned int&);

private:
  double get_syst_weight_(const double&, const double&, const double&, const double&);

};


//_______________________________________________________
//                       Constructor
AnalysisBase::AnalysisBase() { }


//_______________________________________________________
//                       Destructor
AnalysisBase::~AnalysisBase() { }


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


//_______________________________________________________
//                 List of Histograms

TH1D* h_totweight;
TH1D* h_pileup_data;
TH1D* h_pileup_mc;
TH1D* h_pileup_weight;
TH1D* h_pileup_weight_down;
TH1D* h_pileup_weight_up;
TH1D* h_nvtx;
TH1D* h_nvtx_rw;

//_______________________________________________________
//              Define Histograms here
void
AnalysisBase::declare_common_histos()
{
  // total weight
  h_totweight     = new TH1D("totweight",     "MC;;Total (generator) event weight", 1,0,1);
  // pileup
  h_pileup_data        = new TH1D("pileup_data",        "Pile-up distribution - Data (Nominal);Pile-up", 100,0,100);
  h_pileup_mc          = new TH1D("pileup_mc",          "Pile-up distribution - MC;Pile-up",   100,0,100);
  h_pileup_weight      = new TH1D("pileup_weight",      "Pile-up weights - Nominal MB X-sec (69 mb);Pile-up;Weight",    100,0,100);
  h_pileup_weight_down = new TH1D("pileup_weight_down", "Pile-up weights - MB X-sec up 5% (72.45 mb);Pile-up;Weight",   100,0,100);
  h_pileup_weight_up   = new TH1D("pileup_weight_up",   "Pile-up weights - MB X-sec down 5% (65.55 mb);Pile-up;Weight", 100,0,100);
  h_nvtx          = new TH1D("nvtx",   "Number of vertices - Nominal;N_{Vertices}",    100,0,100);
  h_nvtx_rw       = new TH1D("nvtx_rw","Number of vertices - Pile-up reweighted (MC only);N_{Vertices}", 100,0,100);
}


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
  for (auto filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    h_totweight->Add((TH1D*)f->Get(histoname.c_str()));
    f->Close();
  }
  return h_totweight->GetBinContent(1);
}


//_______________________________________________________
//             Load pile-up reweighting infos
void
AnalysisBase::init_pileup_reweightin(const std::string& pileupDir, const std::string& mcPileupHistoName, const std::vector<std::string>& filenames)
{
  // Get data histogram (generated by pileupCalc.py script)
  TFile* f_pileup_data = TFile::Open((pileupDir+"data_pileup.root").c_str());
  h_pileup_data->Add((TH1D*)f_pileup_data->Get("pileup"));
  f_pileup_data->Close();
  // Get mc histogram saved inside the ntuple (unfiltered pileup distribution)
  for (auto filename : filenames) {
    TFile* f_pileup_mc = TFile::Open(filename.c_str());
    h_pileup_mc->Add((TH1D*)f_pileup_mc->Get(mcPileupHistoName.c_str()));
    f_pileup_mc->Close();
  }
  // Divide normalized data histo by normalized mc histo to get pileup weights for each bin
  h_pileup_weight->Divide(h_pileup_data, h_pileup_mc, 1/h_pileup_data->Integral(), 1/h_pileup_mc->Integral());
  //for (int bin=1; bin<=100; ++bin) cout << bin << " " << h_pileup_weight->GetBinContent(bin) << endl;
  // Also get systematic weights
  TFile* f_pileup_data_down = TFile::Open((pileupDir+"data_pileup_down.root").c_str());
  TH1D* h_pileup_data_down = (TH1D*)f_pileup_data_down->Get("pileup");
  h_pileup_weight_down->Divide(h_pileup_data_down, h_pileup_mc, 1/h_pileup_data_down->Integral(), 1/h_pileup_mc->Integral());    
  f_pileup_data_down->Close();
  TFile* f_pileup_data_up = TFile::Open((pileupDir+"data_pileup_up.root").c_str());
  TH1D* h_pileup_data_up = (TH1D*)f_pileup_data_up->Get("pileup");
  h_pileup_weight_up->Divide(h_pileup_data_up, h_pileup_mc, 1/h_pileup_data_up->Integral(), 1/h_pileup_mc->Integral());    
  f_pileup_data_up->Close();
}


//_______________________________________________________
//        private function to get scaled weight
double
AnalysisBase::get_syst_weight_(const double& weight_nominal, const double& weight_up, const double& weight_down, const double& nSigma)
{
  // Compute the weight according to the systematic variation considered
  // Use difference between nominal and up/down as 1 sigma variation 
  double dw_up = weight_up - weight_nominal;
  double dw_down = weight_nominal - weight_down;
  double w = weight_nominal;
  if (nSigma >= 0.) {
    w += nSigma*dw_up; 
  } else {
    w += nSigma*dw_down;
  }
  return w; 
}


//_______________________________________________________
//                  Get pile-up weight
double
AnalysisBase::get_pileup_weight(const int& NtrueInt, const bool& doSystematics, const double& nSigmaPU)
{
  int pu_bin = NtrueInt+1; // eg. pileup 0, is filled in bin 1
  double w_pileup = h_pileup_weight->GetBinContent(pu_bin);
  if ( doSystematics ) {
    double w_pileup_up = h_pileup_weight_up->GetBinContent(pu_bin);
    double w_pileup_down = h_pileup_weight_down->GetBinContent(pu_bin);
    w_pileup = get_syst_weight_(w_pileup, w_pileup_up, w_pileup_down, nSigmaPU);
  }
  return w_pileup;
}


//_______________________________________________________
//                  Get alpha_s weight
double
AnalysisBase::get_alphas_weight(const std::vector<float>& alphas_Weights, const double& nSigmaAlphaS, const int& LHA_PDF_ID)
{
  // A set of two weights corresponding to 
  // Powheg:  alpha_s = 0.118 -+ 0.002 
  // aMC@NLO: alpha_s = 0.118 -+ 0.001
  // Recommendation is to use +- 0.0015 --> rescale difference by 0.75 or 1.5
  // Treat weight as usual, gaussian, rescale to desired nSigma
  double w_alphas = 1;
  double w_alphas_up   = alphas_Weights[1];
  double w_alphas_down = alphas_Weights[0];
  double nSigma_0_0015 = nSigmaAlphaS;
  if (LHA_PDF_ID==260000||LHA_PDF_ID==260400) {
    // Powheg samples have -+ 0.001
    nSigma_0_0015 *= 1.5;
  } else {
    // aMC@NLO samples have -+ 0.002
    nSigma_0_0015 *= 0.75;
  }
  w_alphas = get_syst_weight_(w_alphas, w_alphas_up, w_alphas_down, 0.75*nSigmaAlphaS);
  return w_alphas;
}


//_______________________________________________________
//                  Get scale weight
double
AnalysisBase::get_scale_weight(const std::vector<float>& scale_Weights, const double& nSigmaScale, const unsigned int& numScale)
{
  /*
    Typical LHE run info:
    <weightgroup combine="envelope" type="Central scale variation">
      <weight id="1"> mur=1 muf=1 </weight>
      <weight id="2"> mur=1 muf=2 </weight>     --> save [0]
      <weight id="3"> mur=1 muf=0.5 </weight>   --> save [1]
      <weight id="4"> mur=2 muf=1 </weight>     --> save [2]
      <weight id="5"> mur=2 muf=2 </weight>     --> save [3]
      <weight id="6"> mur=2 muf=0.5 </weight>
      <weight id="7"> mur=0.5 muf=1 </weight>   --> save [4]
      <weight id="8"> mur=0.5 muf=2 </weight>
      <weight id="9"> mur=0.5 muf=0.5 </weight> --> save [5]
    </weightgroup>

    https://github.com/jkarancs/B2GTTrees/blob/master/plugins/B2GEdmExtraVarProducer.cc#L195-L202
    We save only ids: 2,3,4,5,7,9 (in this order)

    The idea here is to randomly choose to vary mu_f or mu_r or both simulataneously
    and rescale weight difference the usual way by desired nSigma
  */
  double w_scale = 1;
  double w_scale_up = 1;
  double w_scale_down = 1;
  if (numScale==1) {
    // fix mu_r = 1.0, vary mu_f = 2.0, 0.5
    w_scale_up   = scale_Weights[0];
    w_scale_down = scale_Weights[1];
  } else if (numScale==2) {
    // fix mu_f = 1.0, vary mu_r = 2.0, 0.5
    w_scale_up   = scale_Weights[2];
    w_scale_down = scale_Weights[4];
  } else if (numScale==3) {
    // vary simulataneously mu_r = mu_f = 2.0, 0.5
    w_scale_up   = scale_Weights[3];
    w_scale_down = scale_Weights[5];
  }
  w_scale = get_syst_weight_(w_scale, w_scale_up, w_scale_down, nSigmaScale);
  return w_scale;
}
