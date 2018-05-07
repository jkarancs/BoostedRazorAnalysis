#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
#include "TH3.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStopwatch.h"

class Postfixes {
public:
  Postfixes() {}
  ~Postfixes() {}
  typedef struct Postfix { std::function<size_t()> sel; std::vector<std::string> vec; std::vector<std::string> leg; std::string colz; } Postfix;
  
private:
  std::map<std::string, Postfix> pf_map_;
  
  double get_dbl_(std::string& str) {
    std::stringstream ss(str);
    double d; ss>>d;
    std::string rest; ss>>rest;
    str.erase(0, str.size()-rest.size());
    return d;
  }
  
  std::vector<double> str_to_vec_dbl_(std::string s) {
    std::vector<double> v_dbl;
    std::stringstream ss(s);
    std::string sub;
    while(std::getline(ss, sub, ',')) {
      if (sub.find("to")!=std::string::npos) {
	double i,j,k=1;
	i = get_dbl_(sub);
	sub.erase(0,2);
	j = get_dbl_(sub);
	if (sub.size()) {
	  sub.erase(0,1);
	  k = get_dbl_(sub);
	} else k=(i<j)?1:-1;
	if (k>=0) for (double ii=i; ii<=j; ii+=k) v_dbl.push_back(ii);
	else for (double ii=i; ii>=j; ii+=k) v_dbl.push_back(ii);
      } else v_dbl.push_back(get_dbl_(sub));
    }
    return v_dbl;
  }
  
  std::vector<std::string> interpret_substring_(std::string s, bool format=0) {
    std::vector<std::vector<std::string> > v_parts;
    std::stringstream ss(s+"[");
    std::string sub;
    while(std::getline(ss, sub, '[')) {
      size_t f_beg = sub.find("[");
      size_t f_end = sub.find("]");
      if (f_beg==std::string::npos&&f_end==std::string::npos) {
        std::vector<std::string> part;
        part.push_back(sub);
        if (v_parts.size()) part[part.size()-1].insert(0,"[");
        v_parts.push_back(part);
      } else if (f_end!=std::string::npos) {
        std::vector<std::string> part1;
        for (auto dbl : str_to_vec_dbl_(sub.substr(0,f_end))) {
	  std::stringstream ss2; ss2<<dbl;
	  part1.push_back(ss2.str());
        }
        v_parts.push_back(part1);
        std::vector<std::string> part2;
        part2.push_back(sub.substr(f_end+1,sub.size()-f_end-1));
        v_parts.push_back(part2);
      }
    }
    size_t n_max=0;
    for (std::vector<std::string>& part : v_parts) if (part.size()>n_max) n_max=part.size();
    std::vector<std::string> v_str;
    for (size_t i=0; i<n_max; ++i) v_str.push_back("");
    for (std::vector<std::string>& part : v_parts) for (size_t i=0; i<n_max; ++i) 
      v_str[i] += part[std::min(i,part.size()-1)];
    if (format) for (size_t i=0, n=v_str.size(); i<n; ++i) {
      size_t f = v_str[i].find(".");
      if (f!=std::string::npos) v_str[i].replace(f,1,"p");
    }
    return v_str;
  }
  
  std::vector<std::string> interpret_string_(std::string s, bool format=0) {
    std::vector<std::string> v_str;
    std::stringstream ss(s+";");
    std::string sub;
    while(std::getline(ss, sub, ';')) {
      std::vector<std::string> subs = interpret_substring_(sub,format);
      for (size_t i=0; i<subs.size(); ++i) v_str.push_back(subs[i]);
    }
    return v_str;
  }
  
public:
  void AddNew(std::string name, std::function<size_t()> sel, std::string pf, std::string leg, std::string colz) { 
    pf_map_.insert(std::pair<std::string, Postfix>(name, { .sel=sel, .vec=interpret_string_(pf,1), .leg=interpret_string_(leg), .colz=colz })); 
  }
  
  Postfix GetPostfix(std::string name) { 
    size_t count = pf_map_.count(name);
    if (!count) std::cout<<"!!! ERROR: Postfixes::GetPostfix: Postfix with name = "<<name<<" was not found."<<std::endl;
    return (count) ? pf_map_[name] : Postfix({0,std::vector<std::string>(),std::vector<std::string>(),""});
  }
};

class Cut {
public:
  Cut(std::function<bool()> func) { func_=func; Reset(); }
  ~Cut() {}
  
  void Reset() { cut_=-1; }
  bool Eval() { if (cut_==-1) cut_ = func_(); return cut_; }
  
private:
  int cut_;
  std::function<bool()> func_;
};

class Cuts {
public:
  Cuts() {}
  ~Cuts() {}
  
private:
  std::map<std::string, Cut*> cut_map_;
  
public:
  void AddNew(std::string name, std::function<bool()> cut) { cut_map_.insert(std::pair<std::string, Cut*>(name, new Cut(cut))); }
  
  Cut* GetCut(std::string name) {
    size_t count = cut_map_.count(name);
    if (!count) std::cout<<"!!! ERROR: Cuts::GetCut: Cut with name = "<<name<<" was not found."<<std::endl;
    return (count) ? cut_map_[name] : new Cut([](){ return 0; });
  }

  void ResetAllCut() { for (auto cut : cut_map_) cut.second->Reset(); }
};


class SmartHisto {
  
public:
  // constructors, destructor
  SmartHisto(std::string name, std::vector<std::string>& pf_names, std::vector<Postfixes::Postfix> pfs,
	     std::vector<std::function<double()> >& ffs, std::vector<std::function<double()> >& weights, std::vector<Cut*> & cuts, 
	     std::string draw, std::string opt, std::vector<double> ranges, std::vector<std::map<int, std::string> >& bin_labels,
	     std::vector<std::vector<std::string> >& spec, std::vector<std::vector<std::string> >& spec2) { 
    name_=name;
    pf_names_=pf_names;
    pfs_ = pfs;
    npf_=pfs_.size();
    ndim_=ffs.size();
    if (ndim_>0) fill_1d_ = ffs[0];
    if (ndim_>1) fill_2d_ = ffs[1];
    if (ndim_>2) fill_3d_ = ffs[2];
    nweight_=weights.size();
    if (nweight_>0) weight1_ = weights[0];
    if (nweight_>1) weight2_ = weights[1];
    if (nweight_>2) weight3_ = weights[2];
    ncut_=cuts.size();
    if (ncut_>0) cut1_ = cuts[0];
    if (ncut_>1) cut2_ = cuts[1];
    if (ncut_>2) cut3_ = cuts[2];
    if (ncut_>3) cut4_ = cuts[3];
    if (ncut_>4) cut5_ = cuts[4];
    if (ncut_>5) cut6_ = cuts[5];
    draw_=draw;
    norm_  = opt.find("Norm")!=std::string::npos;
    logx_  = opt.find("logX")!=std::string::npos;
    log_  = opt.find("Log")!=std::string::npos;
    keep_color_  = opt.find("KeepColor")!=std::string::npos;
    keep_order_  = opt.find("KeepOrder")!=std::string::npos;
    stat_  = opt.find("Stat")!=std::string::npos ? 1110 : 0;
    sumw2_ = opt.find("Sumw2")!=std::string::npos;
    stack_  = opt.find("Stack")!=std::string::npos;
    ratio_  = opt.find("AddRatio")!=std::string::npos;
    approval_ = opt.find("Approval")!=std::string::npos;
    months_ = opt.find("Months")!=std::string::npos;
    weeks_ = opt.find("Weeks")!=std::string::npos;
    dates_ = opt.find("Dates")!=std::string::npos;
    twocol_ = opt.find("TwoCol")!=std::string::npos;
    addint_ = opt.find("AddInt")!=std::string::npos;
    syst_ = name.find("Counts_vs_")!=std::string::npos;
    doublex_ = opt.find("DoubleX")!=std::string::npos;
    n_nostack_ = 0;
    if (stack_)    { std::stringstream ss; ss<<opt.substr(opt.find("Stack")   +5,1); ss>>n_nostack_; }
    if (twocol_)   { std::stringstream ss; ss<<opt.substr(opt.find("TwoCol")  +6,2); ss>>twocol_; }
    if (approval_) { std::stringstream ss; ss<<opt.substr(opt.find("Approval")+8,2); ss>>approval_; }
    ranges_=ranges;
    bin_labels_=bin_labels;
    if (npf_>5) std::cout<<"!!! ERROR: SmartHisto::constructor: Fixme! - More than 5 postfixes, only use max 4, or redefine functions!\n";
    if (ndim_>3) std::cout<<"!!! ERROR: SmartHisto::constructor: More than 3 dimension, define a maximum of 3!\n";
    if (ncut_>6) std::cout<<"!!! ERROR: SmartHisto::constructor: Fixme! - More than 5 cuts specified, please add new variables to store them!\n";
    if (nweight_>3) std::cout<<"!!! ERROR: SmartHisto::constructor: Fixme! - More than 3 weights specified, please add new variables to store them!\n";
    spec_ =spec;
    spec2_=spec2;
    init_();
  }
  ~SmartHisto() {}
  
private:
  std::string name_;
  std::vector<std::string> pf_names_;
  std::vector<Postfixes::Postfix> pfs_;
  size_t npf_;
  
  // Fill functions specify how to fill histos
  size_t ndim_;
  std::function<double()> fill_1d_;
  std::function<double()> fill_2d_;
  std::function<double()> fill_3d_;
  
  // Weights: pointers to floats
  size_t nweight_;
  std::function<double()> weight1_;
  std::function<double()> weight2_;
  std::function<double()> weight3_;
  
  // Cuts: pointers to bools
  size_t ncut_;
  Cut* cut1_;
  Cut* cut2_;
  Cut* cut3_;
  Cut* cut4_;
  Cut* cut5_;
  Cut* cut6_;
  
  // Draw arguments and plot options
  // options:
  // NORM - normalize
  // LOG - log scale
  // STAT - draw stat boxes)
  std::string draw_;
  bool norm_; // Normalize histo
  bool logx_; // X log scale
  bool log_; // Y/Z log scale
  bool keep_color_; // Keep colors/markers in order
  bool keep_order_; // Keep histos in original order on stacked plots
  Int_t stat_; // Draw Stat boxes
  bool sumw2_;
  bool stack_; // Create Stacked plot
  bool ratio_; // Add ratio - For stacked: below plot, otherwise Draw ratio on top
  int approval_; // Add sqrt(s) and CMS (Preliminary/Work in progess) text
  bool months_; // Use X-axis Time format: dd/mm
  bool weeks_;  // Use X-axis Time format: ww
  bool dates_;  // Use X-axis Time format: Jan Feb Mar ....
  int  twocol_;  // Split legend into two columns, number specify how to split,
  //                eg. 57 = first 5 entries in the first colum, 7 in the second
  bool addint_;  // Add plot integrals to a separate column in the legend
  size_t n_nostack_; // Do no stack first n plots
  bool plot_asymm_err_; // Decide automatically if histo should be plotted with asymmetric errors (Using TGraphAE)
  bool syst_;
  bool doublex_; // Double the X-size of the plot area (1000x500)
  
  // axis ranges: xlow, xhigh, ylow, yhigh, zlow, zhigh
  // if low==high -> do not set
  std::vector<double> ranges_;
  
  // Title for each bin for the 3 axis (pair<int, label>)
  std::vector<std::map<int, std::string> > bin_labels_;
  
  // histo containers
  TH1D* h1d_0p_;
  TH2D* h2d_0p_;
  TH3D* h3d_0p_;
  std::vector<TH1D*> h1d_1p_;
  std::vector<TH2D*> h2d_1p_;
  std::vector<TH3D*> h3d_1p_;
  std::vector<std::vector<TH1D*> > h1d_2p_;
  std::vector<std::vector<TH2D*> > h2d_2p_;
  std::vector<std::vector<TH3D*> > h3d_2p_;
  std::vector<std::vector<std::vector<TH1D*> > > h1d_3p_;
  std::vector<std::vector<std::vector<TH2D*> > > h2d_3p_;
  std::vector<std::vector<std::vector<TH3D*> > > h3d_3p_;
  std::vector<std::vector<std::vector<std::vector<TH1D*> > > > h1d_4p_;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > h2d_4p_;
  std::vector<std::vector<std::vector<std::vector<TH3D*> > > > h3d_4p_;
  std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*> > > > > h1d_5p_;
  std::vector<std::vector<std::vector<std::vector<std::vector<TH2D*> > > > > h2d_5p_;
  std::vector<std::vector<std::vector<std::vector<std::vector<TH3D*> > > > > h3d_5p_;
  
  // --------------------------------------------------------------------------
  //                      Special Histogram Calculations:
  
  // Define plots that are calculated from +1 dimensional objects
  // Eg.: MPV from cluster charge distribution
  std::map<TH1D*, TH2D*> mother_2d_;
  std::map<TH2D*, TH3D*> mother_3d_;
  std::vector<std::vector<std::string> > spec_;
  std::vector<std::vector<std::string> > spec2_;
  TF1* ring_fit_[3][4];
  void init_() {
    // Get Hit Eff vs DCol Eff functions (previously measured)
    //std::string fname[3] = {"hiteff_vs_dcol_l1", "hiteff_vs_dcol_l2", "hiteff_vs_dcol_l3"};
    //std::string ringname[4] = {"_ring1", "_ring2", "_ring3", "_ring4"};
    //TFile *f_func = TFile::Open("/data/jkarancs/CMSSW/TimingStudy/CMSSW_7_1_0_pre1/src/DPGAnalysis/PixelTimingStudy/test/DynIneff_scale_factors/HitEffvsDColFunctions.root");
    //if (f_func) {
    //  for (int lay=0; lay<3; ++lay) for (int ring=0; ring<4; ring++) {
    //    TF1* f = (TF1*)f_func->Get((fname[lay]+ringname[ring]).c_str());
    //    ring_fit_[lay][ring] = (TF1*)f->Clone();
    //  }
    //  f_func->Close();
    //} else std::cout<<"SmartHisto::init_() - File not found: "<<f_func->GetName()<<std::endl;
    
    // Initialize some variables
    h1d_0p_ = 0;
    h2d_0p_ = 0;
    h3d_0p_ = 0;
    
    plot_asymm_err_ = 0;
  }
  
  size_t find_spec_(std::string name) {
    size_t find = -1;
    for (size_t s=0; s<spec_.size(); ++s) if (name.find(spec_[s][0])!=std::string::npos) find = s;
    return find;
  }
  size_t find_spec2_(std::string name) {
    size_t find = -1;
    for (size_t s=0; s<spec2_.size(); ++s) 
      if ((s==0)? name.find(spec2_[s][0])==0 : name.find(spec2_[s][0])!=std::string::npos)
	find = s;
    return find;
  }
  
  // Define which functions to use for each special histo
  void calc_spec_1d_(TH1D* h1d, TH2D* h2d, bool savemother=1) {
    size_t s = find_spec2_(name_);
    if (s==0) calc_eff_1d_(h1d, h2d, 1, savemother); // Average (Use ProfileX)
    else if (s==1) calc_mpv_1d_(h1d, h2d, savemother); // Fit landau + gaus, extract maximum
    else if (name_.find("DColEfficiency")!=std::string::npos) calc_dcol_1d_(h1d, h2d, savemother);
    else if (name_.find("Counts")!=std::string::npos) calc_syst_1d_(h1d, h2d, savemother);
    else if (name_.find("SignalSignificance")!=std::string::npos) calc_soverb_1d_(h1d, h2d, savemother);
    else calc_eff_1d_(h1d, h2d, 1, savemother); // Efficiencies/Fractions/Rates
  }
  
  void calc_spec_2d_(TH2D* h2d, TH3D* h3d, bool savemother=0) {
    // Only 2D S/B needs special treatment due to background bin (eg. MLSP vs MGluino)
    if (name_.find("SignalSignificance")!=std::string::npos) calc_soverb_2d_(h2d, h3d);
    else {
      for (int i=1, ni=h3d->GetNbinsX(); i<=ni; ++i) {
	TH1D h1d_temp("h1d_temp","", h3d->GetNbinsY(), h3d->GetYaxis()->GetXmin(), h3d->GetYaxis()->GetXmax());
	TH2D h2d_temp("h2d_temp","", h3d->GetNbinsY(), h3d->GetYaxis()->GetXmin(), h3d->GetYaxis()->GetXmax(), 
		      h3d->GetNbinsZ(), h3d->GetZaxis()->GetXmin(), h3d->GetZaxis()->GetXmax());
	double slice_entries = 0;
	for (int j=1, nj=h3d->GetNbinsY(); j<=nj; ++j) for (int k=1, nk=h3d->GetNbinsZ(); k<=nk; ++k) {
	  slice_entries += h3d->GetBinContent(i, j, k);
	  h2d_temp.SetBinContent(j, k, h3d->GetBinContent(i, j, k));
	  h2d_temp.SetBinError  (j, k, h3d->GetBinError  (i, j, k));
	}
	h2d->SetEntries(slice_entries);
	calc_spec_1d_(&h1d_temp, &h2d_temp, savemother);
	for (int j=1, nj=h3d->GetNbinsY(); j<=nj; ++j) {
	  h2d->SetBinContent(i, j, h1d_temp.GetBinContent(j));
	  h2d->SetBinError  (i, j, h1d_temp.GetBinError  (j));
	}
      }
      h2d->SetEntries(h3d->GetEntries());
      if (savemother) mother_3d_[h2d]=h3d;
    }
  }
  //void calc_spec_2d_(TH2D* h2d, TH3D* h3d) {
  //  if (find_spec2_(name_)==0) calc_eff_2d_(h2d, h3d); // Average (Use 3DProfile)
  //  else if (name_.find("DColEfficiency")!=std::string::npos) /* calc_dcol_2d_(h2d, h3d) */;
  //  else calc_eff_2d_(h2d, h3d); // Efficiencies/Fractions/Rates
  //}
  
  // Define functions below: (eg: Efficiency, MPV etc...)
  // ********************  Efficiencies *********************
  void calc_eff_1d_(TH1D* h1d, TH2D* h2d, int err_type = 1, bool savemother=1) {
    if (err_type==1) { // Use TProfile - Normal Approximation interval
      TProfile* p = h2d->ProfileX();
      for (int i=1; i<=p->GetNbinsX(); ++i) {
	h1d->SetBinContent(i,p->GetBinContent(i));
	h1d->SetBinError(i,p->GetBinError(i));
      }
      delete p;
    } else if (err_type==2) { // Wilson Score Interval
      double z = 1; // N Sigma confidence
      for (int i=1; i<=h2d->GetNbinsX(); ++i) {
	double val = h2d->GetBinContent(i,2), mis = h2d->GetBinContent(i,1);
	if (val+mis>0) {
	  double eff = val / (val + mis);
	  double n = val + mis;
	  double cen = (eff+(z*z/(2*n))) / (1.0 + (z*z/n));
	  double halfwidth = z*sqrt( eff*(1.0-eff)/n + (z*z/(4*n*n)) ) / (1.0 + (z*z/n));
	  double err = halfwidth + fabs(cen-eff);
	  // Assymmetric error -> Choose larger for a conservative error estimate
	  h1d->SetBinContent(i,eff);
	  h1d->SetBinError(i,err);
	}
      }
    }
    h1d->SetEntries(h2d->GetEntries());
    if (savemother) mother_2d_[h1d]=h2d;
    if (h2d->GetNbinsY()==2&&h2d->GetYaxis()->GetBinCenter(1)==0&&h2d->GetYaxis()->GetBinCenter(2)==1) plot_asymm_err_=1;
  }
  
  //void calc_eff_2d_(TH2D* h2d, TH3D* h3d) {
  //  TProfile2D* p = h3d->Project3DProfile("yx");
  //  for (int i=1; i<=p->GetNbinsX(); ++i) for (int j=1; j<=p->GetNbinsY(); ++j) {
  //    h2d->SetBinContent(i,j,p->GetBinContent(i,j));
  //    h2d->SetBinError(i,j,p->GetBinError(i,j));
  //  }
  //  delete p;
  //  h2d->SetEntries(h3d->GetEntries());
  //  if (savemother) mother_3d_[h2d]=h3d;
  //}
  
  // ****************** Most Probable Value ******************
  void calc_mpv_1d_(TH1D* h1d, TH2D* h2d, bool savemother=1) {
    if (h2d->GetEntries()) {
      int nbiny = h2d->GetNbinsY();
      double ylow  = h2d->GetYaxis()->GetXmin();
      double yhigh = h2d->GetYaxis()->GetXmax();
      for (int i=1; i<=h2d->GetNbinsX(); i++) {
	// Get Y-slice in bin of X
	double entries = 0;
	TH1D *slice_x = new TH1D("slice_x","slice_x",nbiny,ylow,yhigh);
	for (int j=1; j<=nbiny; j++) {
	  double cont = h2d->GetBinContent(i,j);
	  entries += cont;
	  slice_x->SetBinContent(j,cont);
	}
	if (entries>100) {
	  double mpv, error;
	  fit_landau_plus_gaus_(slice_x, mpv, error);
	  h1d->SetBinContent(i,mpv);
	  h1d->SetBinError(i,error);
	}
	delete slice_x;
      }
    }
    h1d->SetEntries(h2d->GetEntries());
    if (savemother) mother_2d_[h1d]=h2d;
  }
  
  void fit_landau_plus_gaus_(TH1D* h, double &mpv, double &error) {
    double xmax = h->GetBinCenter(h->GetMaximumBin());
    double ymax = h->GetMaximum();
    // First fit with landau to get some parameter estimates
    TF1 *f = new TF1((std::string(h->GetName())+"_fit").c_str(),"landau", (xmax<25 ? 0 : xmax-25), xmax+25);
    f->SetParameter(0, 5.5*ymax);
    f->SetParLimits(0, 4*ymax,7*ymax);
    f->SetParameter(1, xmax);
    f->SetParLimits(1, xmax<5 ? 0: xmax-5, xmax + 5);
    f->SetParameter(2, 2.3);
    f->SetParLimits(2, 0.1, 5);
    h->Fit((std::string(h->GetName())+"_fit").c_str(),"RBQ0");
    //double Chi2NDoF = f->GetChisquare()/f->GetNDF();
    double Const = f->GetParameter(0);
    //double MPV = f->GetParameter(1);
    double Sigma = f->GetParameter(2);
    //double YMax = f->GetMaximum();
    double XMax = f->GetMaximumX();
    //std::cout<<printf( "Chi2/NDoF: %10.1f  Const: %2.1f  MPV: %3.1f  Sigma: %2.2f", Chi2NDoF, Const/ymax, MPV, Sigma)<<std::endl;
    delete f;
    // Fit again with added gaus
    f = new TF1((std::string(h->GetName())+"_fit").c_str(),"landau(0)+gaus(3)", (XMax<25 ? 0 : XMax-25), XMax+25);
    f->SetParameter(0, 1.4*Const);
    f->SetParLimits(0, 2*ymax,10*ymax);
    f->SetParameter(1, XMax);
    f->SetParLimits(1, XMax<6.5 ? 1.5: XMax-5, XMax + 5);
    f->SetParameter(2, Sigma);
    f->SetParLimits(2, 0.1,10);
    f->SetParameter(3, -0.75*ymax);
    f->SetParLimits(3, -1.5*ymax,ymax);
    f->SetParameter(4, XMax);
    f->SetParLimits(4, 1, 45);
    f->SetParameter(5, Sigma);
    f->SetParLimits(5, 0, 8);
    h->Fit((std::string(h->GetName())+"_fit").c_str(),"RBQ0+");
    //Chi2NDoF = f->GetChisquare()/f->GetNDF();
    //double old_Const = Const;
    //Const = f->GetParameter(0);
    //MPV = f->GetParameter(1);
    Sigma = f->GetParameter(2);
    //double Const2 = f->GetParameter(3);
    //double Mean = f->GetParameter(4);
    //double Sigma2 = f->GetParameter(5);
    //YMax = f->GetMaximum();
    XMax = f->GetMaximumX();
    //std::cout<<printf( "Chi2/NDoF: %4.1f  C1: %2.1f  MPV: %3.1f  S1: %2.2f  C2: %5.2f  M: %3.1f  S2: %2.2f", Chi2NDoF, Const/ymax, MPV, Sigma, Const2/ymax, Mean, Sigma2)<<std::endl;
    mpv = XMax;
    error = f->GetParError(1);
  }
  
  // ******************** DColEfficiency *********************
  void calc_dcol_1d_(TH1D* h1d, TH2D* h2d, bool savemother=1) {
    // Old Phase 0 way: with hiteff to dcol function
    //++  std::string name =h1d->GetName();
    //++  // Check which Layer/Ring is plotted (or if it is a module plot)
    //++  bool mod_plot = (std::string(h1d->GetName()).find("Modules")!=std::string::npos);
    //++  int Lay = 0; // Assume layer 1 if no such postfix is given
    //++  for (int lay=0; lay<3; ++lay) {
    //++    std::stringstream ss; ss<<"Lay"<<lay+1;
    //++    if (name.find(ss.str())!=std::string::npos) Lay = lay;
    //++  }
    //++  int Ring = -1;
    //++  for (int ring=0; ring<4; ring++) {
    //++    std::stringstream ss; ss<<"Mod"<<ring+1;
    //++    if (name.find(ss.str())!=std::string::npos) Ring = ring;
    //++  }
    //++  // Get HitEfficiency and use hiteff_vs_dcol functions to get DCol Efficiency
    //++  calc_eff_1d_(h1d, h2d, 1, savemother);
    //++  for (int bin=1; bin<=h1d->GetNbinsX(); ++bin) {
    //++    if (mod_plot) Ring = abs(bin-5)-1;
    //++    if (Lay==0&&Ring==3) Ring = 2;
    //++    double hiteff = h1d->GetBinContent(bin);
    //++    double dcoleff = (Ring!=-1) ? ( hiteff>ring_fit_[Lay][Ring]->Eval(0.75) ? ring_fit_[Lay][Ring]->GetX(hiteff) : 0 ): 0;
    //++    h1d->SetBinContent(bin, dcoleff);
    //++    h1d->SetBinError(bin, 0);
    //++  }
    //++  h1d->SetEntries(h2d->GetEntries());
    // Direct measurement in Phase 1
    // Wilson Score Interval
    double z = 2; // N Sigma confidence
    for (int i=1; i<=h2d->GetNbinsX(); ++i) {
      double even = h2d->GetBinContent(i,1), odd = h2d->GetBinContent(i,2);
      if (even+odd>0) {
	double eff = std::min(even>0 ? odd/even : 0, 1.0);
	double n = std::max(even,1.0);
	double cen = (eff+(z*z/(2*n))) / (1.0 + (z*z/n));
	double halfwidth = z*sqrt( eff*(1.0-eff)/n + (z*z/(4*n*n)) ) / (1.0 + (z*z/n));
	double err = halfwidth + fabs(cen-eff);
	//double err = z * std::sqrt(1.0/(odd+even)*eff*(1.0-eff)); // Normal approx interval
	// Assymmetric error -> Choose larger for a conservative error estimate
	h1d->SetBinContent(i,eff);
	h1d->SetBinError(i,err);
      }
    }
    h1d->SetEntries(h2d->GetEntries());
    if (savemother) mother_2d_[h1d]=h2d;
    //if (h2d->GetNbinsY()==2&&h2d->GetYaxis()->GetBinCenter(1)==0&&h2d->GetYaxis()->GetBinCenter(2)==1) plot_asymm_err_=1;
    h1d->SetEntries(h2d->GetEntries());
  }

  // ************** Systematic Unc. Calculation **************
  void calc_syst_1d_(TH1D* h1d, TH2D* h2d, bool savemother=1) {
    std::vector<double> v_up, v_down;
    for (int binx=1; binx<=h2d->GetNbinsX(); ++binx) {
      h1d->SetBinContent(binx, h2d->GetBinContent(binx, 1));
      h1d->SetBinError  (binx, h2d->GetBinError(binx, 1));
    }
    h1d->SetEntries(h2d->GetEntries());
    if (savemother) mother_2d_[h1d]=h2d;
  }
  // ******************* Signal Significance  *******************
  void calc_soverb_1d_(TH1D* h1d, TH2D* h2d, bool savemother=0) {
    for (int binx=1; binx<=h2d->GetNbinsX(); ++binx) {
      double bkg = h2d->GetBinContent(binx, 1);
      double sig = h2d->GetBinContent(binx, 2);
      double soverb = sig+bkg>0 ? sig / std::sqrt(sig+bkg) : 0; 
      h1d->SetBinContent(binx, soverb);
    }
    h1d->SetEntries(h2d->GetEntries());
    if (savemother) mother_2d_[h1d]=h2d;
  }
  void calc_soverb_2d_(TH2D* h2d, TH3D* h3d, bool savemother=0) {
    double bkg_underflow = h3d->GetBinContent(0, 0, 1); // MLSP/MGluino/MStop goes to underflow bin
    for (int binx=1; binx<=h3d->GetNbinsX(); ++binx) for (int biny=1; biny<=h3d->GetNbinsY(); ++biny) {
      double bkg = h3d->GetBinContent(binx, biny, 1);
      double sig = h3d->GetBinContent(binx, biny, 2);
      if (bkg_underflow>0) bkg = bkg_underflow;
      double soverb = sig+bkg>0 ? sig / std::sqrt(sig+bkg) : 0; 
      h2d->SetBinContent(binx, biny, soverb);
    }
    h2d->SetEntries(h3d->GetEntries());
    if (savemother) mother_3d_[h2d]=h3d;
  }
  
  //TF1* get_eff_vs_dcol_func_(TH1* h) {
  //  // Check which Layer, Ring is plotted
  //  std::string name =h->GetName();
  //  int Lay = 0; // Assume layer 1 if no such postfix is given
  //  for (int lay=0; lay<3; ++lay) {
  //    std::stringstream ss; ss<<"Lay"<<lay+1;
  //    if (name.find(ss.str())!=std::string::npos) Lay = lay;
  //  }
  //  int Ring = -1;
  //  for (int ring=0; ring<4; ring++) {
  //    std::stringstream ss; ss<<"Ring"<<ring+1;
  //    if (name.find(ss.str())!=std::string::npos) Ring = ring;
  //  }
  //  if (Lay==0&&Ring==3) Ring = 2;
  //  if (Ring!=-1) {
  //    // Get Hit Eff vs DCol Eff function (previously measured)
  //    std::string name[3] = {"hiteff_vs_dcol_l1", "hiteff_vs_dcol_l2", "hiteff_vs_dcol_l3"};
  //    std::string ringname[4] = {"_mod1", "_mod2", "_mod3", "_mod4"};
  //    TFile *f_func = TFile::Open("/data/jkarancs/CMSSW/TimingStudy/CMSSW_7_1_0_pre1/src/DPGAnalysis/PixelTimingStudy/test/DynIneff_scale_factors/HitEffvsDColFunctions.root");
  //    TF1* ring_fit = (TF1*)f_func->Get((name[Lay]+ringname[Ring]).c_str());
  //    f_func->Close();
  //    std::cout<<name<<" "<<name[Lay]+ringname[Ring]<<std::endl;
  //    return ring_fit;
  //  } else return NULL;
  //}
  //
  //void calc_dcol_2d_(TH2D* h2d, TH3D* h3d) {
  //  TF1* ring_fit = get_eff_vs_dcol_func_(h2d);
  //  if (ring_fit) {
  //    // Get HitEfficiency and use hiteff(dcoleff) to get DCol Efficiency
  //    calc_eff_2d_(h2d, h3d);
  //    //for (int binx=1; binx<=h2d->GetNbinsX(); ++binx) for (int biny=1; biny<=h2d->GetNbinsY(); ++biny) {
  //    //  double hiteff = h2d->GetBinContent(binx, biny);
  //    //  double dcoleff = ring_fit->GetX(hiteff);
  //    //  h2d->SetBinContent(binx, biny, dcoleff);
  //    //}
  //  }
  //}
  
  void calc_specials_() {
    if (find_spec_(name_)!=(size_t)-1||find_spec2_(name_)!=(size_t)-1) {
      switch (npf_) {
      case 0:
	switch (ndim_) {
	case 2: calc_spec_1d_(h1d_0p_, h2d_0p_); break;
	case 3: calc_spec_2d_(h2d_0p_, h3d_0p_); break;
	} break;
      case 1:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_1p_.size(); ++i)
	    calc_spec_1d_(h1d_1p_[i], h2d_1p_[i]); break;
	case 3: 
	  for (size_t i=0; i<h3d_1p_.size(); ++i)
	    calc_spec_2d_(h2d_1p_[i], h3d_1p_[i]); break;
	} break;
      case 2:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j)
	    calc_spec_1d_(h1d_2p_[i][j], h2d_2p_[i][j]); break;
	case 3:
	  for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j)
	    calc_spec_2d_(h2d_2p_[i][j], h3d_2p_[i][j]); break;
	} break;
      case 3:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
	    for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) calc_spec_1d_(h1d_3p_[i][j][k], h2d_3p_[i][j][k]); break;
	case 3:
	  for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j)
	    for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) calc_spec_2d_(h2d_3p_[i][j][k], h3d_3p_[i][j][k]); break;
	} break;
      case 4:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) for (size_t k=0; k<h2d_4p_[i][j].size(); ++k)
	    for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) calc_spec_1d_(h1d_4p_[i][j][k][l], h2d_4p_[i][j][k][l]); break;
	case 3:
	  for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) for (size_t k=0; k<h3d_4p_[i][j].size(); ++k)
	    for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) calc_spec_2d_(h2d_4p_[i][j][k][l], h3d_4p_[i][j][k][l]); break;
	} break;
      case 5:
	switch (ndim_) {
	case 2:
	  for (size_t i=0; i<h2d_5p_.size(); ++i) for (size_t j=0; j<h2d_5p_[i].size(); ++j) for (size_t k=0; k<h2d_5p_[i][j].size(); ++k)
	    for (size_t l=0; l<h2d_5p_[i][j][k].size(); ++l) for (size_t m=0; m<h2d_5p_[i][j][k][l].size(); ++m)
	      calc_spec_1d_(h1d_5p_[i][j][k][l][m], h2d_5p_[i][j][k][l][m]); break;
	case 3:
	  for (size_t i=0; i<h3d_5p_.size(); ++i) for (size_t j=0; j<h3d_5p_[i].size(); ++j) for (size_t k=0; k<h3d_5p_[i][j].size(); ++k)
	    for (size_t l=0; l<h3d_5p_[i][j][k].size(); ++l) for (size_t m=0; m<h3d_5p_[i][j][k][l].size(); ++m)
	      calc_spec_2d_(h2d_5p_[i][j][k][l][m], h3d_5p_[i][j][k][l][m]); break;
	} break;
      }
    }
  }
  
  std::string spec_dirname_(TObject* obj) {
    bool is_spec = false;
    std::string dirname = name_;
    std::string hname = obj->GetName();
    if (hname.find(name_)==std::string::npos) for (size_t s=0; s<spec_.size(); ++s) if (hname.find(spec_[s][1])!=std::string::npos) {
      is_spec = true;
      dirname.replace(dirname.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
    }
    if (!is_spec) dirname = "";
    else dirname += "/";
    return dirname;
  }
  
  std::string name_only_pf_(TObject* obj) {
    std::string name = obj->GetName();
    if (npf_>0) {
      std::string del = spec_dirname_(obj);
      if (del.size()) del.replace(del.size()-1,1,"_");
      else del = name_ + "_";
      if (name.find(del)!=std::string::npos)
	name.erase(name.find(del),del.size());
    }
    return name;
  }
  
  void write_(TObject* obj) {
    if (obj) {
      std::string main_dir = name_;
      std::string spec_dir = spec_dirname_(obj);
      bool is_spec = (spec_dir.size()!=0);
      if (!gDirectory->GetKey(main_dir.c_str())) gDirectory->mkdir(main_dir.c_str());
      gDirectory->cd(main_dir.c_str());
      if (is_spec) {
        if (!gDirectory->GetKey(spec_dir.c_str())) gDirectory->mkdir(spec_dir.c_str()); 
        gDirectory->cd(spec_dir.c_str());
      }
      obj->Write(name_only_pf_(obj).c_str()); 
      gDirectory->cd(is_spec?"../..":"..");
    }
  }
  
  void load_(TFile* f, TH1D*& h, bool add = 0) {
    if (h) {
      std::string name = name_ + "/" + spec_dirname_(h) + name_only_pf_(h);
      TH1D* h_temp = (TH1D*)f->Get(name.c_str());
      if (h_temp!=0&&h_temp->GetEntries()>0) { 
	if (add) {
	  double nentries = h->GetEntries();
	  nentries += h_temp->GetEntries();
	  h->Add(h_temp);
	  //for (int i=1; i<=h->GetNbinsX(); ++i) {
	  //  double xval = h_temp->GetBinCenter(i);
	  //  double cont = h_temp->GetBinContent(i);
	  //  double err  = h_temp->GetBinError(i);
	  //  int bin = h->FindBin(xval); 
	  //  h->SetBinContent(bin,h->GetBinContent(bin)+cont);
	  //  h->SetBinError(bin,sqrt(h->GetBinError(bin)*h->GetBinError(bin)+err*err));
	  //}
	  h->SetEntries(nentries);
	}
	else { delete h; h = h_temp; h->SetDirectory(0); }
      }
    }
  }
  
  void load_(TFile* f, TH2D*& h, bool add = 0) {
    if (h) {
      std::string name = name_ + "/" + spec_dirname_(h) + name_only_pf_(h);
      TH2D* h_temp = (TH2D*)f->Get(name.c_str());
      if (h_temp!=0&&h_temp->GetEntries()>0) { 
	if (add) { h->Add(h_temp); h->SetEntries(h->GetEntries()+h_temp->GetEntries()); }
	else { delete h; h = h_temp; h->SetDirectory(0); }
      }
    }
  }
  
  void load_(TFile* f, TH3D*& h, bool add = 0) {
    if (h) {
      std::string name = name_ + "/" + spec_dirname_(h) + name_only_pf_(h);
      TH3D* h_temp = (TH3D*)f->Get(name.c_str());
      if (h_temp!=0&&h_temp->GetEntries()>0) { 
	if (add) { h->Add(h_temp); h->SetEntries(h->GetEntries()+h_temp->GetEntries()); }
	else { delete h; h = h_temp; h->SetDirectory(0); }
      }
    }
  }
  
  // Load/Add all histos from a file
  void load_all_(TFile* f, bool add = 0) {
    if (npf_==0) {
      load_(f,h1d_0p_,add);
      load_(f,h2d_0p_,add);
      load_(f,h3d_0p_,add);
    } else if (npf_==1) {
      for (size_t i=0; i<h1d_1p_.size(); ++i) load_(f,h1d_1p_[i],add);
      for (size_t i=0; i<h2d_1p_.size(); ++i) load_(f,h2d_1p_[i],add);
      for (size_t i=0; i<h3d_1p_.size(); ++i) load_(f,h3d_1p_[i],add);
    } else if (npf_==2) {
      for (size_t i=0; i<h1d_2p_.size(); ++i) for (size_t j=0; j<h1d_2p_[i].size(); ++j) load_(f,h1d_2p_[i][j],add);
      for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j) load_(f,h2d_2p_[i][j],add);
      for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j) load_(f,h3d_2p_[i][j],add);
    } else if (npf_==3) {
      for (size_t i=0; i<h1d_3p_.size(); ++i) for (size_t j=0; j<h1d_3p_[i].size(); ++j)
        for (size_t k=0; k<h1d_3p_[i][j].size(); ++k) load_(f,h1d_3p_[i][j][k],add);
      for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) load_(f,h2d_3p_[i][j][k],add);
      for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) load_(f,h3d_3p_[i][j][k],add);
    } else if (npf_==4) {
      for (size_t i=0; i<h1d_4p_.size(); ++i) for (size_t j=0; j<h1d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_4p_[i][j].size(); ++k) for (size_t l=0; l<h1d_4p_[i][j][k].size(); ++l) load_(f,h1d_4p_[i][j][k][l],add);
      for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_4p_[i][j].size(); ++k) for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) load_(f,h2d_4p_[i][j][k][l],add);
      for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_4p_[i][j].size(); ++k) for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) load_(f,h3d_4p_[i][j][k][l],add);
    } else if (npf_==5) {
      for (size_t i=0; i<h1d_5p_.size(); ++i) for (size_t j=0; j<h1d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_5p_[i][j].size(); ++k) for (size_t l=0; l<h1d_5p_[i][j][k].size(); ++l) 
	  for (size_t m=0; m<h1d_5p_[i][j][k][l].size(); ++m) load_(f,h1d_5p_[i][j][k][l][m],add);
      for (size_t i=0; i<h2d_5p_.size(); ++i) for (size_t j=0; j<h2d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_5p_[i][j].size(); ++k) for (size_t l=0; l<h2d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h2d_5p_[i][j][k][l].size(); ++m) load_(f,h2d_5p_[i][j][k][l][m],add);
      for (size_t i=0; i<h3d_5p_.size(); ++i) for (size_t j=0; j<h3d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_5p_[i][j].size(); ++k) for (size_t l=0; l<h3d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h3d_5p_[i][j][k][l].size(); ++m) load_(f,h3d_5p_[i][j][k][l][m],add);
    }
  }
  
  bool pass_cuts_() {
    switch (ncut_) {
    case 0: return 1; break;
    case 1: return (cut1_->Eval()); break;
    case 2: return (cut1_->Eval())&&(cut2_->Eval()); break;
    case 3: return (cut1_->Eval())&&(cut2_->Eval())&&(cut3_->Eval()); break;
    case 4: return (cut1_->Eval())&&(cut2_->Eval())&&(cut3_->Eval())&&(cut4_->Eval()); break;
    case 5: return (cut1_->Eval())&&(cut2_->Eval())&&(cut3_->Eval())&&(cut4_->Eval())&&(cut5_->Eval()); break;
    case 6: return (cut1_->Eval())&&(cut2_->Eval())&&(cut3_->Eval())&&(cut4_->Eval())&&(cut5_->Eval())&&(cut6_->Eval()); break;
    default: return 1; // Warning in constructor
    }
  }
  
  void set_histo_options_(TH1* h) {
    h->SetDirectory(0);
    if (sumw2_) h->Sumw2();
    if (months_||weeks_||dates_) {
      h->GetXaxis()->SetTimeDisplay(1);
      // Set suitable offset (correct for leap days or monday shift for weeks)
      if (months_||dates_) {
	if (TString(h->GetName()).Contains("2012")||TString(h->GetName()).Contains("2016")||TString(h->GetName()).Contains("2020"))
	  h->GetXaxis()->SetTimeOffset(TDatime(2012,1,1,0,0,0).Convert());
	else
	  h->GetXaxis()->SetTimeOffset(TDatime(2018,1,1,0,0,0).Convert());
      } else h->GetXaxis()->SetTimeOffset(TDatime(2018,1,1,0,0,0).Convert()); // 2018 starts on Monday, also not leap year
      // Set time format
      if (months_) h->GetXaxis()->SetTimeFormat("%b");
      if (weeks_)  h->GetXaxis()->SetTimeFormat("%V");
      if (dates_)  h->GetXaxis()->SetTimeFormat("%d/%m");
    }
    for (size_t iaxis=0; iaxis<ndim_; ++iaxis) {
      TAxis* axis = iaxis==0 ? h->GetXaxis() : iaxis==1 ? h->GetYaxis() : h->GetZaxis();
      int set_opt = 0;
      for (auto pair : bin_labels_[iaxis]) { 
	axis->SetBinLabel(pair.first, pair.second.c_str());
	if (pair.second.size()>8) set_opt = 1;
      }
      //if (set_opt==1) axis->LabelsOption("d"); // down 
      if (set_opt==1) axis->LabelsOption("v"); // vertical
      axis->SetLabelSize(0.05);
    }
  }
  
public:
  const std::string& GetName() { return name_; }
  const std::vector<std::string>& GetPFNames() { return pf_names_; }
  const std::vector<std::vector<std::string> >& GetSpec2() { return spec2_; }
  
  // Add New SmartHisto to the container vector
  // and set Filling properties, Postfixes, titles etc...
  // 2D/3D: If one of the Variable is special, also create
  // -1D object for the result
  // Eg: 2D: ClusterCharge vs. Variable
  // --> 1D: MPV vs. Variable
  //
  // 1D
  void AddNew(std::string name, std::string title,
	      size_t nbinsx, std::vector<double> binsx) {
    Double_t xbins[nbinsx+1]; for (size_t i=0; i<binsx.size(); ++i) if (i<nbinsx+1) xbins[i] = binsx[i];
    if (binsx.size()==2&&nbinsx>1) for (size_t i=0; i<=nbinsx; ++i) xbins[i] = binsx[0] + i*(binsx[1]-binsx[0])/nbinsx;
    if (npf_==0) {
      if (binsx.size()==2) h1d_0p_ = new TH1D(name.c_str(), title.c_str(), nbinsx, binsx[0], binsx[1]);
      else h1d_0p_ = new TH1D(name.c_str(), title.c_str(), nbinsx, xbins);
      set_histo_options_(h1d_0p_);
    } else {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	if (npf_==1) {
	  if (binsx.size()==2) h1d_1p_.push_back(new TH1D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(), nbinsx, binsx[0], binsx[1]));
	  else h1d_1p_.push_back(new TH1D((name+"_"+pfs_[0].vec[i]).c_str(), title.c_str(), nbinsx, xbins));
	  set_histo_options_(h1d_1p_[i]);
	} else {
	  if (npf_==2) h1d_2p_.push_back(std::vector<TH1D*>());
	  else if (npf_==3) h1d_3p_.push_back(std::vector<std::vector<TH1D*> > ());
          else if (npf_==4) h1d_4p_.push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	  else if (npf_==5) h1d_5p_.push_back(std::vector<std::vector<std::vector<std::vector<TH1D*> > > >());
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    if (npf_==2) {
	      if (binsx.size()==2) h1d_2p_[i].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
								 title.c_str(), nbinsx, binsx[0], binsx[1]));
	      else h1d_2p_[i].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
								 title.c_str(), nbinsx, xbins));
	      set_histo_options_(h1d_2p_[i][j]);
	    } else {
	      if (npf_==3) h1d_3p_[i].push_back(std::vector<TH1D*>());
	      else if (npf_==4) h1d_4p_[i].push_back(std::vector<std::vector<TH1D*> >());
	      else if (npf_==5) h1d_5p_[i].push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	      for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
		if (npf_==3) {
		  if (binsx.size()==2) h1d_3p_[i][j].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
									title.c_str(), nbinsx, binsx[0], binsx[1]));
		  else h1d_3p_[i][j].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
									title.c_str(), nbinsx, xbins));
		  set_histo_options_(h1d_3p_[i][j][k]);
		} else {
		  if (npf_==4) h1d_4p_[i][j].push_back(std::vector<TH1D*>());
		  else if (npf_==5) h1d_5p_[i][j].push_back(std::vector<std::vector<TH1D*> >());
		  for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
		    if (npf_==4) {
		      if (binsx.size()==2) h1d_4p_[i][j][k].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
									       title.c_str(), nbinsx, binsx[0], binsx[1]));
		      else h1d_4p_[i][j][k].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
									       title.c_str(), nbinsx, xbins));
		      set_histo_options_(h1d_4p_[i][j][k][l]);
		    } else {
		      if (npf_==5) h1d_5p_[i][j][k].push_back(std::vector<TH1D*>());
		      for (size_t m=0; m<pfs_[4].vec.size(); ++m) {
			if (npf_==5) {
			  if (binsx.size()==2) h1d_5p_[i][j][k][l].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
										      title.c_str(), nbinsx, binsx[0], binsx[1]));
			  else h1d_5p_[i][j][k][l].push_back(new TH1D((name+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
										      title.c_str(), nbinsx, xbins));
			  set_histo_options_(h1d_5p_[i][j][k][l][m]);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  // 2D
  void AddNew(std::string name_1d, std::string title_1d,
	      std::string name_2d, std::string title_2d,
	      size_t nbinsx, std::vector<double> binsx,
	      size_t nbinsy, std::vector<double> binsy) {
    bool spec = name_1d!=name_2d || title_1d != title_2d;
    Double_t xbins[nbinsx+1]; for (size_t i=0; i<binsx.size(); ++i) if (i<nbinsx+1) xbins[i] = binsx[i];
    Double_t ybins[nbinsy+1]; for (size_t i=0; i<binsy.size(); ++i) if (i<nbinsy+1) ybins[i] = binsy[i];
    if (binsx.size()==2&&nbinsx>1) for (size_t i=0; i<=nbinsx; ++i) xbins[i] = binsx[0] + i*(binsx[1]-binsx[0])/nbinsx;
    if (binsy.size()==2&&nbinsy>1) for (size_t i=0; i<=nbinsy; ++i) ybins[i] = binsy[0] + i*(binsy[1]-binsy[0])/nbinsy;
    if (npf_==0) {
      if (binsx.size()==2&&binsy.size()==2) {
	h2d_0p_ = new TH2D(name_2d.c_str(), title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]);
	if (spec) h1d_0p_ = new TH1D(name_1d.c_str(), title_1d.c_str(), nbinsx, binsx[0], binsx[1]);
      } else {
	h2d_0p_ = new TH2D(name_2d.c_str(), title_2d.c_str(), nbinsx, xbins, nbinsy, ybins);
	if (spec) h1d_0p_ = new TH1D(name_1d.c_str(), title_1d.c_str(), nbinsx, xbins);
      }
      set_histo_options_(h2d_0p_);
      if (spec) set_histo_options_(h1d_0p_);
    } else {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	if (npf_==1) {
	  if (binsx.size()==2&&binsy.size()==2) {
	    h2d_1p_.push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]).c_str(), title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
	    if (spec) h1d_1p_.push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]).c_str(), title_1d.c_str(), nbinsx, binsx[0], binsx[1]));
	  } else {
	    h2d_1p_.push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]).c_str(), title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
	    if (spec) h1d_1p_.push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]).c_str(), title_1d.c_str(), nbinsx, xbins));
	  }
	  set_histo_options_(h2d_1p_[i]);
	  if (spec) set_histo_options_(h1d_1p_[i]);
	} else {
	  if (npf_==2) h2d_2p_.push_back(std::vector<TH2D*>());
	  else if (npf_==3) h2d_3p_.push_back(std::vector<std::vector<TH2D*> > ());
	  else if (npf_==4) h2d_4p_.push_back(std::vector<std::vector<std::vector<TH2D*> > >());
	  else if (npf_==5) h2d_5p_.push_back(std::vector<std::vector<std::vector<std::vector<TH2D*> > > >());
	  if (spec) {
	    if (npf_==2) h1d_2p_.push_back(std::vector<TH1D*>());
	    else if (npf_==3) h1d_3p_.push_back(std::vector<std::vector<TH1D*> > ());
	    else if (npf_==4) h1d_4p_.push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	    else if (npf_==5) h1d_5p_.push_back(std::vector<std::vector<std::vector<std::vector<TH1D*> > > >());
	  }
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    if (npf_==2) {
	      if (binsx.size()==2&&binsy.size()==2) {
		h2d_2p_[i].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
					      title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
		if (spec) h1d_2p_[i].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), 
							title_1d.c_str(), nbinsx, binsx[0], binsx[1]));
	      } else {
		h2d_2p_[i].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
					      title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
		if (spec) h1d_2p_[i].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), 
							title_1d.c_str(), nbinsx, xbins));
	      }
	      set_histo_options_(h2d_2p_[i][j]);
	      if (spec) set_histo_options_(h1d_2p_[i][j]);
	    } else {
	      if (npf_==3) h2d_3p_[i].push_back(std::vector<TH2D*>());
	      else if (npf_==4) h2d_4p_[i].push_back(std::vector<std::vector<TH2D*> >());
	      else if (npf_==5) h2d_5p_[i].push_back(std::vector<std::vector<std::vector<TH2D*> > >());
	      if (spec) {
		if (npf_==3) h1d_3p_[i].push_back(std::vector<TH1D*>());
		else if (npf_==4) h1d_4p_[i].push_back(std::vector<std::vector<TH1D*> >());
		else if (npf_==5) h1d_5p_[i].push_back(std::vector<std::vector<std::vector<TH1D*> > >());
	      }
	      for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
		if (npf_==3) {
		  if (binsx.size()==2&&binsy.size()==2) {
		    h2d_3p_[i][j].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
						     title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
		    if (spec) h1d_3p_[i][j].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
							       title_1d.c_str(), nbinsx, binsx[0], binsx[1]));
		  } else {
		    h2d_3p_[i][j].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
						     title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
		    if (spec) h1d_3p_[i][j].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
							       title_1d.c_str(), nbinsx, xbins));
		  }
		  set_histo_options_(h2d_3p_[i][j][k]);
		  if (spec) set_histo_options_(h1d_3p_[i][j][k]);
		} else {
		  if (npf_==4) h2d_4p_[i][j].push_back(std::vector<TH2D*>());
		  else if (npf_==5) h2d_5p_[i][j].push_back(std::vector<std::vector<TH2D*> >());
		  if (spec) {
		    if (npf_==4) h1d_4p_[i][j].push_back(std::vector<TH1D*>());
		    else if (npf_==5) h1d_5p_[i][j].push_back(std::vector<std::vector<TH1D*> >());
		  }
		  for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
		    if (npf_==4) {
		      if (binsx.size()==2&&binsy.size()==2) {
			h2d_4p_[i][j][k].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
							    title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
			if (spec) h1d_4p_[i][j][k].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
								      title_1d.c_str(), nbinsx, binsx[0], binsx[1]));
		      } else {
			h2d_4p_[i][j][k].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
							    title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
			if (spec) h1d_4p_[i][j][k].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
								      title_1d.c_str(), nbinsx, xbins));
		      }
		      set_histo_options_(h2d_4p_[i][j][k][l]);
		      if (spec) set_histo_options_(h1d_4p_[i][j][k][l]);
		    } else {
		      if (npf_==5) h2d_5p_[i][j][k].push_back(std::vector<TH2D*>());
		      if (spec) { if (npf_==5) h1d_5p_[i][j][k].push_back(std::vector<TH1D*>()); }
		      for (size_t m=0; m<pfs_[4].vec.size(); ++m) {
			if (npf_==5) {
			  if (binsx.size()==2&&binsy.size()==2) {
			    h2d_5p_[i][j][k][l].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
								   title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
			    if (spec) h1d_5p_[i][j][k][l].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
									     title_1d.c_str(), nbinsx, binsx[0], binsx[1]));
			  } else {
			    h2d_5p_[i][j][k][l].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
								   title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
			    if (spec) h1d_5p_[i][j][k][l].push_back(new TH1D((name_1d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
									     title_1d.c_str(), nbinsx, xbins));
			  }
			  set_histo_options_(h2d_5p_[i][j][k][l][m]);
			  if (spec) set_histo_options_(h1d_5p_[i][j][k][l][m]);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // 3D
  void AddNew(std::string name_2d, std::string title_2d,
	      std::string name_3d, std::string title_3d,
	      size_t nbinsx, std::vector<double> binsx,
	      size_t nbinsy, std::vector<double> binsy,
	      size_t nbinsz, std::vector<double> binsz) {
    bool spec = name_2d!=name_3d || title_2d != title_3d;
    Double_t xbins[nbinsx+1]; for (size_t i=0; i<binsx.size(); ++i) if (i<nbinsx+1) xbins[i] = binsx[i];
    Double_t ybins[nbinsy+1]; for (size_t i=0; i<binsy.size(); ++i) if (i<nbinsy+1) ybins[i] = binsy[i];
    Double_t zbins[nbinsz+1]; for (size_t i=0; i<binsz.size(); ++i) if (i<nbinsz+1) zbins[i] = binsz[i];
    if (binsx.size()==2&&nbinsx>1) for (size_t i=0; i<=nbinsx; ++i) xbins[i] = binsx[0] + i*(binsx[1]-binsx[0])/nbinsx;
    if (binsy.size()==2&&nbinsy>1) for (size_t i=0; i<=nbinsy; ++i) ybins[i] = binsy[0] + i*(binsy[1]-binsy[0])/nbinsy;
    if (binsz.size()==2&&nbinsz>1) for (size_t i=0; i<=nbinsz; ++i) zbins[i] = binsz[0] + i*(binsz[1]-binsz[0])/nbinsz;
    if (npf_==0) {
      if (binsx.size()==2&&binsy.size()==2&&binsz.size()==2) {
	h3d_0p_ = new TH3D(name_3d.c_str(), title_3d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1], nbinsz, binsz[0], binsz[1]);
	if (spec) h2d_0p_ = new TH2D(name_2d.c_str(), title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]);
      } else {
	h3d_0p_ = new TH3D(name_3d.c_str(), title_3d.c_str(), nbinsx, xbins, nbinsy, ybins, nbinsz, zbins);
	if (spec) h2d_0p_ = new TH2D(name_2d.c_str(), title_2d.c_str(), nbinsx, xbins, nbinsy, ybins);
      }
      set_histo_options_(h3d_0p_);
      if (spec) set_histo_options_(h2d_0p_);
    } else {
      for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	if (npf_==1) {
	  if (binsx.size()==2&&binsy.size()==2&&binsz.size()==2) {
	    h3d_1p_.push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]).c_str(), title_3d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1], nbinsz, binsz[0], binsz[1]));
	    if (spec) h2d_1p_.push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]).c_str(), title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
	  } else {
	    h3d_1p_.push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]).c_str(), title_3d.c_str(), nbinsx, xbins, nbinsy, ybins, nbinsz, zbins));
	    if (spec) h2d_1p_.push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]).c_str(), title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
	  }
	  set_histo_options_(h3d_1p_[i]);
	  if (spec) set_histo_options_(h2d_1p_[i]);
	} else {
	  if (npf_==2) h3d_2p_.push_back(std::vector<TH3D*>());
	  else if (npf_==3) h3d_3p_.push_back(std::vector<std::vector<TH3D*> > ());
	  else if (npf_==4) h3d_4p_.push_back(std::vector<std::vector<std::vector<TH3D*> > >());
	  else if (npf_==5) h3d_5p_.push_back(std::vector<std::vector<std::vector<std::vector<TH3D*> > > >());
	  if (spec) {
	    if (npf_==2) h2d_2p_.push_back(std::vector<TH2D*>());
	    else if (npf_==3) h2d_3p_.push_back(std::vector<std::vector<TH2D*> > ());
	    else if (npf_==4) h2d_4p_.push_back(std::vector<std::vector<std::vector<TH2D*> > >());
	    else if (npf_==5) h2d_5p_.push_back(std::vector<std::vector<std::vector<std::vector<TH2D*> > > >());
	  }
	  for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	    if (npf_==2) {
	      if (binsx.size()==2&&binsy.size()==2&&binsz.size()==2) {
		h3d_2p_[i].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
					      title_3d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1], nbinsz, binsz[0], binsz[1]));
		if (spec) h2d_2p_[i].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), 
							title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
	      } else {
		h3d_2p_[i].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(),
					      title_3d.c_str(), nbinsx, xbins, nbinsy, ybins, nbinsz, zbins));
		if (spec) h2d_2p_[i].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]).c_str(), 
							title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
	      }
	      set_histo_options_(h3d_2p_[i][j]);
	      if (spec) set_histo_options_(h2d_2p_[i][j]);
	    } else {
	      if (npf_==3) h3d_3p_[i].push_back(std::vector<TH3D*>());
	      else if (npf_==4) h3d_4p_[i].push_back(std::vector<std::vector<TH3D*> >());
	      else if (npf_==5) h3d_5p_[i].push_back(std::vector<std::vector<std::vector<TH3D*> > >());
	      if (spec) {
		if (npf_==3) h2d_3p_[i].push_back(std::vector<TH2D*>());
		else if (npf_==4) h2d_4p_[i].push_back(std::vector<std::vector<TH2D*> >());
		else if (npf_==5) h2d_5p_[i].push_back(std::vector<std::vector<std::vector<TH2D*> > >());
	      }
	      for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
		if (npf_==3) {
		  if (binsx.size()==2&&binsy.size()==2&&binsz.size()==2) {
		    h3d_3p_[i][j].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
						     title_3d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1], nbinsz, binsz[0], binsz[1]));
		    if (spec) h2d_3p_[i][j].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
							       title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
		  } else {
		    h3d_3p_[i][j].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
						     title_3d.c_str(), nbinsx, xbins, nbinsy, ybins, nbinsz, zbins));
		    if (spec) h2d_3p_[i][j].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]).c_str(),
							       title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
		  }
		  set_histo_options_(h3d_3p_[i][j][k]);
		  if (spec) set_histo_options_(h2d_3p_[i][j][k]);
		} else {
		  if (npf_==4) h3d_4p_[i][j].push_back(std::vector<TH3D*>());
		  else if (npf_==5) h3d_5p_[i][j].push_back(std::vector<std::vector<TH3D*> >());
		  if (spec) {
		    if (npf_==4) h2d_4p_[i][j].push_back(std::vector<TH2D*>());
		    else if (npf_==5) h2d_5p_[i][j].push_back(std::vector<std::vector<TH2D*> >());
		  }
		  for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
		    if (npf_==4) {
		      if (binsx.size()==2&&binsy.size()==2&&binsz.size()==2) {
			h3d_4p_[i][j][k].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
							    title_3d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1], nbinsz, binsz[0], binsz[1]));
			if (spec) h2d_4p_[i][j][k].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
								      title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
		      } else {
			h3d_4p_[i][j][k].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
							    title_3d.c_str(), nbinsx, xbins, nbinsy, ybins, nbinsz, zbins));
			if (spec) h2d_4p_[i][j][k].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]).c_str(),
								      title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
		      }
		      set_histo_options_(h3d_4p_[i][j][k][l]);
		      if (spec) set_histo_options_(h2d_4p_[i][j][k][l]);
		    } else {
		      if (npf_==5) h3d_5p_[i][j][k].push_back(std::vector<TH3D*>());
		      if (spec) { if (npf_==5) h2d_5p_[i][j][k].push_back(std::vector<TH2D*>()); }
		      for (size_t m=0; m<pfs_[4].vec.size(); ++m) {
			if (npf_==5) {
			  if (binsx.size()==2&&binsy.size()==2&&binsz.size()==2) {
			    h3d_5p_[i][j][k][l].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
								   title_3d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1], nbinsz, binsz[0], binsz[1]));
			    if (spec) h2d_5p_[i][j][k][l].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
									     title_2d.c_str(), nbinsx, binsx[0], binsx[1], nbinsy, binsy[0], binsy[1]));
			  } else {
			    h3d_5p_[i][j][k][l].push_back(new TH3D((name_3d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
								   title_3d.c_str(), nbinsx, xbins, nbinsy, ybins, nbinsz, zbins));
			    if (spec) h2d_5p_[i][j][k][l].push_back(new TH2D((name_2d+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]).c_str(),
									     title_2d.c_str(), nbinsx, xbins, nbinsy, ybins));
			  }
			  set_histo_options_(h3d_5p_[i][j][k][l][m]);
			  if (spec) set_histo_options_(h2d_5p_[i][j][k][l][m]);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  // Fill Histograms using the std::function<double()>
  void Fill(const bool debug = 0) {
    if (debug) {
    //if (name_=="HitEfficiency_vs_InstLumi") {
    //if (name_=="HitEfficiency_vs_LayersDisks"&&npf_==0) {
    //if (name_=="Counts_vs_NVtx"&&npf_==1) {
      std::cout<<name_<<" ";
      for (size_t i=0, n=pf_names_.size(); i<n; ++i) std::cout<<pf_names_[i]<<" ";
      std::cout<<" - pass_cuts: "<<pass_cuts_();
      if (npf_) {
	std::cout<<" pfs: "<<pfs_[0].sel();
	if (npf_>1) std::cout<<", "<<pfs_[1].sel();
	if (npf_>2) std::cout<<", "<<pfs_[2].sel();
	if (npf_>3) std::cout<<", "<<pfs_[3].sel();
	if (npf_>4) std::cout<<", "<<pfs_[4].sel();
      }
      std::cout<<" fill: "<<fill_1d_();
      if (ndim_>1) std::cout<<", "<<fill_2d_();
      if (ndim_>2) std::cout<<", "<<fill_3d_();
      std::cout<<" weight: ";
      if (nweight_==0) std::cout<<1;
      if (nweight_==1) std::cout<<weight1_();
      if (nweight_==2) std::cout<<weight1_()<<" * "<<weight2_();
      if (nweight_==3) std::cout<<weight1_()<<" * "<<weight2_()<<" * "<<weight3_();
      std::cout<<"\n";
    //}
    }
    //std::cout<<pfs_[0].sel()<<std::endl;
    //std::cout<<fill_1d_()<<std::endl;
    //std::cout<<weight1_()<<std::endl;
    //if (name_=="OnCluChargeNorm"&&npf_==2&&pfs_[1].sel()==0&&pass_cuts_()==1) std::cout<<pfs_[0].sel()<<" "<<pfs_[1].sel()<<" "<<fill_1d_()<<" "<<weight1_()<<std::endl;
    //if (name_=="HitEfficiency_vs_LayersDisks"&&npf_==1&&pass_cuts_()==1) std::cout<<pfs_[0].sel()<<" "<<fill_1d_()<<" "<<weight1_()<<std::endl;
    //if (name_=="JetMatchedGenTopPtCoarse"&&npf_==2&&pass_cuts_()==1) std::cout<<pfs_[0].sel()<<" "<<pfs_[1].sel()<<" "<<fill_1d_()<<" "<<weight1_()<<std::endl;
    //if (name_=="NVerticesReweighted") std::cout<<pfs_[0].sel()<<" "<<pfs_[1].sel()<<" "<<fill_1d_()<<" "<<weight1_()<<std::endl;
    if (pass_cuts_()) {
      double weight = 0;
      switch (nweight_) {
      case 0: weight = 1; break;
      case 1: weight = (weight1_()); break;
      case 2: weight = (weight1_())*(weight2_()); break;
      case 3: weight = (weight1_())*(weight2_())*(weight3_()); break;
      default: break;
      }
      switch (npf_) {
      case 0:
        switch (ndim_) {
        case 1: h1d_0p_->Fill(fill_1d_(),weight); break;
        case 2: h2d_0p_->Fill(fill_1d_(),fill_2d_(),weight); break;
        case 3: h3d_0p_->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
        } break;
      case 1:
        if (pfs_[0].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_1p_[pfs_[0].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_1p_[pfs_[0].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_1p_[pfs_[0].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 2:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_2p_[pfs_[0].sel()][pfs_[1].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_2p_[pfs_[0].sel()][pfs_[1].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_2p_[pfs_[0].sel()][pfs_[1].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 3:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1&&pfs_[2].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_3p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_3p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_3p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 4:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1&&pfs_[2].sel()!=(size_t)-1&&pfs_[3].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_4p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_4p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_4p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      case 5:
        if (pfs_[0].sel()!=(size_t)-1&&pfs_[1].sel()!=(size_t)-1&&pfs_[2].sel()!=(size_t)-1&&pfs_[3].sel()!=(size_t)-1&&pfs_[4].sel()!=(size_t)-1) { switch (ndim_) {
          case 1: h1d_5p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()][pfs_[4].sel()]->Fill(fill_1d_(),weight); break;
          case 2: h2d_5p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()][pfs_[4].sel()]->Fill(fill_1d_(),fill_2d_(),weight); break;
          case 3: h3d_5p_[pfs_[0].sel()][pfs_[1].sel()][pfs_[2].sel()][pfs_[3].sel()][pfs_[4].sel()]->Fill(fill_1d_(),fill_2d_(),fill_3d_(),weight); break;
          }
        } break;
      }
    }
  }
  
  void Load(TFile* f) { load_all_(f); calc_specials_(); }
  
  void Add(TFile* f) { load_all_(f, 1); calc_specials_(); }
  
  void CalcSpecials() { calc_specials_(); }

  // Write histos in a file
  void Write(const bool debug = 0) {
    if (debug) {
      std::cout<<name_<<" ";
      for (size_t i=0, n=pf_names_.size(); i<n; ++i) std::cout<<pf_names_[i]<<" ";
      std::cout<<std::endl;
    }
    if (npf_==0) {
      if (h1d_0p_) write_(h1d_0p_);
      if (h2d_0p_) write_(h2d_0p_);
      if (h3d_0p_) write_(h3d_0p_);
    } else if (npf_==1) {
      for (size_t i=0; i<h1d_1p_.size(); ++i) if (h1d_1p_[i]->GetEntries()) write_(h1d_1p_[i]);
      for (size_t i=0; i<h2d_1p_.size(); ++i) if (h2d_1p_[i]->GetEntries()) write_(h2d_1p_[i]);
      for (size_t i=0; i<h3d_1p_.size(); ++i) if (h3d_1p_[i]->GetEntries()) write_(h3d_1p_[i]);
    } else if (npf_==2) {
      for (size_t i=0; i<h1d_2p_.size(); ++i) for (size_t j=0; j<h1d_2p_[i].size(); ++j) if (h1d_2p_[i][j]->GetEntries()) write_(h1d_2p_[i][j]);
      for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j) if (h2d_2p_[i][j]->GetEntries()) write_(h2d_2p_[i][j]);
      for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j) if (h3d_2p_[i][j]->GetEntries()) write_(h3d_2p_[i][j]);
    } else if (npf_==3) {
      for (size_t i=0; i<h1d_3p_.size(); ++i) for (size_t j=0; j<h1d_3p_[i].size(); ++j)
        for (size_t k=0; k<h1d_3p_[i][j].size(); ++k) if (h1d_3p_[i][j][k]->GetEntries()) write_(h1d_3p_[i][j][k]);
      for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) if (h2d_3p_[i][j][k]->GetEntries()) write_(h2d_3p_[i][j][k]);
      for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) if (h3d_3p_[i][j][k]->GetEntries()) write_(h3d_3p_[i][j][k]);
    } else if (npf_==4) {
      for (size_t i=0; i<h1d_4p_.size(); ++i) for (size_t j=0; j<h1d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_4p_[i][j].size(); ++k) for (size_t l=0; l<h1d_4p_[i][j][k].size(); ++l) 
	  if (h1d_4p_[i][j][k][l]->GetEntries()) write_(h1d_4p_[i][j][k][l]);
      for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_4p_[i][j].size(); ++k) for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) 
	  if (h2d_4p_[i][j][k][l]->GetEntries()) write_(h2d_4p_[i][j][k][l]);
      for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_4p_[i][j].size(); ++k) for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) 
	  if (h3d_4p_[i][j][k][l]->GetEntries()) write_(h3d_4p_[i][j][k][l]);
    } else if (npf_==5) {
      for (size_t i=0; i<h1d_5p_.size(); ++i) for (size_t j=0; j<h1d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_5p_[i][j].size(); ++k) for (size_t l=0; l<h1d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h1d_5p_[i][j][k][l].size(); ++m) if (h1d_5p_[i][j][k][l][m]->GetEntries()) write_(h1d_5p_[i][j][k][l][m]);
      for (size_t i=0; i<h2d_5p_.size(); ++i) for (size_t j=0; j<h2d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_5p_[i][j].size(); ++k) for (size_t l=0; l<h2d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h2d_5p_[i][j][k][l].size(); ++m) if (h2d_5p_[i][j][k][l][m]->GetEntries()) write_(h2d_5p_[i][j][k][l][m]);
      for (size_t i=0; i<h3d_5p_.size(); ++i) for (size_t j=0; j<h3d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_5p_[i][j].size(); ++k) for (size_t l=0; l<h3d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h3d_5p_[i][j][k][l].size(); ++m) if (h3d_5p_[i][j][k][l][m]->GetEntries()) write_(h3d_5p_[i][j][k][l][m]);
    }
  }

  int GetTotalNCells() {
    int total = 0;
    if (npf_==0) {
      if (h1d_0p_) total += h1d_0p_->GetNcells();
      if (h2d_0p_) total += h2d_0p_->GetNcells();
      if (h3d_0p_) total += h3d_0p_->GetNcells();
    } else if (npf_==1) {
      for (size_t i=0; i<h1d_1p_.size(); ++i) total += h1d_1p_[i]->GetNcells();
      for (size_t i=0; i<h2d_1p_.size(); ++i) total += h2d_1p_[i]->GetNcells();
      for (size_t i=0; i<h3d_1p_.size(); ++i) total += h3d_1p_[i]->GetNcells();
    } else if (npf_==2) {
      for (size_t i=0; i<h1d_2p_.size(); ++i) for (size_t j=0; j<h1d_2p_[i].size(); ++j) total += h1d_2p_[i][j]->GetNcells();
      for (size_t i=0; i<h2d_2p_.size(); ++i) for (size_t j=0; j<h2d_2p_[i].size(); ++j) total += h2d_2p_[i][j]->GetNcells();
      for (size_t i=0; i<h3d_2p_.size(); ++i) for (size_t j=0; j<h3d_2p_[i].size(); ++j) total += h3d_2p_[i][j]->GetNcells();
    } else if (npf_==3) {
      for (size_t i=0; i<h1d_3p_.size(); ++i) for (size_t j=0; j<h1d_3p_[i].size(); ++j)
        for (size_t k=0; k<h1d_3p_[i][j].size(); ++k) total += h1d_3p_[i][j][k]->GetNcells();
      for (size_t i=0; i<h2d_3p_.size(); ++i) for (size_t j=0; j<h2d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_3p_[i][j].size(); ++k) total += h2d_3p_[i][j][k]->GetNcells();
      for (size_t i=0; i<h3d_3p_.size(); ++i) for (size_t j=0; j<h3d_3p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_3p_[i][j].size(); ++k) total += h3d_3p_[i][j][k]->GetNcells();
    } else if (npf_==4) {
      for (size_t i=0; i<h1d_4p_.size(); ++i) for (size_t j=0; j<h1d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_4p_[i][j].size(); ++k) for (size_t l=0; l<h1d_4p_[i][j][k].size(); ++l) 
	  total += h1d_4p_[i][j][k][l]->GetNcells();
      for (size_t i=0; i<h2d_4p_.size(); ++i) for (size_t j=0; j<h2d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_4p_[i][j].size(); ++k) for (size_t l=0; l<h2d_4p_[i][j][k].size(); ++l) 
	  total += h2d_4p_[i][j][k][l]->GetNcells();
      for (size_t i=0; i<h3d_4p_.size(); ++i) for (size_t j=0; j<h3d_4p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_4p_[i][j].size(); ++k) for (size_t l=0; l<h3d_4p_[i][j][k].size(); ++l) 
	  total += h3d_4p_[i][j][k][l]->GetNcells();
    } else if (npf_==5) {
      for (size_t i=0; i<h1d_5p_.size(); ++i) for (size_t j=0; j<h1d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h1d_5p_[i][j].size(); ++k) for (size_t l=0; l<h1d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h1d_5p_[i][j][k][l].size(); ++m) total += h1d_5p_[i][j][k][l][m]->GetNcells();
      for (size_t i=0; i<h2d_5p_.size(); ++i) for (size_t j=0; j<h2d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h2d_5p_[i][j].size(); ++k) for (size_t l=0; l<h2d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h2d_5p_[i][j][k][l].size(); ++m) total += h2d_5p_[i][j][k][l][m]->GetNcells();
      for (size_t i=0; i<h3d_5p_.size(); ++i) for (size_t j=0; j<h3d_5p_[i].size(); ++j) 
        for (size_t k=0; k<h3d_5p_[i][j].size(); ++k) for (size_t l=0; l<h3d_5p_[i][j][k].size(); ++l)
	  for (size_t m=0; m<h3d_5p_[i][j][k][l].size(); ++m) total += h3d_5p_[i][j][k][l][m]->GetNcells();
    }
    return total;
  }

private:
  //______________________________________________________________________
  //                       Multidraw functions
  
  TCanvas* custom_can_(TH1D* h, std::string canname, int gx = 1, int gy = 1,
        	       int histosize_x = 500, int histosize_y = 500,
        	       int mar_left = 90, int mar_right = 20, int mar_top = 20, int mar_bottom = 60,
		       int title_align = 33, float title_x = 0.99, float title_y = 0.99) {
    if (std::string(h->GetTitle()).size()>0||approval_) mar_top += 25;
    // Increase bottom margin to accomodate bin labels
    bool set_bin_labels = 0;
    if (!bin_labels_[0].empty()) {
      size_t maxlength=0;
      for (const auto& label : bin_labels_[0]) if (label.second.size()>maxlength) maxlength = label.second.size();
      if (maxlength>6) mar_bottom = maxlength*10;
      set_bin_labels = 1;
    }
    float titlefontsize = 32;
    float labelfontsize = 20;
    float yoffset_x = mar_left - titlefontsize - 4;
    float xoffset_y = mar_bottom - titlefontsize - 4;
    float zoffset_x = mar_right - titlefontsize - 4;
    float padsize_x = histosize_x + mar_left + mar_right;
    float padsize_y = histosize_y + mar_top + mar_bottom;
    float padsize = std::min(padsize_x, padsize_y);
    float padratio_yx = padsize_y/padsize_x > 1 ? 1 : padsize_y/padsize_x;
    float padratio_xy = padsize_x/padsize_y > 1 ? 1 : padsize_x/padsize_y;
    float xoffset = (xoffset_y/titlefontsize+0.5) * padratio_xy /1.6;
    float yoffset = (yoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    float zoffset = (zoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    float titlesize = titlefontsize/padsize;
    float labelsize = labelfontsize/padsize;
    if (std::string(h->GetTitle()).size()) {
      gStyle->SetOptTitle(1);
      gStyle->SetTitleH(titlefontsize/padsize);
      gStyle->SetTitleFontSize(titlesize*0.8);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleAlign(title_align);
      gStyle->SetTitleX(title_x);
      gStyle->SetTitleY(title_y);
    }
    h->SetTitleFont(42,"xyz");
    h->SetLabelFont(42,"xyz");
    h->SetTitleSize(titlesize,"xyz");
    h->SetLabelSize(labelsize,"xyz");
    if (set_bin_labels) h->GetXaxis()->SetLabelSize(1.5*labelsize);
    h->GetXaxis()->SetTitleOffset(xoffset);
    h->GetYaxis()->SetTitleOffset(yoffset);
    h->GetZaxis()->SetTitleOffset(zoffset);
    h->GetYaxis()->SetDecimals(1);
    h->GetZaxis()->SetDecimals(1);
    //gStyle->SetPadLeftMargin(mar_left/padsize_x);
    //gStyle->SetPadRightMargin(mar_right/padsize_x);
    //gStyle->SetPadTopMargin(mar_top/padsize_y);
    //gStyle->SetPadBottomMargin(mar_bottom/padsize_y);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleH(titlefontsize/padsize);
    gStyle->SetTitleFontSize(titlesize);
    TCanvas* canvas = new TCanvas(canname.c_str(), h->GetTitle(), padsize_x + 4, padsize_y + 26);
    TVirtualPad* pad = canvas->cd(1);
    pad->SetLeftMargin((float)mar_left/padsize_x);
    pad->SetRightMargin((float)mar_right/padsize_x);
    pad->SetTopMargin((float)mar_top/padsize_y);
    pad->SetBottomMargin((float)mar_bottom/padsize_y);
    canvas->SetGrid(gx,gy);
    if (logx_) canvas->SetLogx(1);
    if (log_) canvas->SetLogy(1);
    return canvas;
  }
  
  TCanvas* custom_can_(TH2D* h, std::string canname, int gx = 0, int gy = 0, 
		       int histosize_x = 500, int histosize_y = 500,
		       int mar_left = 80, int mar_right = 120, int mar_top = 20, int mar_bottom = 60,
		       int title_align = 33, float title_x = 0.99, float title_y = 0.99) {
    if (std::string(h->GetTitle()).size()>0||approval_) mar_top += 25;
    int titlefontsize = 32;
    int labelfontsize = 20;
    int pal_offset_x = 5;
    int pal_width_x = 25;
    int xoffset_y = mar_bottom - titlefontsize - 4;
    int yoffset_x = mar_left - titlefontsize - 4;
    int zoffset_x = mar_right - pal_offset_x - pal_width_x - titlefontsize;
    int padsize_x = histosize_x + mar_left + mar_right;
    int padsize_y = histosize_y + mar_top + mar_bottom;
    int padsize = ((padsize_x<=padsize_y) ? padsize_x : padsize_y);
    float padratio_yx = (float)padsize_y/padsize_x > 1 ? 1 : (float)padsize_y/padsize_x;
    float padratio_xy = (float)padsize_x/padsize_y > 1 ? 1 : (float)padsize_x/padsize_y;
    float xoffset = ((float)xoffset_y/titlefontsize+0.5) * padratio_xy /1.6;
    float yoffset = ((float)yoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    float zoffset = ((float)zoffset_x/titlefontsize+0.5) * padratio_yx /1.6;
    float titlesize = (float)titlefontsize/padsize;
    float labelsize = (float)labelfontsize/padsize;
    h->SetTitleFont(42,"xyz");
    h->SetLabelFont(42,"xyz");
    h->SetTitleSize(titlesize,"xyz");
    h->SetLabelSize(labelsize,"xyz");
    h->GetXaxis()->SetTitleOffset(xoffset);
    h->GetYaxis()->SetTitleOffset(yoffset);
    h->GetZaxis()->SetTitleOffset(zoffset);
    h->GetZaxis()->RotateTitle(1);
    h->GetYaxis()->SetDecimals(1);
    h->GetZaxis()->SetDecimals(1);
    if (histosize_y<250) h->GetZaxis()->SetNdivisions(505);
    if (std::string(h->GetTitle()).size()) {
      gStyle->SetOptTitle(1);
      gStyle->SetTitleH(titlefontsize/padsize);
      gStyle->SetTitleFontSize(titlesize*0.8);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleAlign(title_align);
      gStyle->SetTitleX(title_x);
      gStyle->SetTitleY(title_y);
    }
    //gStyle->SetPadLeftMargin((float)mar_left/padsize_x);
    //gStyle->SetPadRightMargin((float)mar_right/padsize_x);
    //gStyle->SetPadTopMargin((float)mar_top/padsize_y);
    //gStyle->SetPadBottomMargin((float)mar_bottom/padsize_y);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleH(titlefontsize/padsize);
    gStyle->SetTitleFontSize(titlesize);
    TCanvas* canvas = new TCanvas(canname.c_str(), h->GetTitle(), padsize_x + 4, padsize_y + 26);
    TVirtualPad* pad = canvas->cd(1);
    pad->SetLeftMargin((float)mar_left/padsize_x);
    pad->SetRightMargin((float)mar_right/padsize_x);
    pad->SetTopMargin((float)mar_top/padsize_y);
    pad->SetBottomMargin((float)mar_bottom/padsize_y);
    canvas->SetGrid(gx,gy);
    if (norm_&&h->Integral()>0) h = (TH2D*)h->DrawNormalized(draw_.c_str());
    else h->Draw(draw_.c_str());
    if (h->Integral()>0&&draw_=="COLZ") {
      gPad->Update();
      TPaletteAxis* palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
      if (palette) {
	palette->SetX1NDC(1 - (float)(mar_right - pal_offset_x)/padsize_x);
	palette->SetX2NDC(1 - (float)(mar_right - pal_offset_x - pal_width_x)/padsize_x);
	palette->SetY1NDC((float)mar_bottom/padsize_y);
	palette->SetY2NDC(1 - (float)mar_top/padsize_y);
      }
    }
    gStyle->SetOptTitle(0);
    if (logx_) canvas->SetLogx(1);
    if (log_) canvas->SetLogz(1);
    return canvas;
  }

  void add_labels_(TCanvas* c) {
    double xmin = ((TFrame*)c->GetListOfPrimitives()->At(0))->GetX1();
    double xmax = ((TFrame*)c->GetListOfPrimitives()->At(0))->GetX2();
    double ymin = ((TFrame*)c->GetListOfPrimitives()->At(0))->GetY1();
    double ymax = ((TFrame*)c->GetListOfPrimitives()->At(0))->GetY2();
    era_and_prelim_lat_(xmin, xmax, ymin, ymax);
  }
  void add_labels_(TH1D* h) {
    double xmin = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
    double xmax = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast());
    double ymin = h->GetMinimum();
    double ymax = h->GetMaximum();
    era_and_prelim_lat_(xmin, xmax, ymin, ymax);
  }
  void add_labels_(TH2D* h) {
    double xmin = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
    double xmax = h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast());
    double ymin = h->GetYaxis()->GetBinLowEdge(h->GetYaxis()->GetFirst());
    double ymax = h->GetYaxis()->GetBinUpEdge(h->GetYaxis()->GetLast());
    era_and_prelim_lat_(xmin, xmax, ymin, ymax);
  }
  void era_and_prelim_lat_(double xmin, double xmax, double ymin, double ymax, bool in=0) {
    int app = approval_/10;
    if (app) {
      // Latex example: #font[22]{Times bold} and #font[12]{Times Italic}
      std::string text = "";
      if (app==1) text = "CMS #scale[0.7]{#font[52]{Work in progress}}";
      if (app==2) text = "CMS #scale[0.7]{#font[52]{Preliminary}}";
      if (app==3) text = "CMS";
      if (app==4) text = "#scale[0.8]{CMS Simulation }#scale[0.6]{#font[52]{Work in progress}}";
      if (app==5) text = "#scale[0.8]{CMS Simulation }#scale[0.6]{#font[52]{Preliminary}}";
      if (app==6) text = "CMS Simulation";
      if (app==7) text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
      if (app==8) text = "CMS #scale[0.7]{#font[52]{Preliminary 2016}}";
      if (app==9) text = "CMS #scale[0.7]{#font[52]{Preliminary 2017}}";
      double x = in ? xmin+(xmax-xmin)/20.0 : xmin;
      double y = in ? ymax-(ymax-ymin)/10.0 : ymax+(ymax-ymin)/40.0;
      if (log_&&ymin>0) y = in ? std::exp(std::log(ymax)-(std::log(ymax)-std::log(ymin))/10.0) : 
	std::exp(std::log(ymax)+(std::log(ymax)-std::log(ymin))/40);
      TLatex* cms_lat = new TLatex(x, y, text.c_str()); 
      cms_lat->SetLineWidth(2);
      cms_lat->Draw();
    }
    int era = approval_%10;
    if (era) {
      std::string text = "";
      if (era==1) text = "#sqrt{s}=7 TeV";
      if (era==2) text = "#sqrt{s}=8 TeV";
      if (era==3) text = "#sqrt{s}=13 TeV";
      if (era==4) text = "Run 2, #sqrt{s}=13 TeV";
      if (era==5) text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
      if (era==6) text = "#scale[0.65]{W ana, 35.9 fb^{-1} (13 TeV)}";
      if (era==7) text = "#scale[0.65]{Top ana, 35.9 fb^{-1} (13 TeV)}";
      double y = log_&&ymin>0 ? std::exp(std::log(ymax)+(std::log(ymax)-std::log(ymin))/25.0) : ymax+(ymax-ymin)/25.0;
      TLatex* era_lat = new TLatex(xmax, y, text.c_str());
      era_lat->SetTextAlign(32);
      era_lat->SetTextSize(0.04);
      era_lat->SetTextFont(42);
      era_lat->SetLineWidth(2);
      era_lat->Draw();
    }
    gPad->Update();
  }

  void draw_mr_bins(TH1D *h, double ymin, double ymax, bool combine_bins=true,
		    std::vector<double> mrbins = { 0.8, 1.0, 1.2, 1.6, 2.0, 4.0 },
		    std::vector<double> r2bins = { 0.08, 0.12, 0.16, 0.24, 0.4, 2.0 }) {
    for (size_t i=0; i<mrbins.size(); ++i) {
      int bin1 = i*5+1, bin2 = i*5+5;
      if (combine_bins) {
	if (i==3) {
	  bin1 = 16;
	  bin2 = 19;
	} else if (i==4) {
	  bin1 = 20;
	  bin2 = 22;
	}
      }
      double maxcont = -9999;
      for (int binx=bin1; binx<=bin2; ++binx) {
	if (h->GetBinContent(binx)>maxcont) 
	  maxcont = h->GetBinContent(binx);
      }
      std::vector<double> y2 = {maxcont+(ymax-ymin)*0.05, maxcont+(ymax-ymin)*0.125};
      if (ymin!=0) y2 = {maxcont*std::pow(ymax/ymin,0.05), maxcont*std::pow(ymax/ymin,0.125)};
      double x = i*5;
      if (i==4 && combine_bins) x = x-1;
      if (i!=0) {
	TLine *line = new TLine(x,ymin,x,y2[0]);
	line->SetLineStyle(2);
	line->Draw();
      }
      x = 2.0 + i*5;
      if (i==3 && combine_bins) x = x-0.5;
      if (i==4 && combine_bins) x = x-2;
      if (i==0) {
	TLatex* bin_lat = new TLatex(x, y2[1], "M_{R} (TeV)");
	bin_lat->SetTextAlign(22);
	bin_lat->SetTextSize(0.04);
	bin_lat->Draw();
      }
      std::stringstream ss;
      ss<<"["<<std::setprecision(2)<<mrbins[i]<<", "<<mrbins[i+1]<<"]";
      TLatex* lat = new TLatex(x, y2[0], ss.str().c_str());
      lat->SetTextAlign(22);
      lat->SetTextSize(0.04);
      lat->Draw();
    }
  }

  void draw_ak8pt_bins(TH1D *h, bool high,
		       std::vector<double> HTB = {400, 500, 600, 700, 750, 800, 850, 900, 950, 1000, 1500, 10000},
		       std::vector<double> PtB = {200, 300, 400, 450, 500, 550, 600, 1000, 10000}) {
    size_t ilow  = high ? PtB.size()/2 : 0;
    size_t ihigh = high ? PtB.size()-1 : PtB.size()/2;
    for (size_t i=ilow; i<ihigh; ++i) {
      double x = (i-ilow)*(HTB.size()-1);
      if (i!=ilow) {
	TLine *line = new TLine(x,0,x,1);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->Draw();
      }
      x = (HTB.size()-1)/4.0 + x + (i>5)*0.5+(i>6)*0.5;
      if (i==ilow) {
	TLatex* bin_lat = new TLatex(x, 0.55, "p_{T} (GeV)");
	bin_lat->SetTextAlign(22);
	bin_lat->SetTextSize(0.04);
	bin_lat->Draw();
      }
      std::stringstream ss;
      ss<<"["<<PtB[i]<<", "<<PtB[i+1]<<"]";
      TLatex* lat = new TLatex(x, 0.475, ss.str().c_str());
      lat->SetTextAlign(22);
      lat->SetTextSize(0.04);
      lat->Draw();
    }
  }
  
  std::vector<int> string_to_vector_(std::string val) {
    std::vector<int> vec;
    // Reading values from string separated by commas, (ranges x-y can be specified too)
    // eg. "3-6,20,22,24-27,11," - add comma at the end of the string
    std::stringstream ss(val);
    while (ss) {
      while (ss && !isdigit(ss.peek ())) ss.get ();
      if (! ss)	break;
      int i = 0;
      int j = 0;
      ss >> i;
      while (ss && !ispunct(ss.peek ())) ss.get ();
      char a; ss>>a;
      if (a=='-') {
        while (ss && !isdigit(ss.peek ())) ss.get ();
        ss >> j;
      } else j=i;
      for (int n=i; n<=j; n++) vec.push_back(n);
    }
    return vec;
  }
  
  void set_stat_(TH1D* h, Color_t col, int i) {
    if (i<7) {
      gPad->Update();
      TPaveStats* stats = (TPaveStats*)h->FindObject("stats");
      stats->SetLineColor(col);
      stats->SetTextColor(col);
      stats->SetX1NDC(0.74);
      stats->SetX2NDC(0.94);
      stats->SetY1NDC(0.8-i*0.11);
      stats->SetY2NDC(0.9-i*0.11);
    }
  }

  void add_integ_(std::vector<TLegendEntry*>& legentries, TH1D* h) {
    if (addint_) {
      char integ[100]; sprintf(integ, "#font[82]{%7.1f}", h->Integral());
      legentries.push_back(new TLegendEntry((TObject*)0, integ, ""));
    }
  }
  
  void multidraw_with_legend_(size_t skip, std::vector<TH1D*>& hvec, std::vector<std::string> pf, std::string colz,
			      std::string legtitle="", float x1=0.15, float y2=0.9) {
    // Draw multiple histograms, set their marker/line color/style
    // Then Draw legend for all histo with titles from a postfix
    std::vector<int> col = string_to_vector_(colz);
    // full  : circle, square, tri-up, tri-down, star
    // hollow: circle, square, tri-up, diam, star
    std::vector<Int_t> marker = { 20, 21, 22, 23, 29, 24, 25, 26, 27, 30}; 
    size_t nrow = 0, ncol = (1+(twocol_>0))*(1+addint_);
    // Set styles for histos
    for (size_t i=skip; i<hvec.size(); i++) if (hvec[i]->GetEntries()>0) {
      if (!stack_) {
	if (draw_.find("P")!=std::string::npos) {
	  hvec[i]->SetLineColor((Color_t)col[i-(keep_color_?skip:0)]); 
	  hvec[i]->SetMarkerColor((Color_t)col[i-(keep_color_?skip:0)]); 
	  hvec[i]->SetMarkerStyle(marker[(i-(keep_color_?skip:0))%10]); 
	} else {
	  hvec[i]->SetLineColor((Color_t)col[i-(keep_color_?skip:0)]);
	  hvec[i]->SetLineWidth(2);
	}
      }
      if (stat_) hvec[i]->SetStats(1);
      ++nrow;
    }
    if (twocol_) nrow = std::max(twocol_/10, twocol_%10);
    if (legtitle.size()>0) ++nrow;
    std::vector<TLegendEntry*> legentries;
    std::string same = (stat_ ? "SAMES" : "SAME") + draw_;
    std::vector<TH1D*> vh;
    std::vector<std::string> vlegtext;
    std::vector<TGraphAsymmErrors*> graphs;
    bool draw_axis=1;
    size_t n_nonstack_drawn = 0;
    for (size_t i=skip; i<hvec.size() ; ++i) if (hvec[i]->GetEntries()>0) {
      std::stringstream colored_text;
      colored_text<<"#color["<<(Color_t)col[i-(keep_color_?skip:0)]<<"]{"<<pf[i]<<"}";
      if (stack_) {
	// Set special styles for Stack histos
	if (i<n_nostack_) {
	  if (i==0) {
	    // Data
	    hvec[i]->SetFillColor(0);
	    //hvec[i]->SetLineColor(0);
	    hvec[i]->SetMarkerColor((Color_t)col[i-(keep_color_?skip:0)]); 
	    hvec[i]->SetMarkerStyle(marker[(i-(keep_color_?skip:0))%10]); 
	    hvec[i]->Draw("PE1");
	    legentries.push_back(new TLegendEntry(hvec[i], colored_text.str().c_str(), "P"));
	    n_nonstack_drawn++;
	  } else {
	    // Signal
	    hvec[i]->SetLineStyle(2);
	    hvec[i]->SetLineWidth(2);
	    hvec[i]->SetLineColor((Color_t)col[i-(keep_color_?skip:0)]);
	    if (skip!=0&&i==skip) hvec[i]->Draw("HIST");
	    legentries.push_back(new TLegendEntry(hvec[i], colored_text.str().c_str(), "L"));
	    n_nonstack_drawn++;
	  }
	  add_integ_(legentries,hvec[i]);
	} else {
	  hvec[i]->SetLineColor((Color_t)col[i-(keep_color_?skip:0)]);
	  hvec[i]->SetFillColor((Color_t)col[i-(keep_color_?skip:0)]);
	  vh.push_back(hvec[i]);
	  vlegtext.push_back(colored_text.str());
	}
      } else {
	// Draw ordinary histos and legend
	if (norm_&&hvec[i]->Integral()>0) {
	  hvec[i]->DrawNormalized((i==skip) ? draw_.c_str() : same.c_str());
	  legentries.push_back(new TLegendEntry(hvec[i], colored_text.str().c_str(), draw_.find("P")!=std::string::npos ? "P" : "L"));
	  add_integ_(legentries, hvec[i]);
	} else if (plot_asymm_err_) {
	  hvec[i]->SetLineWidth(2);
	  //if (name_=="HitEfficiency_vs_BiasVoltage") std::cout<<hvec[i]->GetName()<<std::endl;
	  TGraphAsymmErrors* tgae = asym_(hvec[i], mother_2d_[hvec[i]]->ProjectionX());
	  if (tgae->GetN()) {
	    if (draw_axis) {
	      tgae->Draw("AP");
	      // X-axis ranges are messed up by default, so need to draw first
	      //  tgae->GetXaxis()->SetRangeUser(hvec[i]->GetXaxis()->GetBinLowEdge(hvec[i]->GetXaxis()->GetFirst()),
	      //    			     hvec[i]->GetXaxis()->GetBinUpEdge(hvec[i]->GetXaxis()->GetLast()));
	      //  tgae->Draw("AP");	      
	      draw_axis=0;
	    } else tgae->Draw("SAMEP");
	    legentries.push_back(new TLegendEntry(hvec[i], colored_text.str().c_str(), draw_.find("P")!=std::string::npos ? "P" : "L"));
	    add_integ_(legentries, hvec[i]);
	  }
	  //tgae->Draw((i==skip) ? "AP" : "SAMEP");
	  if (i==skip) asym_labels_(hvec[i], tgae, 0);
	  TGraphAsymmErrors* tgae_allbins = asym_(hvec[i], mother_2d_[hvec[i]]->ProjectionX(), true);
	  graphs.push_back(tgae_allbins);
	  if (ratio_&&graphs.size()==2) {
	    TGraphAsymmErrors* ratio = get_tgae_ratio_(graphs[0], graphs[1]);
	    ratio->SetMarkerColor(417);
	    ratio->SetMarkerStyle(22);
	    ratio->SetLineColor(417);
	    ratio->Draw("SAMEP");
	    legentries.push_back(new TLegendEntry(ratio, "#color[417]{Ratio}", "P"));
	    if (addint_) legentries.push_back(new TLegendEntry((TObject*)0, "", ""));
	  }
	} else {
	  hvec[i]->Draw((i==skip) ? draw_.c_str() : same.c_str());
	  legentries.push_back(new TLegendEntry(hvec[i], colored_text.str().c_str(), draw_.find("P")!=std::string::npos ? "P" : "L"));
	  add_integ_(legentries, hvec[i]);
	}
      }
      if (stat_) set_stat_(hvec[i], (Color_t)col[i-(keep_color_?skip:0)], i-skip);
    }
    if (stack_) {
      // Draw stacked histos and legend
      THStack *s = new THStack("s","");
      std::vector<TH1D*> vh_max;
      // Reorder histos based on integral (smallest drawn on bottom)
      if (!keep_order_) {
	while (vh.size()) {
	  double max_int = -1; size_t imax=-1;
          for (size_t i=0; i<vh.size(); ++i) {
            if (vh[i]->Integral()>max_int) {
              max_int = vh[i]->Integral();
              imax = i;
            }
          }
          if (imax==(size_t)-1) imax = 0;
          vh_max.push_back(vh[imax]);
          legentries.push_back(new TLegendEntry(vh[imax], vlegtext[imax].c_str(), "F"));
	  add_integ_(legentries, vh[imax]);
          vh.erase(vh.begin()+imax);
          vlegtext.erase(vlegtext.begin()+imax);
        }
      } else {
	// or keep original order
	for (int i=0, n=vh.size(); i<n; ++i) {
	  vh_max.push_back(vh[i]);
	  legentries.push_back(new TLegendEntry(vh[i], vlegtext[i].c_str(), "F"));
	  add_integ_(legentries, vh[i]);
	}
      }
      for (int i=(int)vh_max.size()-1; i>=0; --i) s->Add(vh_max[i]);
      s->SetTitle(std::to_string(vh_max.size()).c_str());
      s->Draw(n_nostack_ ? "SAMEHIST" : "HIST");
      for (int i=n_nostack_-1; i>=(int)skip; --i)
	hvec[i]->Draw(i==0 ? "SAMEPE1" : "SAMEHIST");
    } else if (ratio_&&!plot_asymm_err_) {
      // Additionally draw ratio of first 2 plots
      TH1D* ratio = (TH1D*)hvec[0]->Clone();
      TH1D* den_err = (TH1D*)hvec[1]->Clone();
      // Instead of Divide(), scale the error of num, and plot error of den around 1
      //ratio->Divide(hvec[1]);
      for (int bin=1; bin<=ratio->GetNbinsX(); ++bin) {
	ratio  ->SetBinContent(bin, hvec[0]->GetBinContent(bin)/hvec[1]->GetBinContent(bin));
	ratio  ->SetBinError  (bin, hvec[0]->GetBinError(bin)  /hvec[1]->GetBinContent(bin));
	den_err->SetBinContent(bin, 1);
	den_err->SetBinError  (bin, hvec[1]->GetBinError(bin)  /hvec[1]->GetBinContent(bin));
      }
      ratio->SetMarkerColor(417);
      ratio->SetMarkerStyle(22);
      ratio->SetLineColor(417);
      den_err->Draw("SAME E2");
      ratio->Draw(same.c_str());
      legentries.push_back(new TLegendEntry(ratio, "#color[417]{Ratio}", "P"));
      if (addint_) legentries.push_back(new TLegendEntry((TObject*)0, "", ""));
    }
    // Legend
    float textsize = ratio_ ? 0.028 : 0.04;
    float x2 = x1 + 0.2 * (1+(twocol_>0))*(1+addint_*0.3);
    float y1 = y2 - nrow * (textsize*1.2); // 0.1 + textsize + 0.1 (fEntrySeparation = 0.1)
    TLegend *leg = new TLegend(x1,y1,x2,y2,legtitle.c_str());
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(textsize);
    leg->SetNColumns(ncol);
    // Reorganize legend entries into two columns before drawing
    if (twocol_) {
      size_t nrow1 = twocol_/10, nrow2 = twocol_%10;
      if (n_nonstack_drawn>0) nrow1 = std::min(n_nonstack_drawn, nrow1);
      for (size_t irow=0, n = std::max(nrow1,nrow2); irow<n; ++irow) {
	// add left column(s) first
	size_t i = irow * (1+addint_);
	if (i+addint_<legentries.size() && irow<nrow1) {
	  leg->AddEntry(legentries[i]->GetObject(), legentries[i]->GetLabel(), legentries[i]->GetOption());
	  if (addint_) leg->AddEntry(legentries[i+1]->GetObject(), legentries[i+1]->GetLabel(), legentries[i+1]->GetOption());
	} else {
	  leg->AddEntry((TObject*)0,"","");
	  if (addint_) leg->AddEntry((TObject*)0,"","");
	}
	// then right column(s)
	i = (nrow1 + irow) * (1+addint_);
	if (i+addint_<legentries.size() && irow<nrow2) {
	  leg->AddEntry(legentries[i]->GetObject(), legentries[i]->GetLabel(), legentries[i]->GetOption());
	  if (addint_) leg->AddEntry(legentries[i+1]->GetObject(), legentries[i+1]->GetLabel(), legentries[i+1]->GetOption());
	} else {
	  leg->AddEntry((TObject*)0,"","");
	  if (addint_) leg->AddEntry((TObject*)0,"","");
	}
      }
    } else for (auto e : legentries) leg->AddEntry(e->GetObject(), e->GetLabel(), e->GetOption());
    leg->Draw("SAME");
  }
  
  void add_stack_ratio_plot_(TCanvas*& c, double xmin, double xmax , bool remove=false, bool debug = 0) {
    // Canvas division sizes
    float mar_top    = 45;
    float y1         = 365;
    float mid2       = 10;
    float y2         = 115;
    float mar_bottom = 60;
    // Allow larger space for TLatex bin labels
    bool set_bin_labels = false;
    if (!bin_labels_[0].empty()) {
      size_t maxlength=0;
      for (const auto& label : bin_labels_[0]) if (label.second.size()>maxlength) maxlength = label.second.size();
      if (maxlength>6) mar_bottom = maxlength*10;
      set_bin_labels = true;
    }
    float mar_left   = 90;
    float x          = 500;
    float mar_right  = 20;
    float x_can    = mar_left+x+mar_right;
    float y_can    = mar_top+y1+mid2*2+y2+mar_bottom;
    float padsize1 = std::min(mar_top+y1+mid2,    x_can);
    float padsize2 = std::min(mar_bottom+y2+mid2, x_can);
    float labelfontsize = 20;
    float titlefontsize = 32;
    float leg_y2 = 0.87;
    bool ok = 0;
    if (debug) std::cout<<"Start debugging: "<<c->GetName()<<std::endl;
    if (debug) std::cout<<"ok"<<std::endl;
    if (c->GetListOfPrimitives()->GetEntries()>2) {
      // Histos
      if (debug) std::cout<<"ok1"<<std::endl;
      TH1D* Data = (TH1D*)c->GetListOfPrimitives()->At(1);
      if (debug) std::cout<<"ok1"<<std::endl;
      THStack* MCstack = (THStack*)c->GetListOfPrimitives()->At(2);
      if (debug) std::cout<<"ok1"<<std::endl;
      if (std::string(MCstack->GetTitle())!="0") {
        TH1D* ratio = (TH1D*)Data->Clone();
	if (debug) std::cout<<"ok2"<<std::endl;
        TH1D* mc_sum = (TH1D*)MCstack->GetHists()->At(0)->Clone();
        TH2D* mc_sum_syst = 0;
	if (syst_) mc_sum_syst = (TH2D*)mother_2d_[(TH1D*)MCstack->GetHists()->At(0)]->Clone();
	if (debug) std::cout<<"ok2"<<std::endl;
        for (int iStack=1; iStack<MCstack->GetHists()->GetEntries(); ++iStack) {
	  TH1D* h = (TH1D*)MCstack->GetHists()->At(iStack);
	  mc_sum->Add((TH1D*)h->Clone());
	  if (syst_) mc_sum_syst->Add(mother_2d_[h]);
	}
	if (debug) std::cout<<"ok2"<<std::endl;
	TH1D* den_stat_err = (TH1D*)mc_sum->Clone("den_stat_err");
	TH1D* den_all_up_err = (TH1D*)mc_sum->Clone("den_all_up_err");
	TH1D* den_all_down_err = (TH1D*)mc_sum->Clone("den_all_down_err");
	if (debug) std::cout<<"ok2"<<std::endl;
	// Instead of Divide(), scale the error of num, and plot error of den around 1
        ratio->Divide(mc_sum);
	if (debug) std::cout<<"ok2"<<std::endl;
	for (int bin=1; bin<=ratio->GetNbinsX(); ++bin) {
	  if (mc_sum->GetBinContent(bin)!=0) {
	    ratio  ->SetBinContent(bin, Data->GetBinContent(bin)/mc_sum->GetBinContent(bin));
	    ratio  ->SetBinError  (bin, Data->GetBinError(bin)  /mc_sum->GetBinContent(bin));
	    den_stat_err->SetBinContent(bin, 1);
	    den_stat_err->SetBinError  (bin, mc_sum->GetBinError(bin)  /mc_sum->GetBinContent(bin));
	    if (syst_) {
	      // Loop on systematics bins and look for std dev of up/down variations
	      double var_up = 0, var_down = 0, n_up = 0, n_down = 0;
	      double mean = mc_sum_syst->GetBinContent(bin,1);
	      for (int biny=2; biny<=mc_sum_syst->GetNbinsY(); ++biny) /*if (biny<12||biny>13)*/ {
		double cont = mc_sum_syst->GetBinContent(bin,biny);
		if (cont>mean) {
		  var_up += (cont-mean)*(cont-mean);
		  n_up++;
		} else if (cont<mean) {
		  var_down += (cont-mean)*(cont-mean);
		  n_down++;
		}
	      }
	      double up_err = n_up>0 ? std::sqrt(var_up) : 0;
	      double down_err = n_down>0 ? std::sqrt(var_down) : 0;
	      // add stat error
	      up_err   = std::sqrt(up_err*up_err + mc_sum->GetBinError(bin)*mc_sum->GetBinError(bin));
	      down_err = std::sqrt(down_err*down_err + mc_sum->GetBinError(bin)*mc_sum->GetBinError(bin));
	      up_err   /= mc_sum->GetBinContent(bin);
	      down_err /= mc_sum->GetBinContent(bin);
	      // Maximize the error to the visible range*2
	      up_err   = std::min(up_err,  1.999);
	      down_err = std::min(down_err,1.999);
	      den_all_up_err  ->SetBinContent(bin, 1+up_err/2);
	      den_all_down_err->SetBinContent(bin, 1-down_err/2);
	      den_all_up_err  ->SetBinError(bin, up_err/2);
	      den_all_down_err->SetBinError(bin, down_err/2);
	    }
	  } else {
	    ratio  ->SetBinContent(bin, 0);
	    ratio  ->SetBinError  (bin, 0);
	    den_stat_err->SetBinContent(bin, 1);
	    den_stat_err->SetBinError  (bin, 0);
	  }
	}
	if (debug) std::cout<<"ok2"<<std::endl;
	den_stat_err->SetFillColor(1);
	den_stat_err->SetFillStyle(3004);
	if (syst_) {
	  den_all_up_err->SetFillColor(kGray);
	  den_all_up_err->SetFillStyle(1001);
	  den_all_down_err->SetFillColor(kGray);
	  den_all_down_err->SetFillStyle(1001);
	}
        // Legend
        TLegend* leg = 0;
        // Remove Non-Data non-stack plots (eg. signal)
        std::vector<TH1D*> rest;
        // indices:
        // 0: TFrame, 1: Data, 2: stack, 3: Data again, 4+: (signals), 4+nsig: Legend
	if (debug) std::cout<<"ok2"<<std::endl;
        for (int i=3; i<c->GetListOfPrimitives()->GetEntries(); ++i) {
          std::string prim_name = c->GetListOfPrimitives()->At(i)->GetName();
          if (prim_name=="TPave") {
	    leg = (TLegend*)c->GetListOfPrimitives()->At(i);
          } else if (!remove&&prim_name!=Data->GetName()&&prim_name.size()>0)
            rest.push_back((TH1D*)c->GetListOfPrimitives()->At(i));
        }
	if (debug) std::cout<<"ok2"<<std::endl;
        if (remove) {
          for (int i=leg->GetListOfPrimitives()->GetEntries(); i>0; --i) {
            bool match = 0;
            std::string name = ((TLegendEntry*)leg->GetListOfPrimitives()->At(i))->GetObject()->GetName();
            for (int j=0; j<MCstack->GetHists()->GetEntries(); ++j) if (name==MCstack->GetHists()->At(j)->GetName()) match = 1;
            if (!match) leg->GetListOfPrimitives()->RemoveAt(i);
          }
        }
	if (debug) std::cout<<"ok2"<<std::endl;
        // Styles
	float heightratio1 = padsize1/y_can;
        Data ->SetTitleSize  (Data->GetYaxis()->GetTitleSize()  /heightratio1,"y");
        Data ->SetTitleOffset(Data->GetYaxis()->GetTitleOffset()*heightratio1,"y");
        Data ->SetLabelSize(labelfontsize/padsize1,"xyz");
        ratio->SetLabelSize(labelfontsize/padsize2,"xyz");
	if (set_bin_labels) ratio->SetLabelSize(labelfontsize*1.5/padsize2,"x");
        ratio->SetTitleSize(titlefontsize/padsize2,"xyz");
        ratio->GetYaxis()->SetRangeUser(0,2);
        ratio->GetYaxis()->SetNdivisions(305);
        ratio->GetYaxis()->SetTitle("Data/MC");
	if (debug) std::cout<<"ok2"<<std::endl;
	float heightratio2 = padsize2/y_can;
        ratio->SetTitleOffset(ratio->GetYaxis()->GetTitleOffset()*heightratio2,"y");
        ratio->SetTitle("");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerColor(1);
        ratio->SetLineColor(1);
	if (debug) std::cout<<"ok2"<<std::endl;
        // New Canvas
	float left_mar = c->GetLeftMargin(), right_mar = c->GetRightMargin();
        bool logScale = c->GetLogy();
        c = new TCanvas((std::string(c->GetName())+"_Ratio").c_str(), c->GetTitle(), x_can+4,y_can+26); // 600, 600
        c->Divide(1,2);
	if (debug) std::cout<<"ok2"<<std::endl;
        // Pad 1 (x: 90+500+20 x y: 45+350+10)
        TVirtualPad* p = c->cd(1);
	if (debug) std::cout<<"ok2"<<std::endl;
        p->SetPad(0,padsize2/y_can,1,1);
        p->SetTopMargin(mar_top/(mar_top+y1+mid2));
        p->SetBottomMargin(0);
	p->SetLeftMargin(left_mar);
	p->SetRightMargin(right_mar);
	if (debug) std::cout<<"ok2"<<std::endl;
        if (logScale) p->SetLogy(1);
        Data->Draw("PE1");
        MCstack->Draw("SAMEHIST");
	if (debug) std::cout<<"ok2"<<std::endl;
        for (size_t i=0; i<rest.size(); ++i) rest[i]->Draw("SAMEHIST");
        Data->Draw("SAMEPE1");
	if (debug) std::cout<<"ok3"<<std::endl;
	// Move/Resize legend to fit divided canvas better
	leg->SetTextSize(0.04);
	leg->SetY1(leg_y2 - leg->GetNRows()*0.04);
	leg->SetY2(leg_y2);
        leg->Draw("SAME");
	if (debug) std::cout<<"ok3"<<std::endl;
	add_labels_(Data);
	if (debug) std::cout<<"ok3"<<std::endl;
        gPad->Update();
	if (debug) std::cout<<"ok3"<<std::endl;
        // Pad 2 (x: 90+500+20 x y: 60+150+10)
        p = c->cd(2);
        p->SetGrid(0,1);
        p->SetPad(0,0,1,padsize2/y_can);
        p->SetTopMargin(mid2/padsize2);
        p->SetBottomMargin(mar_bottom/padsize2);
	p->SetLeftMargin(left_mar);
	p->SetRightMargin(right_mar);
	if (debug) std::cout<<"ok3"<<std::endl;
        ratio->Draw("PE1");
	if (syst_) {
	  den_all_up_err->Draw("SAME E2");
	  den_all_down_err->Draw("SAME E2");
	}
	den_stat_err->Draw("SAME E2");
        ratio->Draw("SAME PE1");
	if (debug) std::cout<<"ok3"<<std::endl;
	if (xmin==xmax) {
	  xmin = ratio->GetYaxis()->GetXmin();
	  xmax = ratio->GetYaxis()->GetXmax();
	}
	if (debug) std::cout<<"ok3"<<std::endl;
	TLine* l = new TLine(xmin, 1, xmax, 1);
	l->SetLineWidth(2);
	//l->SetLineColor(2);
	l->SetLineStyle(2);
	l->Draw();
	if (debug) std::cout<<"ok3"<<std::endl;
	ok = 1;
	gPad->Update();
      }
    }
    if (ok) {
      if (debug) std::cout<<"ok4"<<std::endl;
      TPad* pad = (TPad*)c->GetListOfPrimitives()->At(0);
      // Primitives order in Divided TPad: TFrame, data, stack, signals (n_nostack_-1), data, legend, cms, era
      if (debug) std::cout<<"ok5"<<std::endl;
      if (pad->GetListOfPrimitives()->GetEntries()>(int)n_nostack_+4) {
        // Resize CMS label and era text
        int index = n_nostack_+4;
        if (debug) std::cout<<"ok5"<<std::endl;
        if (approval_/10>0) {
          TLatex* cms_lat = (TLatex*)pad->GetListOfPrimitives()->At(index++);
          cms_lat->SetTextSize(cms_lat->GetTextSize()*(y1+y2)/y1);
        }
        if (debug) std::cout<<"ok5"<<std::endl;
        if (approval_%10>0) {
          TLatex* era_lat = (TLatex*)pad->GetListOfPrimitives()->At(index);
          era_lat->SetTextSize(era_lat->GetTextSize()*(y1+y2)/y1);
        }
        if (debug) std::cout<<"ok5"<<std::endl;
      }
    }
  }
  
  // __________________________________________________________
  //     Draw Asymmetric errors for Efficiencies/Fractions
  
  void asym1d_draw_(TH1D *h, std::string opt = "AP", int hor_vert_up_down = 0) { 
    TGraphAsymmErrors* tgae = asym_(h, mother_2d_[h]->ProjectionX());
    tgae->Draw(opt.c_str());
    if (opt == "AP") asym_labels_(h, tgae, hor_vert_up_down);
  }
  
  std::vector<double> bincoordx_;
  std::vector<std::string> binlabels_;
  std::map<const char*, TH1F*> tgae_hist_;
  
  TGraphAsymmErrors* asym_(TH1D* eff, TH1D* den, bool allbins = false) {
    bincoordx_.clear();
    binlabels_.clear();
    //TGraphAsymmErrors* tgae = new TGraphAsymmErrors(eff);
    // Get ranges
    double xlow  = eff->GetXaxis()->GetBinLowEdge(eff->GetXaxis()->GetFirst());
    double xhigh = eff->GetXaxis()->GetBinUpEdge(eff->GetXaxis()->GetLast());
    double ylow  = eff->GetMinimum();
    double yhigh = eff->GetMaximum();
    // Check if we need to add dummy points in order to correctly set axis ranges
    bool add_low = 0, add_high = 0;
    if (!allbins) {
      int first = -1, last = -1;
      for (int i=1; i<=den->GetNbinsX(); ++i) if (den->GetBinContent(i)) {
	if (den->GetXaxis()->GetBinLowEdge(i)>=xlow&&den->GetXaxis()->GetBinUpEdge(i)<xlow) first = i;
	if (den->GetXaxis()->GetBinLowEdge(i)>xhigh&&den->GetXaxis()->GetBinUpEdge(i)<=xhigh) last = i;	
      }
      if (first!=-1) add_low  = 0;
      if (last !=-1) add_high = 0;
    }
    int n = 0;
    if (allbins) n = den->GetNbinsX();
    else for (int i=0; i<den->GetNbinsX(); ++i) if (den->GetBinContent(i+1)>0) n++;
    TGraphAsymmErrors* tgae = new TGraphAsymmErrors(n+add_low+add_high); // add 2 dummy points (see below)
    int m = add_low;
    const int Method = 1;
    // Method 1
    // Calculate the asymmetric Wilson Score Interval
    if (Method==1) {
      double z = 1; // 1 Sigma confidence
      for (int i=0; i<den->GetNbinsX(); ++i) {
        double x = eff->GetBinCenter(i+1);
        double y = eff->GetBinContent(i+1);
        double denum = den->GetBinContent(i+1);
        double y_center = denum>0 ? (y+(z*z/(2*denum))) / (1.0 + (z*z/denum)) : 0;
        double y_halfwidth = denum>0 ? z*sqrt( y*(1.0-y)/denum + (z*z/(4*denum*denum)) ) / (1.0 + (z*z/denum)) : 0;
        double x_halfwidth = eff->GetXaxis()->GetBinWidth(i+1)/2;
        if (allbins||denum>0) {
          tgae->SetPoint(m,x,y);
          tgae->SetPointEXlow (m,x_halfwidth);
          tgae->SetPointEXhigh(m,x_halfwidth);
          tgae->SetPointEYlow (m,y_halfwidth+y-y_center);
          tgae->SetPointEYhigh(m,y_halfwidth-y+y_center);
          m++;
        }
      }
    } else if (Method==2) {
      // Method 2
      // Simply using BayesDivide
      TH1D* pass = (TH1D*)eff->Clone();
      pass->Multiply(den);
      tgae->BayesDivide(pass, den);
    }
    // Add a 2 dummy points (not visible) to make sure to be able to set X axis range
    if (add_low) {
      tgae->SetPoint(0,xlow,ylow-(yhigh-ylow)*10);
      tgae->SetPointEXlow (0,0);
      tgae->SetPointEXhigh(0,0);
      tgae->SetPointEYlow (0,0);
      tgae->SetPointEYhigh(0,0);
    }
    if (add_high) {
      tgae->SetPoint(n+1,xhigh,ylow-(yhigh-ylow)*10);
      tgae->SetPointEXlow (n+1,0);
      tgae->SetPointEXhigh(n+1,0);
      tgae->SetPointEYlow (n+1,0);
      tgae->SetPointEYhigh(n+1,0);
    }
    tgae->GetXaxis()->SetRangeUser(xlow,xhigh);
    tgae->GetYaxis()->SetRangeUser(ylow,yhigh);
    tgae->SetMinimum(ylow);
    tgae->SetMaximum(yhigh);
    tgae->SetTitle(eff->GetTitle());
    tgae->GetXaxis()->SetTitle(eff->GetXaxis()->GetTitle());
    tgae->GetYaxis()->SetTitle(eff->GetYaxis()->GetTitle());
    tgae->GetXaxis()->SetTitleSize(eff->GetXaxis()->GetTitleSize());
    tgae->GetYaxis()->SetTitleSize(eff->GetYaxis()->GetTitleSize());
    tgae->GetXaxis()->SetTitleOffset(eff->GetXaxis()->GetTitleOffset());
    tgae->GetYaxis()->SetTitleOffset(eff->GetYaxis()->GetTitleOffset());
    tgae->GetXaxis()->SetTitleFont(eff->GetXaxis()->GetTitleFont());
    tgae->GetYaxis()->SetTitleFont(eff->GetYaxis()->GetTitleFont());
    tgae->GetXaxis()->SetNdivisions(eff->GetNdivisions("X"));
    tgae->GetYaxis()->SetNdivisions(eff->GetNdivisions("Y"));
    tgae->GetYaxis()->SetDecimals(1);
    tgae->GetXaxis()->SetTimeDisplay(eff->GetXaxis()->GetTimeDisplay());
    tgae->GetXaxis()->SetTimeFormat(eff->GetXaxis()->GetTimeFormat());
    tgae->SetLineColor(eff->GetLineColor());
    tgae->SetLineWidth(eff->GetLineWidth());
    tgae->SetMarkerStyle(eff->GetMarkerStyle());
    tgae->SetMarkerColor(eff->GetMarkerColor());
    tgae->SetMarkerSize(eff->GetMarkerSize());
    double minbinwidth = 99999999, min = eff->GetXaxis()->GetXmin(), max = eff->GetXaxis()->GetXmax();
    for (int bin=1; bin<=eff->GetNbinsX(); ++bin) if (eff->GetBinWidth(bin)<minbinwidth) minbinwidth = eff->GetBinWidth(bin);
    int nbin = (max-min)/minbinwidth;
    // TGraph X-axis ranges are messed up a lot of times by badly chosen automatic bins
    // Make a histogram used to draw the plots and correctly set the ranges
    TH1F* hist;
    if (tgae_hist_.count(eff->GetName())==0) {
      hist = new TH1F((std::string(eff->GetName())+"_forrange").c_str(),"",nbin,min,max);
      tgae_hist_[eff->GetName()] = hist;
    } else {
      hist = tgae_hist_[eff->GetName()];
    }
    hist->GetXaxis()->SetRangeUser(xlow,xhigh);
    hist->GetYaxis()->SetRangeUser(ylow,yhigh);
    hist->SetMinimum(ylow);
    hist->SetMaximum(yhigh);
    hist->SetTitle(eff->GetTitle());
    hist->GetXaxis()->SetTitle(eff->GetXaxis()->GetTitle());
    hist->GetYaxis()->SetTitle(eff->GetYaxis()->GetTitle());
    hist->GetXaxis()->SetTitleSize(eff->GetXaxis()->GetTitleSize());
    hist->GetYaxis()->SetTitleSize(eff->GetYaxis()->GetTitleSize());
    hist->GetXaxis()->SetTitleOffset(eff->GetXaxis()->GetTitleOffset());
    hist->GetYaxis()->SetTitleOffset(eff->GetYaxis()->GetTitleOffset());
    hist->GetXaxis()->SetTitleFont(eff->GetXaxis()->GetTitleFont());
    hist->GetYaxis()->SetTitleFont(eff->GetYaxis()->GetTitleFont());
    hist->GetXaxis()->SetNdivisions(eff->GetNdivisions("X"));
    hist->GetYaxis()->SetNdivisions(eff->GetNdivisions("Y"));
    hist->GetYaxis()->SetDecimals(1);
    hist->GetXaxis()->SetTimeDisplay(eff->GetXaxis()->GetTimeDisplay());
    hist->GetXaxis()->SetTimeFormat(eff->GetXaxis()->GetTimeFormat());
    hist->SetLineColor(eff->GetLineColor());
    hist->SetLineWidth(eff->GetLineWidth());
    hist->SetMarkerStyle(eff->GetMarkerStyle());
    hist->SetMarkerColor(eff->GetMarkerColor());
    hist->SetMarkerSize(eff->GetMarkerSize());
    tgae->SetHistogram(hist);
    // If Bin Labels exist remove current X axis labels and draw similar bin labels but with TLatex
    for (int bin=1; bin<=eff->GetNbinsX(); ++bin) {
      if (std::string(eff->GetXaxis()->GetBinLabel(bin)).size()) {
	bincoordx_.push_back(eff->GetXaxis()->GetBinCenter(bin));
	binlabels_.push_back(eff->GetXaxis()->GetBinLabel(bin));
      }
    }
    return tgae;
  }
  
  void asym_labels_(TH1D *orig, TGraphAsymmErrors* tgae, int hor_vert_up_down = 1) {
    if (binlabels_.size()>10) hor_vert_up_down = 1;
    int angle = hor_vert_up_down==0 ? 0 : hor_vert_up_down==1 ? 90 : hor_vert_up_down==2 ? 20 : -20 ;
    int align = hor_vert_up_down==0 ? 23 : hor_vert_up_down==1 ? 32 : hor_vert_up_down==2 ? 33 : 13;
    if (binlabels_.size()>0) tgae->GetXaxis()->SetLabelColor(0);
    double labelsize = orig->GetXaxis()->GetLabelSize()/1.5;
    double min = orig->GetMinimum(), max = orig->GetMaximum();
    double offset = (max-min) * orig->GetXaxis()->GetLabelOffset() * 5;
    for (size_t i=0; i<binlabels_.size(); ++i) {
      TLatex *lat = new TLatex(bincoordx_[i], min-offset, binlabels_[i].c_str());
      lat->SetTextAlign(align);
      lat->SetTextAngle(angle);
      lat->SetTextFont(orig->GetXaxis()->GetLabelFont());
      lat->SetTextSize(labelsize);
      lat->Draw();
    }
  }
  
  TGraphAsymmErrors* get_tgae_ratio_(TGraphAsymmErrors* data, TGraphAsymmErrors* mc) {
    TGraphAsymmErrors* ratio = (TGraphAsymmErrors*)data->Clone();
    for (int i=0; i<ratio->GetN(); i++) {
      Double_t x_data, y_data;
      data->GetPoint(i, x_data, y_data);
      Double_t xh_data = data->GetErrorXhigh(i);
      Double_t xl_data = data->GetErrorXlow(i);
      Double_t yh_data = data->GetErrorYhigh(i);
      Double_t yl_data = data->GetErrorYlow(i);
      Double_t x_mc, y_mc;
      mc->GetPoint(i, x_mc, y_mc);
      Double_t yh_mc = mc->GetErrorYhigh(i);
      Double_t yl_mc = mc->GetErrorYlow(i);
      if (y_mc==0.) continue;
      
      Double_t y_ratio = y_data/y_mc;
      Double_t yh = sqrt((yh_data*yh_data)/(y_data*y_data) + (yh_mc*yh_mc)/(y_mc*y_mc))*y_ratio;
      //yh = TMath::Min(yh, 1.-y_ratio);
      if (y_data==0.) yh = yh_mc;
      Double_t yl = sqrt((yl_data*yl_data)/(y_data*y_data) + (yl_mc*yl_mc)/(y_mc*y_mc))*y_ratio;
      yl = TMath::Min(yl, y_ratio);
      
      ratio->SetPoint(i, x_data, y_ratio);
      ratio->SetPointError(i, xl_data, xh_data, yl, yh);
    }
    return ratio;
  }
  
  //-----------------------------------------------------------
  
  void draw_one_(TH1D* h) {
    if ((draw_.find("P")!=std::string::npos)) { 
      h->SetMarkerColor(1); 
      h->SetMarkerStyle(20); 
    } else { 
      h->SetLineColor(1);
      h->SetLineWidth(2);
    }
    if (stat_) h->SetStats(1);
    if (norm_&&h->Integral()>0) h->DrawNormalized(draw_.c_str());
    else h->Draw(draw_.c_str());
    if (stat_) set_stat_(h,1,0);
  }
  
  typedef struct DrawParams1D { std::vector<TH1D*> hvec; std::string canname; std::string legtitle;  } DrawParams1D;
  std::vector<DrawParams1D> get_dps_1stpf_1d_() {
    std::vector<DrawParams1D> dps;
    bool is_spec = (find_spec_(name_)!=(size_t)-1)||(find_spec2_(name_)!=(size_t)-1);
    //std::string appr = approval_ ? "_ForApproval" : "";
    std::string appr = "";
    if (ndim_==1||(ndim_==2&&is_spec)) {
      if (npf_==0) {
	dps.push_back({ .hvec={h1d_0p_}, .canname=name_+appr, .legtitle="" });
      } else if (npf_==1) {
	dps.push_back({ .hvec={}, .canname=name_+"_"+pf_names_[0]+appr, .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_1p_[i]);
      } else if (npf_==2) for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_2p_[i][j]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j]+appr;
	dps[dps.size()-1].legtitle=pfs_[1].leg[j];
      } else if (npf_==3) for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_3p_[i][j][k]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+appr;
	dps[dps.size()-1].legtitle=pfs_[1].leg[j];
	if (pfs_[2].leg[k].size()) { if (dps[dps.size()-1].legtitle.size()) dps[dps.size()-1].legtitle+=", "; dps[dps.size()-1].legtitle+=pfs_[2].leg[k]; }
      } else if (npf_==4) for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k) for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_4p_[i][j][k][l]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+appr;
	dps[dps.size()-1].legtitle=pfs_[1].leg[j];
	if (pfs_[2].leg[k].size()) { if (dps[dps.size()-1].legtitle.size()) dps[dps.size()-1].legtitle+=", "; dps[dps.size()-1].legtitle+=pfs_[2].leg[k]; }
	if (pfs_[3].leg[l].size()) { if (dps[dps.size()-1].legtitle.size()) dps[dps.size()-1].legtitle+=", "; dps[dps.size()-1].legtitle+=pfs_[3].leg[l]; }
      } else if (npf_==5) for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k) 
	for (size_t l=0; l<pfs_[3].vec.size(); ++l) for (size_t m=0; m<pfs_[4].vec.size(); ++m) {
	dps.push_back({ .hvec={}, .canname="", .legtitle="" });
	for (size_t i=0; i<pfs_[0].vec.size(); ++i) dps[dps.size()-1].hvec.push_back(h1d_5p_[i][j][k][l][m]);
	dps[dps.size()-1].canname=name_+"_"+pf_names_[0]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[m]+appr;
	dps[dps.size()-1].legtitle=pfs_[1].leg[j];
	if (pfs_[2].leg[k].size()) { if (dps[dps.size()-1].legtitle.size()) dps[dps.size()-1].legtitle+=", "; dps[dps.size()-1].legtitle+=pfs_[2].leg[k]; }
	if (pfs_[3].leg[l].size()) { if (dps[dps.size()-1].legtitle.size()) dps[dps.size()-1].legtitle+=", "; dps[dps.size()-1].legtitle+=pfs_[3].leg[l]; }
	if (pfs_[4].leg[m].size()) { if (dps[dps.size()-1].legtitle.size()) dps[dps.size()-1].legtitle+=", "; dps[dps.size()-1].legtitle+=pfs_[4].leg[m]; }
      }
    }
    return dps;
  }
  
  typedef struct DrawParams2D { TH2D* h; std::string canname; std::string label; } DrawParams2D;
  std::vector<DrawParams2D> get_dps_2d_() {
    std::vector<DrawParams2D> dps;
    bool is_spec = (find_spec_(name_)!=(size_t)-1)||(find_spec2_(name_)!=(size_t)-1);
    //std::string appr = approval_ ? "_ForApproval" : "";
    std::string appr = "";
    if ((ndim_==2&&!is_spec)||(ndim_==3&&is_spec)) {
      if (npf_==0)
	dps.push_back({ .h=h2d_0p_, .canname=name_+appr, .label="" });
      else if (npf_==1) for (size_t i=0; i<pfs_[0].vec.size(); ++i) {
	dps.push_back({ .h=h2d_1p_[i], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+appr, .label="" });
	if (pfs_[0].leg[i].size()) dps[dps.size()-1].label = pfs_[0].leg[i];
      }
      else if (npf_==2) for (size_t i=0; i<pfs_[0].vec.size(); ++i) for (size_t j=0; j<pfs_[1].vec.size(); ++j) {
	dps.push_back({ .h=h2d_2p_[i][j], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+appr, .label="" });
	if (pfs_[0].leg[i].size()) dps[dps.size()-1].label = pfs_[0].leg[i];
	if (pfs_[1].leg[j].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[1].leg[j]; }
      }
      else if (npf_==3) for (size_t i=0; i<pfs_[0].vec.size(); ++i)
	for (size_t j=0; j<pfs_[1].vec.size(); ++j) for (size_t k=0; k<pfs_[2].vec.size(); ++k) {
	  dps.push_back({ .h=h2d_3p_[i][j][k], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+appr, .label="" });
	  if (pfs_[0].leg[i].size()) dps[dps.size()-1].label = pfs_[0].leg[i];
	  if (pfs_[1].leg[j].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[1].leg[j]; }
	  if (pfs_[2].leg[k].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[2].leg[k]; }
	}
      else if (npf_==4) for (size_t i=0; i<pfs_[0].vec.size(); ++i) for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	for (size_t k=0; k<pfs_[2].vec.size(); ++k) for (size_t l=0; l<pfs_[3].vec.size(); ++l) {
	  dps.push_back({ .h=h2d_4p_[i][j][k][l], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+appr, .label="" });
	  if (pfs_[0].leg[i].size()) dps[dps.size()-1].label = pfs_[0].leg[i];
	  if (pfs_[1].leg[j].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[1].leg[j]; }
	  if (pfs_[2].leg[k].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[2].leg[k]; }
	  if (pfs_[3].leg[l].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[3].leg[l]; }
	}
      else if (npf_==5) for (size_t i=0; i<pfs_[0].vec.size(); ++i) for (size_t j=0; j<pfs_[1].vec.size(); ++j)
	for (size_t k=0; k<pfs_[2].vec.size(); ++k) for (size_t l=0; l<pfs_[3].vec.size(); ++l) for (size_t m=0; m<pfs_[4].vec.size(); ++m) {
	  dps.push_back({ .h=h2d_5p_[i][j][k][l][m], .canname=name_+"_"+pf_names_[0]+"_"+pfs_[0].vec[i]+"_"+pfs_[1].vec[j]+"_"+pfs_[2].vec[k]+"_"+pfs_[3].vec[l]+"_"+pfs_[4].vec[l]+appr, .label="" });
	  if (pfs_[0].leg[i].size()) dps[dps.size()-1].label = pfs_[0].leg[i];
	  if (pfs_[1].leg[j].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[1].leg[j]; }
	  if (pfs_[2].leg[k].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[2].leg[k]; }
	  if (pfs_[3].leg[l].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[3].leg[l]; }
	  if (pfs_[4].leg[m].size()) { if (dps[dps.size()-1].label.size()) dps[dps.size()-1].label+=", "; dps[dps.size()-1].label+=pfs_[4].leg[m]; }
	}
    }
    return dps;
  }
  
public:
  void DrawPlots(bool debug=0) {
    if (debug) {
      std::cout<<name_<<" ";
      for (size_t i=0, n=pf_names_.size(); i<n; ++i) std::cout<<pf_names_[i]<<" ";
      std::cout<<"\n";
    }
    calc_specials_();
    //gStyle->SetOptStat(stat_);
    // 1D plots
    std::vector<DrawParams1D> dps1d = get_dps_1stpf_1d_();
    for (size_t i=0; i<dps1d.size(); ++i) {
      size_t skip = 0;
      while (dps1d[i].hvec[skip]->GetEntries()==0) { 
	++skip;
	if (skip==dps1d[i].hvec.size()) break; 
      }
      // Do not draw data also if the plot is blinded
      //if (std::string(dps1d[i].hvec[skip]->GetName()).find("Blind")!=std::string::npos&&
      //    (approval_/10==5)&&skip==0) ++skip;
      if (skip<dps1d[i].hvec.size()) {
	TCanvas *c = custom_can_(dps1d[i].hvec[skip], dps1d[i].canname, 1,1, (1+doublex_)*500,500);
	bool y_range_set = 0;
	if (ranges_.size()>=2) if (ranges_[0]!=ranges_[1]) 
	  dps1d[i].hvec[skip]->GetXaxis()->SetRangeUser(ranges_[0],ranges_[1]);
	if (ranges_.size()>=4) if (ranges_[2]!=ranges_[3]) {
	  dps1d[i].hvec[skip]->GetYaxis()->SetRangeUser(ranges_[2]*(norm_?dps1d[i].hvec[skip]->Integral():1),ranges_[3]*(norm_?dps1d[i].hvec[skip]->Integral():1));
	  y_range_set = 1;
	}
	if (npf_)  {
	  if (ranges_.size()==6)
	    multidraw_with_legend_(skip, dps1d[i].hvec, pfs_[0].leg, pfs_[0].colz, dps1d[i].legtitle, ranges_[4], ranges_[5]);
	  else multidraw_with_legend_(skip, dps1d[i].hvec, pfs_[0].leg, pfs_[0].colz, dps1d[i].legtitle);
	  if (TString(name_).Contains("HTJet1AK8PtLow"))  draw_ak8pt_bins(dps1d[i].hvec[skip], 0);
	  if (TString(name_).Contains("HTJet1AK8PtHigh")) draw_ak8pt_bins(dps1d[i].hvec[skip], 1);
	} else draw_one_(dps1d[i].hvec[0]);
	if (!norm_ && y_range_set) add_labels_(dps1d[i].hvec[skip]);
	write_(c);
	if (stack_&&ratio_&&skip==0) {
	  gPad->Update();
	  add_stack_ratio_plot_(c,ranges_[0],ranges_[1]);
	  if (TString(name_).Contains("Razor")) { c->cd(1); draw_mr_bins(dps1d[i].hvec[skip], ranges_[2], ranges_[3]); }
	  write_(c);
	}
      }
    }
    // 2D Plots
    std::vector<DrawParams2D> dps2d = get_dps_2d_();
    for (size_t i=0; i<dps2d.size(); ++i) {
      if (dps2d[i].h->GetEntries()) {
	if (ranges_.size()>=2) if (ranges_[0]!=ranges_[1]) 
	  dps2d[i].h->GetXaxis()->SetRangeUser(ranges_[0],ranges_[1]);
	if (ranges_.size()>=4) if (ranges_[2]!=ranges_[3]) 
	  dps2d[i].h->GetYaxis()->SetRangeUser(ranges_[2],ranges_[3]);
	if (ranges_.size()>=6) if (ranges_[4]!=ranges_[5]) 
	  dps2d[i].h->GetZaxis()->SetRangeUser(ranges_[4],ranges_[5]);
	TCanvas *c = custom_can_(dps2d[i].h, dps2d[i].canname, 0,0, (1+doublex_)*500,500);
	add_labels_(dps2d[i].h);
	// Add more custom labels (similar to TLegend header)
	if (ranges_.size()==8) {
	  double x  = ranges_[0]+(ranges_[1]-ranges_[0])*ranges_[6];
	  double y  = ranges_[2]+(ranges_[3]-ranges_[2])*ranges_[7];
	  TLatex* lat = new TLatex(x, y, dps2d[i].label.c_str()); 
	  lat->SetLineWidth(2);
	  lat->SetTextSize(0.04);
	  lat->SetTextFont(42);
	  lat->Draw();
	}
	write_(c);
      }
    }
  }
  
  // histo containers
  TH1D*& GetH1D() { return h1d_0p_; }
  TH2D*& GetH2D() { return h2d_0p_; }
  TH3D*& GetH3D() { return h3d_0p_; }
  std::vector<TH1D*>& GetH1D1P() { return h1d_1p_; }
  std::vector<TH2D*>& GetH2D1P() { return h2d_1p_; }
  std::vector<TH3D*>& GetH3D1P() { return h3d_1p_; }
  std::vector<std::vector<TH1D*> >& GetH1D2P() { return h1d_2p_; }
  std::vector<std::vector<TH2D*> >& GetH2D2P() { return h2d_2p_; }
  std::vector<std::vector<TH3D*> >& GetH3D2P() { return h3d_2p_; }
  std::vector<std::vector<std::vector<TH1D*> > >& GetH1D3P() { return h1d_3p_; }
  std::vector<std::vector<std::vector<TH2D*> > >& GetH2D3P() { return h2d_3p_; }
  std::vector<std::vector<std::vector<TH3D*> > >& GetH3D3P() { return h3d_3p_; }
  std::vector<std::vector<std::vector<std::vector<TH1D*> > > >& GetH1D4P() { return h1d_4p_; }
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > >& GetH2D4P() { return h2d_4p_; }
  std::vector<std::vector<std::vector<std::vector<TH3D*> > > >& GetH3D4P() { return h3d_4p_; }
  std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*> > > > >& GetH1D5P() { return h1d_5p_; }
  std::vector<std::vector<std::vector<std::vector<std::vector<TH2D*> > > > >& GetH2D5P() { return h2d_5p_; }
  std::vector<std::vector<std::vector<std::vector<std::vector<TH3D*> > > > >& GetH3D5P() { return h3d_5p_; }
  
};

class SmartHistos {
  
public:
  // constructor, destructor
  SmartHistos() { 
    pf_ = new Postfixes();
    cuts_ = new Cuts();
    weights_={}; 
    //                Name Pre/Suffix, Title Pre/Suffix
    spec2_.push_back({"Avg",           "Avg. "});
    spec2_.push_back({"MPV",           " MPV"});
    gErrorIgnoreLevel = kWarning;
  }
  ~SmartHistos() {}
  typedef struct FillParams { size_t nbin; std::vector<double> bins; std::function<double()> fill; std::string axis_title; std::vector<double> def_range; std::map<int, std::string> bin_labels; } FillParams;
  typedef struct FillParams2 { size_t nbin; std::vector<double> bins; std::map<int, std::string> bin_labels; std::function<double()> fill; std::string axis_title; } FillParams2;
  
private:
  // Postfix container
  Postfixes* pf_;
  
  // Cut container
  Cuts* cuts_;
  
  // Histo weights
  std::vector<std::function<double()> > weights_;
  
  // Special calculations
  std::vector<std::vector<std::string> > spec_;
  std::vector<std::vector<std::string> > spec2_;
  
  // SmartHisto containers
  std::map<std::string, std::vector<SmartHisto*> > sh_;
  
  // FillParams container
  std::map<std::string, FillParams> hp_map_;
  
  // Name of the objects
  std::map<std::string, std::string> objects_;
  
  FillParams get_hp_(std::string name) {
    // Check if name has a special pre/suffix (Avg/MPV etc)
    // and remove them
    for (size_t s=0; s<spec2_.size(); ++s) 
      if ((s==0)? name.find(spec2_[s][0])==0 : name.find(spec2_[s][0])!=std::string::npos)
	name.erase(name.find(spec2_[s][0]), spec2_[s][0].size());
    size_t count = hp_map_.count(name);
    if (count) return hp_map_[name];
    else {
      std::cout<<"!!! ERROR: SmartHistos::get_hp_: FillParams with name = "<<name<<" was not found."<<std::endl;
      return FillParams({.nbin=0,.bins={},.fill={ [](){return -9999.9;} },.axis_title="", .def_range={}, .bin_labels=std::map<int, std::string>()});
    }
  }
  
  std::pair<std::vector<std::string>, std::vector<FillParams> > get_hp_vec_(std::string fill) {
    std::vector<FillParams> fillparams;
    std::vector<std::string> names;
    size_t found1 = 0;
    size_t found2 = fill.find("_vs_",found1);
    while (found2 != std::string::npos) {
      fill.erase(found2,4);
      std::string name = fill.substr(found1,found2-found1);
      fillparams.push_back(get_hp_(name));
      names.push_back(name);
      found1=found2;
      found2 = fill.find("_vs_",found1);
    }
    std::string name = fill.substr(found1,fill.length()-found1);
    fillparams.push_back(get_hp_(name));
    names.push_back(name);
    std::pair<std::vector<std::string>, std::vector<FillParams> > pair = std::make_pair(names,fillparams);
    return pair;
  }
  
  void set_default_style_() {
    gStyle->SetPaperSize(20.,20.);
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetErrorX(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineWidth(2);
    //gStyle->SetLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(1);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleX(1);
    gStyle->SetTitleY(1);
    gStyle->SetTitleAlign(33);
  }
  
public:
  typedef struct Special { std::string name; std::string name_plus_1d; std::string axis; std::string axis_plus_1d; } Special;
  void AddSpecial(Special spec) { 
    std::vector<std::string> param;
    param.push_back(spec.name);
    param.push_back(spec.name_plus_1d);
    param.push_back(spec.axis);
    param.push_back(spec.axis_plus_1d);
    spec_.push_back(param); 
  }
  
  void AddNewFillParam(std::string name, FillParams hp) {
    if (hp.nbin>=hp.bins.size()&&hp.bins.size()>2) 
      std::cout<<"SmartHistos::AddNewFillParam() Warning: Smaller number of bins defined, decrease nbin param or add new bins!\n";
    hp_map_.insert(std::pair<std::string, FillParams>(name, hp)); 
  }
  void AddNewFillParam(std::string name, FillParams2 hp) {
    if (hp.nbin<hp.bin_labels.size()) 
      std::cout<<"SmartHistos::AddNewFillParam() Warning: Smaller number of bins defined, increase nbin param or remove bin titles!\n";
    FillParams fp({.nbin=hp.nbin,.bins=hp.bins,.fill=hp.fill,.axis_title=hp.axis_title,.def_range={},.bin_labels=hp.bin_labels});
    hp_map_.insert(std::pair<std::string, FillParams>(name, fp )); 
  }
  
  void AddNewPostfix(std::string name, std::function<size_t()> sel, std::string pf, std::string leg, std::string colz) { pf_->AddNew(name, sel, pf, leg, colz); }
  
  void AddNewCut(std::string name, std::function<bool()> cut) { cuts_->AddNew(name, cut); }
  
  void SetHistoWeights(std::vector<std::function<double()> > weights) { weights_ = weights; }
  
  void AddHistoType(std::string type, std::string objects) { sh_[type] = std::vector<SmartHisto*>(); objects_[type] = objects; }
  
  typedef struct HistoParams    { std::string fill; std::vector<std::string> pfs; std::vector<std::string> cuts; std::string draw; std::string opt; std::vector<double> ranges; } HistoParams;
  typedef struct HistoParamsNew { std::string fill; std::vector<std::string> pfs; std::vector<std::string> cuts; std::string draw; std::string opt; std::vector<double> ranges={}; } HistoParamsNew;
  
  void AddHistos(std::string name, HistoParams hp, bool AddCutsToTitle = true, const bool debug = 0) {
    if (debug) std::cout<<"Start adding: "<<name<<" "<<hp.fill<<std::endl;
    if (sh_.count(name)) {
      std::pair<std::vector<std::string>, std::vector<FillParams> > hp_vec = get_hp_vec_(hp.fill);
      bool valid = 1;
      for (size_t i=0; i<hp_vec.second.size(); ++i) valid = (valid&&(hp_vec.second[i].nbin!=0));
      if (valid) {
        std::vector<Postfixes::Postfix> pfs;
	if (debug) std::cout<<"ok"<<std::endl;
        for (size_t i=0; i<hp.pfs.size(); ++i) pfs.push_back(pf_->GetPostfix(hp.pfs[i]));
	if (debug) std::cout<<"ok"<<std::endl;
	// Determine special histo name and axes
	std::string spec_histo_name = hp.fill, histo_name = hp.fill;
        std::string spec_axis_titles = "", axis_titles = "";
        std::vector<std::function<double()> > fillfuncs;
	std::vector<std::map<int, std::string> > bin_labels;
        std::vector<double> ranges;
	if (debug) std::cout<<"ok"<<std::endl;
	if (AddCutsToTitle) for (size_t icut=0; icut<hp.cuts.size(); ++icut) {
	  if (icut>0) {
	    axis_titles += ", ";
	    spec_axis_titles += ", ";
	  }
	  axis_titles += hp.cuts[icut];
	  spec_axis_titles += hp.cuts[icut];
	}
	if (debug) std::cout<<axis_titles<<std::endl;
	if (debug) std::cout<<spec_axis_titles<<std::endl;
	bool isSpecial = 0;
        for (size_t i=hp_vec.second.size(), i_hp=0; i>0; --i, ++i_hp) {
	  std::string hp_name = hp_vec.first[i-1];
	  if (debug) std::cout<<i<<" "<<hp_name<<std::endl;
	  FillParams fill_params = hp_vec.second[i-1];
	  std::string spec_axis_title = fill_params.axis_title;
	  std::string axis_title = fill_params.axis_title;
	  if (debug) std::cout<<axis_title<<std::endl;
	  if (debug) std::cout<<spec_axis_title<<std::endl;
	  // First find special axes and treat them differently
	  size_t s = -1, s2 = -1;
	  for (size_t j=0; j<spec_.size(); ++j) if (hp_name.find(spec_[j][0])!=std::string::npos) s = j;
	  if (debug) std::cout<<"s = "<<s<<std::endl;
	  if (s!=(size_t)-1) {
	    isSpecial = 1;
	    if (debug) std::cout<<"spec found"<<std::endl;
	    histo_name = std::string(histo_name).replace(histo_name.find(spec_[s][0]),spec_[s][0].size(),spec_[s][1]);
	    if (debug) std::cout<<histo_name<<std::endl;
	    // Replace secondary histo (+1 dimensional) axis title
	    if (debug) std::cout<<"axis_title: "<<axis_title<<" spec: "<<spec_[s][2]<<" find: "<<axis_title.find(spec_[s][2])<<std::endl;
	    axis_title.replace(axis_title.find(spec_[s][2]), spec_[s][2].size(), spec_[s][3]);
	    if (debug) std::cout<<axis_title<<std::endl;
	  }
	  if (debug) std::cout<<"ok "<<i<<std::endl;
	  for (size_t j=0; j<spec2_.size(); ++j) {
	    if (j==0&&hp_name.find(spec2_[j][0])==0) s2 = j;
	    else if (j!=0&&hp_name.find(spec2_[j][0])!=std::string::npos) s2 = j;
	  }
	  if (debug) std::cout<<"s2 = "<<s2<<std::endl;
	  if (s2!=(size_t)-1) {
	    isSpecial = 1;
	    if (debug) std::cout<<"spec2 found"<<std::endl;
	    histo_name = std::string(histo_name).insert(histo_name.find(spec2_[s2][0]),"For");
	    if (debug) std::cout<<histo_name<<std::endl;
	    if (s2==0) /* Add "Avg. " */  spec_axis_title = std::string(spec_axis_title).insert(0,spec2_[s2][1]);
	    else if (s2==1) {
	      /* Add " MPV" before axis unit (if exist) */
	      size_t find1 = spec_axis_title.find(" (");
	      size_t find2 = spec_axis_title.find(" [");
	      if (debug) std::cout<<find1<<std::endl;
	      if (debug) std::cout<<find2<<std::endl;
	      if (find1 != std::string::npos) spec_axis_title.insert(find1, spec2_[s2][1]);
	      else if (find2 != std::string::npos) spec_axis_title.insert(find2, spec2_[s2][1]);
	      else spec_axis_title.insert(spec_axis_title.size(), spec2_[s2][1]);
	    }
	    if (debug) std::cout<<spec_axis_title<<std::endl;
	  }
	  if (debug) std::cout<<"ok "<<i<<std::endl;
          spec_axis_titles += ";" + spec_axis_title;
	  if (debug) std::cout<<spec_axis_titles<<std::endl;
          axis_titles += ";" + axis_title;
	  if (debug) std::cout<<axis_titles<<std::endl;
          fillfuncs.push_back(fill_params.fill);
	  if (debug) std::cout<<"ok "<<i<<std::endl;
	  // Set ranges, first the default from FillParams, then from AddHistos
	  double min = 0, max = 0;
	  if (fill_params.def_range.size()==2) {
	    double def_min = fill_params.def_range[0], def_max = fill_params.def_range[1];
	    if (def_min!=def_max) { min=def_min; max=def_max; }
	  } else if (fill_params.def_range.size()!=0) {
	    std::cout<<"!!! ERROR: SmartHistos::AddNewFillParam: name = "<<hp_name<<", .def_range has too many elements: "<<fill_params.def_range.size()<<std::endl;
	  }
	  if (hp.ranges.size()>=(i_hp+1)*2) {
	    double ran_min = hp.ranges[i_hp*2], ran_max = hp.ranges[i_hp*2+1];
	    if (ran_min!=ran_max) { min=ran_min; max=ran_max; }
	  }
	  ranges.push_back(min);
	  ranges.push_back(max);
	  if (debug) std::cout<<"ok "<<i<<std::endl;
	  bin_labels.push_back(fill_params.bin_labels);
	  if (debug) std::cout<<"ok "<<i<<std::endl;
        }
	if (hp.ranges.size()>ranges.size()) for (size_t i=ranges.size(); i<hp.ranges.size(); ++i) ranges.push_back(hp.ranges[i]);
	if (debug) std::cout<<"ok"<<std::endl;
	// Add object name in y/z axis (only for non-special histos)
	if (!isSpecial) {
	  axis_titles      += ";" + objects_[name];
	  spec_axis_titles += ";" + objects_[name];
	  if (hp_vec.second.size()==1&&hp_vec.second[0].bins.size()==2) {
	    // Add also bin size for 1D plots
	    double binsize = (hp_vec.second[0].bins[1]-hp_vec.second[0].bins[0])/hp_vec.second[0].nbin;
	    // And unit if available
	    //std::string unit = " units";
	    std::string unit = "";
	    size_t beg = hp_vec.second[0].axis_title.find("(");
	    size_t end = hp_vec.second[0].axis_title.find(")");
	    if (beg!=std::string::npos&&end!=std::string::npos) {
	      unit = " "+hp_vec.second[0].axis_title.substr(beg+1,end-beg-1);
	    }
	    std::stringstream ss;
	    ss<<" / "<<binsize<<unit;
	    axis_titles      += ss.str();
	    spec_axis_titles += ss.str();
	  } else if (hp_vec.second.size()==2&&hp_vec.second[0].bins.size()==2) {
	    // or " / bin" for 2D ones
	    axis_titles      += " / bin";
	    spec_axis_titles += " / bin";	  
	  }
	}
	if (debug) std::cout<<"Add object names ok"<<std::endl;
	std::vector<Cut*> cuts;
        for (size_t i=0; i<hp.cuts.size(); ++i) cuts.push_back(cuts_->GetCut(hp.cuts[i]));
	if (debug) std::cout<<"ok"<<std::endl;
        sh_[name].push_back(new SmartHisto(hp.fill.c_str(), hp.pfs, pfs, fillfuncs, weights_, cuts, hp.draw, hp.opt, ranges, bin_labels, spec_, spec2_));
        if (hp_vec.second.size()==1) sh_[name][sh_[name].size()-1]->AddNew(spec_histo_name, spec_axis_titles, hp_vec.second[0].nbin, hp_vec.second[0].bins);
        else if (hp_vec.second.size()==2) sh_[name][sh_[name].size()-1]->AddNew(spec_histo_name, spec_axis_titles,
									 histo_name, axis_titles,
									 hp_vec.second[1].nbin, hp_vec.second[1].bins,
									 hp_vec.second[0].nbin, hp_vec.second[0].bins);
        else if (hp_vec.second.size()==3) sh_[name][sh_[name].size()-1]->AddNew(spec_histo_name, spec_axis_titles,
									 histo_name, axis_titles,
									 hp_vec.second[2].nbin, hp_vec.second[2].bins,
									 hp_vec.second[1].nbin, hp_vec.second[1].bins,
									 hp_vec.second[0].nbin, hp_vec.second[0].bins);
	if (debug) std::cout<<"ok"<<std::endl;
      }
    }
  }
  
  void AddHistos(std::string name, HistoParamsNew hp_new, bool AddCutsToTitle = true, const bool debug = 0) {
    AddHistos(name, HistoParams({.fill=hp_new.fill, .pfs=hp_new.pfs, .cuts=hp_new.cuts, .draw=hp_new.draw, .opt=hp_new.opt, .ranges=hp_new.ranges }));
  }
  
  void PrintNames() { 
    for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it) {
      for (size_t i=0; i<it->second.size(); ++i)  {
	std::cout<<it->second[i]->GetName();
	size_t npf = it->second[i]->GetPFNames().size();
	if (npf>0) {
	  std::cout<<" - ";
	  for (size_t j=0; j<npf; ++j) {
	    std::cout<<it->second[i]->GetPFNames()[j];
	    if (j+1!=npf) std::cout<<", ";
	  }
	}
	std::cout<<"\n";
      }
    }
  }

  void GetTotalNCells() {
    for(const auto& sh : sh_) {
      for (size_t i=0; i<sh.second.size(); ++i)  {
	std::cout<<sh.second[i]->GetTotalNCells()<<" - ";
	std::cout<<sh.second[i]->GetName();
	size_t npf = sh.second[i]->GetPFNames().size();
	if (npf>0) {
	  std::cout<<" - ";
	  for (size_t j=0; j<npf; ++j) {
	    std::cout<<sh.second[i]->GetPFNames()[j];
	    if (j+1!=npf) std::cout<<", ";
	  }
	}
	std::cout<<"\n";
      }
    }
  }
  
  void Load(std::string filename) {
    TFile* f = TFile::Open(filename.c_str()); 
    for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it)
      for (size_t i=0; i<it->second.size(); ++i) it->second[i]->Load(f); 
    f->Close();
  }
  
  void Add(std::string filenames) {
    TChain* fc = new TChain("fc");
    fc->Add(filenames.c_str());
    TObjArray* fileElements=fc->GetListOfFiles();
    TIter next(fileElements); TChainElement* chEl=0;
    while (( chEl=(TChainElement*)next() )) {
      TFile* f = TFile::Open(chEl->GetTitle()); 
      for(std::map<std::string, std::vector<SmartHisto*> >::iterator it = sh_.begin(); it != sh_.end(); ++it)
	for (size_t i=0; i<it->second.size(); ++i) it->second[i]->Add(f); 
      f->Close();
    }
  }
  
  void CalcSpecials() { for (auto vs : sh_) for (auto s : vs.second) s->CalcSpecials(); }
  
  void Fill(std::string name) {
    cuts_->ResetAllCut();
    for (size_t i=0; i<sh_[name].size(); ++i) sh_[name][i]->Fill(); 
  }
  
  void DrawPlots(bool time=0) {
    std::cout<<"Drawing plots ..."<<std::endl;
    TStopwatch sw;
    set_default_style_();
    for(auto h : sh_) {
      sw.Start();
      for (size_t i=0; i<h.second.size(); ++i) h.second[i]->DrawPlots(); 
      sw.Stop();
      if (time) std::cout<<h.first<<" took "<<sw.RealTime()<<"s"<<std::endl;
    }
  }
  
  void Write(std::string name = "", bool time=0) { 
    std::cout<<"Writing histograms ..."<<std::endl;
    TStopwatch sw;
    if (name.size()) for (size_t i=0; i<sh_[name].size(); ++i) sh_[name][i]->Write();
    else for(auto h : sh_) {
      sw.Start();
      for (size_t i=0; i<h.second.size(); ++i) h.second[i]->Write();
      sw.Stop();
      if (time) std::cout<<h.first<<" took "<<sw.RealTime()<<"s"<<std::endl;
    }
  }
  
  SmartHisto* GetHistos(std::string name, size_t i) { return sh_[name][i]; }
  
};
