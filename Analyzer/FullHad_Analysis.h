#include "common/AnalysisBase.h"

class Analysis : public AnalysisBase
{
public:
  Analysis();
  ~Analysis();

  void declare_histograms();

  void define_selections(const DataStruct&);

  virtual bool signal_selection(const DataStruct&);

  void fill_histograms(const DataStruct&, const double&);

  std::vector<Cut> analysis_cuts;

private:
  bool _apply_ncut(size_t);
};


//_______________________________________________________
//                       Constructor
Analysis::Analysis() : AnalysisBase() { }


//_______________________________________________________
//                       Destructor
Analysis::~Analysis() { }


//_______________________________________________________
//         Define Analysis event selection cuts
void
Analysis::define_selections(const DataStruct& data)
{ 
  // jet 1 pt >= 400
  analysis_cuts.push_back({ .name="jet1_pt",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8.size==0) return 0; // for safety
			      if (data.jetsAK8.Pt[0]<400) return 0;
			      return 1;
			    } });

  // jet 2 pt >= 400
  analysis_cuts.push_back({ .name="jet2_pt",   .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8.size<2) return 0; // for safety
			      if (data.jetsAK8.Pt[1]<400) return 0;
			      return 1;
			    } });

  // 110 <= jet 1 mass (softdrop) < 210
  analysis_cuts.push_back({ .name="jet1_mass", .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8.size==0) return 0; // for safety
			      if (data.jetsAK8.softDropMass[0]<110) return 0;
			      if (data.jetsAK8.softDropMass[0]>=210) return 0;
			      return 1;
			    } });

  // 110 <= jet 2 mass (softdrop) < 210
  analysis_cuts.push_back({ .name="jet2_mass", .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8.size<2) return 0; // for safety
			      if (data.jetsAK8.softDropMass[1]<110) return 0;
			      if (data.jetsAK8.softDropMass[1]>=210) return 0;
			      return 1;
			    } });
  
  // Full-hadronic trigger
  analysis_cuts.push_back({ .name="hlt_ak8ht700_mass50", .func = [&data](){
			      // Define cut function here:
			      return data.evt.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50==1; 
			    } });

  // | DeltaPhi | < 2.7
  analysis_cuts.push_back({ .name="delta_phi", .func = [&data](){ 
			      // Define cut function here:
			      if (data.evt.NTopHadPreTag<2) return 0; // for safety
			      if (data.evt.TTHadDPhi>=2.7) return 0;
			      return 1;
			    } });
}

bool
Analysis::_apply_ncut(size_t ncut) {
  if (ncut>analysis_cuts.size()) return 0;
  for (size_t i=0; i<ncut; ++i) if ( ! analysis_cuts[i].func() ) return 0;
  return 1;
}

//_______________________________________________________
//                 Signal Region
//     Must define it, because we blind it in data!
bool
Analysis::signal_selection(const DataStruct& data) {
  if ( data.evt.NTopHad<2 ) return 0;
  if ( data.evt.TTHadDPhi>=2.7) return 0;
  if ( data.evt.R<0.4 ) return 0;
  return 1;
}

//_______________________________________________________
//                 List of Histograms
TH1D* h_nvtx;
TH2D* h_abcd;

//_______________________________________________________
//              Define Histograms here
void
Analysis::declare_histograms()
{
  h_nvtx = new TH1D("h_nvtx",";N_{Vertices}", 100,0,100);
  Double_t R_bins[3] = { 0, 0.4, 2 };
  Double_t Ntop_bins[3] = { 0, 2, 3 };
  h_abcd = new TH2D("h_abcd",";R;N_{top-tag}", 2, R_bins, 2, Ntop_bins );
}


//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_histograms(const DataStruct& data, const double& weight)
{
  h_nvtx->Fill(data.evt.NGoodVtx);
  if ( _apply_ncut( analysis_cuts.size() ) ) h_abcd->Fill(data.evt.R, data.evt.NTopHad, weight);
}
