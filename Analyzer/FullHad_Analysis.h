#include "common/AnalysisBase.h"

class Analysis : public AnalysisBase
{
public:
  Analysis();
  ~Analysis();

  void declare_histograms();

  void define_selections(const DataStruct&);

  void fill_histograms(const DataStruct&);

  std::vector<Cut> analysis_cuts;
  
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
}


//_______________________________________________________
//                 List of Histograms
TH1D* h_nvtx;


//_______________________________________________________
//              Define Histograms here
void
Analysis::declare_histograms()
{
  h_nvtx = new TH1D("h_nvtx",";N_{Vertices}", 100,0,100);
}


//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_histograms(const DataStruct& data)
{
  h_nvtx->Fill(data.evt.NGoodVtx);
}
