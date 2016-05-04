#include "TLorentzVector.h"
#include "common/AnalysisBase.h"

class Analysis : public AnalysisBase
{
public:
  Analysis();
  ~Analysis();

  void calculate_variables(DataStruct&);

  double get_analysis_weight(DataStruct&);

  bool pass_skimming(DataStruct&);

  void define_selections(const DataStruct&);

  virtual bool signal_selection(const DataStruct&);

  void define_histo_options(const double&, const DataStruct&, std::string, bool);

  void init_analysis_histos();

  void fill_analysis_histos(DataStruct&, const double&);

  void load_analysis_histos(std::string);

  void save_analysis_histos(bool);

  std::vector<Cut> analysis_cuts;

private:
  bool _apply_cut(std::string);
  bool _apply_ncut(size_t);
};


//_______________________________________________________
//                       Constructor
Analysis::Analysis() : AnalysisBase() { }


//_______________________________________________________
//                       Destructor
Analysis::~Analysis() { }

//_______________________________________________________
//                  Calculate variables

// AK8 CHS jets
std::vector<int>   pass_Loose_Jet_ID;

// event
unsigned int nLooseJet;
double AK4_Ht;

void
Analysis::calculate_variables(DataStruct& data)
{
  // Jet variables (initialize)
  nLooseJet = 0;
  pass_Loose_Jet_ID.clear();

  // Loop on jets
  while(data.jetsAK8.Loop()) {
    // _______________________________________________________
    //                       Jet ID
    
    double eta = data.jetsAK8.Eta[data.jetsAK8.it];
    double NHF = data.jetsAK8.neutralHadronEnergy[data.jetsAK8.it] / data.jetsAK8.E[data.jetsAK8.it];
    double NEMF = data.jetsAK8.neutralEmEnergy[data.jetsAK8.it] / data.jetsAK8.E[data.jetsAK8.it];
    double CHF = data.jetsAK8.chargedHadronEnergy[data.jetsAK8.it]/data.jetsAK8.E[data.jetsAK8.it];
    double CEMF = data.jetsAK8.chargedEmEnergy[data.jetsAK8.it]/data.jetsAK8.E[data.jetsAK8.it];
    int NumConst = data.jetsAK8.chargedMultiplicity[data.jetsAK8.it] + data.jetsAK8.neutralMultiplicity[data.jetsAK8.it];
    int CHM = data.jetsAK8.chargedMultiplicity[data.jetsAK8.it];
    int NumNeutralParticle = data.jetsAK8.neutralMultiplicity[data.jetsAK8.it];
    bool pass_Loose_ID = ( (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta)>2.4) && fabs(eta)<=3.0 )
      || ( NEMF<0.9 && NumNeutralParticle>10 && fabs(eta)>3.0 );
    pass_Loose_Jet_ID.push_back(pass_Loose_ID);
    if (pass_Loose_ID) nLooseJet++;
    
  } // end of AK8 jet loop
  
  AK4_Ht = 0;
  // Loop on AK4 jets
  while(data.jetsAK4.Loop()) {
    AK4_Ht += data.jetsAK4.Pt[data.jetsAK4.it];
  }
}

//_______________________________________________________
//          Define Analysis specific weights
double
Analysis::get_analysis_weight(DataStruct& data)
{
  double w = 1;
  
  /*
  //____________________________________________________
  //               Jet pT reweighting
  
  // Use linear function calculated by scripts/JetPtScaleFactors.C script
  // reweight event corresponding to the product of SFs for each jet
  // linear function: p0 + p1 * jet pt
  
  const double p0 = 0.970718;
  const double p1 = -0.000331894;
  while(data.jetsAK8.Loop()) {
    w *= p0 + p1 * data.jetsAK8.Pt[data.jetsAK8.it];
  }
  */
  
  return w;
}

//_______________________________________________________
//          Define Analysis event selection cuts

bool
Analysis::pass_skimming(DataStruct& data)
{
  float pt_threshold = 300;
  int N_CHS = 0, N_Puppi = 0;
  while(data.jetsAK8.Loop()) if (data.jetsAK8.Pt[data.jetsAK8.it]>=pt_threshold) ++N_CHS;
  while(data.jetsAK8Puppi.Loop()) if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=pt_threshold) ++N_Puppi;
  return (N_CHS >= 1 || N_Puppi >= 1);
}

void
Analysis::define_selections(const DataStruct& data)
{
  // cut1: njet >= 1
  analysis_cuts.push_back({ .name="1jet",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8.size<1) return 0;
			      return 1;
			    } });

  // cut2: jet 1 pass loose jet id
  analysis_cuts.push_back({ .name="jet1_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8.size<1) return 0; // for safety
			      return pass_Loose_Jet_ID[0];
			    } });

  // cut3: jet 1 eta < 2.4
  analysis_cuts.push_back({ .name="jet1_eta",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8.size<1) return 0; // for safety
			      if (fabs(data.jetsAK8.Eta[0])>=2.4) return 0;
			      return 1;
			    } });

  // cut4: jet 1 pt >= 400
  analysis_cuts.push_back({ .name="jet1_pt",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8.size<1) return 0; // for safety
			      if (data.jetsAK8.Pt[0]<400) return 0;
			      return 1;
			    } });

  // cut5: Full-hadronic trigger
  analysis_cuts.push_back({ .name="hlt_ak8ht700_mass50", .func = [&data](){
			      // Define cut function here:
			      return data.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; 
			    } });

}

bool
Analysis::_apply_cut(std::string cut_name) {
  for (auto cut : analysis_cuts) if (cut_name == cut.name) return cut.func();
  return 0;
}

bool
Analysis::_apply_ncut(size_t ncut) {
  if (ncut>analysis_cuts.size()) return 0;
  for (size_t i=0; i<ncut; ++i) if ( ! analysis_cuts[i].func() ) return 0;
  return 1;
}

//_______________________________________________________
//                 Signal Region
// Must define at some point, because we need to blind this region in data

bool
Analysis::signal_selection(const DataStruct& data) {
  //if ( data.evt.NTopHad<2 ) return 0;
  //if ( data.evt.TTHadDPhi>=2.7) return 0;
  //if ( data.evt.R<0.4 ) return 0;
  return 0;
}

//_______________________________________________________
//                 List of Histograms
TH1D* h_njet;
TH1D* h_ht;
TH1D* h_jet1_pt;

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos()
{
  h_njet    = new TH1D("njet",   ";N_{AK8-jet, loose ID}",  20, 0,  20);
  h_ht      = new TH1D("ht",     ";H_{T, AK4}",            200, 0,2000);
  h_jet1_pt = new TH1D("jet1_pt",";p_{T, jet1}",           200, 0,2000);

  h_njet->Sumw2();
  h_ht->Sumw2();
  h_jet1_pt->Sumw2();
}


//_______________________________________________________
//               Fill Histograms here
void
Analysis::fill_analysis_histos(DataStruct& data, const double& weight)
{
  h_njet->Fill(nLooseJet, weight);
  
  h_ht->Fill(AK4_Ht, weight);
  
  if (_apply_ncut(2)) {
    h_jet1_pt->Fill(data.jetsAK8.Pt[0], weight);
  }
}

// Methods used by SmartHistos (Plotter)
// Can leave them empty
void
Analysis::define_histo_options(const double& weight, const DataStruct& d, std::string dirname, bool runOnSkim=false)
{
}

void
Analysis::load_analysis_histos(std::string inputfile)
{
}

void
Analysis::save_analysis_histos(bool draw=0)
{
}
