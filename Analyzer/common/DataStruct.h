
#ifndef DataStruct_h
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

  class EventData {
  public:
    EventData() { init(); };
    
    int npv;
    double vx;
    double vy;
    double vz;
    unsigned int RunNumber;
    unsigned int LumiBlock;
    long EventNumber;
    double rho;
    int NGoodVtx;
    int NLep;
    int NTopHad;
    int NTopHadPreTag;
    int NTopLep;
    int NTop;
    int LHA_PDF_ID;
    float HtLep;
    float HtTop;
    float Ht;
    float HtAll;
    float HtEx;
    float HtExFr;
    float HtTopFr;
    float TTHadDR;
    float TTHadDPhi;
    float TTHadDEta;
    float TTHadMass;
    float TTHadSumPt;
    float TTHadPz;
    float TTHadHz;
    float TTHadDPz;
    float TTHadMR;
    float TTHadMTR;
    float TTHadR;
    float TTHadR2;
    float MR;
    float MTR;
    float R;
    float R2;
    float AK8Puppi_MR;
    float AK8Puppi_MTR;
    float AK8Puppi_R;
    float AK8Puppi_R2;
    float AK4_MR;
    float AK4_MTR;
    float AK4_R;
    float AK4_R2;
    float XSec;
    float NEvent_Corr;
    float Lumi_Weight;
    float Gen_Weight;
    float Gen_Ht;
    float SUSY_Gluino_Mass;
    float SUSY_LSP_Mass;
    
    void init() {
      npv=NOVAL_I;
      vx=NOVAL_F;
      vy=NOVAL_F;
      vz=NOVAL_F;
      RunNumber=9999;
      LumiBlock=9999;
      EventNumber=9999;
      rho=NOVAL_F;
      NGoodVtx=NOVAL_I;
      NLep=NOVAL_I;
      NTopHad=NOVAL_I;
      NTopHadPreTag=NOVAL_I;
      NTopLep=NOVAL_I;
      NTop=NOVAL_I;
      LHA_PDF_ID=NOVAL_I;
      HtLep=NOVAL_F;
      HtTop=NOVAL_F;
      Ht=NOVAL_F;
      HtAll=NOVAL_F;
      HtEx=NOVAL_F;
      HtExFr=NOVAL_F;
      HtTopFr=NOVAL_F;
      TTHadDR=NOVAL_F;
      TTHadDPhi=NOVAL_F;
      TTHadDEta=NOVAL_F;
      TTHadMass=NOVAL_F;
      TTHadSumPt=NOVAL_F;
      TTHadPz=NOVAL_F;
      TTHadHz=NOVAL_F;
      TTHadDPz=NOVAL_F;
      TTHadMR=NOVAL_F;
      TTHadMTR=NOVAL_F;
      TTHadR=NOVAL_F;
      TTHadR2=NOVAL_F;
      MR=NOVAL_F;
      MTR=NOVAL_F;
      R=NOVAL_F;
      R2=NOVAL_F;
      AK8Puppi_MR=NOVAL_F;
      AK8Puppi_MTR=NOVAL_F;
      AK8Puppi_R=NOVAL_F;
      AK8Puppi_R2=NOVAL_F;
      AK4_MR=NOVAL_F;
      AK4_MTR=NOVAL_F;
      AK4_R=NOVAL_F;
      AK4_R2=NOVAL_F;
      XSec=NOVAL_F;
      NEvent_Corr=NOVAL_F;
      Lumi_Weight=NOVAL_F;
      Gen_Weight=NOVAL_F;
      Gen_Ht=NOVAL_F;
      SUSY_Gluino_Mass=NOVAL_F;
      SUSY_LSP_Mass=NOVAL_F;
    }
    
  } evt;
  
  class METData {
  public:
    METData() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Px;
    std::vector<float> Py;
    std::vector<float> Phi;
    std::vector<float> uncorPt;
    std::vector<float> uncorPhi;
    std::vector<float> uncorSumEt;
    
    void init() {
      size=9999;
      init_vec(Pt);
      init_vec(Px);
      init_vec(Py);
      init_vec(Phi);
      init_vec(uncorPt);
      init_vec(uncorPhi);
      init_vec(uncorSumEt);
    }
    
  } met;
  
  class PileupData {
  public:
    PileupData() { init(); };
    
    int NtrueInt;
    unsigned int size;
    std::vector<int> BX;
    std::vector<int> NInt;
    
    unsigned int it;
    
    void init() {
      it = -1;
      NtrueInt=NOVAL_I;
      size=9999;
      init_vec(BX);
      init_vec(NInt);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } pu;
  
  class VertexData {
  public:
    VertexData() { init(); };
    
    unsigned int size;
    std::vector<int> ndof;
    std::vector<float> z;
    std::vector<float> rho;
    std::vector<float> chi;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(ndof);
      init_vec(z);
      init_vec(rho);
      init_vec(chi);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } vtx;
  
  class SystScaleData {
  public:
    SystScaleData() { init(); };
    
    unsigned int size;
    std::vector<float> Weights;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Weights);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } syst_scale;
  
  class SystPDFData {
  public:
    SystPDFData() { init(); };
    
    unsigned int size;
    std::vector<float> Weights;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Weights);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } syst_pdf;
  
  class SystAlphaSData {
  public:
    SystAlphaSData() { init(); };
    
    unsigned int size;
    std::vector<float> Weights;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Weights);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } syst_alphas;
  
  class FilterData {
  public:
    FilterData() { init(); };
    
    int HBHEIsoNoiseFilterResult;
    int HBHENoiseFilterResult;
    int HBHENoiseFilterResultRun1;
    int HBHENoiseFilterResultRun2Loose;
    int HBHENoiseFilterResultRun2Tight;
    int goodVertices;
    int eeBadScFilter;
    int ecalLaserCorrFilter;
    int EcalDeadCellTriggerPrimitiveFilter;
    int EcalDeadCellBoundaryEnergyFilter;
    int HcalStripHaloFilter;
    int hcalLaserEventFilter;
    int HBHENoiseFilter;
    int HBHENoiseIsoFilter;
    int CSCTightHaloFilter;
    int CSCTightHaloTrkMuUnvetoFilter;
    int CSCTightHalo2015Filter;
    int muonBadTrackFilter;
    int chargedHadronTrackResolutionFilter;
    int trkPOGFilters;
    int trkPOG_manystripclus53X;
    int trkPOG_toomanystripclus53X;
    int trkPOG_logErrorTooManyClusters;
    int METFilters;
    
    void init() {
      HBHEIsoNoiseFilterResult=NOVAL_I;
      HBHENoiseFilterResult=NOVAL_I;
      HBHENoiseFilterResultRun1=NOVAL_I;
      HBHENoiseFilterResultRun2Loose=NOVAL_I;
      HBHENoiseFilterResultRun2Tight=NOVAL_I;
      goodVertices=NOVAL_I;
      eeBadScFilter=NOVAL_I;
      ecalLaserCorrFilter=NOVAL_I;
      EcalDeadCellTriggerPrimitiveFilter=NOVAL_I;
      EcalDeadCellBoundaryEnergyFilter=NOVAL_I;
      HcalStripHaloFilter=NOVAL_I;
      hcalLaserEventFilter=NOVAL_I;
      HBHENoiseFilter=NOVAL_I;
      HBHENoiseIsoFilter=NOVAL_I;
      CSCTightHaloFilter=NOVAL_I;
      CSCTightHaloTrkMuUnvetoFilter=NOVAL_I;
      CSCTightHalo2015Filter=NOVAL_I;
      muonBadTrackFilter=NOVAL_I;
      chargedHadronTrackResolutionFilter=NOVAL_I;
      trkPOGFilters=NOVAL_I;
      trkPOG_manystripclus53X=NOVAL_I;
      trkPOG_toomanystripclus53X=NOVAL_I;
      trkPOG_logErrorTooManyClusters=NOVAL_I;
      METFilters=NOVAL_I;
    }
    
  } filter;
  
  class HLTData {
  public:
    HLTData() { init(); };
    
    int AK8PFJet360_TrimMass30;
    int AK8PFJet360_TrimMass30_prescale;
    int PFJet40;
    int PFJet40_prescale;
    int PFJet60;
    int PFJet60_prescale;
    int PFJet80;
    int PFJet80_prescale;
    int PFJet140;
    int PFJet140_prescale;
    int PFJet200;
    int PFJet200_prescale;
    int PFJet260;
    int PFJet260_prescale;
    int PFJet320;
    int PFJet320_prescale;
    int PFJet400;
    int PFJet400_prescale;
    int PFJet450;
    int PFJet450_prescale;
    int PFJet500;
    int PFJet500_prescale;
    int DiPFJetAve40;
    int DiPFJetAve40_prescale;
    int DiPFJetAve60;
    int DiPFJetAve60_prescale;
    int DiPFJetAve80;
    int DiPFJetAve80_prescale;
    int DiPFJetAve140;
    int DiPFJetAve140_prescale;
    int DiPFJetAve200;
    int DiPFJetAve200_prescale;
    int DiPFJetAve260;
    int DiPFJetAve260_prescale;
    int DiPFJetAve320;
    int DiPFJetAve320_prescale;
    int DiPFJetAve400;
    int DiPFJetAve400_prescale;
    int DiPFJetAve500;
    int DiPFJetAve500_prescale;
    int DiPFJetAve60_HFJEC;
    int DiPFJetAve60_HFJEC_prescale;
    int DiPFJetAve80_HFJEC;
    int DiPFJetAve80_HFJEC_prescale;
    int DiPFJetAve100_HFJEC;
    int DiPFJetAve100_HFJEC_prescale;
    int DiPFJetAve160_HFJEC;
    int DiPFJetAve160_HFJEC_prescale;
    int DiPFJetAve220_HFJEC;
    int DiPFJetAve220_HFJEC_prescale;
    int DiPFJetAve300_HFJEC;
    int DiPFJetAve300_HFJEC_prescale;
    int AK8DiPFJet250_200_TrimMass30_BTagCSV0p45;
    int AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale;
    int AK8DiPFJet280_200_TrimMass30_BTagCSV0p45;
    int AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale;
    int AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45;
    int AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale;
    int AK8PFHT650_TrimR0p1PT0p03Mass50;
    int AK8PFHT650_TrimR0p1PT0p03Mass50_prescale;
    int AK8PFHT700_TrimR0p1PT0p03Mass50;
    int AK8PFHT700_TrimR0p1PT0p03Mass50_prescale;
    int PFHT550_4JetPt50;
    int PFHT550_4JetPt50_prescale;
    int PFHT650_4JetPt50;
    int PFHT650_4JetPt50_prescale;
    int PFHT750_4JetPt50;
    int PFHT750_4JetPt50_prescale;
    int ECALHT800;
    int ECALHT800_prescale;
    int PFHT600;
    int PFHT600_prescale;
    int PFHT650;
    int PFHT650_prescale;
    int PFHT800;
    int PFHT800_prescale;
    int PFHT200;
    int PFHT200_prescale;
    int PFHT250;
    int PFHT250_prescale;
    int PFHT300;
    int PFHT300_prescale;
    int PFHT350;
    int PFHT350_prescale;
    int PFHT400;
    int PFHT400_prescale;
    int PFHT475;
    int PFHT475_prescale;
    int Rsq0p25;
    int Rsq0p25_prescale;
    int Rsq0p30;
    int Rsq0p30_prescale;
    int RsqMR240_Rsq0p09_MR200;
    int RsqMR240_Rsq0p09_MR200_prescale;
    int RsqMR240_Rsq0p09_MR200_4jet;
    int RsqMR240_Rsq0p09_MR200_4jet_prescale;
    int RsqMR270_Rsq0p09_MR200;
    int RsqMR270_Rsq0p09_MR200_prescale;
    int RsqMR270_Rsq0p09_MR200_4jet;
    int RsqMR270_Rsq0p09_MR200_4jet_prescale;
    int Mu30_eta2p1_PFJet150_PFJet50;
    int Mu30_eta2p1_PFJet150_PFJet50_prescale;
    int Mu40_eta2p1_PFJet200_PFJet50;
    int Mu40_eta2p1_PFJet200_PFJet50_prescale;
    int Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;
    int Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale;
    int Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
    int Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale;
    int Ele8_CaloIdM_TrackIdM_PFJet30;
    int Ele8_CaloIdM_TrackIdM_PFJet30_prescale;
    int Ele12_CaloIdM_TrackIdM_PFJet30;
    int Ele12_CaloIdM_TrackIdM_PFJet30_prescale;
    int Ele18_CaloIdM_TrackIdM_PFJet30;
    int Ele18_CaloIdM_TrackIdM_PFJet30_prescale;
    int Ele23_CaloIdM_TrackIdM_PFJet30;
    int Ele23_CaloIdM_TrackIdM_PFJet30_prescale;
    int Ele33_CaloIdM_TrackIdM_PFJet30;
    int Ele33_CaloIdM_TrackIdM_PFJet30_prescale;
    int Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
    int Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;
    int Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
    int Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;
    int Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30;
    int Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;
    int Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07;
    int Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;
    int Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07;
    int Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;
    int Ele27_eta2p1_WPLoose_Gsf_HT200;
    int Ele27_eta2p1_WPLoose_Gsf_HT200_prescale;
    int Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF;
    int Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale;
    int Ele15_IsoVVVL_BTagCSV0p72_PFHT400;
    int Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;
    int Ele15_IsoVVVL_PFHT350_PFMET50;
    int Ele15_IsoVVVL_PFHT350_PFMET50_prescale;
    int Ele15_IsoVVVL_PFHT600;
    int Ele15_IsoVVVL_PFHT600_prescale;
    int Ele15_IsoVVVL_PFHT350;
    int Ele15_IsoVVVL_PFHT350_prescale;
    int Mu3er_PFHT140_PFMET125;
    int Mu3er_PFHT140_PFMET125_prescale;
    int Mu6_PFHT200_PFMET80_BTagCSV0p72;
    int Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale;
    int Mu6_PFHT200_PFMET100;
    int Mu6_PFHT200_PFMET100_prescale;
    int Mu10_CentralPFJet30_BTagCSV0p54PF;
    int Mu10_CentralPFJet30_BTagCSV0p54PF_prescale;
    int Mu15_IsoVVVL_BTagCSV0p72_PFHT400;
    int Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;
    int Mu15_IsoVVVL_PFHT350_PFMET50;
    int Mu15_IsoVVVL_PFHT350_PFMET50_prescale;
    int Mu15_IsoVVVL_PFHT600;
    int Mu15_IsoVVVL_PFHT600_prescale;
    int Mu15_IsoVVVL_PFHT350;
    int Mu15_IsoVVVL_PFHT350_prescale;
    int Mu8;
    int Mu8_prescale;
    int Mu17;
    int Mu17_prescale;
    int Mu20;
    int Mu20_prescale;
    int TkMu20;
    int TkMu20_prescale;
    int Mu24_eta2p1;
    int Mu24_eta2p1_prescale;
    int TkMu24_eta2p1;
    int TkMu24_eta2p1_prescale;
    int Mu27;
    int Mu27_prescale;
    int TkMu27;
    int TkMu27_prescale;
    int Mu45_eta2p1;
    int Mu45_eta2p1_prescale;
    int Mu50;
    int Mu50_prescale;
    int Mu55;
    int Mu55_prescale;
    int Mu300;
    int Mu300_prescale;
    int Mu350;
    int Mu350_prescale;
    int Ele105_CaloIdVT_GsfTrkIdT;
    int Ele105_CaloIdVT_GsfTrkIdT_prescale;
    int Ele115_CaloIdVT_GsfTrkIdT;
    int Ele115_CaloIdVT_GsfTrkIdT_prescale;
    int IsoMu17_eta2p1;
    int IsoMu17_eta2p1_prescale;
    int IsoMu18;
    int IsoMu18_prescale;
    int OldIsoMu18;
    int OldIsoMu18_prescale;
    int IsoTkMu18;
    int IsoTkMu18_prescale;
    int OldIsoTkMu18;
    int OldIsoTkMu18_prescale;
    int IsoMu20;
    int IsoMu20_prescale;
    int IsoMu22;
    int IsoMu22_prescale;
    int IsoMu27;
    int IsoMu27_prescale;
    int IsoTkMu24_eta2p1;
    int IsoTkMu24_eta2p1_prescale;
    int IsoTkMu27;
    int IsoTkMu27_prescale;
    int Mu8_TrkIsoVVL;
    int Mu8_TrkIsoVVL_prescale;
    int Mu17_TrkIsoVVL;
    int Mu17_TrkIsoVVL_prescale;
    int Ele22_eta2p1_WPLoose_Gsf;
    int Ele22_eta2p1_WPLoose_Gsf_prescale;
    int Ele22_eta2p1_WPTight_Gsf;
    int Ele22_eta2p1_WPTight_Gsf_prescale;
    int Ele23_WPLoose_Gsf;
    int Ele23_WPLoose_Gsf_prescale;
    int Ele27_WPLoose_Gsf;
    int Ele27_WPLoose_Gsf_prescale;
    int Ele27_eta2p1_WPLoose_Gsf;
    int Ele27_eta2p1_WPLoose_Gsf_prescale;
    int Ele27_eta2p1_WPTight_Gsf;
    int Ele27_eta2p1_WPTight_Gsf_prescale;
    int Ele32_eta2p1_WPTight_Gsf;
    int Ele32_eta2p1_WPTight_Gsf_prescale;
    int Ele12_CaloIdL_TrackIdL_IsoVL;
    int Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
    int Ele17_CaloIdL_TrackIdL_IsoVL;
    int Ele17_CaloIdL_TrackIdL_IsoVL_prescale;
    int Ele23_CaloIdL_TrackIdL_IsoVL;
    int Ele23_CaloIdL_TrackIdL_IsoVL_prescale;
    
    void init() {
      AK8PFJet360_TrimMass30=NOVAL_I;
      AK8PFJet360_TrimMass30_prescale=NOVAL_I;
      PFJet40=NOVAL_I;
      PFJet40_prescale=NOVAL_I;
      PFJet60=NOVAL_I;
      PFJet60_prescale=NOVAL_I;
      PFJet80=NOVAL_I;
      PFJet80_prescale=NOVAL_I;
      PFJet140=NOVAL_I;
      PFJet140_prescale=NOVAL_I;
      PFJet200=NOVAL_I;
      PFJet200_prescale=NOVAL_I;
      PFJet260=NOVAL_I;
      PFJet260_prescale=NOVAL_I;
      PFJet320=NOVAL_I;
      PFJet320_prescale=NOVAL_I;
      PFJet400=NOVAL_I;
      PFJet400_prescale=NOVAL_I;
      PFJet450=NOVAL_I;
      PFJet450_prescale=NOVAL_I;
      PFJet500=NOVAL_I;
      PFJet500_prescale=NOVAL_I;
      DiPFJetAve40=NOVAL_I;
      DiPFJetAve40_prescale=NOVAL_I;
      DiPFJetAve60=NOVAL_I;
      DiPFJetAve60_prescale=NOVAL_I;
      DiPFJetAve80=NOVAL_I;
      DiPFJetAve80_prescale=NOVAL_I;
      DiPFJetAve140=NOVAL_I;
      DiPFJetAve140_prescale=NOVAL_I;
      DiPFJetAve200=NOVAL_I;
      DiPFJetAve200_prescale=NOVAL_I;
      DiPFJetAve260=NOVAL_I;
      DiPFJetAve260_prescale=NOVAL_I;
      DiPFJetAve320=NOVAL_I;
      DiPFJetAve320_prescale=NOVAL_I;
      DiPFJetAve400=NOVAL_I;
      DiPFJetAve400_prescale=NOVAL_I;
      DiPFJetAve500=NOVAL_I;
      DiPFJetAve500_prescale=NOVAL_I;
      DiPFJetAve60_HFJEC=NOVAL_I;
      DiPFJetAve60_HFJEC_prescale=NOVAL_I;
      DiPFJetAve80_HFJEC=NOVAL_I;
      DiPFJetAve80_HFJEC_prescale=NOVAL_I;
      DiPFJetAve100_HFJEC=NOVAL_I;
      DiPFJetAve100_HFJEC_prescale=NOVAL_I;
      DiPFJetAve160_HFJEC=NOVAL_I;
      DiPFJetAve160_HFJEC_prescale=NOVAL_I;
      DiPFJetAve220_HFJEC=NOVAL_I;
      DiPFJetAve220_HFJEC_prescale=NOVAL_I;
      DiPFJetAve300_HFJEC=NOVAL_I;
      DiPFJetAve300_HFJEC_prescale=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_BTagCSV0p45=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_BTagCSV0p45=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale=NOVAL_I;
      AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45=NOVAL_I;
      AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale=NOVAL_I;
      AK8PFHT650_TrimR0p1PT0p03Mass50=NOVAL_I;
      AK8PFHT650_TrimR0p1PT0p03Mass50_prescale=NOVAL_I;
      AK8PFHT700_TrimR0p1PT0p03Mass50=NOVAL_I;
      AK8PFHT700_TrimR0p1PT0p03Mass50_prescale=NOVAL_I;
      PFHT550_4JetPt50=NOVAL_I;
      PFHT550_4JetPt50_prescale=NOVAL_I;
      PFHT650_4JetPt50=NOVAL_I;
      PFHT650_4JetPt50_prescale=NOVAL_I;
      PFHT750_4JetPt50=NOVAL_I;
      PFHT750_4JetPt50_prescale=NOVAL_I;
      ECALHT800=NOVAL_I;
      ECALHT800_prescale=NOVAL_I;
      PFHT600=NOVAL_I;
      PFHT600_prescale=NOVAL_I;
      PFHT650=NOVAL_I;
      PFHT650_prescale=NOVAL_I;
      PFHT800=NOVAL_I;
      PFHT800_prescale=NOVAL_I;
      PFHT200=NOVAL_I;
      PFHT200_prescale=NOVAL_I;
      PFHT250=NOVAL_I;
      PFHT250_prescale=NOVAL_I;
      PFHT300=NOVAL_I;
      PFHT300_prescale=NOVAL_I;
      PFHT350=NOVAL_I;
      PFHT350_prescale=NOVAL_I;
      PFHT400=NOVAL_I;
      PFHT400_prescale=NOVAL_I;
      PFHT475=NOVAL_I;
      PFHT475_prescale=NOVAL_I;
      Rsq0p25=NOVAL_I;
      Rsq0p25_prescale=NOVAL_I;
      Rsq0p30=NOVAL_I;
      Rsq0p30_prescale=NOVAL_I;
      RsqMR240_Rsq0p09_MR200=NOVAL_I;
      RsqMR240_Rsq0p09_MR200_prescale=NOVAL_I;
      RsqMR240_Rsq0p09_MR200_4jet=NOVAL_I;
      RsqMR240_Rsq0p09_MR200_4jet_prescale=NOVAL_I;
      RsqMR270_Rsq0p09_MR200=NOVAL_I;
      RsqMR270_Rsq0p09_MR200_prescale=NOVAL_I;
      RsqMR270_Rsq0p09_MR200_4jet=NOVAL_I;
      RsqMR270_Rsq0p09_MR200_4jet_prescale=NOVAL_I;
      Mu30_eta2p1_PFJet150_PFJet50=NOVAL_I;
      Mu30_eta2p1_PFJet150_PFJet50_prescale=NOVAL_I;
      Mu40_eta2p1_PFJet200_PFJet50=NOVAL_I;
      Mu40_eta2p1_PFJet200_PFJet50_prescale=NOVAL_I;
      Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50=NOVAL_I;
      Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale=NOVAL_I;
      Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50=NOVAL_I;
      Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale=NOVAL_I;
      Ele8_CaloIdM_TrackIdM_PFJet30=NOVAL_I;
      Ele8_CaloIdM_TrackIdM_PFJet30_prescale=NOVAL_I;
      Ele12_CaloIdM_TrackIdM_PFJet30=NOVAL_I;
      Ele12_CaloIdM_TrackIdM_PFJet30_prescale=NOVAL_I;
      Ele18_CaloIdM_TrackIdM_PFJet30=NOVAL_I;
      Ele18_CaloIdM_TrackIdM_PFJet30_prescale=NOVAL_I;
      Ele23_CaloIdM_TrackIdM_PFJet30=NOVAL_I;
      Ele23_CaloIdM_TrackIdM_PFJet30_prescale=NOVAL_I;
      Ele33_CaloIdM_TrackIdM_PFJet30=NOVAL_I;
      Ele33_CaloIdM_TrackIdM_PFJet30_prescale=NOVAL_I;
      Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30=NOVAL_I;
      Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale=NOVAL_I;
      Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30=NOVAL_I;
      Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale=NOVAL_I;
      Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30=NOVAL_I;
      Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale=NOVAL_I;
      Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07=NOVAL_I;
      Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale=NOVAL_I;
      Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07=NOVAL_I;
      Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf_HT200=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf_HT200_prescale=NOVAL_I;
      Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF=NOVAL_I;
      Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale=NOVAL_I;
      Ele15_IsoVVVL_BTagCSV0p72_PFHT400=NOVAL_I;
      Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT350_PFMET50=NOVAL_I;
      Ele15_IsoVVVL_PFHT350_PFMET50_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT600=NOVAL_I;
      Ele15_IsoVVVL_PFHT600_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT350=NOVAL_I;
      Ele15_IsoVVVL_PFHT350_prescale=NOVAL_I;
      Mu3er_PFHT140_PFMET125=NOVAL_I;
      Mu3er_PFHT140_PFMET125_prescale=NOVAL_I;
      Mu6_PFHT200_PFMET80_BTagCSV0p72=NOVAL_I;
      Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale=NOVAL_I;
      Mu6_PFHT200_PFMET100=NOVAL_I;
      Mu6_PFHT200_PFMET100_prescale=NOVAL_I;
      Mu10_CentralPFJet30_BTagCSV0p54PF=NOVAL_I;
      Mu10_CentralPFJet30_BTagCSV0p54PF_prescale=NOVAL_I;
      Mu15_IsoVVVL_BTagCSV0p72_PFHT400=NOVAL_I;
      Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale=NOVAL_I;
      Mu15_IsoVVVL_PFHT350_PFMET50=NOVAL_I;
      Mu15_IsoVVVL_PFHT350_PFMET50_prescale=NOVAL_I;
      Mu15_IsoVVVL_PFHT600=NOVAL_I;
      Mu15_IsoVVVL_PFHT600_prescale=NOVAL_I;
      Mu15_IsoVVVL_PFHT350=NOVAL_I;
      Mu15_IsoVVVL_PFHT350_prescale=NOVAL_I;
      Mu8=NOVAL_I;
      Mu8_prescale=NOVAL_I;
      Mu17=NOVAL_I;
      Mu17_prescale=NOVAL_I;
      Mu20=NOVAL_I;
      Mu20_prescale=NOVAL_I;
      TkMu20=NOVAL_I;
      TkMu20_prescale=NOVAL_I;
      Mu24_eta2p1=NOVAL_I;
      Mu24_eta2p1_prescale=NOVAL_I;
      TkMu24_eta2p1=NOVAL_I;
      TkMu24_eta2p1_prescale=NOVAL_I;
      Mu27=NOVAL_I;
      Mu27_prescale=NOVAL_I;
      TkMu27=NOVAL_I;
      TkMu27_prescale=NOVAL_I;
      Mu45_eta2p1=NOVAL_I;
      Mu45_eta2p1_prescale=NOVAL_I;
      Mu50=NOVAL_I;
      Mu50_prescale=NOVAL_I;
      Mu55=NOVAL_I;
      Mu55_prescale=NOVAL_I;
      Mu300=NOVAL_I;
      Mu300_prescale=NOVAL_I;
      Mu350=NOVAL_I;
      Mu350_prescale=NOVAL_I;
      Ele105_CaloIdVT_GsfTrkIdT=NOVAL_I;
      Ele105_CaloIdVT_GsfTrkIdT_prescale=NOVAL_I;
      Ele115_CaloIdVT_GsfTrkIdT=NOVAL_I;
      Ele115_CaloIdVT_GsfTrkIdT_prescale=NOVAL_I;
      IsoMu17_eta2p1=NOVAL_I;
      IsoMu17_eta2p1_prescale=NOVAL_I;
      IsoMu18=NOVAL_I;
      IsoMu18_prescale=NOVAL_I;
      OldIsoMu18=NOVAL_I;
      OldIsoMu18_prescale=NOVAL_I;
      IsoTkMu18=NOVAL_I;
      IsoTkMu18_prescale=NOVAL_I;
      OldIsoTkMu18=NOVAL_I;
      OldIsoTkMu18_prescale=NOVAL_I;
      IsoMu20=NOVAL_I;
      IsoMu20_prescale=NOVAL_I;
      IsoMu22=NOVAL_I;
      IsoMu22_prescale=NOVAL_I;
      IsoMu27=NOVAL_I;
      IsoMu27_prescale=NOVAL_I;
      IsoTkMu24_eta2p1=NOVAL_I;
      IsoTkMu24_eta2p1_prescale=NOVAL_I;
      IsoTkMu27=NOVAL_I;
      IsoTkMu27_prescale=NOVAL_I;
      Mu8_TrkIsoVVL=NOVAL_I;
      Mu8_TrkIsoVVL_prescale=NOVAL_I;
      Mu17_TrkIsoVVL=NOVAL_I;
      Mu17_TrkIsoVVL_prescale=NOVAL_I;
      Ele22_eta2p1_WPLoose_Gsf=NOVAL_I;
      Ele22_eta2p1_WPLoose_Gsf_prescale=NOVAL_I;
      Ele22_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele22_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
      Ele23_WPLoose_Gsf=NOVAL_I;
      Ele23_WPLoose_Gsf_prescale=NOVAL_I;
      Ele27_WPLoose_Gsf=NOVAL_I;
      Ele27_WPLoose_Gsf_prescale=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf_prescale=NOVAL_I;
      Ele27_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele27_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
      Ele32_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele32_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
      Ele12_CaloIdL_TrackIdL_IsoVL=NOVAL_I;
      Ele12_CaloIdL_TrackIdL_IsoVL_prescale=NOVAL_I;
      Ele17_CaloIdL_TrackIdL_IsoVL=NOVAL_I;
      Ele17_CaloIdL_TrackIdL_IsoVL_prescale=NOVAL_I;
      Ele23_CaloIdL_TrackIdL_IsoVL=NOVAL_I;
      Ele23_CaloIdL_TrackIdL_IsoVL_prescale=NOVAL_I;
    }
    
  } hlt;
  
  class GenVars {
  public:
    GenVars() { init(); };
    
    unsigned int size;
    std::vector<int> ID;
    std::vector<int> Status;
    std::vector<int> Mom0ID;
    std::vector<int> Mom0Status;
    std::vector<int> Mom1ID;
    std::vector<int> Mom1Status;
    std::vector<int> Dau0ID;
    std::vector<int> Dau0Status;
    std::vector<int> Dau1ID;
    std::vector<int> Dau1Status;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Mass;
    std::vector<float> Charge;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(ID);
      init_vec(Status);
      init_vec(Mom0ID);
      init_vec(Mom0Status);
      init_vec(Mom1ID);
      init_vec(Mom1Status);
      init_vec(Dau0ID);
      init_vec(Dau0Status);
      init_vec(Dau1ID);
      init_vec(Dau1Status);
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Mass);
      init_vec(Charge);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } gen;
  
  class ElectronVars {
  public:
    ElectronVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> Key;
    std::vector<float> Iso03;
    std::vector<float> Iso03db;
    std::vector<float> MiniIso;
    std::vector<float> rho;
    std::vector<float> EA;
    std::vector<float> sumChargedHadronPt;
    std::vector<float> sumNeutralHadronEt;
    std::vector<float> sumPhotonEt;
    std::vector<float> sumPUPt;
    std::vector<float> D0;
    std::vector<float> Dz;
    std::vector<float> dEtaIn;
    std::vector<float> dPhiIn;
    std::vector<float> HoE;
    std::vector<float> full5x5siee;
    std::vector<float> ooEmooP;
    std::vector<float> missHits;
    std::vector<float> hasMatchedConVeto;
    std::vector<float> vidVeto;
    std::vector<float> vidLoose;
    std::vector<float> vidTight;
    std::vector<float> vidMedium;
    std::vector<float> vidHEEP;
    std::vector<float> vidHEEPnoiso;
    std::vector<float> SCEta;
    std::vector<float> SCPhi;
    std::vector<int> IsPartOfNearAK4Jet;
    std::vector<int> IsPartOfNearAK8Jet;
    std::vector<int> IsPartOfNearSubjet;
    std::vector<float> DRNearGenEleFromSLTop;
    std::vector<float> PtNearGenEleFromSLTop;
    std::vector<float> PtNearGenTop;
    std::vector<float> LepAK4JetFrac;
    std::vector<float> LepAK8JetFrac;
    std::vector<float> LepSubjetFrac;
    std::vector<float> LepAK4JetMassDrop;
    std::vector<float> LepAK8JetMassDrop;
    std::vector<float> LepSubjetMassDrop;
    std::vector<float> AK4JetV1DR;
    std::vector<float> AK4JetV2DR;
    std::vector<float> AK4JetV3DR;
    std::vector<float> AK8JetV1DR;
    std::vector<float> AK8JetV2DR;
    std::vector<float> AK8JetV3DR;
    std::vector<float> SubjetV1DR;
    std::vector<float> SubjetV2DR;
    std::vector<float> SubjetV3DR;
    std::vector<float> AK4JetV1PtRel;
    std::vector<float> AK4JetV2PtRel;
    std::vector<float> AK4JetV3PtRel;
    std::vector<float> AK8JetV1PtRel;
    std::vector<float> AK8JetV2PtRel;
    std::vector<float> AK8JetV3PtRel;
    std::vector<float> SubjetV1PtRel;
    std::vector<float> SubjetV2PtRel;
    std::vector<float> SubjetV3PtRel;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(Key);
      init_vec(Iso03);
      init_vec(Iso03db);
      init_vec(MiniIso);
      init_vec(rho);
      init_vec(EA);
      init_vec(sumChargedHadronPt);
      init_vec(sumNeutralHadronEt);
      init_vec(sumPhotonEt);
      init_vec(sumPUPt);
      init_vec(D0);
      init_vec(Dz);
      init_vec(dEtaIn);
      init_vec(dPhiIn);
      init_vec(HoE);
      init_vec(full5x5siee);
      init_vec(ooEmooP);
      init_vec(missHits);
      init_vec(hasMatchedConVeto);
      init_vec(vidVeto);
      init_vec(vidLoose);
      init_vec(vidTight);
      init_vec(vidMedium);
      init_vec(vidHEEP);
      init_vec(vidHEEPnoiso);
      init_vec(SCEta);
      init_vec(SCPhi);
      init_vec(IsPartOfNearAK4Jet);
      init_vec(IsPartOfNearAK8Jet);
      init_vec(IsPartOfNearSubjet);
      init_vec(DRNearGenEleFromSLTop);
      init_vec(PtNearGenEleFromSLTop);
      init_vec(PtNearGenTop);
      init_vec(LepAK4JetFrac);
      init_vec(LepAK8JetFrac);
      init_vec(LepSubjetFrac);
      init_vec(LepAK4JetMassDrop);
      init_vec(LepAK8JetMassDrop);
      init_vec(LepSubjetMassDrop);
      init_vec(AK4JetV1DR);
      init_vec(AK4JetV2DR);
      init_vec(AK4JetV3DR);
      init_vec(AK8JetV1DR);
      init_vec(AK8JetV2DR);
      init_vec(AK8JetV3DR);
      init_vec(SubjetV1DR);
      init_vec(SubjetV2DR);
      init_vec(SubjetV3DR);
      init_vec(AK4JetV1PtRel);
      init_vec(AK4JetV2PtRel);
      init_vec(AK4JetV3PtRel);
      init_vec(AK8JetV1PtRel);
      init_vec(AK8JetV2PtRel);
      init_vec(AK8JetV3PtRel);
      init_vec(SubjetV1PtRel);
      init_vec(SubjetV2PtRel);
      init_vec(SubjetV3PtRel);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } ele;
  
  class MuonVars {
  public:
    MuonVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> Key;
    std::vector<float> Iso04;
    std::vector<float> MiniIso;
    std::vector<float> D0;
    std::vector<float> D0err;
    std::vector<float> Dxy;
    std::vector<float> Dxyerr;
    std::vector<float> Dz;
    std::vector<float> Dzerr;
    std::vector<float> IsSoftMuon;
    std::vector<float> IsLooseMuon;
    std::vector<float> IsMediumMuon;
    std::vector<float> IsTightMuon;
    std::vector<float> IsHighPtMuon;
    std::vector<float> IsPFMuon;
    std::vector<float> IsGlobalMuon;
    std::vector<float> IsTrackerMuon;
    std::vector<float> GlbTrkNormChi2;
    std::vector<float> NumberValidMuonHits;
    std::vector<float> NumberMatchedStations;
    std::vector<float> NumberValidPixelHits;
    std::vector<float> NumberTrackerLayers;
    std::vector<float> NumberOfValidTrackerHits;
    std::vector<float> NumberOfPixelLayers;
    std::vector<float> InTrkNormChi2;
    std::vector<float> SumChargedHadronPt;
    std::vector<float> SumNeutralHadronPt;
    std::vector<float> SumPhotonPt;
    std::vector<float> SumPUPt;
    std::vector<float> GenMuonY;
    std::vector<float> GenMuonEta;
    std::vector<float> GenMuonPhi;
    std::vector<float> GenMuonPt;
    std::vector<float> GenMuonE;
    std::vector<float> GenMuonCharge;
    std::vector<int> IsPartOfNearAK4Jet;
    std::vector<int> IsPartOfNearAK8Jet;
    std::vector<int> IsPartOfNearSubjet;
    std::vector<float> DRNearGenMuFromSLTop;
    std::vector<float> PtNearGenMuFromSLTop;
    std::vector<float> PtNearGenTop;
    std::vector<float> LepAK4JetFrac;
    std::vector<float> LepAK8JetFrac;
    std::vector<float> LepSubjetFrac;
    std::vector<float> LepAK4JetMassDrop;
    std::vector<float> LepAK8JetMassDrop;
    std::vector<float> LepSubjetMassDrop;
    std::vector<float> AK4JetV1DR;
    std::vector<float> AK4JetV2DR;
    std::vector<float> AK4JetV3DR;
    std::vector<float> AK8JetV1DR;
    std::vector<float> AK8JetV2DR;
    std::vector<float> AK8JetV3DR;
    std::vector<float> SubjetV1DR;
    std::vector<float> SubjetV2DR;
    std::vector<float> SubjetV3DR;
    std::vector<float> AK4JetV1PtRel;
    std::vector<float> AK4JetV2PtRel;
    std::vector<float> AK4JetV3PtRel;
    std::vector<float> AK8JetV1PtRel;
    std::vector<float> AK8JetV2PtRel;
    std::vector<float> AK8JetV3PtRel;
    std::vector<float> SubjetV1PtRel;
    std::vector<float> SubjetV2PtRel;
    std::vector<float> SubjetV3PtRel;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(Key);
      init_vec(Iso04);
      init_vec(MiniIso);
      init_vec(D0);
      init_vec(D0err);
      init_vec(Dxy);
      init_vec(Dxyerr);
      init_vec(Dz);
      init_vec(Dzerr);
      init_vec(IsSoftMuon);
      init_vec(IsLooseMuon);
      init_vec(IsMediumMuon);
      init_vec(IsTightMuon);
      init_vec(IsHighPtMuon);
      init_vec(IsPFMuon);
      init_vec(IsGlobalMuon);
      init_vec(IsTrackerMuon);
      init_vec(GlbTrkNormChi2);
      init_vec(NumberValidMuonHits);
      init_vec(NumberMatchedStations);
      init_vec(NumberValidPixelHits);
      init_vec(NumberTrackerLayers);
      init_vec(NumberOfValidTrackerHits);
      init_vec(NumberOfPixelLayers);
      init_vec(InTrkNormChi2);
      init_vec(SumChargedHadronPt);
      init_vec(SumNeutralHadronPt);
      init_vec(SumPhotonPt);
      init_vec(SumPUPt);
      init_vec(GenMuonY);
      init_vec(GenMuonEta);
      init_vec(GenMuonPhi);
      init_vec(GenMuonPt);
      init_vec(GenMuonE);
      init_vec(GenMuonCharge);
      init_vec(IsPartOfNearAK4Jet);
      init_vec(IsPartOfNearAK8Jet);
      init_vec(IsPartOfNearSubjet);
      init_vec(DRNearGenMuFromSLTop);
      init_vec(PtNearGenMuFromSLTop);
      init_vec(PtNearGenTop);
      init_vec(LepAK4JetFrac);
      init_vec(LepAK8JetFrac);
      init_vec(LepSubjetFrac);
      init_vec(LepAK4JetMassDrop);
      init_vec(LepAK8JetMassDrop);
      init_vec(LepSubjetMassDrop);
      init_vec(AK4JetV1DR);
      init_vec(AK4JetV2DR);
      init_vec(AK4JetV3DR);
      init_vec(AK8JetV1DR);
      init_vec(AK8JetV2DR);
      init_vec(AK8JetV3DR);
      init_vec(SubjetV1DR);
      init_vec(SubjetV2DR);
      init_vec(SubjetV3DR);
      init_vec(AK4JetV1PtRel);
      init_vec(AK4JetV2PtRel);
      init_vec(AK4JetV3PtRel);
      init_vec(AK8JetV1PtRel);
      init_vec(AK8JetV2PtRel);
      init_vec(AK8JetV3PtRel);
      init_vec(SubjetV1PtRel);
      init_vec(SubjetV2PtRel);
      init_vec(SubjetV3PtRel);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } mu;
  
  class AK4CHSJetVars {
  public:
    AK4CHSJetVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> CSVv2;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> CMVA;
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> neutralHadronMultiplicity;
    std::vector<float> neutralHadronEnergy;
    std::vector<float> neutralEmEnergy;
    std::vector<float> chargedEmEnergy;
    std::vector<float> chargedHadronEnergy;
    std::vector<float> photonMultiplicity;
    std::vector<float> electronMultiplicity;
    std::vector<float> HFHadronMultiplicity;
    std::vector<float> HFEMMultiplicity;
    std::vector<float> ChargeMuEnergy;
    std::vector<float> neutralMultiplicity;
    std::vector<float> neutralHadronEnergyFrac;
    std::vector<float> neutralEmEnergyFrac;
    std::vector<float> chargedHadronEnergyFrac;
    std::vector<float> muonEnergyFrac;
    std::vector<float> chargedEmEnergyFrac;
    std::vector<float> chargedMultiplicity;
    std::vector<float> NumConstituents;
    std::vector<float> jecFactor0;
    std::vector<float> jecFactorL3Absolute;
    std::vector<float> jetArea;
    std::vector<float> nSV;
    std::vector<float> SV0mass;
    std::vector<float> SV1mass;
    std::vector<float> QGL;
    std::vector<float> jecUncertainty;
    std::vector<float> PtResolution;
    std::vector<float> JERSF;
    std::vector<float> JERSFUp;
    std::vector<float> JERSFDown;
    std::vector<float> SmearedPt;
    std::vector<float> SmearedPEta;
    std::vector<float> SmearedPhi;
    std::vector<float> SmearedE;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSVv2);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(CMVA);
      init_vec(GenPartonY);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetY);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(muonMultiplicity);
      init_vec(PhotonEnergy);
      init_vec(ElectronEnergy);
      init_vec(MuonEnergy);
      init_vec(HFHadronEnergy);
      init_vec(HFEMEnergy);
      init_vec(ChargedHadronMultiplicity);
      init_vec(numberOfDaughters);
      init_vec(neutralHadronMultiplicity);
      init_vec(neutralHadronEnergy);
      init_vec(neutralEmEnergy);
      init_vec(chargedEmEnergy);
      init_vec(chargedHadronEnergy);
      init_vec(photonMultiplicity);
      init_vec(electronMultiplicity);
      init_vec(HFHadronMultiplicity);
      init_vec(HFEMMultiplicity);
      init_vec(ChargeMuEnergy);
      init_vec(neutralMultiplicity);
      init_vec(neutralHadronEnergyFrac);
      init_vec(neutralEmEnergyFrac);
      init_vec(chargedHadronEnergyFrac);
      init_vec(muonEnergyFrac);
      init_vec(chargedEmEnergyFrac);
      init_vec(chargedMultiplicity);
      init_vec(NumConstituents);
      init_vec(jecFactor0);
      init_vec(jecFactorL3Absolute);
      init_vec(jetArea);
      init_vec(nSV);
      init_vec(SV0mass);
      init_vec(SV1mass);
      init_vec(QGL);
      init_vec(jecUncertainty);
      init_vec(PtResolution);
      init_vec(JERSF);
      init_vec(JERSFUp);
      init_vec(JERSFDown);
      init_vec(SmearedPt);
      init_vec(SmearedPEta);
      init_vec(SmearedPhi);
      init_vec(SmearedE);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } jetsAK4;
  
  class AK4PuppiJetVars {
  public:
    AK4PuppiJetVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> CSVv2;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> CMVA;
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> neutralHadronMultiplicity;
    std::vector<float> neutralHadronEnergy;
    std::vector<float> neutralEmEnergy;
    std::vector<float> chargedEmEnergy;
    std::vector<float> chargedHadronEnergy;
    std::vector<float> photonMultiplicity;
    std::vector<float> electronMultiplicity;
    std::vector<float> HFHadronMultiplicity;
    std::vector<float> HFEMMultiplicity;
    std::vector<float> ChargeMuEnergy;
    std::vector<float> neutralMultiplicity;
    std::vector<float> neutralHadronEnergyFrac;
    std::vector<float> neutralEmEnergyFrac;
    std::vector<float> chargedHadronEnergyFrac;
    std::vector<float> muonEnergyFrac;
    std::vector<float> chargedEmEnergyFrac;
    std::vector<float> chargedMultiplicity;
    std::vector<float> NumConstituents;
    std::vector<float> jecFactor0;
    std::vector<float> jecFactorL3Absolute;
    std::vector<float> jetArea;
    std::vector<float> nSV;
    std::vector<float> SV0mass;
    std::vector<float> SV1mass;
    std::vector<float> jecUncertainty;
    std::vector<float> PtResolution;
    std::vector<float> JERSF;
    std::vector<float> JERSFUp;
    std::vector<float> JERSFDown;
    std::vector<float> SmearedPt;
    std::vector<float> SmearedPEta;
    std::vector<float> SmearedPhi;
    std::vector<float> SmearedE;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSVv2);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(CMVA);
      init_vec(GenPartonY);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetY);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(muonMultiplicity);
      init_vec(PhotonEnergy);
      init_vec(ElectronEnergy);
      init_vec(MuonEnergy);
      init_vec(HFHadronEnergy);
      init_vec(HFEMEnergy);
      init_vec(ChargedHadronMultiplicity);
      init_vec(numberOfDaughters);
      init_vec(neutralHadronMultiplicity);
      init_vec(neutralHadronEnergy);
      init_vec(neutralEmEnergy);
      init_vec(chargedEmEnergy);
      init_vec(chargedHadronEnergy);
      init_vec(photonMultiplicity);
      init_vec(electronMultiplicity);
      init_vec(HFHadronMultiplicity);
      init_vec(HFEMMultiplicity);
      init_vec(ChargeMuEnergy);
      init_vec(neutralMultiplicity);
      init_vec(neutralHadronEnergyFrac);
      init_vec(neutralEmEnergyFrac);
      init_vec(chargedHadronEnergyFrac);
      init_vec(muonEnergyFrac);
      init_vec(chargedEmEnergyFrac);
      init_vec(chargedMultiplicity);
      init_vec(NumConstituents);
      init_vec(jecFactor0);
      init_vec(jecFactorL3Absolute);
      init_vec(jetArea);
      init_vec(nSV);
      init_vec(SV0mass);
      init_vec(SV1mass);
      init_vec(jecUncertainty);
      init_vec(PtResolution);
      init_vec(JERSF);
      init_vec(JERSFUp);
      init_vec(JERSFDown);
      init_vec(SmearedPt);
      init_vec(SmearedPEta);
      init_vec(SmearedPhi);
      init_vec(SmearedE);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } jetsAK4Puppi;
  
  class AK8CHSJetVars {
  public:
    AK8CHSJetVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> CSVv2;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> CMVA;
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> neutralHadronMultiplicity;
    std::vector<float> neutralHadronEnergy;
    std::vector<float> neutralEmEnergy;
    std::vector<float> chargedEmEnergy;
    std::vector<float> chargedHadronEnergy;
    std::vector<float> photonMultiplicity;
    std::vector<float> electronMultiplicity;
    std::vector<float> HFHadronMultiplicity;
    std::vector<float> HFEMMultiplicity;
    std::vector<float> ChargeMuEnergy;
    std::vector<float> neutralMultiplicity;
    std::vector<float> neutralHadronEnergyFrac;
    std::vector<float> neutralEmEnergyFrac;
    std::vector<float> chargedHadronEnergyFrac;
    std::vector<float> muonEnergyFrac;
    std::vector<float> chargedEmEnergyFrac;
    std::vector<float> chargedMultiplicity;
    std::vector<float> NumConstituents;
    std::vector<float> jecFactor0;
    std::vector<float> jecFactorL3Absolute;
    std::vector<float> jetArea;
    std::vector<float> nSV;
    std::vector<float> SV0mass;
    std::vector<float> SV1mass;
    std::vector<float> jecUncertainty;
    std::vector<float> PtResolution;
    std::vector<float> JERSF;
    std::vector<float> JERSFUp;
    std::vector<float> JERSFDown;
    std::vector<float> SmearedPt;
    std::vector<float> SmearedPEta;
    std::vector<float> SmearedPhi;
    std::vector<float> SmearedE;
    std::vector<float> vSubjetIndex0;
    std::vector<float> vSubjetIndex1;
    std::vector<float> tau1;
    std::vector<float> tau2;
    std::vector<float> tau3;
    std::vector<float> softDropMass;
    std::vector<float> trimmedMass;
    std::vector<float> prunedMass;
    std::vector<float> filteredMass;
    std::vector<int> HasNearGenTop;
    std::vector<int> NearGenTopIsHadronic;
    std::vector<int> NearGenWIsHadronic;
    std::vector<int> NearGenWToENu;
    std::vector<int> NearGenWToMuNu;
    std::vector<int> NearGenWToTauNu;
    std::vector<int> PassTopTag;
    std::vector<float> DRNearGenTop;
    std::vector<float> DRNearGenWFromTop;
    std::vector<float> DRNearGenBFromTop;
    std::vector<float> DRNearGenLepFromSLTop;
    std::vector<float> DRNearGenNuFromSLTop;
    std::vector<float> PtNearGenTop;
    std::vector<float> PtNearGenBFromTop;
    std::vector<float> PtNearGenWFromTop;
    std::vector<float> PtNearGenLepFromSLTop;
    std::vector<float> PtNearGenNuFromSLTop;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSVv2);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(CMVA);
      init_vec(GenPartonY);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetY);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(muonMultiplicity);
      init_vec(PhotonEnergy);
      init_vec(ElectronEnergy);
      init_vec(MuonEnergy);
      init_vec(HFHadronEnergy);
      init_vec(HFEMEnergy);
      init_vec(ChargedHadronMultiplicity);
      init_vec(numberOfDaughters);
      init_vec(neutralHadronMultiplicity);
      init_vec(neutralHadronEnergy);
      init_vec(neutralEmEnergy);
      init_vec(chargedEmEnergy);
      init_vec(chargedHadronEnergy);
      init_vec(photonMultiplicity);
      init_vec(electronMultiplicity);
      init_vec(HFHadronMultiplicity);
      init_vec(HFEMMultiplicity);
      init_vec(ChargeMuEnergy);
      init_vec(neutralMultiplicity);
      init_vec(neutralHadronEnergyFrac);
      init_vec(neutralEmEnergyFrac);
      init_vec(chargedHadronEnergyFrac);
      init_vec(muonEnergyFrac);
      init_vec(chargedEmEnergyFrac);
      init_vec(chargedMultiplicity);
      init_vec(NumConstituents);
      init_vec(jecFactor0);
      init_vec(jecFactorL3Absolute);
      init_vec(jetArea);
      init_vec(nSV);
      init_vec(SV0mass);
      init_vec(SV1mass);
      init_vec(jecUncertainty);
      init_vec(PtResolution);
      init_vec(JERSF);
      init_vec(JERSFUp);
      init_vec(JERSFDown);
      init_vec(SmearedPt);
      init_vec(SmearedPEta);
      init_vec(SmearedPhi);
      init_vec(SmearedE);
      init_vec(vSubjetIndex0);
      init_vec(vSubjetIndex1);
      init_vec(tau1);
      init_vec(tau2);
      init_vec(tau3);
      init_vec(softDropMass);
      init_vec(trimmedMass);
      init_vec(prunedMass);
      init_vec(filteredMass);
      init_vec(HasNearGenTop);
      init_vec(NearGenTopIsHadronic);
      init_vec(NearGenWIsHadronic);
      init_vec(NearGenWToENu);
      init_vec(NearGenWToMuNu);
      init_vec(NearGenWToTauNu);
      init_vec(PassTopTag);
      init_vec(DRNearGenTop);
      init_vec(DRNearGenWFromTop);
      init_vec(DRNearGenBFromTop);
      init_vec(DRNearGenLepFromSLTop);
      init_vec(DRNearGenNuFromSLTop);
      init_vec(PtNearGenTop);
      init_vec(PtNearGenBFromTop);
      init_vec(PtNearGenWFromTop);
      init_vec(PtNearGenLepFromSLTop);
      init_vec(PtNearGenNuFromSLTop);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } jetsAK8;
  
  class AK8PuppiJetVars {
  public:
    AK8PuppiJetVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> CSVv2;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> CMVA;
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> neutralHadronMultiplicity;
    std::vector<float> neutralHadronEnergy;
    std::vector<float> neutralEmEnergy;
    std::vector<float> chargedEmEnergy;
    std::vector<float> chargedHadronEnergy;
    std::vector<float> photonMultiplicity;
    std::vector<float> electronMultiplicity;
    std::vector<float> HFHadronMultiplicity;
    std::vector<float> HFEMMultiplicity;
    std::vector<float> ChargeMuEnergy;
    std::vector<float> neutralMultiplicity;
    std::vector<float> neutralHadronEnergyFrac;
    std::vector<float> neutralEmEnergyFrac;
    std::vector<float> chargedHadronEnergyFrac;
    std::vector<float> muonEnergyFrac;
    std::vector<float> chargedEmEnergyFrac;
    std::vector<float> chargedMultiplicity;
    std::vector<float> NumConstituents;
    std::vector<float> jecFactor0;
    std::vector<float> jecFactorL3Absolute;
    std::vector<float> jetArea;
    std::vector<float> nSV;
    std::vector<float> SV0mass;
    std::vector<float> SV1mass;
    std::vector<float> jecUncertainty;
    std::vector<float> PtResolution;
    std::vector<float> JERSF;
    std::vector<float> JERSFUp;
    std::vector<float> JERSFDown;
    std::vector<float> SmearedPt;
    std::vector<float> SmearedPEta;
    std::vector<float> SmearedPhi;
    std::vector<float> SmearedE;
    std::vector<float> vSubjetIndex0;
    std::vector<float> vSubjetIndex1;
    std::vector<float> tau1;
    std::vector<float> tau2;
    std::vector<float> tau3;
    std::vector<float> softDropMass;
    std::vector<float> trimmedMass;
    std::vector<float> prunedMass;
    std::vector<float> filteredMass;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSVv2);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(CMVA);
      init_vec(GenPartonY);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetY);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(muonMultiplicity);
      init_vec(PhotonEnergy);
      init_vec(ElectronEnergy);
      init_vec(MuonEnergy);
      init_vec(HFHadronEnergy);
      init_vec(HFEMEnergy);
      init_vec(ChargedHadronMultiplicity);
      init_vec(numberOfDaughters);
      init_vec(neutralHadronMultiplicity);
      init_vec(neutralHadronEnergy);
      init_vec(neutralEmEnergy);
      init_vec(chargedEmEnergy);
      init_vec(chargedHadronEnergy);
      init_vec(photonMultiplicity);
      init_vec(electronMultiplicity);
      init_vec(HFHadronMultiplicity);
      init_vec(HFEMMultiplicity);
      init_vec(ChargeMuEnergy);
      init_vec(neutralMultiplicity);
      init_vec(neutralHadronEnergyFrac);
      init_vec(neutralEmEnergyFrac);
      init_vec(chargedHadronEnergyFrac);
      init_vec(muonEnergyFrac);
      init_vec(chargedEmEnergyFrac);
      init_vec(chargedMultiplicity);
      init_vec(NumConstituents);
      init_vec(jecFactor0);
      init_vec(jecFactorL3Absolute);
      init_vec(jetArea);
      init_vec(nSV);
      init_vec(SV0mass);
      init_vec(SV1mass);
      init_vec(jecUncertainty);
      init_vec(PtResolution);
      init_vec(JERSF);
      init_vec(JERSFUp);
      init_vec(JERSFDown);
      init_vec(SmearedPt);
      init_vec(SmearedPEta);
      init_vec(SmearedPhi);
      init_vec(SmearedE);
      init_vec(vSubjetIndex0);
      init_vec(vSubjetIndex1);
      init_vec(tau1);
      init_vec(tau2);
      init_vec(tau3);
      init_vec(softDropMass);
      init_vec(trimmedMass);
      init_vec(prunedMass);
      init_vec(filteredMass);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } jetsAK8Puppi;
  
  class AK8CHSSubjetVars {
  public:
    AK8CHSSubjetVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> CSVv2;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> CMVA;
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> neutralHadronMultiplicity;
    std::vector<float> neutralHadronEnergy;
    std::vector<float> neutralEmEnergy;
    std::vector<float> chargedEmEnergy;
    std::vector<float> chargedHadronEnergy;
    std::vector<float> photonMultiplicity;
    std::vector<float> electronMultiplicity;
    std::vector<float> HFHadronMultiplicity;
    std::vector<float> HFEMMultiplicity;
    std::vector<float> ChargeMuEnergy;
    std::vector<float> neutralMultiplicity;
    std::vector<float> neutralHadronEnergyFrac;
    std::vector<float> neutralEmEnergyFrac;
    std::vector<float> chargedHadronEnergyFrac;
    std::vector<float> muonEnergyFrac;
    std::vector<float> chargedEmEnergyFrac;
    std::vector<float> chargedMultiplicity;
    std::vector<float> NumConstituents;
    std::vector<float> jecFactor0;
    std::vector<float> jecFactorL3Absolute;
    std::vector<float> jetArea;
    std::vector<float> nSV;
    std::vector<float> SV0mass;
    std::vector<float> SV1mass;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSVv2);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(CMVA);
      init_vec(GenPartonY);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetY);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(muonMultiplicity);
      init_vec(PhotonEnergy);
      init_vec(ElectronEnergy);
      init_vec(MuonEnergy);
      init_vec(HFHadronEnergy);
      init_vec(HFEMEnergy);
      init_vec(ChargedHadronMultiplicity);
      init_vec(numberOfDaughters);
      init_vec(neutralHadronMultiplicity);
      init_vec(neutralHadronEnergy);
      init_vec(neutralEmEnergy);
      init_vec(chargedEmEnergy);
      init_vec(chargedHadronEnergy);
      init_vec(photonMultiplicity);
      init_vec(electronMultiplicity);
      init_vec(HFHadronMultiplicity);
      init_vec(HFEMMultiplicity);
      init_vec(ChargeMuEnergy);
      init_vec(neutralMultiplicity);
      init_vec(neutralHadronEnergyFrac);
      init_vec(neutralEmEnergyFrac);
      init_vec(chargedHadronEnergyFrac);
      init_vec(muonEnergyFrac);
      init_vec(chargedEmEnergyFrac);
      init_vec(chargedMultiplicity);
      init_vec(NumConstituents);
      init_vec(jecFactor0);
      init_vec(jecFactorL3Absolute);
      init_vec(jetArea);
      init_vec(nSV);
      init_vec(SV0mass);
      init_vec(SV1mass);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } subjetsAK8;
  
  class AK8PuppiSubjetVars {
  public:
    AK8PuppiSubjetVars() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> CSVv2;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> CMVA;
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> neutralHadronMultiplicity;
    std::vector<float> neutralHadronEnergy;
    std::vector<float> neutralEmEnergy;
    std::vector<float> chargedEmEnergy;
    std::vector<float> chargedHadronEnergy;
    std::vector<float> photonMultiplicity;
    std::vector<float> electronMultiplicity;
    std::vector<float> HFHadronMultiplicity;
    std::vector<float> HFEMMultiplicity;
    std::vector<float> ChargeMuEnergy;
    std::vector<float> neutralMultiplicity;
    std::vector<float> neutralHadronEnergyFrac;
    std::vector<float> neutralEmEnergyFrac;
    std::vector<float> chargedHadronEnergyFrac;
    std::vector<float> muonEnergyFrac;
    std::vector<float> chargedEmEnergyFrac;
    std::vector<float> chargedMultiplicity;
    std::vector<float> NumConstituents;
    std::vector<float> jecFactor0;
    std::vector<float> jecFactorL3Absolute;
    std::vector<float> jetArea;
    std::vector<float> nSV;
    std::vector<float> SV0mass;
    std::vector<float> SV1mass;
    
    unsigned int it;
    
    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSVv2);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(CMVA);
      init_vec(GenPartonY);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetY);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(muonMultiplicity);
      init_vec(PhotonEnergy);
      init_vec(ElectronEnergy);
      init_vec(MuonEnergy);
      init_vec(HFHadronEnergy);
      init_vec(HFEMEnergy);
      init_vec(ChargedHadronMultiplicity);
      init_vec(numberOfDaughters);
      init_vec(neutralHadronMultiplicity);
      init_vec(neutralHadronEnergy);
      init_vec(neutralEmEnergy);
      init_vec(chargedEmEnergy);
      init_vec(chargedHadronEnergy);
      init_vec(photonMultiplicity);
      init_vec(electronMultiplicity);
      init_vec(HFHadronMultiplicity);
      init_vec(HFEMMultiplicity);
      init_vec(ChargeMuEnergy);
      init_vec(neutralMultiplicity);
      init_vec(neutralHadronEnergyFrac);
      init_vec(neutralEmEnergyFrac);
      init_vec(chargedHadronEnergyFrac);
      init_vec(muonEnergyFrac);
      init_vec(chargedEmEnergyFrac);
      init_vec(chargedMultiplicity);
      init_vec(NumConstituents);
      init_vec(jecFactor0);
      init_vec(jecFactorL3Absolute);
      init_vec(jetArea);
      init_vec(nSV);
      init_vec(SV0mass);
      init_vec(SV1mass);
    }
    
    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }
    
  } subjetsAK8Puppi;

  class AK8GenJetVars {
  public:
    AK8GenJetVars() { init(); };

    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;

    unsigned int it;

    void init() {
      it = -1;
      size=9999;
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
    }

    bool Loop() {
      ++it;
      if (it<size) {
        return 1;
      } else {
        it=-1;
        return 0;
      }
    }

  } genjetsAK8;

};

#endif

