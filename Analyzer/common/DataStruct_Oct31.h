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
    
    unsigned int RunNumber;
    unsigned int LumiBlock;
    long EventNumber;
    int NGoodVtx;
    int LHA_PDF_ID;
    float MR;
    float MTR;
    float R;
    float R2;
    float AK8_MR;
    float AK8_MTR;
    float AK8_R;
    float AK8_R2;
    float XSec;
    float Gen_Weight;
    float Gen_Ht;
    float SUSY_Stop_Mass;
    float SUSY_Gluino_Mass;
    float SUSY_LSP_Mass;
    
    void init() {
      RunNumber=9999;
      LumiBlock=9999;
      EventNumber=9999;
      NGoodVtx=NOVAL_I;
      LHA_PDF_ID=NOVAL_I;
      MR=NOVAL_F;
      MTR=NOVAL_F;
      R=NOVAL_F;
      R2=NOVAL_F;
      AK8_MR=NOVAL_F;
      AK8_MTR=NOVAL_F;
      AK8_R=NOVAL_F;
      AK8_R2=NOVAL_F;
      XSec=NOVAL_F;
      Gen_Weight=NOVAL_F;
      Gen_Ht=NOVAL_F;
      SUSY_Stop_Mass=NOVAL_F;
      SUSY_Gluino_Mass=NOVAL_F;
      SUSY_LSP_Mass=NOVAL_F;
    }
    
  } evt;
  
  class METData {
  public:
    METData() { init(); };
    
    unsigned int size;
    std::vector<float> Pt;
    std::vector<float> Phi;
    
    void init() {
      size=9999;
      init_vec(Pt);
      init_vec(Phi);
    }
    
  } met;
  
  class PileupData {
  public:
    PileupData() { init(); };
    
    int NtrueInt;
    
    void init() {
      NtrueInt=NOVAL_I;
    }
    
  } pu;
  
  class VertexData {
  public:
    VertexData() { init(); };
    
    
    void init() {
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
    
    int BadPFMuonFilter;
    int BadChargedCandidateFilter;
    int HBHENoiseFilter;
    int HBHENoiseIsoFilter;
    int CSCTightHaloFilter;
    int CSCTightHaloTrkMuUnvetoFilter;
    int CSCTightHalo2015Filter;
    int globalTightHalo2016Filter;
    int globalSuperTightHalo2016Filter;
    int HcalStripHaloFilter;
    int hcalLaserEventFilter;
    int EcalDeadCellTriggerPrimitiveFilter;
    int EcalDeadCellBoundaryEnergyFilter;
    int goodVertices;
    int eeBadScFilter;
    int ecalLaserCorrFilter;
    int trkPOGFilters;
    int chargedHadronTrackResolutionFilter;
    int muonBadTrackFilter;
    int trkPOG_manystripclus53X;
    int trkPOG_toomanystripclus53X;
    int trkPOG_logErrorTooManyClusters;
    int METFilters;
    
    void init() {
      BadPFMuonFilter=NOVAL_I;
      BadChargedCandidateFilter=NOVAL_I;
      HBHENoiseFilter=NOVAL_I;
      HBHENoiseIsoFilter=NOVAL_I;
      CSCTightHaloFilter=NOVAL_I;
      CSCTightHaloTrkMuUnvetoFilter=NOVAL_I;
      CSCTightHalo2015Filter=NOVAL_I;
      globalTightHalo2016Filter=NOVAL_I;
      globalSuperTightHalo2016Filter=NOVAL_I;
      HcalStripHaloFilter=NOVAL_I;
      hcalLaserEventFilter=NOVAL_I;
      EcalDeadCellTriggerPrimitiveFilter=NOVAL_I;
      EcalDeadCellBoundaryEnergyFilter=NOVAL_I;
      goodVertices=NOVAL_I;
      eeBadScFilter=NOVAL_I;
      ecalLaserCorrFilter=NOVAL_I;
      trkPOGFilters=NOVAL_I;
      chargedHadronTrackResolutionFilter=NOVAL_I;
      muonBadTrackFilter=NOVAL_I;
      trkPOG_manystripclus53X=NOVAL_I;
      trkPOG_toomanystripclus53X=NOVAL_I;
      trkPOG_logErrorTooManyClusters=NOVAL_I;
      METFilters=NOVAL_I;
    }
    
  } filter;
  
  class HLTData {
  public:
    HLTData() { init(); };
    
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
    int AK8PFJet360_TrimMass30;
    int AK8PFJet360_TrimMass30_prescale;
    int AK8PFJet40;
    int AK8PFJet40_prescale;
    int AK8PFJet60;
    int AK8PFJet60_prescale;
    int AK8PFJet80;
    int AK8PFJet80_prescale;
    int AK8PFJet140;
    int AK8PFJet140_prescale;
    int AK8PFJet200;
    int AK8PFJet200_prescale;
    int AK8PFJet260;
    int AK8PFJet260_prescale;
    int AK8PFJet320;
    int AK8PFJet320_prescale;
    int AK8PFJet400;
    int AK8PFJet400_prescale;
    int AK8PFJet450;
    int AK8PFJet450_prescale;
    int AK8PFJet500;
    int AK8PFJet500_prescale;
    int AK8DiPFJet250_200_TrimMass30;
    int AK8DiPFJet250_200_TrimMass30_prescale;
    int AK8DiPFJet280_200_TrimMass30;
    int AK8DiPFJet280_200_TrimMass30_prescale;
    int AK8DiPFJet250_200_TrimMass30_BTagCSV_p20;
    int AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_prescale;
    int AK8DiPFJet280_200_TrimMass30_BTagCSV_p20;
    int AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_prescale;
    int AK8DiPFJet250_200_TrimMass30_BTagCSV0p45;
    int AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale;
    int AK8DiPFJet280_200_TrimMass30_BTagCSV0p45;
    int AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale;
    int AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v2;
    int AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v2_prescale;
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
    int PFHT600;
    int PFHT600_prescale;
    int PFHT650;
    int PFHT650_prescale;
    int PFHT800;
    int PFHT800_prescale;
    int PFHT900;
    int PFHT900_prescale;
    int MET200;
    int MET200_prescale;
    int MET250;
    int MET250_prescale;
    int MET300;
    int MET300_prescale;
    int MET600;
    int MET600_prescale;
    int MET700;
    int MET700_prescale;
    int PFMET300;
    int PFMET300_prescale;
    int PFMET400;
    int PFMET400_prescale;
    int PFMET500;
    int PFMET500_prescale;
    int PFMET600;
    int PFMET600_prescale;
    int PFHT300_PFMET100;
    int PFHT300_PFMET100_prescale;
    int PFHT300_PFMET110;
    int PFHT300_PFMET110_prescale;
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
    int Ele50_CaloIdVT_GsfTrkIdT_PFJet140;
    int Ele50_CaloIdVT_GsfTrkIdT_PFJet140_prescale;
    int Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
    int Ele50_CaloIdVT_GsfTrkIdT_PFJet165_prescale;
    int Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF;
    int Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale;
    int Ele15_IsoVVVL_BTagCSV0p72_PFHT400;
    int Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;
    int Ele15_IsoVVVL_PFHT350_PFMET50;
    int Ele15_IsoVVVL_PFHT350_PFMET50_prescale;
    int Ele15_IsoVVVL_PFHT400_PFMET50;
    int Ele15_IsoVVVL_PFHT400_PFMET50_prescale;
    int Ele15_IsoVVVL_PFHT600;
    int Ele15_IsoVVVL_PFHT600_prescale;
    int Ele15_IsoVVVL_PFHT350;
    int Ele15_IsoVVVL_PFHT350_prescale;
    int Ele15_IsoVVVL_PFHT400;
    int Ele15_IsoVVVL_PFHT400_prescale;
    int Ele50_IsoVVVL_PFHT400;
    int Ele50_IsoVVVL_PFHT400_prescale;
    int Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07;
    int Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;
    int Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07;
    int Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;
    int Ele27_eta2p1_WPLoose_Gsf_HT200;
    int Ele27_eta2p1_WPLoose_Gsf_HT200_prescale;
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
    int Mu15_IsoVVVL_PFHT400_PFMET50;
    int Mu15_IsoVVVL_PFHT400_PFMET50_prescale;
    int Mu15_IsoVVVL_PFHT600;
    int Mu15_IsoVVVL_PFHT600_prescale;
    int Mu15_IsoVVVL_PFHT350;
    int Mu15_IsoVVVL_PFHT350_prescale;
    int Mu15_IsoVVVL_PFHT400;
    int Mu15_IsoVVVL_PFHT400_prescale;
    int Mu50_IsoVVVL_PFHT400;
    int Mu50_IsoVVVL_PFHT400_prescale;
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
    int IsoMu18;
    int IsoMu18_prescale;
    int IsoTkMu18;
    int IsoTkMu18_prescale;
    int IsoMu20;
    int IsoMu20_prescale;
    int IsoTkMu20;
    int IsoTkMu20_prescale;
    int IsoMu22;
    int IsoMu22_prescale;
    int IsoTkMu22;
    int IsoTkMu22_prescale;
    int IsoMu24;
    int IsoMu24_prescale;
    int IsoTkMu24;
    int IsoTkMu24_prescale;
    int IsoTkMu24_eta2p1;
    int IsoTkMu24_eta2p1_prescale;
    int IsoMu27;
    int IsoMu27_prescale;
    int IsoTkMu27;
    int IsoTkMu27_prescale;
    int Ele22_eta2p1_WPLoose_Gsf;
    int Ele22_eta2p1_WPLoose_Gsf_prescale;
    int Ele22_eta2p1_WPTight_Gsf;
    int Ele22_eta2p1_WPTight_Gsf_prescale;
    int Ele23_WPLoose_Gsf;
    int Ele23_WPLoose_Gsf_prescale;
    int Ele24_WPLoose_Gsf;
    int Ele24_WPLoose_Gsf_prescale;
    int Ele24_eta2p1_WPLoose_Gsf;
    int Ele24_eta2p1_WPLoose_Gsf_prescale;
    int Ele25_eta2p1_WPLoose_Gsf;
    int Ele25_eta2p1_WPLoose_Gsf_prescale;
    int Ele25_eta2p1_WPTight_Gsf;
    int Ele25_eta2p1_WPTight_Gsf_prescale;
    int Ele27_WPLoose_Gsf;
    int Ele27_WPLoose_Gsf_prescale;
    int Ele27_WPTight_Gsf;
    int Ele27_WPTight_Gsf_prescale;
    int Ele27_eta2p1_WPLoose_Gsf;
    int Ele27_eta2p1_WPLoose_Gsf_prescale;
    int Ele27_eta2p1_WPTight_Gsf;
    int Ele27_eta2p1_WPTight_Gsf_prescale;
    int Ele32_eta2p1_WPTight_Gsf;
    int Ele32_eta2p1_WPTight_Gsf_prescale;
    
    void init() {
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
      AK8PFJet360_TrimMass30=NOVAL_I;
      AK8PFJet360_TrimMass30_prescale=NOVAL_I;
      AK8PFJet40=NOVAL_I;
      AK8PFJet40_prescale=NOVAL_I;
      AK8PFJet60=NOVAL_I;
      AK8PFJet60_prescale=NOVAL_I;
      AK8PFJet80=NOVAL_I;
      AK8PFJet80_prescale=NOVAL_I;
      AK8PFJet140=NOVAL_I;
      AK8PFJet140_prescale=NOVAL_I;
      AK8PFJet200=NOVAL_I;
      AK8PFJet200_prescale=NOVAL_I;
      AK8PFJet260=NOVAL_I;
      AK8PFJet260_prescale=NOVAL_I;
      AK8PFJet320=NOVAL_I;
      AK8PFJet320_prescale=NOVAL_I;
      AK8PFJet400=NOVAL_I;
      AK8PFJet400_prescale=NOVAL_I;
      AK8PFJet450=NOVAL_I;
      AK8PFJet450_prescale=NOVAL_I;
      AK8PFJet500=NOVAL_I;
      AK8PFJet500_prescale=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_prescale=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_prescale=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_BTagCSV_p20=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_prescale=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_BTagCSV_p20=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_prescale=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_BTagCSV0p45=NOVAL_I;
      AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_BTagCSV0p45=NOVAL_I;
      AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale=NOVAL_I;
      AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v2=NOVAL_I;
      AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v2_prescale=NOVAL_I;
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
      PFHT600=NOVAL_I;
      PFHT600_prescale=NOVAL_I;
      PFHT650=NOVAL_I;
      PFHT650_prescale=NOVAL_I;
      PFHT800=NOVAL_I;
      PFHT800_prescale=NOVAL_I;
      PFHT900=NOVAL_I;
      PFHT900_prescale=NOVAL_I;
      MET200=NOVAL_I;
      MET200_prescale=NOVAL_I;
      MET250=NOVAL_I;
      MET250_prescale=NOVAL_I;
      MET300=NOVAL_I;
      MET300_prescale=NOVAL_I;
      MET600=NOVAL_I;
      MET600_prescale=NOVAL_I;
      MET700=NOVAL_I;
      MET700_prescale=NOVAL_I;
      PFMET300=NOVAL_I;
      PFMET300_prescale=NOVAL_I;
      PFMET400=NOVAL_I;
      PFMET400_prescale=NOVAL_I;
      PFMET500=NOVAL_I;
      PFMET500_prescale=NOVAL_I;
      PFMET600=NOVAL_I;
      PFMET600_prescale=NOVAL_I;
      PFHT300_PFMET100=NOVAL_I;
      PFHT300_PFMET100_prescale=NOVAL_I;
      PFHT300_PFMET110=NOVAL_I;
      PFHT300_PFMET110_prescale=NOVAL_I;
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
      Ele50_CaloIdVT_GsfTrkIdT_PFJet140=NOVAL_I;
      Ele50_CaloIdVT_GsfTrkIdT_PFJet140_prescale=NOVAL_I;
      Ele50_CaloIdVT_GsfTrkIdT_PFJet165=NOVAL_I;
      Ele50_CaloIdVT_GsfTrkIdT_PFJet165_prescale=NOVAL_I;
      Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF=NOVAL_I;
      Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale=NOVAL_I;
      Ele15_IsoVVVL_BTagCSV0p72_PFHT400=NOVAL_I;
      Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT350_PFMET50=NOVAL_I;
      Ele15_IsoVVVL_PFHT350_PFMET50_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT400_PFMET50=NOVAL_I;
      Ele15_IsoVVVL_PFHT400_PFMET50_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT600=NOVAL_I;
      Ele15_IsoVVVL_PFHT600_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT350=NOVAL_I;
      Ele15_IsoVVVL_PFHT350_prescale=NOVAL_I;
      Ele15_IsoVVVL_PFHT400=NOVAL_I;
      Ele15_IsoVVVL_PFHT400_prescale=NOVAL_I;
      Ele50_IsoVVVL_PFHT400=NOVAL_I;
      Ele50_IsoVVVL_PFHT400_prescale=NOVAL_I;
      Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07=NOVAL_I;
      Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale=NOVAL_I;
      Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07=NOVAL_I;
      Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf_HT200=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf_HT200_prescale=NOVAL_I;
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
      Mu15_IsoVVVL_PFHT400_PFMET50=NOVAL_I;
      Mu15_IsoVVVL_PFHT400_PFMET50_prescale=NOVAL_I;
      Mu15_IsoVVVL_PFHT600=NOVAL_I;
      Mu15_IsoVVVL_PFHT600_prescale=NOVAL_I;
      Mu15_IsoVVVL_PFHT350=NOVAL_I;
      Mu15_IsoVVVL_PFHT350_prescale=NOVAL_I;
      Mu15_IsoVVVL_PFHT400=NOVAL_I;
      Mu15_IsoVVVL_PFHT400_prescale=NOVAL_I;
      Mu50_IsoVVVL_PFHT400=NOVAL_I;
      Mu50_IsoVVVL_PFHT400_prescale=NOVAL_I;
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
      IsoMu18=NOVAL_I;
      IsoMu18_prescale=NOVAL_I;
      IsoTkMu18=NOVAL_I;
      IsoTkMu18_prescale=NOVAL_I;
      IsoMu20=NOVAL_I;
      IsoMu20_prescale=NOVAL_I;
      IsoTkMu20=NOVAL_I;
      IsoTkMu20_prescale=NOVAL_I;
      IsoMu22=NOVAL_I;
      IsoMu22_prescale=NOVAL_I;
      IsoTkMu22=NOVAL_I;
      IsoTkMu22_prescale=NOVAL_I;
      IsoMu24=NOVAL_I;
      IsoMu24_prescale=NOVAL_I;
      IsoTkMu24=NOVAL_I;
      IsoTkMu24_prescale=NOVAL_I;
      IsoTkMu24_eta2p1=NOVAL_I;
      IsoTkMu24_eta2p1_prescale=NOVAL_I;
      IsoMu27=NOVAL_I;
      IsoMu27_prescale=NOVAL_I;
      IsoTkMu27=NOVAL_I;
      IsoTkMu27_prescale=NOVAL_I;
      Ele22_eta2p1_WPLoose_Gsf=NOVAL_I;
      Ele22_eta2p1_WPLoose_Gsf_prescale=NOVAL_I;
      Ele22_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele22_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
      Ele23_WPLoose_Gsf=NOVAL_I;
      Ele23_WPLoose_Gsf_prescale=NOVAL_I;
      Ele24_WPLoose_Gsf=NOVAL_I;
      Ele24_WPLoose_Gsf_prescale=NOVAL_I;
      Ele24_eta2p1_WPLoose_Gsf=NOVAL_I;
      Ele24_eta2p1_WPLoose_Gsf_prescale=NOVAL_I;
      Ele25_eta2p1_WPLoose_Gsf=NOVAL_I;
      Ele25_eta2p1_WPLoose_Gsf_prescale=NOVAL_I;
      Ele25_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele25_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
      Ele27_WPLoose_Gsf=NOVAL_I;
      Ele27_WPLoose_Gsf_prescale=NOVAL_I;
      Ele27_WPTight_Gsf=NOVAL_I;
      Ele27_WPTight_Gsf_prescale=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf=NOVAL_I;
      Ele27_eta2p1_WPLoose_Gsf_prescale=NOVAL_I;
      Ele27_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele27_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
      Ele32_eta2p1_WPTight_Gsf=NOVAL_I;
      Ele32_eta2p1_WPTight_Gsf_prescale=NOVAL_I;
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
    std::vector<float> Dxy;
    std::vector<float> Dz;
    std::vector<float> DB;
    std::vector<float> DBerr;
    std::vector<float> vidVeto;
    std::vector<float> vidLoose;
    std::vector<float> vidMedium;
    std::vector<float> vidTight;
    std::vector<float> vidHEEP;
    std::vector<float> vidVetonoiso;
    std::vector<float> vidLoosenoiso;
    std::vector<float> vidMediumnoiso;
    std::vector<float> vidTightnoiso;
    std::vector<float> vidHEEPnoiso;
    std::vector<int> IsPartOfNearAK4Jet;
    std::vector<int> IsPartOfNearAK8Jet;
    std::vector<int> IsPartOfNearSubjet;
    std::vector<int> IsoVeto;
    std::vector<int> IsoLoose;
    std::vector<int> IsoMedium;
    std::vector<int> IsoTight;
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
      init_vec(Dxy);
      init_vec(Dz);
      init_vec(DB);
      init_vec(DBerr);
      init_vec(vidVeto);
      init_vec(vidLoose);
      init_vec(vidMedium);
      init_vec(vidTight);
      init_vec(vidHEEP);
      init_vec(vidVetonoiso);
      init_vec(vidLoosenoiso);
      init_vec(vidMediumnoiso);
      init_vec(vidTightnoiso);
      init_vec(vidHEEPnoiso);
      init_vec(IsPartOfNearAK4Jet);
      init_vec(IsPartOfNearAK8Jet);
      init_vec(IsPartOfNearSubjet);
      init_vec(IsoVeto);
      init_vec(IsoLoose);
      init_vec(IsoMedium);
      init_vec(IsoTight);
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
    std::vector<float> Dxy;
    std::vector<float> Dz;
    std::vector<float> DB;
    std::vector<float> DBerr;
    std::vector<float> IsSoftMuon;
    std::vector<float> IsLooseMuon;
    std::vector<float> IsMediumMuon;
    std::vector<float> IsTightMuon;
    std::vector<float> IsHighPtMuon;
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
      init_vec(Dxy);
      init_vec(Dz);
      init_vec(DB);
      init_vec(DBerr);
      init_vec(IsSoftMuon);
      init_vec(IsLooseMuon);
      init_vec(IsMediumMuon);
      init_vec(IsTightMuon);
      init_vec(IsHighPtMuon);
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
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> jecFactor0;
    std::vector<float> jecUncertainty;
    std::vector<float> JERSF;
    std::vector<float> JERSFUp;
    std::vector<float> JERSFDown;
    std::vector<float> SmearedPt;
    std::vector<std::vector<int> > Keys;
    std::vector<int> looseJetID;
    std::vector<int> tightJetID;
    std::vector<int> tightLepVetoJetID;
    
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
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(jecFactor0);
      init_vec(jecUncertainty);
      init_vec(JERSF);
      init_vec(JERSFUp);
      init_vec(JERSFDown);
      init_vec(SmearedPt);
      init_vec(Keys);
      init_vec(looseJetID);
      init_vec(tightJetID);
      init_vec(tightLepVetoJetID);
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
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> jecFactor0;
    std::vector<float> jecUncertainty;
    std::vector<float> JERSF;
    std::vector<float> JERSFUp;
    std::vector<float> JERSFDown;
    std::vector<float> SmearedPt;
    std::vector<float> DoubleBAK8;
    std::vector<float> DoubleBCA15;
    std::vector<float> vSubjetIndex0;
    std::vector<float> vSubjetIndex1;
    std::vector<float> tau1;
    std::vector<float> tau2;
    std::vector<float> tau3;
    std::vector<float> softDropMass;
    std::vector<std::vector<int> > Keys;
    std::vector<int> HasNearGenTop;
    std::vector<int> NearGenTopIsHadronic;
    std::vector<int> NearGenWIsHadronic;
    std::vector<int> NearGenWToENu;
    std::vector<int> NearGenWToMuNu;
    std::vector<int> NearGenWToTauNu;
    std::vector<int> looseJetID;
    std::vector<int> tightJetID;
    std::vector<int> tightLepVetoJetID;
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
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(jecFactor0);
      init_vec(jecUncertainty);
      init_vec(JERSF);
      init_vec(JERSFUp);
      init_vec(JERSFDown);
      init_vec(SmearedPt);
      init_vec(DoubleBAK8);
      init_vec(DoubleBCA15);
      init_vec(vSubjetIndex0);
      init_vec(vSubjetIndex1);
      init_vec(tau1);
      init_vec(tau2);
      init_vec(tau3);
      init_vec(softDropMass);
      init_vec(Keys);
      init_vec(HasNearGenTop);
      init_vec(NearGenTopIsHadronic);
      init_vec(NearGenWIsHadronic);
      init_vec(NearGenWToENu);
      init_vec(NearGenWToMuNu);
      init_vec(NearGenWToTauNu);
      init_vec(looseJetID);
      init_vec(tightJetID);
      init_vec(tightLepVetoJetID);
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
    std::vector<float> CMVAv2;
    std::vector<float> CvsL;
    std::vector<float> CvsB;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    std::vector<float> jecFactor0;
    std::vector<std::vector<int> > Keys;
    std::vector<int> looseJetID;
    std::vector<int> tightJetID;
    std::vector<int> tightLepVetoJetID;
    
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
      init_vec(CMVAv2);
      init_vec(CvsL);
      init_vec(CvsB);
      init_vec(GenPartonEta);
      init_vec(GenPartonPhi);
      init_vec(GenPartonPt);
      init_vec(GenPartonE);
      init_vec(GenPartonCharge);
      init_vec(PartonFlavour);
      init_vec(HadronFlavour);
      init_vec(GenJetEta);
      init_vec(GenJetPhi);
      init_vec(GenJetPt);
      init_vec(GenJetE);
      init_vec(GenJetCharge);
      init_vec(jecFactor0);
      init_vec(Keys);
      init_vec(looseJetID);
      init_vec(tightJetID);
      init_vec(tightLepVetoJetID);
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

