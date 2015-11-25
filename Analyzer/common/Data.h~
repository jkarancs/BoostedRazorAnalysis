#ifndef Data_h
#define Data_h

#define NOVAL_I -9999
#define NOVAL_F -9999.0

#define NLEP 60
#define NJET 100
#define NGEN 500

#include <cassert>
#include <map>
#include <vector>
#include <iostream>

#include "TLorentzVector.h"

using namespace std;

class Data {
public:
  Data() {}
  ~Data() {}
  
  class GenVars {
  public:
    GenVars() { init(); }
    ~GenVars() {}
    
    // Basic
    float Mass[NGEN];
    float Pt[NGEN];
    float Eta[NGEN];
    float Y[NGEN];
    float Phi[NGEN];
    float E[NGEN];
    float Charge[NGEN];
    float _ID[NGEN];
    float _Status[NGEN];
    float _MomID[NGEN];
    int ID[NGEN];
    int Status[NGEN];
    int MomID[NGEN];
    
    size_t it;
    int size;
    
    void init() {
      for (size_t i=0; i<NGEN; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        _ID[i]=NOVAL_F;
	_Status[i]=NOVAL_F;
        _MomID[i]=NOVAL_F;
        ID[i]=NOVAL_F;
	Status[i]=NOVAL_F;
        MomID[i]=NOVAL_F;
      }
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
#if PHYS14 == 1
	ID[it]=_ID[it];
	Status[it]=_Status[it];
        MomID[it]=_MomID[it];
#endif
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
    
  } gen;
  
  class ElectronVars {
  public:
    ElectronVars() { init(); }
    ~ElectronVars() {}
    
    // Basic
    float Mass[NLEP];
    float Pt[NLEP];
    float Eta[NLEP];
    float Y[NLEP];
    float Phi[NLEP];
    float E[NLEP];
    float Charge[NLEP];
    // ElectronVars
    float Iso03[NLEP];
    float D0[NLEP];
    float Dz[NLEP];
    float dEtaIn[NLEP];
    float dPhiIn[NLEP];
    float HoE[NLEP];
    float full5x5siee[NLEP];
    float ooEmooP[NLEP];
    float missHits[NLEP];
    float hasMatchedConVeto[NLEP];
    float isEB[NLEP];
    float isVeto[NLEP];
    float isLoose[NLEP];
    float isTight[NLEP];
    float isMedium[NLEP];
    float scEta[NLEP];
    
    size_t it;
    int size;
    
    float DRNearGenEleFromSLTop[NLEP];
    float PtNearGenEleFromSLTop[NLEP];
    float PtNearGenTop[NLEP];
    int IsPartOfNearAK4Jet[NLEP];
    int IsPartOfNearAK8Jet[NLEP];
    int IsPartOfNearSubjet[NLEP];
    float LepAK4JetFrac[NLEP];
    float LepAK8JetFrac[NLEP];
    float LepSubjetFrac[NLEP];
    float LepAK4JetMassDrop[NLEP];
    float LepAK8JetMassDrop[NLEP];
    float LepSubjetMassDrop[NLEP];
    float AK4JetV1DR[NLEP];
    float AK4JetV2DR[NLEP];
    float AK4JetV3DR[NLEP];
    float AK8JetV1DR[NLEP];
    float AK8JetV2DR[NLEP];
    float AK8JetV3DR[NLEP];
    float SubjetV1DR[NLEP];
    float SubjetV2DR[NLEP];
    float SubjetV3DR[NLEP];
    float AK4JetV1PtRel[NLEP];
    float AK4JetV2PtRel[NLEP];
    float AK4JetV3PtRel[NLEP];
    float AK8JetV1PtRel[NLEP];
    float AK8JetV2PtRel[NLEP];
    float AK8JetV3PtRel[NLEP];
    float SubjetV1PtRel[NLEP];
    float SubjetV2PtRel[NLEP];
    float SubjetV3PtRel[NLEP];
    
    void init() {
      for (size_t i=0; i<NLEP; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        Iso03[i]=NOVAL_F;
        D0[i]=NOVAL_F;
        Dz[i]=NOVAL_F;
        dEtaIn[i]=NOVAL_F;
        dPhiIn[i]=NOVAL_F;
        HoE[i]=NOVAL_F;
        full5x5siee[i]=NOVAL_F;
        ooEmooP[i]=NOVAL_F;
        missHits[i]=NOVAL_F;
        hasMatchedConVeto[i]=NOVAL_F;
        isEB[i]=NOVAL_F;
        isVeto[i]=NOVAL_F;
        isLoose[i]=NOVAL_F;
        isTight[i]=NOVAL_F;
        isMedium[i]=NOVAL_F;
        scEta[i]=NOVAL_F;
      }
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
    
  } ele;

  class MuonVars {
  public:
    MuonVars() { init(); }
    ~MuonVars() {}
    
    // Basic
    float Mass[NLEP];
    float Pt[NLEP];
    float Eta[NLEP];
    float Y[NLEP];
    float Phi[NLEP];
    float E[NLEP];
    float Charge[NLEP];
    // MuonVars
    float Iso04[NLEP];
    float D0[NLEP];
    float D0err[NLEP];
    float Dxy[NLEP];
    float Dxyerr[NLEP];
    float Dz[NLEP];
    float Dzerr[NLEP];
    float IsLooseMuon[NLEP];
    float IsSoftMuon[NLEP];
    float IsTightMuon[NLEP];
    float IsPFMuon[NLEP];
    float IsGlobalMuon[NLEP];
    float IsTrackerMuon[NLEP];
    float GlbTrkNormChi2[NLEP];
    float NumberValidMuonHits[NLEP];
    float NumberMatchedStations[NLEP];
    float NumberValidPixelHits[NLEP];
    float NumberTrackerLayers[NLEP];
    float NumberOfValidTrackerHits[NLEP];
    float NumberOfPixelLayers[NLEP];
    float InTrkNormChi2[NLEP];
    float SumChargedHadronPt[NLEP];
    float SumNeutralHadronPt[NLEP];
    float SumPhotonPt[NLEP];
    float SumPUPt[NLEP];
    float GenMuonY[NLEP];
    float GenMuonEta[NLEP];
    float GenMuonPhi[NLEP];
    float GenMuonPt[NLEP];
    float GenMuonE[NLEP];
    float GenMuonCharge[NLEP];
    
    size_t it;
    int size;

    float DRNearGenMuFromSLTop[NLEP];
    float PtNearGenMuFromSLTop[NLEP];
    float PtNearGenTop[NLEP];
    int IsPartOfNearAK4Jet[NLEP];
    int IsPartOfNearAK8Jet[NLEP];
    int IsPartOfNearSubjet[NLEP];
    float LepAK4JetFrac[NLEP];
    float LepAK8JetFrac[NLEP];
    float LepSubjetFrac[NLEP];
    float LepAK4JetMassDrop[NLEP];
    float LepAK8JetMassDrop[NLEP];
    float LepSubjetMassDrop[NLEP];
    float AK4JetV1DR[NLEP];
    float AK4JetV2DR[NLEP];
    float AK4JetV3DR[NLEP];
    float AK8JetV1DR[NLEP];
    float AK8JetV2DR[NLEP];
    float AK8JetV3DR[NLEP];
    float SubjetV1DR[NLEP];
    float SubjetV2DR[NLEP];
    float SubjetV3DR[NLEP];
    float AK4JetV1PtRel[NLEP];
    float AK4JetV2PtRel[NLEP];
    float AK4JetV3PtRel[NLEP];
    float AK8JetV1PtRel[NLEP];
    float AK8JetV2PtRel[NLEP];
    float AK8JetV3PtRel[NLEP];
    float SubjetV1PtRel[NLEP];
    float SubjetV2PtRel[NLEP];
    float SubjetV3PtRel[NLEP];

    void init() {
      for (size_t i=0; i<NLEP; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        Iso04[i]=NOVAL_F;
        D0[i]=NOVAL_F;
        D0err[i]=NOVAL_F;
        Dxy[i]=NOVAL_F;
        Dxyerr[i]=NOVAL_F;
        Dz[i]=NOVAL_F;
        Dzerr[i]=NOVAL_F;
        IsLooseMuon[i]=NOVAL_F;
        IsSoftMuon[i]=NOVAL_F;
        IsTightMuon[i]=NOVAL_F;
        IsPFMuon[i]=NOVAL_F;
        IsGlobalMuon[i]=NOVAL_F;
        IsTrackerMuon[i]=NOVAL_F;
        GlbTrkNormChi2[i]=NOVAL_F;
        NumberValidMuonHits[i]=NOVAL_F;
        NumberMatchedStations[i]=NOVAL_F;
        NumberValidPixelHits[i]=NOVAL_F;
        NumberTrackerLayers[i]=NOVAL_F;
        NumberOfValidTrackerHits[i]=NOVAL_F;
        NumberOfPixelLayers[i]=NOVAL_F;
        InTrkNormChi2[i]=NOVAL_F;
        SumChargedHadronPt[i]=NOVAL_F;
        SumNeutralHadronPt[i]=NOVAL_F;
        SumPhotonPt[i]=NOVAL_F;
        SumPUPt[i]=NOVAL_F;
        GenMuonY[i]=NOVAL_F;
        GenMuonEta[i]=NOVAL_F;
        GenMuonPhi[i]=NOVAL_F;
        GenMuonPt[i]=NOVAL_F;
        GenMuonE[i]=NOVAL_F;
        GenMuonCharge[i]=NOVAL_F;
      }
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
    
  } mu;
  
  class JetVars {
  public:
    // Basic
    float Mass[NJET];
    float Pt[NJET];
    float Eta[NJET];
    float Y[NJET];
    float Phi[NJET];
    float E[NJET];
    float Charge[NJET];
    // B-TAGGING
    float CSV[NJET];
    float CSVV1[NJET];
    // GEN PARTON
    float GenPartonY[NJET];
    float GenPartonEta[NJET];
    float GenPartonPhi[NJET];
    float GenPartonPt[NJET];
    float GenPartonE[NJET];
    float GenPartonCharge[NJET];
    float PartonFlavour[NJET];
    float HadronFlavour[NJET];
    // GEN JET
    float GenJetY[NJET];
    float GenJetEta[NJET];
    float GenJetPhi[NJET];
    float GenJetPt[NJET];
    float GenJetE[NJET];
    float GenJetCharge[NJET];
    // CONSTITUENTS
    float muonMultiplicity[NJET];
    float PhotonEnergy[NJET];
    float ElectronEnergy[NJET];
    float MuonEnergy[NJET];
    float HFHadronEnergy[NJET];
    float HFEMEnergy[NJET];
    float ChargedHadronMultiplicity[NJET];
    float numberOfDaughters[NJET];
    float chargedMultiplicity[NJET];
    float neutralHadronMultiplicity[NJET];
    float neutralHadronEnergy[NJET];
    float neutralEmEnergy[NJET];
    float chargedEmEnergy[NJET];
    float chargedHadronEnergy[NJET];
    float photonMultiplicity[NJET];
    float electronMultiplicity[NJET];
    float HFHadronMultiplicity[NJET];
    float HFEMMultiplicity[NJET];
    float ChargeMuEnergy[NJET];
    float neutralMultiplicity[NJET];
    //FOR JEC
    float jecFactor0[NJET];
    float jetArea[NJET];
    // FOR SYSTEMATICS
    float SmearedPt[NJET];
    float SmearedPEta[NJET];
    float SmearedPhi[NJET];
    float SmearedE[NJET];
    float JERup[NJET];
    float JERdown[NJET];
    
    int HasNearGenTop[NJET];
    int NearGenTopIsHadronic[NJET];
    int NearGenWIsHadronic[NJET];
    int NearGenWToENu[NJET];
    int NearGenWToMuNu[NJET];
    int NearGenWToTauNu[NJET];
    int PassTopTag[NJET];
    float DRNearGenTop[NJET];
    float DRNearGenWFromTop[NJET];
    float DRNearGenBFromTop[NJET];
    float DRNearGenLepFromSLTop[NJET];
    float DRNearGenNuFromSLTop[NJET];
    float PtNearGenTop[NJET];
    float PtNearGenBFromTop[NJET];
    float PtNearGenWFromTop[NJET];
    float PtNearGenLepFromSLTop[NJET];
    float PtNearGenNuFromSLTop[NJET];
    
    void init() {
      for (size_t i=0; i<NJET; ++i) {
        Mass[i]=NOVAL_F;
        Pt[i]=NOVAL_F;
        Eta[i]=NOVAL_F;
        Y[i]=NOVAL_F;
        Phi[i]=NOVAL_F;
        E[i]=NOVAL_F;
        Charge[i]=NOVAL_F;
        CSV[i]=NOVAL_F;
        CSVV1[i]=NOVAL_F;
        GenPartonY[i]=NOVAL_F;
        GenPartonEta[i]=NOVAL_F;
        GenPartonPhi[i]=NOVAL_F;
        GenPartonPt[i]=NOVAL_F;
        GenPartonE[i]=NOVAL_F;
        GenPartonCharge[i]=NOVAL_F;
        PartonFlavour[i]=NOVAL_F;
        HadronFlavour[i]=NOVAL_F;
        GenJetY[i]=NOVAL_F;
        GenJetEta[i]=NOVAL_F;
        GenJetPhi[i]=NOVAL_F;
        GenJetPt[i]=NOVAL_F;
        GenJetE[i]=NOVAL_F;
        GenJetCharge[i]=NOVAL_F;
        muonMultiplicity[i]=NOVAL_F;
        PhotonEnergy[i]=NOVAL_F;
        ElectronEnergy[i]=NOVAL_F;
        MuonEnergy[i]=NOVAL_F;
        HFHadronEnergy[i]=NOVAL_F;
        HFEMEnergy[i]=NOVAL_F;
        ChargedHadronMultiplicity[i]=NOVAL_F;
        numberOfDaughters[i]=NOVAL_F;
        chargedMultiplicity[i]=NOVAL_F;
        neutralHadronMultiplicity[i]=NOVAL_F;
        neutralHadronEnergy[i]=NOVAL_F;
        neutralEmEnergy[i]=NOVAL_F;
        chargedEmEnergy[i]=NOVAL_F;
        chargedHadronEnergy[i]=NOVAL_F;
        photonMultiplicity[i]=NOVAL_F;
        electronMultiplicity[i]=NOVAL_F;
        HFHadronMultiplicity[i]=NOVAL_F;
        HFEMMultiplicity[i]=NOVAL_F;
        ChargeMuEnergy[i]=NOVAL_F;
        neutralMultiplicity[i]=NOVAL_F;
	jecFactor0[i]=NOVAL_F;
	jetArea[i]=NOVAL_F;
        SmearedPt[i]=NOVAL_F;
        SmearedPEta[i]=NOVAL_F;
        SmearedPhi[i]=NOVAL_F;
        SmearedE[i]=NOVAL_F;
        JERup[i]=NOVAL_F;
        JERdown[i]=NOVAL_F;
      }
    }
    
  } jet;
  
  class AK8Vars {
  public:
    float vSubjetIndex0[NJET];
    float vSubjetIndex1[NJET];
    float topSubjetIndex0[NJET];
    float topSubjetIndex1[NJET];
    float topSubjetIndex2[NJET];
    float topSubjetIndex3[NJET];
    float tau1[NJET];
    float tau2[NJET];
    float tau3[NJET];
    float softDropMass[NJET];
    float trimmedMass[NJET];
    float prunedMass[NJET];
    float filteredMass[NJET];
    float topMass[NJET];
    float wMass[NJET];
    float nSubJets[NJET];
    float minmass[NJET];
    
    void init() {
      for (size_t it=0; it<NJET; ++it) {
	vSubjetIndex0[it]=NOVAL_F;
	vSubjetIndex1[it]=NOVAL_F;
	topSubjetIndex0[it]=NOVAL_F;
	topSubjetIndex1[it]=NOVAL_F;
	topSubjetIndex2[it]=NOVAL_F;
	topSubjetIndex3[it]=NOVAL_F;
	tau1[it]=NOVAL_F;
	tau2[it]=NOVAL_F;
	tau3[it]=NOVAL_F;
	softDropMass[it]=NOVAL_F;
	trimmedMass[it]=NOVAL_F;
	prunedMass[it]=NOVAL_F;
	filteredMass[it]=NOVAL_F;
	topMass[it]=NOVAL_F;
	wMass[it]=NOVAL_F;
	nSubJets[it]=NOVAL_F;
	minmass[it]=NOVAL_F;
      }
    }
    
  } AK8;
  
  class AK4JetVars : public JetVars {
  public:
    AK4JetVars() { init(); }
    ~AK4JetVars() {}
    
    size_t it;
    int size;
    
    // Razor variables
    float MR;
    float MTR;
    float R;
    float R2;
    
    void init() {
      JetVars::init();
      
      it = -1;
      size = 0;
      
      MR = NOVAL_F;
      MTR = NOVAL_F;
      R = NOVAL_F;
      R2 = NOVAL_F;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } jetsAK4;
  
  class AK8JetVars : public JetVars, public AK8Vars {
  public:
    AK8JetVars() { init(); }
    ~AK8JetVars() {}
    
    size_t it;
    int size;
    
    // Razor variables
    float MR;
    float MTR;
    float R;
    float R2;
    
    void init() {
      JetVars::init();
      AK8Vars::init();
      
      it = -1;
      size = 0;
      
      MR = NOVAL_F;
      MTR = NOVAL_F;
      R = NOVAL_F;
      R2 = NOVAL_F;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } jetsAK8;
  
  class AK8SubJetVars : public JetVars {
  public:
    AK8SubJetVars() { init(); }
    ~AK8SubJetVars() {}
    
    float subjetCSV[NJET];
    
    size_t it;
    int size;
    
    void init() {
      JetVars::init();
      for (size_t it=0; it<NJET; ++it) subjetCSV[it]=NOVAL_F;
      
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	subjetCSV[NJET] = subjetCSV[it];
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } subjetsAK8;
  
  class CmsTopTagSubJetVars : public JetVars {
  public:
    CmsTopTagSubJetVars() { init(); }
    ~CmsTopTagSubJetVars() {}
    
    size_t it;
    int size;
    
    void init() {
      JetVars::init();
      
      it = -1;
      size = 0;
    }
    
    bool Loop() {
      ++it;
      if (it<(size_t)size) {
	return 1;
      } else {
	it=-1;
	return 0;
      }
    }
  } subjetsCmsTopTag;
  
  class MetData {
  public:
    MetData() { init(); };
    
    float Pt;
    float Phi;
    float Px;
    float Py;
    
    void init() {
      Pt=NOVAL_F;
      Phi=NOVAL_F;
      Px=NOVAL_F;
      Py=NOVAL_F;
    }
    
  } met;
  
  class EventData {
  public:
    EventData() { init(); };
    
    int NLep;
    int NTopHad;
    int NTopHadPreTag; // New Oct 16
    int NTopLep;
    int NTop;
    float HtLep;
    float HtTop;
    float Ht;
    float HtAll;
    float HtEx;
    float HtExFr;
    float HtTopFr;
    float TTHadDR;
    float TTHadDPhi;
    float TTHadSumPt; // New Oct 18
    float TTHadDEta;
    float TTHadMass;
    float TTHadPz;
    float TTHadHz;
    float TTHadDPz;
    float TTHadMR;
    float TTHadMTR;
    float TTHadR;
    float TTHadR2;
    float AK4_MR;
    float AK4_MTR;
    float AK4_R;
    float AK4_R2;
    float MR;
    float MTR;
    float R;
    float R2;
    float weight;
    int npv;
    int NGoodVtx; // New Oct 22
    
    // Other Variables
    float HT;
    float HTall;
    float HTtt;
    float HTlep;
    float HTex;
    float HTttFraction;
    float HTexFraction;
    
    // Development 04 March
    int nmu;
    int nele;
    int neletight;
    int nmuveto;
    int neleveto;
    float DRJetLep[NJET];
    float EleDRJet[NLEP];
    float MuDRJet[NLEP];
    float RelPtJetLep[NJET];
    float EleRelPtJet[NLEP];
    float MuRelPtJet[NLEP];
    float EleJetCombMass[NLEP];
    float MuJetCombMass[NLEP];

    
    // Development 20 April
    int JetGenTruth[NJET];
    bool JetHasMatchedGenTop[NJET];
    int JetMatchedGenTopType[NJET];
    bool JetMatchedGenTopIsMerged[NJET];
    float JetMatchedGenTopPt[NJET];
    float JetMatchedGenTopJetDR[NJET];
    float GenBJetDR[NJET];
    float GenWJetDR[NJET];
    float GenWGenBDR[NJET];
    float GenLepJetDR[NJET];
    float GenLepGenBDR[NJET];
    int NGenLepFromTop;
    bool IsGenTop[NGEN];
    int GenTopType[NGEN];
    bool GenTopHasMatchedJet[NGEN];
    bool GenTopHasMatchedTopTagJet[NGEN];
    bool JetIsHadTopTagged[NJET];
    float maxSubjetCSV[NJET];
    
    // Triggers
    // Hadronic
    bool HLT_AK8PFJet360_TrimMass30;
    bool HLT_PFJet450;
    bool HLT_PFJet500;
    bool HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;
    bool HLT_PFHT750_4Jet;
    bool HLT_PFHT750_4JetPt50;
    bool HLT_ECALHT800;
    bool HLT_PFHT800;
    bool HLT_PFHT900;
    // Hadronic - Prescaled Auxilary
    bool HLT_PFHT350;
    bool HLT_PFHT400;
    bool HLT_PFHT475;
    bool HLT_PFHT600;
    bool HLT_PFHT650;
    bool HLT_PFHT550_4Jet;
    bool HLT_PFHT650_4Jet;
    // Razor
    bool HLT_Rsq0p25;
    bool HLT_Rsq0p30;
    bool HLT_RsqMR240_Rsq0p09_MR200_4jet;
    bool HLT_RsqMR240_Rsq0p09_MR200;
    bool HLT_RsqMR270_Rsq0p09_MR200_4jet;
    bool HLT_RsqMR270_Rsq0p09_MR200;
    // Lepton + B-tag
    bool HLT_Mu10_CentralPFJet30_BTagCSV0p5PF;
    bool HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF;
    bool HLT_Ele15_IsoVVVL_BTagtop8CSV07_PFHT400;
    bool HLT_Ele15_IsoVVVL_PFHT600;
    bool HLT_Ele15_PFHT300;
    bool HLT_Mu15_IsoVVVL_BTagCSV07_PFHT400;
    bool HLT_Mu15_IsoVVVL_PFHT600;
    bool HLT_Mu15_PFHT300;
    // Lepton - Non-isolated
    bool HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
    bool HLT_Mu40_eta2p1_PFJet200_PFJet50;
    bool HLT_Mu45_eta2p1;
    bool HLT_Mu50;
    // Lepton - Isolated
    bool HLT_Ele32_eta2p1_WPLoose_Gsf;
    bool HLT_Ele32_eta2p1_WPTight_Gsf;
    bool HLT_IsoMu24_eta2p1;
    bool HLT_IsoMu27;
    bool HLT_IsoTkMu24_eta2p1;
    bool HLT_IsoTkMu27;
    
    // Event filters (these are automatically picked up)
    bool Flag_HBHEIsoNoiseFilterResult;
    bool Flag_HBHENoiseFilterResult;
    bool Flag_HBHENoiseFilterResultRun1;
    bool Flag_HBHENoiseFilterResultRun2Loose;
    bool Flag_HBHENoiseFilterResultRun2Tight;
    bool Flag_trackingFailureFilter;
    bool Flag_goodVertices;
    bool Flag_CSCTightHaloFilter;
    bool Flag_trkPOGFilters;
    bool Flag_trkPOG_logErrorTooManyClusters;
    bool Flag_EcalDeadCellTriggerPrimitiveFilter;
    bool Flag_ecalLaserCorrFilter;
    bool Flag_trkPOG_manystripclus53X;
    bool Flag_eeBadScFilter;
    bool Flag_METFilters;
    bool Flag_HBHENoiseFilter;
    bool Flag_trkPOG_toomanystripclus53X;
    bool Flag_hcalLaserEventFilter;

    // Vertices
    int vtx_size;
    float vtx_z[100];
    float vtx_rho[100];
    float vtx_chi[100];
    int vtx_ndof[100];
    
    unsigned int RunNumber;
    unsigned int LumiBlock;
    unsigned long long EventNumber;
    
    void init() {
      NLep=NOVAL_I;
      NTopHad=NOVAL_I;
      NTopLep=NOVAL_I;
      NTop=NOVAL_I;
      HtLep=NOVAL_F;
      HtTop=NOVAL_F;
      Ht=NOVAL_F;
      HtAll=NOVAL_F;
      HtEx=NOVAL_F;
      HtExFr=NOVAL_F;
      HtTopFr=NOVAL_F;
      TTHadDR=NOVAL_F;
      TTHadDPhi=NOVAL_F;
      TTHadSumPt=NOVAL_F;
      TTHadDEta=NOVAL_F;
      TTHadMass=NOVAL_F;
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
      
      // Other Variables
      HT=NOVAL_F;
      HTall=NOVAL_F;
      HTtt=NOVAL_F;
      HTex=NOVAL_F;
      HTttFraction=NOVAL_F;
      HTexFraction=NOVAL_F;

      // Development 04 March
      nmu=NOVAL_I;
      nele=NOVAL_I;
      neletight=NOVAL_I;
      nmuveto=NOVAL_I;
      neleveto=NOVAL_I;
      for (size_t i=0; i<NJET; ++i) {
	DRJetLep[i]=NOVAL_F;
	RelPtJetLep[i]=NOVAL_F;
      }
      for (size_t i=0; i<NLEP; ++i) {
	EleDRJet[i]=NOVAL_F;
	MuDRJet[i]=NOVAL_F;
	EleRelPtJet[i]=NOVAL_F;
	MuRelPtJet[i]=NOVAL_F;
	EleJetCombMass[i]=NOVAL_F;
	MuJetCombMass[i]=NOVAL_F;
      }
      
      // Development 20 April
      for (size_t i=0; i<NJET; ++i) {
	JetGenTruth[i]=NOVAL_I;
	JetHasMatchedGenTop[i]=0;
	JetMatchedGenTopType[i]=NOVAL_I;
	JetMatchedGenTopIsMerged[i]=0;
	JetMatchedGenTopPt[i]=NOVAL_F;
        JetMatchedGenTopJetDR[i]=NOVAL_F;
        GenBJetDR[i]=NOVAL_F;
        GenWJetDR[i]=NOVAL_F;
        GenWGenBDR[i]=NOVAL_F;
        GenLepJetDR[i]=NOVAL_F;
        GenLepGenBDR[i]=NOVAL_F;
	JetIsHadTopTagged[i]=0;
	maxSubjetCSV[i]=NOVAL_F;
      }
      NGenLepFromTop=NOVAL_I;
      for (size_t i=0; i<NGEN; ++i) {
	IsGenTop[i]=0;
	GenTopType[i]=NOVAL_I;
	GenTopHasMatchedJet[i]=0;
	GenTopHasMatchedTopTagJet[i]=0;
      }
      
      // Triggers
      // Hadronic
      HLT_AK8PFJet360_TrimMass30=NOVAL_I;
      HLT_PFJet450=NOVAL_I;
      HLT_PFJet500=NOVAL_I;
      HLT_AK8PFHT700_TrimR0p1PT0p03Mass50=NOVAL_I;
      HLT_PFHT750_4Jet=NOVAL_I;
      HLT_PFHT750_4JetPt50=NOVAL_I;
      HLT_ECALHT800=NOVAL_I;
      HLT_PFHT800=NOVAL_I;
      HLT_PFHT900=NOVAL_I;
      // Hadronic - Prescaled Auxilary
      HLT_PFHT350=NOVAL_I;
      HLT_PFHT400=NOVAL_I;
      HLT_PFHT475=NOVAL_I;
      HLT_PFHT600=NOVAL_I;
      HLT_PFHT650=NOVAL_I;
      HLT_PFHT550_4Jet=NOVAL_I;
      HLT_PFHT650_4Jet=NOVAL_I;
      // Razor
      HLT_Rsq0p25=NOVAL_I;
      HLT_Rsq0p30=NOVAL_I;
      HLT_RsqMR240_Rsq0p09_MR200_4jet=NOVAL_I;
      HLT_RsqMR240_Rsq0p09_MR200=NOVAL_I;
      HLT_RsqMR270_Rsq0p09_MR200_4jet=NOVAL_I;
      HLT_RsqMR270_Rsq0p09_MR200=NOVAL_I;
      // Lepton + B-tag
      HLT_Mu10_CentralPFJet30_BTagCSV0p5PF=NOVAL_I;
      HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF=NOVAL_I;
      HLT_Ele15_IsoVVVL_BTagtop8CSV07_PFHT400=NOVAL_I;
      HLT_Ele15_IsoVVVL_PFHT600=NOVAL_I;
      HLT_Ele15_PFHT300=NOVAL_I;
      HLT_Mu15_IsoVVVL_BTagCSV07_PFHT400=NOVAL_I;
      HLT_Mu15_IsoVVVL_PFHT600=NOVAL_I;
      HLT_Mu15_PFHT300=NOVAL_I;
      // Lepton - Non-isolated
      HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50=NOVAL_I;
      HLT_Mu40_eta2p1_PFJet200_PFJet50=NOVAL_I;
      HLT_Mu45_eta2p1=NOVAL_I;
      HLT_Mu50=NOVAL_I;
      // Lepton - Isolated
      HLT_Ele32_eta2p1_WPLoose_Gsf=NOVAL_I;
      HLT_Ele32_eta2p1_WPTight_Gsf=NOVAL_I;
      HLT_IsoMu24_eta2p1=NOVAL_I;
      HLT_IsoMu27=NOVAL_I;
      HLT_IsoTkMu24_eta2p1=NOVAL_I;
      HLT_IsoTkMu27=NOVAL_I;
    }
  } evt;
  
  void CalculateAllVariables() {
    //--------------------------------------------------------------------------
    //                                Vertices
    //--------------------------------------------------------------------------
    
    evt.NGoodVtx = 0;
    for (int iVtx=0; iVtx<evt.vtx_size; ++iVtx)
      if (evt.vtx_ndof[iVtx]>=4&&fabs(evt.vtx_z[iVtx])<24&&fabs(evt.vtx_rho[iVtx])<2) ++evt.NGoodVtx;
    
    //--------------------------------------------------------------------------
    //                           Generator Particles
    //--------------------------------------------------------------------------
    
    // Make a list of Generator level objects and save them to vectors
    std::vector<TLorentzVector> gen_top;
    std::vector<size_t > gen_top_index;
    std::vector<int> gen_top_ID;
    std::vector<TLorentzVector> gen_W_from_top;
    std::vector<TLorentzVector> gen_b_from_top;
    std::vector<TLorentzVector> gen_lep_from_W;
    std::vector<TLorentzVector> gen_neu_from_W;
    while(gen.Loop()) {
      evt.IsGenTop[gen.it]=0;
      evt.GenTopType[gen.it] = NOVAL_I;
      evt.GenTopHasMatchedJet[gen.it]=0;
      evt.GenTopHasMatchedTopTagJet[gen.it]=0;
      if (gen.Pt[gen.it]>0) {
	TLorentzVector genp; genp.SetPtEtaPhiE(gen.Pt[gen.it], gen.Eta[gen.it], gen.Phi[gen.it], gen.E[gen.it]);
        if (gen.ID[gen.it]!=gen.MomID[gen.it]) {
          if (abs(gen.ID[gen.it])==6) { 
	    evt.IsGenTop[gen.it]=1; 
	    gen_top.push_back(genp); 
	    gen_top_index.push_back(gen.it); 
	    gen_top_ID.push_back(gen.ID[gen.it]);
	    evt.GenTopType[gen.it] = 0;
	  }
          if (abs(gen.ID[gen.it])==5&&abs(gen.MomID[gen.it])==6) { gen_b_from_top.push_back(genp); }
          if (abs(gen.ID[gen.it])==24&&abs(gen.MomID[gen.it])==6) { gen_W_from_top.push_back(genp); }
          if ((abs(gen.ID[gen.it])==11||abs(gen.ID[gen.it])==13||abs(gen.ID[gen.it])==15)&&(abs(gen.MomID[gen.it])==24)) gen_lep_from_W.push_back(genp);
          if ((abs(gen.ID[gen.it])==12||abs(gen.ID[gen.it])==14||abs(gen.ID[gen.it])==16)&&(abs(gen.MomID[gen.it])==24)) gen_neu_from_W.push_back(genp);
        } else if (gen.ID[gen.it]==gen.MomID[gen.it]) {
	  // tops emit particles and we have to match consecutive tops to the original one
	  if (abs(gen.ID[gen.it])==6) {
	    size_t i=0, i_m_dR = -1, i_m_dE = -1;
	    double min_dE = 9999, min_dR = 9999;
	    while(i<gen_top.size()) {
	      if (gen_top_ID[i]==gen.MomID[gen.it]) {
		double dE = gen_top[i].E()-genp.E();
		double dR = gen_top[i].DeltaR(genp);
		if (fabs(dE)<fabs(min_dE)) {
		  min_dE = dE;
		  i_m_dE = i;
		}
		if (dR<min_dR) {
		  min_dR = dR;
		  i_m_dR = i;
		}
	      }
	      ++i;
	    }
	    //std::cout<<"match: dE: "<<gen_top[imatch].E()-genp.E()<<" Ep: "<<gen_top[imatch].E()<<" Ec: "<<genp.E()<<" dR: "<<genp.DeltaR(gen_top[imatch])<<std::endl;
	    //if (i_m_dE!=i_m_dR) printf("match: dE: %03.1f iE: %d   Ep: %04.1f Ec: %04.1f dR: %1.3f      min_dR: %1.3f Ec2: %4.1f iR: %d\n", double(min_dE), int(i_m_dE), double(gen_top[i_m_dE].E()), double(genp.E()), double(gen_top[i_m_dE].DeltaR(genp)), double(min_dR), double(gen_top[i_m_dR].E()), int(i_m_dR));
	    size_t imatch = (i_m_dE==i_m_dR) ? i_m_dE : ( (fabs(min_dE)/gen_top[i_m_dE].E()<0.1) ? i_m_dE : i_m_dR );
	    evt.IsGenTop[gen_top_index[imatch]]=0;
	    evt.IsGenTop[gen.it]=1;
	    gen_top[imatch]=genp;
	    gen_top_index[imatch]=gen.it;
	  }
	}
      }
    }
    
    // std::cout<<"Start looping on Gen\n";
    // //while(gen.Loop()) if (abs(gen.ID[gen.it])==1000021&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found ~g, it="<<gen.it<<" ID=~"<<(abs(gen.ID[gen.it])-1e6)<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // //while(gen.Loop()) if (abs(gen.ID[gen.it])==1000006&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found ~t, it="<<gen.it<<" ID=~"<<(abs(gen.ID[gen.it])-1e6)<<" MomID=~"<<(abs(gen.MomID[gen.it])-1e6)<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // //while(gen.Loop()) if (abs(gen.ID[gen.it])==1000024&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found ~chi+, it="<<gen.it<<" ID=~"<<(abs(gen.ID[gen.it])-1e6)<<" MomID=~"<<(abs(gen.MomID[gen.it])-1e6)<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // //while(gen.Loop()) if (abs(gen.ID[gen.it])==1000022&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found ~chi0, it="<<gen.it<<" ID=~"<<(abs(gen.ID[gen.it])-1e6)<<" MomID=~"<<(abs(gen.MomID[gen.it])-1e6)<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // while(gen.Loop()) if (abs(gen.ID[gen.it])==6&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found t, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // while(gen.Loop()) if (abs(gen.ID[gen.it])==5&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found b, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // while(gen.Loop()) if (abs(gen.ID[gen.it])==24&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found W, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // while(gen.Loop()) if ((abs(gen.ID[gen.it])==11||abs(gen.ID[gen.it])==13)&&(abs(gen.MomID[gen.it])==24||abs(gen.MomID[gen.it])>1e6)) std::cout<<"Found mu/e, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    // std::cout<<"\n";
    //if ((gen_top[0].Pt()>400&&gen_top[1].Pt()>400)&&gen_lep_from_W.size()>0) {
    //  std::cout<<"Start looping on Gen\n";
    //  while(gen.Loop()) if (abs(gen.ID[gen.it])==6&&gen.ID[gen.it]!=gen.MomID[gen.it]) std::cout<<"Found t, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //  while(gen.Loop()) if (abs(gen.ID[gen.it])==5&&abs(gen.MomID[gen.it])==6) std::cout<<"Found b, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //  while(gen.Loop()) if (abs(gen.ID[gen.it])==24&&abs(gen.MomID[gen.it])==6) std::cout<<"Found W, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //  while(gen.Loop()) if ((abs(gen.ID[gen.it])==11||abs(gen.ID[gen.it])==13)&&abs(gen.MomID[gen.it])==24) std::cout<<"Found mu/e from W mother, it="<<gen.it<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //  if (good_W_matches&&nlep_from_top>0) for (size_t i=0; i<gen_top_matched_W_matched_lep.size(); ++i) {
    //    TLorentzVector b = gen_top_matched_b[top_parent[i]];
    //    TLorentzVector W = gen_top_matched_W[top_parent[i]];
    //    TLorentzVector lep = gen_top_matched_W_matched_lep[i];
    //    double DR_lep_to_b = lep.DeltaR(b);
    //    std::cout<<"DR(lep, b)="<<DR_lep_to_b<<std::endl;
    //  }
    //  std::cout<<"\n";
    //}
    
    // Find bs and Ws
    // Method: bs and Ws with top parent are combined
    // Best pair with lowest combined mass and DR difference is selected
    std::vector<TLorentzVector> gen_top_matched_b;
    std::vector<TLorentzVector> gen_top_matched_W;
    std::vector<bool> W_is_leptonic;
    bool good_W_matches = true;
    for (size_t i=0; i<gen_top.size(); ++i) {
      // Match b and W to t
      size_t j_b = -1, k_W = -1;
      double min_DM = 9999, min_DR = 9999;
      if (gen_b_from_top.size()<gen_top.size()||gen_W_from_top.size()<gen_top.size()) {
	//std::cout<<"Not enough b/W with top parent"<<std::endl;
	good_W_matches = false;
      } else {
        for (size_t j=0; j<gen_b_from_top.size(); ++j) {
          for (size_t k=0; k<gen_W_from_top.size(); ++k) {
            TLorentzVector bW_comb = gen_b_from_top[j]+gen_W_from_top[k];
            double DR = gen_top[i].DeltaR(bW_comb);
            double DM = fabs(gen_top[i].M()-bW_comb.M());
            if (DR<0.8) {
              if (DM<min_DM) {
                min_DM = DM;
		min_DR = DR;
                j_b = j;
                k_W = k;
              }
            }
          }
        }
	//printf("W/b to top match: %.6f %.6f\n", min_DR, min_DM);
	if (min_DR<0.8&&min_DM<1) {
	  gen_top_matched_b.push_back(gen_b_from_top[j_b]);
	  gen_top_matched_W.push_back(gen_W_from_top[k_W]);
	} else {
	  good_W_matches = false;
	}
      }
    }
    //if (!good_W_matches) {
    //  while(gen.Loop()) std::cout<<"it="<<gen.it<<" ID=   "<<gen.ID[gen.it]<<" MomID=   "<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<" E="<<gen.E[gen.it]<<std::endl;
    //  std::cout<<""<<std::endl;
    //}
    //while(gen.Loop()) if (gen.ID[gen.it]==6) std::cout<<"Found t it="<<gen.it<<" ID="<<gen.ID[gen.it]<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //while(gen.Loop()) if (gen.ID[gen.it]==-6) std::cout<<"Found t it="<<gen.it<<" ID="<<gen.ID[gen.it]<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //while(gen.Loop()) if (gen.MomID[gen.it]==6) std::cout<<"Found child it="<<gen.it<<" ID="<<gen.ID[gen.it]<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    //while(gen.Loop()) if (gen.MomID[gen.it]==-6) std::cout<<"Found child it="<<gen.it<<" ID="<<gen.ID[gen.it]<<" MomID="<<gen.MomID[gen.it]<<" Pt="<<gen.Pt[gen.it]<<" Eta="<<gen.Eta[gen.it]<<" Phi="<<gen.Phi[gen.it]<<std::endl;
    
    // If we have lepton from W, find parent
    // Do as above with tops, but use neutrino and lepton instead to find W parent
    // In the end associate with top already found
    evt.NGenLepFromTop = 0;
    std::vector<TLorentzVector> gen_top_matched_W_matched_lep;
    std::vector<TLorentzVector> gen_top_matched_W_matched_neu;
    for (size_t i=0; i<gen_top_matched_W.size(); ++i) {
      TLorentzVector lep, neu;
      // Match lep and neutrino to W
      size_t j_l = -1, k_n = -1;
      double min_DM = 9999, min_DR = 9999;
      for (size_t j=0; j<gen_lep_from_W.size(); ++j) {
	for (size_t k=0; k<gen_neu_from_W.size(); ++k) {
	  TLorentzVector ln_comb = gen_lep_from_W[j]+gen_neu_from_W[k];
	  double DR = gen_top_matched_W[i].DeltaR(ln_comb);
	  double DM = fabs(gen_top_matched_W[i].M()-ln_comb.M());
	  if (DR<0.8) {
	    if (DM<min_DM) {
	      min_DM = DM;
	      min_DR = DR;
	      j_l = j;
	      k_n = k;
	    }
	  }
	}
      }
      bool lep_found = (min_DR<0.8&&min_DM<1);
      W_is_leptonic.push_back(lep_found);
      evt.GenTopType[gen_top_index[i]] = 1;
      //printf("l/v to W match: %.6f %.6f\n", min_DR, min_DM);
      if (lep_found) {
	lep = gen_lep_from_W[j_l];
	neu = gen_neu_from_W[k_n];
	++evt.NGenLepFromTop;
      }
      gen_top_matched_W_matched_lep.push_back(lep);
      gen_top_matched_W_matched_neu.push_back(neu);
    }
    
    // Match jets to tops (find the closest jet to top, sort by distance from gen)
    // Could also do genjet matching (but this is done in B2G)
    std::vector<TLorentzVector> temp = gen_top;
    std::vector<size_t > temp_it = gen_top_index;
    std::map<size_t, size_t > jet_gentop_it;
    const bool verbose = 0;
    while (temp.size()) {
      // find gentop  - jet pair with lowest DR (associate them and remove from temp top colelction)
      double min_DR = 9999, matched_DR = 9999;
      size_t top_m_it = -1, top_closest_it = -1, jet_m_it = -1;
      for (size_t i=0; i<temp.size(); ++i) {
	TLorentzVector top = temp[i];
	while(jetsAK8.Loop()) {
	  TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it], jetsAK8.Phi[jetsAK8.it], jetsAK8.E[jetsAK8.it]);
	  double DR = jet.DeltaR(top);
	  if (DR<min_DR) {
	    min_DR = DR;
	    top_closest_it = i;
	    if (!jet_gentop_it.count(jetsAK8.it)) {
	      matched_DR = DR;
	      top_m_it = i;
	      jet_m_it = jetsAK8.it;
	    }
	  }
	}
      }
      if (matched_DR<0.8) {
	if (verbose) std::cout<<"Top-jet match found, top(gen) it="<<temp_it[top_m_it]<<" jet it="<<jet_m_it<<" dR="<<matched_DR<<std::endl;
	jet_gentop_it[jet_m_it] = top_m_it;
	evt.GenTopHasMatchedJet[temp_it[top_m_it]] = 1;
#if NEW_TOP_DEF == 1
	evt.GenTopHasMatchedTopTagJet[temp_it[top_m_it]] = min_DR<0.8 && (jetsAK8.tau3[jet_m_it]/jetsAK8.tau2[jet_m_it])<TAU32_CUT_NEW && jetsAK8.Pt[jet_m_it]>PT_CUT_NEW && 
	  jetsAK8.softDropMass[jet_m_it] > SD_MASS_CUT_LOW && jetsAK8.softDropMass[jet_m_it] < SD_MASS_CUT_HIGH;
#else
	evt.GenTopHasMatchedTopTagJet[temp_it[top_m_it]] = min_DR<0.8 && jetsAK8.tau3[jet_m_it]/jetsAK8.tau2[jet_m_it]<TAU32_CUT_OLD && jetsAK8.Pt[jet_m_it]>PT_CUT_OLD && jetsAK8.prunedMass[jet_m_it]>PRUNED_MASS_CUT_LOW;
#endif
	temp.erase(temp.begin()+top_m_it);
	temp_it.erase(temp_it.begin()+top_m_it);
      } else if (jetsAK8.size) {
	if (verbose) {
	  std::cout<<"No match  found, possible pairs:"<<std::endl;
	  for (size_t i=0; i<temp.size(); ++i) {
	    TLorentzVector top = temp[i];
	    while(jetsAK8.Loop()) {
	      TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it], jetsAK8.Phi[jetsAK8.it], jetsAK8.E[jetsAK8.it]);
	      double DR = jet.DeltaR(top);
	      std::cout<<"  top(gen) it="<<temp_it[i]<<" jet it="<<jetsAK8.it<<" dR="<<DR<<(jet_gentop_it.count(jetsAK8.it)?" (Already found)":"")<<std::endl;
	    }
	  }
	}
	temp.erase(temp.begin()+top_closest_it);
	temp_it.erase(temp_it.begin()+top_closest_it);
      } else {
	if (verbose) std::cout<<"No jets in event!!!!\n";
	temp.clear();
	temp_it.clear();
      }
    }
    if (verbose) std::cout<<std::endl;
    
    //--------------------------------------------------------------------------
    //                               Leptons
    //--------------------------------------------------------------------------
    
    // find good leptons (for letponic tops)
    int ngoodleptons = 0;
    std::vector<TLorentzVector> goodleps;
    evt.HTlep = 0;
    evt.nmu = 0;
    evt.nele = 0;
    evt.neletight = 0;
    evt.nmuveto = 0;
    evt.neleveto = 0;
    while(ele.Loop()) {
      TLorentzVector el; el.SetPtEtaPhiE(ele.Pt[ele.it], ele.Eta[ele.it], ele.Phi[ele.it], ele.E[ele.it]);
      evt.EleDRJet[ele.it] = 9999;
      evt.EleRelPtJet[ele.it] = 9999;
      evt.EleJetCombMass[ele.it] = 9999;
      while(jetsAK8.Loop()) {
        TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it], jetsAK8.Phi[jetsAK8.it], jetsAK8.E[jetsAK8.it]);
        double dR = el.DeltaR(jet);
	if (dR<evt.EleDRJet[jetsAK8.it]) {
          evt.EleDRJet[ele.it] = dR;
          evt.EleRelPtJet[ele.it] = el.Perp(jet.Vect());
          evt.EleJetCombMass[ele.it] = (jet+el).M();
        }
      }
      if (ele.Pt[ele.it] > 35 && fabs(ele.Eta[ele.it]) < 2.5) {
	++ngoodleptons;
	++evt.nele;
	goodleps.push_back(el);
	evt.HTlep += ele.Pt[ele.it];
	if (ele.isVeto[ele.it]>0) ++evt.neleveto;
	if (ele.isTight[ele.it]>0) { // New 05 March
	  ++evt.neletight;
	}
      }
    }
    while(mu.Loop()) {
      TLorentzVector muon; muon.SetPtEtaPhiE(mu.Pt[mu.it], mu.Eta[mu.it], mu.Phi[mu.it], mu.E[mu.it]);
      evt.MuDRJet[mu.it] = 9999;
      evt.MuRelPtJet[mu.it] = 9999;
      evt.MuJetCombMass[mu.it] = 9999;
      while(jetsAK8.Loop()) {
        TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it], jetsAK8.Phi[jetsAK8.it], jetsAK8.E[jetsAK8.it]);
        double dR = muon.DeltaR(jet);
        if (dR<evt.MuDRJet[jetsAK8.it]) {
          evt.MuDRJet[mu.it] = dR;
          evt.MuRelPtJet[mu.it] = muon.Perp(jet.Vect());
          evt.MuJetCombMass[mu.it] = (jet+muon).M();
        }
      }
      if (mu.Pt[mu.it] > 45 && fabs(mu.Eta[mu.it]) < 2.1) {
	if (mu.IsTightMuon[mu.it]>0) {
	  ++ngoodleptons;
	  ++evt.nmu;
	  goodleps.push_back(muon);
	  evt.HTlep += mu.Pt[mu.it];
	}
	if (mu.IsSoftMuon[mu.it]>0) ++evt.nmuveto;
      }
    }
    
    //--------------------------------------------------------------------------
    //                                 Jets
    //--------------------------------------------------------------------------
    
    std::vector<TLorentzVector> hadtop;
    std::vector<TLorentzVector> hadtop_pre;
    std::vector<TLorentzVector> leptop;
    evt.HT = 0;
    evt.NTopHadPreTag = 0;
#if (NEW_TOP_DEF == 0) || (HIGHEST_PT_JETS_ONLY == 1)
    evt.NTopHad = 0;
#endif
    while(jetsAK8.Loop()) {
      bool is_top = false;
      TLorentzVector jet; jet.SetPtEtaPhiE(jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it], jetsAK8.Phi[jetsAK8.it], jetsAK8.E[jetsAK8.it]);
      // _______________________________________________________
      //                  Gen Particle Truth
      evt.JetGenTruth[jetsAK8.it] = gen_top.size()>0;
      evt.JetHasMatchedGenTop[jetsAK8.it] = 0;
      evt.JetMatchedGenTopType[jetsAK8.it] = NOVAL_I;
      evt.JetMatchedGenTopIsMerged[jetsAK8.it] = false;
      evt.JetMatchedGenTopPt[jetsAK8.it] = NOVAL_F;
      evt.JetMatchedGenTopJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenBJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenWJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenWGenBDR[jetsAK8.it] = NOVAL_F;
      evt.GenLepJetDR[jetsAK8.it] = NOVAL_F;
      evt.GenLepGenBDR[jetsAK8.it] = NOVAL_F;
      if (jet_gentop_it.count(jetsAK8.it)) {
	size_t top_it = jet_gentop_it[jetsAK8.it];
	evt.JetGenTruth[jetsAK8.it] = 2;
	evt.JetHasMatchedGenTop[jetsAK8.it] = 1;
	evt.JetMatchedGenTopType[jetsAK8.it] = 0;
	evt.JetMatchedGenTopPt[jetsAK8.it] = gen_top[top_it].Pt();
	evt.JetMatchedGenTopJetDR[jetsAK8.it] = gen_top[top_it].DeltaR(jet);
	// If W matching was successful, more information is available
	if (good_W_matches) {
	  evt.JetGenTruth[jetsAK8.it] = 2;
	  evt.JetMatchedGenTopType[jetsAK8.it] = W_is_leptonic[top_it];
	  // Both b and Whad/lepton within jet cone
	  if (jet.DeltaR(gen_top_matched_b[top_it])<0.7 && jet.DeltaR(W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it] : gen_top_matched_W[top_it])<0.7) {
	    evt.JetMatchedGenTopIsMerged[jetsAK8.it] = true;
	    evt.JetGenTruth[jetsAK8.it] = 3+W_is_leptonic[top_it];
	  }
	  evt.GenBJetDR[jetsAK8.it] = gen_top_matched_b[top_it].DeltaR(jet);
	  evt.GenWJetDR[jetsAK8.it] = gen_top_matched_W[top_it].DeltaR(jet);
	  evt.GenWGenBDR[jetsAK8.it] = gen_top_matched_W[top_it].DeltaR(gen_top_matched_b[top_it]);
	  evt.GenLepJetDR[jetsAK8.it] = W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it].DeltaR(jet) : NOVAL_F;
	  evt.GenLepGenBDR[jetsAK8.it] = W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it].DeltaR(gen_top_matched_b[top_it]) : NOVAL_F;
	} else {
	  evt.JetGenTruth[jetsAK8.it] = 5;
	}
      }
      
      // Subjets
      evt.maxSubjetCSV[jetsAK8.it] = 0;
      for (size_t i=0; i<jetsAK8.nSubJets[jetsAK8.it]; ++i) if (i<2) {
	size_t subjet_it = i==0 ? jetsAK8.vSubjetIndex0[jetsAK8.it] : jetsAK8.vSubjetIndex1[jetsAK8.it];
	float sjCSV = subjetsAK8.subjetCSV[subjet_it];
	if (sjCSV > evt.maxSubjetCSV[jetsAK8.it]) evt.maxSubjetCSV[jetsAK8.it] = sjCSV;
      }
      //printf("Subjets of jet it: %d  pt: %3.1f  eta: %1.2f  phi: %1.2f minSjCSV: %1.2f\n", int(jetsAK8.it), jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it], jetsAK8.Phi[jetsAK8.it], evt.maxSubjetCSV[jetsAK8.it]);
      //while(subjetsAK8.Loop()) {
      //  TLorentzVector subjet; subjet.SetPtEtaPhiE(subjetsAK8.Pt[jetsAK8.it], subjetsAK8.Eta[jetsAK8.it], subjetsAK8.Phi[jetsAK8.it], subjetsAK8.E[jetsAK8.it]);
      //  if (jet.DeltaR(subjet)<0.8) {
      //    printf("    subjet it: %d  pt: %3.1f  eta: %1.2f  phi: %1.2f CSV: %12f\n", int(subjetsAK8.it), subjetsAK8.Pt[jetsAK8.it], subjetsAK8.Eta[jetsAK8.it], subjetsAK8.Phi[jetsAK8.it], subjetsAK8.subjetCSV[jetsAK8.it]);
      //  }
      //}
      
      // _______________________________________________________
      //                  Hadronic Top Tag definition
      
      evt.JetIsHadTopTagged[jetsAK8.it] = 0;
#if HIGHEST_PT_JETS_ONLY == 1
      if (jetsAK8.it<2) { // New Oct
#endif
#if NEW_TOP_DEF == 0
	if ( (jetsAK8.tau2[jetsAK8.it]>0 && jetsAK8.tau3[jetsAK8.it]>0 ? jetsAK8.tau3[jetsAK8.it]/jetsAK8.tau2[jetsAK8.it] < TAU32_CUT_OLD : 0 ) &&
	     jetsAK8.Pt[jetsAK8.it] > PT_CUT_OLD &&
	     jetsAK8.prunedMass[jetsAK8.it] > PRUNED_MASS_CUT_LOW) { // Latest
	  evt.JetIsHadTopTagged[jetsAK8.it] = 1;
	  is_top = true;
	  ++evt.NTopHadPreTag;
	  ++evt.NTopHad;
	  hadtop.push_back(jet);
	} else if (jetsAK8.Pt[jetsAK8.it] > PT_CUT_OLD && jetsAK8.prunedMass[jetsAK8.it] > PRUNED_MASS_CUT_LOW) { // Top like jets for sideband fitting
	  ++evt.NTopHadPreTag;
	  hadtop_pre.push_back(jet);
	}
#else
	// New hadronic top tag
	if ((jetsAK8.tau2[jetsAK8.it]>0 && jetsAK8.tau3[jetsAK8.it]>0 ? jetsAK8.tau3[jetsAK8.it]/jetsAK8.tau2[jetsAK8.it] < TAU32_CUT_NEW : 0 ) &&
	    jetsAK8.Pt[jetsAK8.it] > PT_CUT_NEW &&
	    jetsAK8.softDropMass[jetsAK8.it] > SD_MASS_CUT_LOW && jetsAK8.softDropMass[jetsAK8.it] < SD_MASS_CUT_HIGH) {
	  evt.JetIsHadTopTagged[jetsAK8.it] = 1;
	  is_top = true;
	  ++evt.NTopHadPreTag;
#if HIGHEST_PT_JETS_ONLY == 1
	  ++evt.NTopHad;
#endif
	  hadtop.push_back(jet);
	  // Top like jets for sideband fitting
	} else if (jetsAK8.Pt[jetsAK8.it] > PT_CUT_NEW && jetsAK8.softDropMass[jetsAK8.it] > SD_MASS_CUT_LOW && jetsAK8.softDropMass[jetsAK8.it] < SD_MASS_CUT_HIGH) {
	  ++evt.NTopHadPreTag;
	  hadtop_pre.push_back(jet);
	}
#endif
#if HIGHEST_PT_JETS_ONLY == 1
      }
#endif
      
      // semi-leptonic tops
      TLorentzVector lep;
      evt.DRJetLep[jetsAK8.it] = 9999;
      evt.RelPtJetLep[jetsAK8.it] = 9999;
      for (size_t i=0; i<goodleps.size(); ++i) {
        if (goodleps[i].DeltaR(jet)< evt.DRJetLep[jetsAK8.it]) {
          evt.DRJetLep[jetsAK8.it] = goodleps[i].DeltaR(jet);
          evt.RelPtJetLep[jetsAK8.it] = goodleps[i].Perp(jet.Vect());
          lep = goodleps[i];
        }
      }
      if (evt.DRJetLep[jetsAK8.it]<1.0) {
        is_top = true;
        TLorentzVector lepjet = lep + jet;
        leptop.push_back(lepjet);
      }
      // Extra - all except above hadronic/leptonic tops
      evt.HT += jetsAK8.Pt[jetsAK8.it];
      if (is_top) evt.HTtt += jetsAK8.Pt[jetsAK8.it];
    } // End jet loop
    evt.HTall = evt.HT + met.Pt + evt.HTlep;
    //std::cout<<evt.HTall<<" "<<evt.HT<<" "<<met.Pt<<" "<<evt.HTlep<<std::endl;
    
    // _______________________________________________________
    //                  Top-pair variables
    
    TLorentzVector top1;
    TLorentzVector top2;
    if (evt.NTopHadPreTag >= 2) {
      if (evt.NTopHad == 2) {
        top1 = hadtop[0];
        top2 = hadtop[1];
      } else if (evt.NTopHad == 1) {
        if (hadtop[0].Pt() > hadtop_pre[0].Pt()) {
          top1 = hadtop[0];
          top2 = hadtop_pre[0];
        } else {
          top1 = hadtop_pre[0];
          top2 = hadtop[0];
        }
      } else if (evt.NTopHad == 0) {
        top1 = hadtop_pre[0];
        top2 = hadtop_pre[1];
      }
    }
    
    evt.HTtt=NOVAL_F;
    evt.HTex=NOVAL_F;
    evt.HTttFraction=NOVAL_F;
    evt.HTexFraction=NOVAL_F;
    if (evt.NTopHadPreTag>=2) {
      evt.HTtt = top1.Pt() + top2.Pt();
      evt.HTex = evt.HTall - evt.HTtt;
      evt.HTttFraction = evt.HTtt / evt.HTall;
      evt.HTexFraction = evt.HTex / evt.HTall;
    }
    evt.TTHadDR=NOVAL_F;
    evt.TTHadDPhi=NOVAL_F;
    evt.TTHadSumPt=NOVAL_F;
    evt.TTHadDEta=NOVAL_F;
    evt.TTHadMass=NOVAL_F;
    evt.TTHadPz=NOVAL_F;
    evt.TTHadHz=NOVAL_F;
    evt.TTHadDPz=NOVAL_F;
    evt.TTHadMR=NOVAL_F;
    evt.TTHadMTR=NOVAL_F;
    evt.TTHadR=NOVAL_F;
    evt.TTHadR2=NOVAL_F;
    if (evt.NTopHadPreTag>=2) {
      evt.TTHadDR = top1.DeltaR(top2);
      evt.TTHadDPhi = top1.DeltaPhi(top2);
      evt.TTHadSumPt = top1.Pt()+top2.Pt();
      evt.TTHadDEta = fabs(top1.Eta() - top2.Eta());
      evt.TTHadMass = (top1 + top2).M();
      evt.TTHadPz = (top1 + top2).Pz();
      evt.TTHadHz = top1.Pz() + top2.Pz();
      evt.TTHadDPz = fabs(top1.Pz() - top2.Pz());
      // Razor for hadronic top pair
      TVector3 metl;
      metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
      evt.TTHadMR = CalcMR_(top1, top2);
      evt.TTHadMTR = CalcMTR_(top1, top2, metl);
      evt.TTHadR = evt.TTHadMTR / evt.TTHadMR;
      evt.TTHadR2 = pow(evt.TTHadR, 2);
    }

    //___________________________________________________________________
    //                       RAZOR Variables
    //calcRazorAK4_();
    //calcRazorAK8_();
    //calcRazorCmsTopTag_();
#if R_TYPE == 2
    // OK in Ntuple
    evt.R=evt.AK4_R;
    evt.MR=evt.AK4_MR;
    evt.MTR=evt.AK4_MTR;
#elif R_TYPE == 3
    // TTHad quantities not OK for Pre-tag
    evt.R=evt.TTHadR;
    evt.MR=evt.TTHadMR;
    evt.MTR=evt.TTHadMTR;
#endif
    evt.R2 = evt.R * evt.R;
  }
    
  private:
  // Razor recipe taken from the RazorBoost gurus: N. Strobbe, S. Sekmen
  //   https://github.com/nstrobbe/RazorBoost/blob/master/analyzer/utils.h
  
  // Hemispheres:
  vector<TLorentzVector> CombineJets_(vector<TLorentzVector> myjets) {
    //std::cout<<"Start CombineJets with "<<myjets.size()<<" jets\n";
    vector<TLorentzVector> mynewjets;
    TLorentzVector j1, j2;
    //bool foundGood = false;
    int N_comb = 1;
    for(unsigned int i = 0; i < myjets.size(); i++){
      N_comb *= 2;
    }
    //std::cout<<"N_comb = "<<N_comb<<std::endl;
    double M_min = 9999999999.0;
    int j_count;
    for(int i = 1; i < N_comb-1; i++){
      TLorentzVector j_temp1, j_temp2;
      int itemp = i;
      j_count = N_comb/2;
      int count = 0;
      while(j_count > 0){
        if(itemp/j_count == 1){
          j_temp1 += myjets[count];
	  //std::cout<<"  1 <- "<<count<<" M2="<<myjets[count].M2()<<std::endl;
        } else {
          j_temp2 += myjets[count];
	  //std::cout<<"  2 <- "<<count<<" M2="<<myjets[count].M2()<<std::endl;
        }
        itemp -= j_count*(itemp/j_count);
        j_count /= 2;
        count++;
      }
      double M_temp = j_temp1.M2()+j_temp2.M2();
      //std::cout<<"  --> M_temp "<<j_temp1.M2()<<" + "<<j_temp2.M2()<<" = "<<M_temp<<std::endl;
      // smallest mass
      if(M_temp < M_min){
        M_min = M_temp;
        j1 = j_temp1;
        j2 = j_temp2;
      }
      //std::cout<<" M_min = "<<M_min<<std::endl;
    }
    if(j2.Pt() > j1.Pt()){
      TLorentzVector temp = j1;
      j1 = j2;
      j2 = temp;
    }
    //std::cout<<"Result: Jet1 pT = "<<j1.Pt()<<" Jet2 pT = "<<j2.Pt()<<"\n\n";
    mynewjets.push_back(j1);
    mynewjets.push_back(j2);
    return mynewjets;
  }
  
  // MR
  double CalcMR_(TLorentzVector ja, TLorentzVector jb){
    double A = ja.P();
    double B = jb.P();
    double az = ja.Pz();
    double bz = jb.Pz();
    TVector3 jaT, jbT;
    jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
    jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
    double ATBT = (jaT+jbT).Mag2();
    double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
      		 (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
    double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
    double mygamma = 1./sqrt(1.-mybeta*mybeta);
    //gamma times MRstar
    temp *= mygamma;
    return temp;
  }
  
  // MTR
  double CalcMTR_(TLorentzVector ja, TLorentzVector jb, TVector3 met){
    double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
    temp /= 2.;
    temp = sqrt(temp);
    return temp;
  }
  
  // MT
  double CalcMT_(TLorentzVector lepton, TLorentzVector pfmet){
    return sqrt( 2 * lepton.Pt() * pfmet.Pt() * ( 1 - cos( pfmet.Phi() - lepton.Phi() ) ) );
  }
  
  // Calculate Razor variables here
  void calcRazorAK4_() {

    // Select the best pair of jets (AK4, pt>40, |eta| < 3.0)
    std::vector<TLorentzVector> sjetl;
    while(jetsAK4.Loop()) {
      if (!(jetsAK4.Pt[jetsAK4.it] > 40) ) continue;
      if (!(fabs(jetsAK4.Eta[jetsAK4.it]) < 3) ) continue;
      TLorentzVector jl;
      jl.SetPtEtaPhiE(jetsAK4.Pt[jetsAK4.it], jetsAK4.Eta[jetsAK4.it],
      		jetsAK4.Phi[jetsAK4.it], jetsAK4.E[jetsAK4.it]);
      sjetl.push_back(jl);
    }
    std::vector<TLorentzVector> hemis = CombineJets_(sjetl);
    
    // ---------------------
    // -- Razor variables --
    // ---------------------
    
    TVector3 metl;
    metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
    
    if (hemis.size() == 2) {
      jetsAK4.MR = CalcMR_(hemis[0], hemis[1]);
      jetsAK4.MTR = CalcMTR_(hemis[0], hemis[1], metl);
      jetsAK4.R = jetsAK4.MTR / jetsAK4.MR;
      jetsAK4.R2 = pow(jetsAK4.R, 2);
    }
  }
  
  void calcRazorAK8_() {
    
    // Select the best pair of jets (AK8, pt>40, |eta| < 3.0)
    std::vector<TLorentzVector> sjetl;
    while(jetsAK8.Loop()) {
      if (!(jetsAK8.Pt[jetsAK8.it] > 40) ) continue;
      if (!(fabs(jetsAK8.Eta[jetsAK8.it]) < 3) ) continue;
      TLorentzVector jl;
      jl.SetPtEtaPhiE(jetsAK8.Pt[jetsAK8.it], jetsAK8.Eta[jetsAK8.it],
      		jetsAK8.Phi[jetsAK8.it], jetsAK8.E[jetsAK8.it]);
      sjetl.push_back(jl);
    }
    std::vector<TLorentzVector> hemis = CombineJets_(sjetl);
    
    // ---------------------
    // -- Razor variables --
    // ---------------------
    
    TVector3 metl;
    metl.SetPtEtaPhi(met.Pt, 0, met.Phi);
    
    if (hemis.size() == 2) {
      jetsAK8.MR = CalcMR_(hemis[0], hemis[1]);
      jetsAK8.MTR = CalcMTR_(hemis[0], hemis[1], metl);
      jetsAK8.R = jetsAK8.MTR / jetsAK8.MR;
      jetsAK8.R2 = pow(jetsAK8.R, 2);
    }
    //if (evt.NTopHad==2)
    //  std::cout<<"AK8: pt1="<<hemis[0].Pt()<<" pt2="<<hemis[1].Pt()<<" MET="<<met.Pt<<" MR="<<jetsAK8.MR<<" MTR="<<jetsAK8.MTR<<" R="<<jetsAK8.R<<std::endl;
  }
  
};

#endif
