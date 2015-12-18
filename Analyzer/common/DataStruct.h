#ifndef DataStruct_h
#define DataStruct_h

#define NOVAL_I -9999
#define NOVAL_F -9999.0

#include <cassert>
#include <map>
#include <vector>
#include <iostream>

#include "TLorentzVector.h"

#define NEW_TOP_DEF 1

/* 
   Top Tagging working points:
   https://twiki.cern.ch/twiki/bin/view/CMS/JetTopTagging#13_TeV_working_points
   e(B) = 3% SD WP2 	3% 	  	SD(beta=0,z=0.1, R=0.8) [110,210] 	tau32 < 0.75, max SD subjet b-discriminant > 0.39 

   --> Use WP2
   --> Use SoftDrop with PUPPI

*/
#define SD_MASS_CUT_LOW 110
#define SD_MASS_CUT_HIGH 210
#define TAU32_CUT_NEW 0.75
#define PT_CUT_NEW 400
#define R_CUT_NEW 0.4
#define R_CUT_LOW 0.2
#define DPHI_CUT_NEW 2.7

#define PRUNED_MASS_CUT_LOW 140
#define TAU32_CUT_OLD 0.75
#define PT_CUT_OLD 400
#define R_CUT_OLD 0.4
#define DPHI_CUT_OLD 2.8

inline void init_vec(std::vector<bool>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=0; }
inline void init_vec(std::vector<int>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=-9999; }
inline void init_vec(std::vector<float>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=-9999; }
inline void init_vec(std::vector<double>& vec) { vec.resize(500); for (int i=0; i<500; ++i) vec[i]=-9999; }

class DataStruct {
public:
  DataStruct() {};
  ~DataStruct() {};
  
  class GenVars {
  public:
    GenVars() { init(); };
    ~GenVars() {};
    
    // Basic
    std::vector<float> Mass;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Y;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    std::vector<float> _ID;
    std::vector<float> _Status;
    std::vector<float> _MomID;
    std::vector<int> ID;
    std::vector<int> Status;
    std::vector<int> MomID;
    
    size_t it;
    int size;
    
    void init() {
      it = -1;
      size = 0;
      init_vec(Mass);
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Y);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(_ID);
      init_vec(_Status);
      init_vec(_MomID);
      init_vec(ID);
      init_vec(Status);
      init_vec(MomID);
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
    
  } gen;
  
  class ElectronVars {
  public:
    ElectronVars() { init(); }
    ~ElectronVars() {}
    
    // Basic
    std::vector<float> Mass;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Y;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    // ElectronVars
    std::vector<float> Iso03;
    std::vector<float> D0;
    std::vector<float> Dz;
    std::vector<float> dEtaIn;
    std::vector<float> dPhiIn;
    std::vector<float> HoE;
    std::vector<float> full5x5siee;
    std::vector<float> ooEmooP;
    std::vector<float> missHits;
    std::vector<float> hasMatchedConVeto;
    std::vector<float> isEB;
    std::vector<float> isVeto;
    std::vector<float> isLoose;
    std::vector<float> isTight;
    std::vector<float> isMedium;
    std::vector<float> scEta;
    
    size_t it;
    int size;
    
    std::vector<float> DRNearGenEleFromSLTop;
    std::vector<float> PtNearGenEleFromSLTop;
    std::vector<float> PtNearGenTop;
    std::vector<int> IsPartOfNearAK4Jet;
    std::vector<int> IsPartOfNearAK8Jet;
    std::vector<int> IsPartOfNearSubjet;
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
    
    void init() {
      it = -1;
      size = 0;
      init_vec(Mass);
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Y);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(Iso03);
      init_vec(D0);
      init_vec(Dz);
      init_vec(dEtaIn);
      init_vec(dPhiIn);
      init_vec(HoE);
      init_vec(full5x5siee);
      init_vec(ooEmooP);
      init_vec(missHits);
      init_vec(hasMatchedConVeto);
      init_vec(isEB);
      init_vec(isVeto);
      init_vec(isLoose);
      init_vec(isTight);
      init_vec(isMedium);
      init_vec(scEta);
      init_vec(DRNearGenEleFromSLTop);
      init_vec(PtNearGenEleFromSLTop);
      init_vec(PtNearGenTop);
      init_vec(IsPartOfNearAK4Jet);
      init_vec(IsPartOfNearAK8Jet);
      init_vec(IsPartOfNearSubjet);
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
    std::vector<float> Mass;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Y;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    // MuonVars
    std::vector<float> Iso04;
    std::vector<float> D0;
    std::vector<float> D0err;
    std::vector<float> Dxy;
    std::vector<float> Dxyerr;
    std::vector<float> Dz;
    std::vector<float> Dzerr;
    std::vector<float> IsLooseMuon;
    std::vector<float> IsSoftMuon;
    std::vector<float> IsTightMuon;
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
    
    size_t it;
    int size;
  
    std::vector<float> DRNearGenMuFromSLTop;
    std::vector<float> PtNearGenMuFromSLTop;
    std::vector<float> PtNearGenTop;
    std::vector<int> IsPartOfNearAK4Jet;
    std::vector<int> IsPartOfNearAK8Jet;
    std::vector<int> IsPartOfNearSubjet;
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
  
    void init() {
      it = -1;
      size = 0;
      init_vec(Mass);
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Y);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(Iso04);
      init_vec(D0);
      init_vec(D0err);
      init_vec(Dxy);
      init_vec(Dxyerr);
      init_vec(Dz);
      init_vec(Dzerr);
      init_vec(IsLooseMuon);
      init_vec(IsSoftMuon);
      init_vec(IsTightMuon);
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
      init_vec(DRNearGenMuFromSLTop);
      init_vec(PtNearGenMuFromSLTop);
      init_vec(PtNearGenTop);
      init_vec(IsPartOfNearAK4Jet);
      init_vec(IsPartOfNearAK8Jet);
      init_vec(IsPartOfNearSubjet);
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
    std::vector<float> Mass;
    std::vector<float> Pt;
    std::vector<float> Eta;
    std::vector<float> Y;
    std::vector<float> Phi;
    std::vector<float> E;
    std::vector<float> Charge;
    // B-TAGGING
    std::vector<float> CSV;
    std::vector<float> CSVV1;
    // GEN PARTON
    std::vector<float> GenPartonY;
    std::vector<float> GenPartonEta;
    std::vector<float> GenPartonPhi;
    std::vector<float> GenPartonPt;
    std::vector<float> GenPartonE;
    std::vector<float> GenPartonCharge;
    std::vector<float> PartonFlavour;
    std::vector<float> HadronFlavour;
    // GEN JET
    std::vector<float> GenJetY;
    std::vector<float> GenJetEta;
    std::vector<float> GenJetPhi;
    std::vector<float> GenJetPt;
    std::vector<float> GenJetE;
    std::vector<float> GenJetCharge;
    // CONSTITUENTS
    std::vector<float> muonMultiplicity;
    std::vector<float> PhotonEnergy;
    std::vector<float> ElectronEnergy;
    std::vector<float> MuonEnergy;
    std::vector<float> HFHadronEnergy;
    std::vector<float> HFEMEnergy;
    std::vector<float> ChargedHadronMultiplicity;
    std::vector<float> numberOfDaughters;
    std::vector<float> chargedMultiplicity;
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
    //FOR JEC
    std::vector<float> jecFactor0;
    std::vector<float> jetArea;
    // FOR SYSTEMATICS
    std::vector<float> SmearedPt;
    std::vector<float> SmearedPEta;
    std::vector<float> SmearedPhi;
    std::vector<float> SmearedE;
    std::vector<float> JERup;
    std::vector<float> JERdown;
    
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
    
    void init() {
      init_vec(Mass);
      init_vec(Pt);
      init_vec(Eta);
      init_vec(Y);
      init_vec(Phi);
      init_vec(E);
      init_vec(Charge);
      init_vec(CSV);
      init_vec(CSVV1);
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
      init_vec(chargedMultiplicity);
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
      init_vec(jecFactor0);
      init_vec(jetArea);
      init_vec(SmearedPt);
      init_vec(SmearedPEta);
      init_vec(SmearedPhi);
      init_vec(SmearedE);
      init_vec(JERup);
      init_vec(JERdown);
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
    
  } jet;
  
  class AK8Vars {
  public:
    std::vector<float> vSubjetIndex0;
    std::vector<float> vSubjetIndex1;
    std::vector<float> topSubjetIndex0;
    std::vector<float> topSubjetIndex1;
    std::vector<float> topSubjetIndex2;
    std::vector<float> topSubjetIndex3;
    std::vector<float> tau1;
    std::vector<float> tau2;
    std::vector<float> tau3;
    std::vector<float> softDropMass;
    std::vector<float> trimmedMass;
    std::vector<float> prunedMass;
    std::vector<float> filteredMass;
    std::vector<float> topMass;
    std::vector<float> wMass;
    std::vector<float> nSubJets;
    std::vector<float> minmass;
    
    void init() {
      init_vec(vSubjetIndex0);
      init_vec(vSubjetIndex1);
      init_vec(topSubjetIndex0);
      init_vec(topSubjetIndex1);
      init_vec(topSubjetIndex2);
      init_vec(topSubjetIndex3);
      init_vec(tau1);
      init_vec(tau2);
      init_vec(tau3);
      init_vec(softDropMass);
      init_vec(trimmedMass);
      init_vec(prunedMass);
      init_vec(filteredMass);
      init_vec(topMass);
      init_vec(wMass);
      init_vec(nSubJets);
      init_vec(minmass);
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
    
    std::vector<float> subjetCSV;
    
    size_t it;
    int size;
    
    void init() {
      JetVars::init();
      
      it = -1;
      size = 0;
      init_vec(subjetCSV);
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
    
    std::vector<float> Pt;
    std::vector<float> Phi;
    std::vector<float> Px;
    std::vector<float> Py;
    
    void init() {
      init_vec(Pt);
      init_vec(Phi);
      init_vec(Px);
      init_vec(Py);
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

    // Generator/LHE variables
    float XSec;
    float NEvent_Corr;
    float Lumi_Weight;
    float Gen_Weight;
    std::vector<float> scale_Weights;
    std::vector<float> pdf_Weights;
    std::vector<float> alphas_Weights;
    int LHA_PDF_ID;
    float SUSY_Gluino_Mass;
    float SUSY_LSP_Mass;

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
    std::vector<float> DRJetLep;
    std::vector<float> EleDRJet;
    std::vector<float> MuDRJet;
    std::vector<float> RelPtJetLep;
    std::vector<float> EleRelPtJet;
    std::vector<float> MuRelPtJet;
    std::vector<float> EleJetCombMass;
    std::vector<float> MuJetCombMass;
  
    
    // Development 20 April
    std::vector<int> JetGenTruth;
    std::vector<bool>  JetHasMatchedGenTop;
    std::vector<int> JetMatchedGenTopType;
    std::vector<bool>  JetMatchedGenTopIsMerged;
    std::vector<float> JetMatchedGenTopPt;
    std::vector<float> JetMatchedGenTopJetDR;
    std::vector<float> GenBJetDR;
    std::vector<float> GenWJetDR;
    std::vector<float> GenWGenBDR;
    std::vector<float> GenLepJetDR;
    std::vector<float> GenLepGenBDR;
    int NGenLepFromTop;
    std::vector<bool>  IsGenTop;
    std::vector<int> GenTopType;
    std::vector<bool>  GenTopHasMatchedJet;
    std::vector<bool>  GenTopHasMatchedTopTagJet;
    std::vector<bool>  JetIsHadTopTagged;
    std::vector<float> maxSubjetCSV;
    
    // Triggers
    // Hadronic
    int HLT_AK8PFJet360_TrimMass30;
    int HLT_PFJet450;
    int HLT_PFJet500;
    int HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;
    int HLT_PFHT750_4Jet;
    int HLT_PFHT750_4JetPt50;
    int HLT_ECALHT800;
    int HLT_PFHT800;
    int HLT_PFHT900;
    // Hadronic - Prescaled Auxilary
    int HLT_PFHT350;
    int HLT_PFHT400;
    int HLT_PFHT475;
    int HLT_PFHT600;
    int HLT_PFHT650;
    int HLT_PFHT550_4Jet;
    int HLT_PFHT650_4Jet;
    // Razor
    int HLT_Rsq0p25;
    int HLT_Rsq0p30;
    int HLT_RsqMR240_Rsq0p09_MR200_4jet;
    int HLT_RsqMR240_Rsq0p09_MR200;
    int HLT_RsqMR270_Rsq0p09_MR200_4jet;
    int HLT_RsqMR270_Rsq0p09_MR200;
    // Lepton + B-tag
    int HLT_Mu10_CentralPFJet30_BTagCSV0p5PF;
    int HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF;
    int HLT_Ele15_IsoVVVL_BTagtop8CSV07_PFHT400;
    int HLT_Ele15_IsoVVVL_PFHT600;
    int HLT_Ele15_PFHT300;
    int HLT_Mu15_IsoVVVL_BTagCSV07_PFHT400;
    int HLT_Mu15_IsoVVVL_PFHT600;
    int HLT_Mu15_PFHT300;
    // Lepton - Non-isolated
    int HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
    int HLT_Mu40_eta2p1_PFJet200_PFJet50;
    int HLT_Mu45_eta2p1;
    int HLT_Mu50;
    // Lepton - Isolated
    int HLT_Ele32_eta2p1_WPLoose_Gsf;
    int HLT_Ele32_eta2p1_WPTight_Gsf;
    int HLT_IsoMu24_eta2p1;
    int HLT_IsoMu27;
    int HLT_IsoTkMu24_eta2p1;
    int HLT_IsoTkMu27;
    
    // Event filters (these are automatically picked up)
    int Flag_HBHEIsoNoiseFilterResult;
    int Flag_HBHENoiseFilterResult;
    int Flag_HBHENoiseFilterResultRun1;
    int Flag_HBHENoiseFilterResultRun2Loose;
    int Flag_HBHENoiseFilterResultRun2Tight;
    int Flag_trackingFailureFilter;
    int Flag_goodVertices;
    int Flag_CSCTightHaloFilter;
    int Flag_trkPOGFilters;
    int Flag_trkPOG_logErrorTooManyClusters;
    int Flag_EcalDeadCellTriggerPrimitiveFilter;
    int Flag_ecalLaserCorrFilter;
    int Flag_trkPOG_manystripclus53X;
    int Flag_eeBadScFilter;
    int Flag_METFilters;
    int Flag_HBHENoiseFilter;
    int Flag_trkPOG_toomanystripclus53X;
    int Flag_hcalLaserEventFilter;
  
    // Vertices
    int vtx_size;
    std::vector<float> vtx_z;
    std::vector<float> vtx_rho;
    std::vector<float> vtx_chi;
    std::vector<int> vtx_ndof;
    
    // Pileup
    int pu_NtrueInt;
    std::vector<int> pu_NInt;
    std::vector<int> pu_BX;
    
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

      // Generator/LHE variables
      XSec=NOVAL_F;
      NEvent_Corr=NOVAL_F;
      Lumi_Weight=NOVAL_F;
      Gen_Weight=NOVAL_F;
      init_vec(scale_Weights);
      init_vec(pdf_Weights);
      init_vec(alphas_Weights);
      LHA_PDF_ID=NOVAL_I;
      SUSY_Gluino_Mass=NOVAL_F;
      SUSY_LSP_Mass=NOVAL_F;
      
      npv=NOVAL_I;
      NGoodVtx=NOVAL_I;
      
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
      
      // Development 20 April
      NGenLepFromTop=NOVAL_I;
      
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
      
      init_vec(DRJetLep);
      init_vec(EleDRJet);
      init_vec(MuDRJet);
      init_vec(RelPtJetLep);
      init_vec(EleRelPtJet);
      init_vec(MuRelPtJet);
      init_vec(EleJetCombMass);
      init_vec(MuJetCombMass);
      init_vec(JetGenTruth);
      init_vec(JetHasMatchedGenTop);
      init_vec(JetMatchedGenTopType);
      init_vec(JetMatchedGenTopIsMerged);
      init_vec(JetMatchedGenTopPt);
      init_vec(JetMatchedGenTopJetDR);
      init_vec(GenBJetDR);
      init_vec(GenWJetDR);
      init_vec(GenWGenBDR);
      init_vec(GenLepJetDR);
      init_vec(GenLepGenBDR);
      init_vec(IsGenTop);
      init_vec(GenTopType);
      init_vec(GenTopHasMatchedJet);
      init_vec(GenTopHasMatchedTopTagJet);
      init_vec(JetIsHadTopTagged);
      init_vec(maxSubjetCSV);
      init_vec(vtx_z);
      init_vec(vtx_rho);
      init_vec(vtx_chi);
      init_vec(vtx_ndof);
      
      pu_NtrueInt = NOVAL_I;
      init_vec(pu_NInt);
      init_vec(pu_BX);
    }
  } evt;
};

#endif
