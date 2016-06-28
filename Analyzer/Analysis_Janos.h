#include "TLorentzVector.h"
#include "common/AnalysisBase.h"
#include "common/SmartHistos.h"

SmartHistos sh;

//_______________________________________________________
//                       Constructor
Analysis::Analysis() : AnalysisBase() { }


//_______________________________________________________
//                       Destructor
Analysis::~Analysis() { }

//_______________________________________________________
//                  Calculate variables

// gen particles
std::vector<bool> is_Gen_Top;
std::vector<int>  gen_Top_Type;
std::vector<bool> has_Matched_Jet;
std::vector<bool> has_Matched_Tagged_Jet;

// AK8 Puppi jets
std::vector<int>   top_Children_Within_Cone;
std::vector<int>  jetsAK8Puppi_HasNearGenTop;
std::vector<int>  jetsAK8Puppi_NearGenTopIsHadronic;
std::vector<float> jetsAK8Puppi_PtNearGenTop;
std::vector<float> jetsAK8Puppi_DRNearGenTop;
std::vector<float> jetsAK8Puppi_DRNearGenBFromTop;
std::vector<float> jetsAK8Puppi_DRNearGenWFromTop;
std::vector<float> jetsAK8Puppi_DRNearGenLepFromSLTop;

// Event
double dPhi;
double dEta;
double dR;
int passTau32Cuts;

void
Analysis::calculate_variables(DataStruct& data, const unsigned int& syst_index)
{
  // It only makes sense to calculate certain variables only once if they don't depend on jet energy
  if (syst_index == 0) {

    is_Gen_Top.clear();
    gen_Top_Type.clear();
    has_Matched_Jet.clear();
    has_Matched_Tagged_Jet.clear();

    //      GEN particles

    // Make a list of Generator level objects and save them to vectors
    std::vector<TLorentzVector> gen_top;
    std::vector<size_t > gen_top_index;
    std::vector<int> gen_top_ID;
    std::vector<TLorentzVector> gen_W_from_top;
    std::vector<TLorentzVector> gen_b_from_top;
    std::vector<TLorentzVector> gen_lep_from_W;
    std::vector<TLorentzVector> gen_neu_from_W;
    std::vector<TLorentzVector> gen_lep_from_top;
    while(data.gen.Loop()) {
      is_Gen_Top.push_back(0);
      gen_Top_Type.push_back(NOVAL_I);
      has_Matched_Jet.push_back(0);
      has_Matched_Tagged_Jet.push_back(0);
      if (data.gen.Pt[data.gen.it]>0) {
        TLorentzVector genp; genp.SetPtEtaPhiE(data.gen.Pt[data.gen.it], data.gen.Eta[data.gen.it], data.gen.Phi[data.gen.it], data.gen.E[data.gen.it]);
        if (data.gen.ID[data.gen.it]!=data.gen.Mom0ID[data.gen.it]) {
          if (abs(data.gen.ID[data.gen.it])==6) { 
            is_Gen_Top[data.gen.it]=1; 
            gen_top.push_back(genp); 
            gen_top_index.push_back(data.gen.it); 
            gen_top_ID.push_back(data.gen.ID[data.gen.it]);
            gen_Top_Type[data.gen.it] = 0;
          }
          if (abs(data.gen.ID[data.gen.it])==5&&abs(data.gen.Mom0ID[data.gen.it])==6) { gen_b_from_top.push_back(genp); }
          if (abs(data.gen.ID[data.gen.it])==24&&abs(data.gen.Mom0ID[data.gen.it])==6) { gen_W_from_top.push_back(genp); }
          if ((abs(data.gen.ID[data.gen.it])==11||abs(data.gen.ID[data.gen.it])==13||abs(data.gen.ID[data.gen.it])==15)&&(abs(data.gen.Mom0ID[data.gen.it])==24)) gen_lep_from_W.push_back(genp);
          if ((abs(data.gen.ID[data.gen.it])==12||abs(data.gen.ID[data.gen.it])==14||abs(data.gen.ID[data.gen.it])==16)&&(abs(data.gen.Mom0ID[data.gen.it])==24)) gen_neu_from_W.push_back(genp);
        } else if (data.gen.ID[data.gen.it]==data.gen.Mom0ID[data.gen.it]) {
          // tops emit particles and we have to match consecutive tops to the original one
          if (abs(data.gen.ID[data.gen.it])==6) {
            unsigned int i=0, i_m_dR = -1, i_m_dE = -1;
            double min_dE = 9999, min_dR = 9999;
            while(i<gen_top.size()) {
              if (gen_top_ID[i]==data.gen.Mom0ID[data.gen.it]) {
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
            unsigned int imatch = (i_m_dE==i_m_dR) ? i_m_dE : ( (fabs(min_dE)/gen_top[i_m_dE].E()<0.1) ? i_m_dE : i_m_dR );
            is_Gen_Top[gen_top_index[imatch]]=0;
            is_Gen_Top[data.gen.it]=1;
            gen_top[imatch]=genp;
            gen_top_index[imatch]=data.gen.it;
          }
        }
      }
    }

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

    // If we have lepton from W, find parent
    // Do as above with tops, but use neutrino and lepton instead to find W parent
    // In the end associate with top already found
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
      gen_Top_Type[gen_top_index[i]] = 1;
      //printf("l/v to W match: %.6f %.6f\n", min_DR, min_DM);
      if (lep_found) {
        lep = gen_lep_from_W[j_l];
        neu = gen_neu_from_W[k_n];
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
        while(data.jetsAK8Puppi.Loop()) {
          TLorentzVector jet; jet.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it], data.jetsAK8Puppi.Eta[data.jetsAK8Puppi.it], data.jetsAK8Puppi.Phi[data.jetsAK8Puppi.it], data.jetsAK8Puppi.E[data.jetsAK8Puppi.it]);
          double DR = jet.DeltaR(top);
          if (DR<min_DR) {
            min_DR = DR;
            top_closest_it = i;
            if (!jet_gentop_it.count(data.jetsAK8Puppi.it)) {
              matched_DR = DR;
              top_m_it = i;
              jet_m_it = data.jetsAK8Puppi.it;
            }
          }
        }
      }
      if (matched_DR<0.8) {
        if (verbose) std::cout<<"Top-jet match found, top(gen) it="<<temp_it[top_m_it]<<" jet it="<<jet_m_it<<" dR="<<matched_DR<<std::endl;
        jet_gentop_it[jet_m_it] = top_m_it;
        has_Matched_Jet[temp_it[top_m_it]] = 1;
        has_Matched_Tagged_Jet[temp_it[top_m_it]] = min_DR<0.8 && passHadTopTag[jet_m_it];
        temp.erase(temp.begin()+top_m_it);
        temp_it.erase(temp_it.begin()+top_m_it);
      } else if (data.jetsAK8Puppi.size) {
        if (verbose) {
          std::cout<<"No match  found, possible pairs:"<<std::endl;
          for (size_t i=0; i<temp.size(); ++i) {
            TLorentzVector top = temp[i];
            while(data.jetsAK8Puppi.Loop()) {
              TLorentzVector jet; jet.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it], data.jetsAK8Puppi.Eta[data.jetsAK8Puppi.it], data.jetsAK8Puppi.Phi[data.jetsAK8Puppi.it], data.jetsAK8Puppi.E[data.jetsAK8Puppi.it]);
              double DR = jet.DeltaR(top);
              std::cout<<"  top(gen) it="<<temp_it[i]<<" jet it="<<data.jetsAK8Puppi.it<<" dR="<<DR<<(jet_gentop_it.count(data.jetsAK8Puppi.it)?" (Already found)":"")<<std::endl;
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

    //      JETs

    top_Children_Within_Cone.clear();
    jetsAK8Puppi_HasNearGenTop.clear();
    jetsAK8Puppi_NearGenTopIsHadronic.clear();
    jetsAK8Puppi_PtNearGenTop.clear();
    jetsAK8Puppi_DRNearGenTop.clear();
    jetsAK8Puppi_DRNearGenBFromTop.clear();
    jetsAK8Puppi_DRNearGenWFromTop.clear();
    jetsAK8Puppi_DRNearGenLepFromSLTop.clear();

    // Loop on jets
    while(data.jetsAK8Puppi.Loop()) {

      TLorentzVector jet; jet.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it], data.jetsAK8Puppi.Eta[data.jetsAK8Puppi.it], data.jetsAK8Puppi.Phi[data.jetsAK8Puppi.it], data.jetsAK8Puppi.E[data.jetsAK8Puppi.it]);
      // _______________________________________________________
      //                  Gen Particle Truth
      top_Children_Within_Cone.push_back(0);
      if (jet_gentop_it.count(data.jetsAK8Puppi.it)) {
        size_t top_it = jet_gentop_it[data.jetsAK8Puppi.it];
        // If W matching was successful, more information is available
        if (good_W_matches) {
          // Both b and Whad/lepton within jet cone
          if (jet.DeltaR(gen_top_matched_b[top_it])<0.8 && jet.DeltaR(W_is_leptonic[top_it] ? gen_top_matched_W_matched_lep[top_it] : gen_top_matched_W[top_it])<0.8) {
            top_Children_Within_Cone[data.jetsAK8Puppi.it] = 1;
          }
        }
      }

      // Match AK8 Puppi jets to CHS jets and get their gen particle properties
      double match_dR = 9999;
      unsigned int match_it = 9999;
      while(data.jetsAK8.Loop()) {
        TLorentzVector jetCHS; jetCHS.SetPtEtaPhiE(data.jetsAK8.Pt[data.jetsAK8.it], data.jetsAK8.Eta[data.jetsAK8.it], data.jetsAK8.Phi[data.jetsAK8.it], data.jetsAK8.E[data.jetsAK8.it]);
        double dR = jet.DeltaR(jetCHS);
        if (dR < match_dR) {
          match_dR = dR;
          match_it = data.jetsAK8.it;
        }
      }
      if (match_dR<0.1) {
        //std::cout<<"Puppi jet["<<data.jetsAK8Puppi.it<<"] match to CHS jet["<<match_it<<"] with dR="<<match_dR<<std::endl;
        jetsAK8Puppi_HasNearGenTop.push_back(data.jetsAK8.HasNearGenTop[match_it]);
        jetsAK8Puppi_NearGenTopIsHadronic.push_back(data.jetsAK8.NearGenTopIsHadronic[match_it]);
        jetsAK8Puppi_PtNearGenTop.push_back(data.jetsAK8.PtNearGenTop[match_it]);
        jetsAK8Puppi_DRNearGenTop.push_back(data.jetsAK8.DRNearGenTop[match_it]);
        jetsAK8Puppi_DRNearGenBFromTop.push_back(data.jetsAK8.DRNearGenBFromTop[match_it]);
        jetsAK8Puppi_DRNearGenWFromTop.push_back(data.jetsAK8.DRNearGenWFromTop[match_it]);
        jetsAK8Puppi_DRNearGenLepFromSLTop.push_back(data.jetsAK8.DRNearGenLepFromSLTop[match_it]);
      } else {
        //if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=400)
        //  std::cout<<"Puppi jet["<<data.jetsAK8Puppi.it<<"] had no matching CHS jet, min dR="<<match_dR<<std::endl;
        jetsAK8Puppi_HasNearGenTop.push_back(NOVAL_I);
        jetsAK8Puppi_NearGenTopIsHadronic.push_back(NOVAL_I);
        jetsAK8Puppi_PtNearGenTop.push_back(NOVAL_F);
        jetsAK8Puppi_DRNearGenTop.push_back(NOVAL_F);
        jetsAK8Puppi_DRNearGenBFromTop.push_back(NOVAL_F);
        jetsAK8Puppi_DRNearGenWFromTop.push_back(NOVAL_F);
        jetsAK8Puppi_DRNearGenLepFromSLTop.push_back(NOVAL_F);
      }

    } // end of AK8 jet loop

    // Delta-phi
    dPhi = -9999;
    dEta = -9999;
    dR = -9999;
    passTau32Cuts = -9999;
    if (data.jetsAK8Puppi.size>=2) {
      TLorentzVector jet1, jet2;
      jet1.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[0], data.jetsAK8Puppi.Eta[0], data.jetsAK8Puppi.Phi[0], data.jetsAK8Puppi.E[0]);
      jet2.SetPtEtaPhiE(data.jetsAK8Puppi.Pt[1], data.jetsAK8Puppi.Eta[1], data.jetsAK8Puppi.Phi[1], data.jetsAK8Puppi.E[1]);
      dPhi = fabs(jet1.DeltaPhi(jet2));
      dEta = fabs(jet1.Eta()-jet2.Eta());
      dR   = fabs(jet1.DeltaR(jet2));
      double jet1_tau32 = data.jetsAK8Puppi.tau3[0];
      if (data.jetsAK8Puppi.tau2[0]!=0) jet1_tau32 /= data.jetsAK8Puppi.tau2[0];
      else jet1_tau32 = 9999;
      double jet2_tau32 = data.jetsAK8Puppi.tau3[1];
      if (data.jetsAK8Puppi.tau2[1]!=0) jet2_tau32 /= data.jetsAK8Puppi.tau2[1];
      else jet2_tau32 = 9999;
      passTau32Cuts = (jet1_tau32 < TOP_TAU32_CUT) && (jet2_tau32 < TOP_TAU32_CUT);
    }
  }
}

//_______________________________________________________
//          Define Analysis specific weights
double
Analysis::get_analysis_weight(DataStruct& data)
{
  double w = 1;

  return w;
}

//_______________________________________________________
//          Define Analysis event selection cuts

bool
Analysis::pass_skimming(DataStruct& data)
{
  float pt_threshold = 300;
  int N_CHS = 0, N_Puppi = 0;
  while(data.jetsAK8Puppi.Loop()) if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=pt_threshold) ++N_CHS;
  while(data.jetsAK8Puppi.Loop()) if (data.jetsAK8Puppi.Pt[data.jetsAK8Puppi.it]>=pt_threshold) ++N_Puppi;
  return (N_CHS >= 1 || N_Puppi >= 1);
}

void
Analysis::define_selections(const DataStruct& data)
{
  // cut1: njet >= 2
  analysis_cuts.push_back({ .name="2jet",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0;
			      return 1;
			    } });

  // cut2: jet 1 pass loose jet id
  analysis_cuts.push_back({ .name="jet1_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      return passLooseJetID[0];
			    } });

  // cut3: jet 2 pass loose jet id
  analysis_cuts.push_back({ .name="jet2_id",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      return passLooseJetID[1];
			    } });

  // cut4: jet 1 eta < 2.4
  analysis_cuts.push_back({ .name="jet1_eta",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (fabs(data.jetsAK8Puppi.Eta[0])>=2.4) return 0;
			      return 1;
			    } });

  // cut5: jet 2 eta < 2.4
  analysis_cuts.push_back({ .name="jet2_eta",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (fabs(data.jetsAK8Puppi.Eta[1])>=2.4) return 0;
			      return 1;
			    } });

  // cut6: jet 1 pt >= 400
  analysis_cuts.push_back({ .name="jet1_pt",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (data.jetsAK8Puppi.Pt[0]<TOP_PT_CUT) return 0;
			      return 1;
			    } });

  // cut7: jet 2 pt >= 400
  analysis_cuts.push_back({ .name="jet2_pt",   .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (data.jetsAK8Puppi.Pt[1]<TOP_PT_CUT) return 0;
			      return 1;
			    } });

  // cut8: 110 <= jet 1 mass (softdrop) < 210
  analysis_cuts.push_back({ .name="jet1_mass", .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      if (data.jetsAK8Puppi.softDropMass[0]<TOP_SD_MASS_CUT_LOW) return 0;
			      if (data.jetsAK8Puppi.softDropMass[0]>=TOP_SD_MASS_CUT_HIGH) return 0;
			      return 1;
			    } });

  // cut9: 110 <= jet 2 mass (softdrop) < 210
  analysis_cuts.push_back({ .name="jet2_mass", .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (data.jetsAK8Puppi.softDropMass[1]<TOP_SD_MASS_CUT_LOW) return 0;
			      if (data.jetsAK8Puppi.softDropMass[1]>=TOP_SD_MASS_CUT_HIGH) return 0;
			      return 1;
			    } });

  // cut10: Full-hadronic trigger
  analysis_cuts.push_back({ .name="hlt_ak8ht700_mass50", .func = [&data](){
			      // Define cut function here:
			      return data.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50==1; 
			    } });

  // cut11: | DeltaPhi | < DPHI_CUT
  analysis_cuts.push_back({ .name="delta_phi", .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      if (dPhi>=DPHI_CUT) return 0;
			      return 1;
			    } });

  // cut12: jet 1 tau32 < TOP_TAU32_CUT
  analysis_cuts.push_back({ .name="jet1_tau32",   .func = [&data](){
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<1) return 0; // for safety
			      double tau32 = data.jetsAK8Puppi.tau3[0];
			      if (data.jetsAK8Puppi.tau2[0]!=0) tau32 /= data.jetsAK8Puppi.tau2[0];
			      else tau32 = 9999;
			      if (tau32>=TOP_TAU32_CUT) return 0;
			      return 1;
			    } });

  // cut13: jet 2 tau32 < TOP_TAU32_CUT
  analysis_cuts.push_back({ .name="jet2_tau32",   .func = [&data](){ 
			      // Define cut function here:
			      if (data.jetsAK8Puppi.size<2) return 0; // for safety
			      double tau32 = data.jetsAK8Puppi.tau3[1];
			      if (data.jetsAK8Puppi.tau2[1]!=0) tau32 /= data.jetsAK8Puppi.tau2[1];
			      else tau32 = 9999;
			      if (tau32>=TOP_TAU32_CUT) return 0;
			      return 1;
			    } });

  // cut14: R >= R_CUT
  analysis_cuts.push_back({ .name="R", .func = [&data](){
			      // Define cut function here:
			      if (data.evt.AK8Puppi_R>=R_CUT) return 0;
			      return 1;
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
//     Must define it, because we blind it in data!
bool
Analysis::signal_selection(const DataStruct& data) {
  return _apply_ncut(analysis_cuts.size());
}

//_______________________________________________________
//      Define Histo options: Filling/Postfixes
Analysis::PostfixOptions
Analysis::get_pf_opts_(std::vector<std::vector<Sample> > lists, std::string dirname) {
  std::vector<Sample> samples;
  for (auto list : lists) samples.insert(samples.end(), list.begin(), list.end());
  PostfixOptions opt{ (size_t)-1, "", "", "" };
  for (size_t i=0; i<samples.size(); ++i) {
    // Find index of matching directory
    for (size_t j=0; j<samples[i].dirs.size(); ++j)
      if (samples[i].dirs[j] == dirname) opt.index = i;
    opt.postfixes += samples[i].postfix;
    opt.legends += samples[i].legend;
    opt.colors += samples[i].color;
    if (i+1!=samples.size()) {
      opt.postfixes +=  ";";
      opt.legends += ";";
      opt.colors += ",";
    }
  }
  return opt;
}

void
Analysis::define_histo_options(const double& weight, const DataStruct& d, const unsigned int& syst_nSyst, const unsigned int& syst_index, std::string dirname, bool runOnSkim=false)
{
  sh.SetHistoWeights({ [&weight](){ return weight; } });

  // --------------------------------------------------------------------
  //                            Colors
  // --------------------------------------------------------------------

  // Common Histo colorings
  // 400 kYellow  800 kOrange
  // 416 kGreen	  820 kSpring
  // 432 kCyan	  840 kTeal
  // 600 kBlue	  860 kAzure
  // 616 kMagenta 880 kViolet
  // 632 kRed     900 kPink
  
  std::string col3_red_to_blue = "633,618,601,"; // red, purple, blue
  std::string col4_cyan_to_red = "434,601,618,633,"; // Cyan, blue, purple, red
  std::string col5_green_to_red = "418,434,601,618,633,"; // green, cyan, blue, purple, red
  std::string col5_red_to_green = "633,618,601,434,418,"; // red, , purple, blue, cyan, green
  std::string col6_rainbow_dark = "601,434,418,402,633,618,"; // blue, cyan, green, yellow, red, purple
  std::string col8 = "1,601,434,418,402,807,633,618,"; // above plus black and orange
  std::string col10 = "4,6,2,800,402,417,433,9,618,633,";
  std::string col12 = "1,4,6,2,800,402,417,433,9,618,633,924,"; // black, blue, magenta, red, orange, darker yellow, darker green, darker cyan, blue-purple, dark purple, dark red
  std::string col12_rainbow = "402,416,433,600,617,632,802,813,833,863,883,892,"; // Go around first bright and then dark colors

  // --------------------------------------------------------------------
  //                            Cuts
  // --------------------------------------------------------------------

  sh.AddNewCut("JetHasMatchedGenTop", [&d](){ return jetsAK8Puppi_HasNearGenTop[d.jetsAK8Puppi.it]==1; });
  sh.AddNewCut("IsGenTop",            [&d](){ return is_Gen_Top[d.gen.it]; });
  //sh.AddNewCut("HLT",               [this](){ return _apply_cut("hlt_ak8ht700_mass50"); });

  //__________________________________
  //            Postfixes
  //__________________________________

  // Postfixes are vector definitions for histograms
  // They attach _<string> after histogram names
  // where <string> is chosen from a vector of strings
  // you need to give a natural number as a vector index
  // for the histogram to choose which histo to fill

  // Notation:
  // AddNewPostfix("Name of Postfix", lambda function returning non-negative integer, "postfix 1;2;3; ...", "legend text 1;2;3; ...", "ROOT color 1,2,3, ...")

  // Systematics variations
  sh.AddNewPostfix("Syst", [&syst_index]() { return syst_index; }, std::string(";Syst[0to")+std::to_string(syst_nSyst)+"]", std::string(";systematics [0to")+std::to_string(syst_nSyst)+"]", "1-999");
  if (syst_nSyst>998) std::cout<<"Warning: Too large number of systematics, define more colors!"<<std::endl;

  // Jets
  //sh.AddNewPostfix("AK4JetsPtOrdered", [&d](){ return d.jetsAK4.it; }, "Jet[1to10]", "1st Jet;2nd Jet;3rd Jet;[4to10]th Jet", "1-10");
  sh.AddNewPostfix("PtOrderedJets", [&d]()
		    { 
		      if (d.jetsAK8Puppi.it>=4) return (size_t)-1;
		      return (size_t)d.jetsAK8Puppi.it;
		    }, "Jet[1to4]", "1st Jet;2nd Jet;3rd Jet;4th Jet", col4_cyan_to_red);

  sh.AddNewPostfix("CutFlow",       [this]()
		    { 
		      for (int i=0; i<13; ++i) if (!analysis_cuts[i].func()) return i;
		      return 13;
		    }, "2Jet;jet1ID;jet2ID;Jet1Eta;Jet2Eta;Jet1Pt;Jet2Pt;Jet1Mass;Jet2Mass;Jet1Tau32;Jet2Tau32;DeltaPhi;R;PassAll", "2 jet;jet1 loose id;jet2 loose id;jet1 eta;jet2 eta;jet1 pt;jet2 pt;jet1 mass;jet2 mass;jet1 tau32;jet2 tau32;#Delta#phi;R;pass all", col10);

  // Cut Postfixes
  sh.AddNewPostfix("Pass0Cuts",       [this]() { return 0;                               }, "Pass0Cuts", "No cuts", "1");
  sh.AddNewPostfix("Pass1Cuts",       [this]() { return _apply_ncut(1) ? 0 : (size_t)-1; }, "Pass1Cuts", "Cuts up to 2jet", "1");
  sh.AddNewPostfix("Pass2Cuts",       [this]() { return _apply_ncut(2) ? 0 : (size_t)-1; }, "Pass2Cuts", "Cuts up to jet1 id", "1");
  sh.AddNewPostfix("Pass3Cuts",       [this]() { return _apply_ncut(3) ? 0 : (size_t)-1; }, "Pass3Cuts", "Cuts up to jet2 id", "1");
  sh.AddNewPostfix("Pass4Cuts",       [this]() { return _apply_ncut(4) ? 0 : (size_t)-1; }, "Pass4Cuts", "Cuts up to jet1 eta", "1");
  sh.AddNewPostfix("Pass5Cuts",       [this]() { return _apply_ncut(5) ? 0 : (size_t)-1; }, "Pass5Cuts", "Cuts up to jet2 eta", "1");
  sh.AddNewPostfix("Pass6Cuts",       [this]() { return _apply_ncut(6) ? 0 : (size_t)-1; }, "Pass6Cuts", "Cuts up to jet1 pt", "1");
  sh.AddNewPostfix("Pass7Cuts",       [this]() { return _apply_ncut(7) ? 0 : (size_t)-1; }, "Pass7Cuts", "Cuts up to jet2 pt", "1");
  sh.AddNewPostfix("Pass8Cuts",       [this]() { return _apply_ncut(8) ? 0 : (size_t)-1; }, "Pass8Cuts", "Cuts up to jet1 mass", "1");
  sh.AddNewPostfix("Pass9Cuts",       [this]() { return _apply_ncut(9) ? 0 : (size_t)-1; }, "Pass9Cuts", "Cuts up to jet2 mass", "1");
  //sh.AddNewPostfix("Pass10Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("jet1_tau32") ? 0 : (size_t)-1; }, "Pass10Cuts", "Cuts up to jet1 tau32", "1");
  //sh.AddNewPostfix("Pass11Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("jet1_tau32") && _apply_cut("jet2_tau32") ? 0 : (size_t)-1; }, "Pass11Cuts", "Cuts up to jet2 tau32", "1");
  //sh.AddNewPostfix("Pass12Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("jet1_tau32") && _apply_cut("jet2_tau32") && _apply_cut("delta_phi") ? 0 : (size_t)-1; }, "Pass12Cuts", "Cuts up to delta phi", "1");
  //sh.AddNewPostfix("Pass13Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("jet1_tau32") && _apply_cut("jet2_tau32") && _apply_cut("delta_phi") && _apply_cut("R") ? 0 : (size_t)-1; }, "Pass13Cuts", "Cuts up to R", "1");
  sh.AddNewPostfix("Pass10Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("delta_phi") ? 0 : (size_t)-1; }, "Pass10Cuts", "Cuts up to delta phi", "1");
  sh.AddNewPostfix("Pass11Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("delta_phi") && _apply_cut("jet1_tau32")  ? 0 : (size_t)-1; }, "Pass11Cuts", "Cuts up to jet1 tau32", "1");
  sh.AddNewPostfix("Pass12Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("delta_phi") && _apply_cut("jet1_tau32") && _apply_cut("jet2_tau32") ? 0 : (size_t)-1; }, "Pass12Cuts", "Cuts up to jet tau32", "1");
  sh.AddNewPostfix("Pass13Cuts",      [this]() { return _apply_ncut(9) && _apply_cut("delta_phi") && _apply_cut("jet1_tau32") && _apply_cut("jet2_tau32") && _apply_cut("R") ? 0 : (size_t)-1; }, "Pass13Cuts", "Cuts up to R", "1");
  sh.AddNewPostfix("PassHLTCut",      [this]() { return _apply_cut("hlt_ak8ht700_mass50") ? 0 : (size_t)-1;                  }, "PassHLT",      "HLT cut", "1");
  sh.AddNewPostfix("PassDPhiCut",     [this]() { return _apply_cut("delta_phi")  ? 0 : (size_t)-1;                           }, "PassDPhi",     "#Delta#phi cut", "1");
  sh.AddNewPostfix("PassJetMassCuts", [this]() { return _apply_cut("jet1_mass") && _apply_cut("jet2_mass") ? 0 : (size_t)-1; }, "PassMassCuts", "Mass cuts", "1");

  //sh.AddNewPostfix("NSubJet",       [&d](){ return (size_t)d.jetsAK8Puppi.nSubJets[d.jetsAK8Puppi.it]; }, "NSubJet[0to4]", "N_{subjet}=[0to4]", "1-5");
  sh.AddNewPostfix("JetPtCut",      [&d](){ return (size_t)(d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it] >= TOP_PT_CUT); }, "PtBelow400;PtAbove400", "p_{T} < 400;p_{T} >= 400", "2,3");
  sh.AddNewPostfix("JetMassCut",    [&d](){ return (size_t)(d.jetsAK8Puppi.softDropMass[d.jetsAK8Puppi.it] >= TOP_SD_MASS_CUT_LOW && d.jetsAK8Puppi.softDropMass[d.jetsAK8Puppi.it] < TOP_SD_MASS_CUT_HIGH); }, "OutMassWindow;InMassWindow", "M_{Soft-drop} < 110 or M_{Soft-drop} >= 210;110 < M_{Soft-drop} < 210", "2,3");
  sh.AddNewPostfix("JetTau32Cut",   [&d](){ return (size_t)(d.jetsAK8Puppi.tau3[d.jetsAK8Puppi.it]/d.jetsAK8Puppi.tau2[d.jetsAK8Puppi.it] < TOP_TAU32_CUT); }, "FailTau32Cut;PassTau32Cut", "Fail #tau_{3}/#tau_{2} cut;Pass #tau_{3}/#tau_{2} cut", "2,3");

  sh.AddNewPostfix("JetGenTruth",          [&d](){ return (size_t)(jetsAK8Puppi_HasNearGenTop[d.jetsAK8Puppi.it]==0 ? 0 : jetsAK8Puppi_NearGenTopIsHadronic[d.jetsAK8Puppi.it]!=1 ? 1 : jetsAK8Puppi_DRNearGenWFromTop[d.jetsAK8Puppi.it]<0.6&&jetsAK8Puppi_DRNearGenBFromTop[d.jetsAK8Puppi.it]<0.6 ? 3 : 2); }, "NotTop;SemiLepTop;NonMergedHadTop;MergedHadTop", "Non-top jet;Semi-leptonic top;Non-Merged hadronic top;Merged hadronic top", "1,601,633,418");
  sh.AddNewPostfix("JetMatchedGenTopType", [&d](){ return (size_t)(jetsAK8Puppi_NearGenTopIsHadronic[d.jetsAK8Puppi.it]!=-9999 ? jetsAK8Puppi_NearGenTopIsHadronic[d.jetsAK8Puppi.it] : -1); }, "MatchedGenTopLep;MatchedGenHadTop", "Semi-leptonic top;Hadronic top", "4,2");
  sh.AddNewPostfix("TopSizeCut",             [&d](){ return (size_t)(jetsAK8Puppi_HasNearGenTop[d.jetsAK8Puppi.it]==1 ? jetsAK8Puppi_DRNearGenWFromTop[d.jetsAK8Puppi.it]<0.6&&jetsAK8Puppi_DRNearGenBFromTop[d.jetsAK8Puppi.it]<0.6 : -1); }, "TopSizeAbove0p6;TopSizeBelow0p6", "Non-merged top;Merged top", "2,4");
  sh.AddNewPostfix("JetPassToptag",      [&d](){ return (size_t)(passHadTopTag[d.jetsAK8Puppi.it]); }, "FailTopTag;PassTopTag", "Non top-tagged AK8 jet;Top-tagged AK8 jet", "2,3");
  sh.AddNewPostfix("JetPassSubjetBTag",  [&d](){ return passSubjetBTag[d.jetsAK8Puppi.it];          }, "FailSJBTag;PassSJBTag", "No subjet B-tag;Pass subjet B-tag", "2,3");
  
  // Event
  sh.AddNewPostfix("RBins",          [&d](){ return (size_t)((d.evt.AK8Puppi_R>=0.1)+(d.evt.AK8Puppi_R>=0.2)+(d.evt.AK8Puppi_R>=0.4)); }, "R0to0p1;R0p1to0p2;R0p2to0p4;R0p4", "0.0<R<0.1;0.1<R<0.2;0.2<R<0.4;R>=0.4", "1,4,418,633");
  sh.AddNewPostfix("AK8ChsRBins",    [&d](){ return (size_t)((d.evt.R>=0.1)+(d.evt.R>=0.2)+(d.evt.R>=0.4)); }, "AK8ChsR0to0p1;AK8ChsR0p1to0p2;AK8ChsR0p2to0p4;AK8ChsR0p4", "0.0<R (AK8 CHS)<0.1;0.1<R (AK8 CHS)<0.2;0.2<R (AK8 CHS)<0.4;R (AK8 CHS)>=0.4", "1,4,418,633");
  sh.AddNewPostfix("AK4ChsRBins",    [&d](){ return (size_t)((d.evt.AK4_R>=0.1)+(d.evt.AK4_R>=0.2)+(d.evt.AK4_R>=0.4)); }, "AK4ChsR0to0p1;AK4ChsR0p1to0p2;AK4ChsR0p2to0p4;AK4ChsR0p4", "0.0<R (AK4 CHS)<0.1;0.1<R (AK4 CHS)<0.2;0.2<R (AK4 CHS)<0.4;R (AK4 CHS)>=0.4", "1,4,418,633");
  sh.AddNewPostfix("TTHadMRBins",    [&d](){ return (size_t)(d.evt.TTHadMR >= 5000 ? -1 : d.evt.TTHadMR/500); }, "MR[250to4750++500]", "MR_{tt} [250to4750++500]#pm250", "1,4,418,401,807,633,618,1,4,418,401,807,633,618");

  sh.AddNewPostfix("RBands",         [&d](){ return d.evt.AK8Puppi_R <= R_CUT_LOW ? -1 : d.evt.AK8Puppi_R >= R_CUT; }, "RBelow0p4;RAbove0p4", "0.2 < R < 0.4;R >= 0.4", "4,2");
  sh.AddNewPostfix("DPhiBands",      [](){ return dPhi >= DPHI_CUT; }, "PassDPhiCut;FailDPhiCut", "Pass #Delta#phi_{t#bar{t}} cut;Fail #Delta#phi_{t#bar{t}} cut", "2,4");
  sh.AddNewPostfix("Tau32Cuts",  [&d](){ return passTau32Cuts==-9999 ? (size_t)-1 : passTau32Cuts; }, "FailAnyTau32Cut;PassTau32Cuts", "Fail any tau_{32} cuts;Pass both tau_{32} cuts", "4,2");
  //sh.AddNewPostfix("NHadTop",     [&d](){ return nHadTopTag; },   "0HadTopTag;1HadTopTag;2HadTopTag;3HadTopTag;4HadTopTag", "N_{top-tag,hadronic}=0;N_{top-tag,hadronic}=1;N_{top-tag,hadronic}=2;N_{top-tag,hadronic}=3;N_{top-tag,hadronic}=4", col5_green_to_red);
  //sh.AddNewPostfix("NHadTopPreTag",     [&d](){ return nHadTopPreTag; },  "0HadTopPreTag;1HadTopPreTag;2HadTopPreTag;3HadTopPreTag;4HadTopPreTag", "N_{top-like,hadronic}=0;N_{top-like,hadronic}=1;N_{top-like,hadronic}=2;N_{top-like,hadronic}=3;N_{top-like,hadronic}=4", col5_green_to_red);
  sh.AddNewPostfix("NSubjetBTag",   [&d](){ return nSubjetBTag>4 ? (size_t)-1 : nSubjetBTag; },   "[0to4]SJBTag", "N_{Subjet B-Tag} = [0to4]", col5_green_to_red);
  sh.AddNewPostfix("ABCD",           [&d](){ return 
		       d.jetsAK8Puppi.size<2 ? (size_t)-1 :
		       d.evt.AK8Puppi_R<R_CUT_LOW ? (size_t)-1 : // Exclude low R region
		       (d.evt.AK8Puppi_R>=R_CUT) + (passTau32Cuts)*2;
		   }, "A;B;C;D", "A: 0.2<R<0.4, Fail any tau_{32} cut;B: R >= 0.4, Fail any tau_{32} cut;C: 0.2<R<0.4, Pass tau_{32} cuts;D (SR): R >= 0.4, Pass tau_{32} cuts", "858,602,628,634");

  sh.AddNewPostfix("CutHtAll",        [&d](){ return d.evt.HtAll >= 1200; },  "HtAllBelow1200;HtAllAbove1200", "H_{T,all} < 1200;H_{T,all} >= 1200", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1450",       [&d](){ return d.evt.HtAll >= 1450; },  "HtAllBelow1450;HtAllAbove1450", "H_{T,all} < 1450;H_{T,all} >= 1450", "4,2"); // Best cut
  sh.AddNewPostfix("HtAll1500",       [&d](){ return d.evt.HtAll >= 1500; },  "HtAllBelow1500;HtAllAbove1500", "H_{T,all} < 1500;H_{T,all} >= 1500", "4,2");
  //IMP sh.AddNewPostfix("NGenLepFromTop",  [&d](){ return (size_t)d.evt.NGenLepFromTop; }, "FullHad;[1to4]LepTop", "[0to4]l (e/#mu, from top)", "1-5");

  // Trigger
  sh.AddNewPostfix("PFHT475",         [&d](){ if (d.hlt.PFHT475==-9999) return (size_t)-1; else return (size_t)d.hlt.PFHT475; }, "NoPassHLT_PFHT475;PassHLT_PFHT475", "Do not pass HLT_PFHT475;Pass HLT_PFHT475", "633;418");

  // Gen Particles
  sh.AddNewPostfix("GenTopType",      [&d](){ return (size_t)(gen_Top_Type[d.gen.it]!=-9999 ? gen_Top_Type[d.gen.it] : -1); }, "GenHadTop;GenTopLep", "Hadronic top;Semi-leptonic top", "2,4");

  // Sample postfixes
  // Determine them from the directory name the input file is in

  // Map directory names to postfix name, legend entry and color
  std::vector<Sample> bkg_ttbars;
  bkg_ttbars.push_back({ .postfix="TTJetsMGHT",    .legend="t#bar{t}",                      .color="634",/*DRed*/   .dirs={ 
			   /*"TTJets_HT-0to600",*/ "TTJets_HT-600to800", "TTJets_HT-800to1200", "TTJets_HT-1200to2500", "TTJets_HT-2500toInf" 
			 } });
  bkg_ttbars.push_back({ .postfix="TTJetsMG",      .legend="t#bar{t} (MadGraph)",           .color="903",/*DPink*/  .dirs={ "TTJets_madgraph" } });
  //bkg_ttbars.push_back({ .postfix="TTJetsMGFS",    .legend="t#bar{t} (Madgraph) - FastSim", .color="901",/*Pink*/   .dirs={ "TTJets_madgraph_FastSim" } });
  bkg_ttbars.push_back({ .postfix="TTJetsNLOFXFX", .legend="t#bar{t} (MC@NLO FxFx)",                 .color="617",/*Magent*/ .dirs={ "TTJets_amcatnloFXFX" } });
  bkg_ttbars.push_back({ .postfix="TTNLO",         .legend="t#bar{t} (MC@NLO, pythia8)",             .color="619",/*DMagen*/ .dirs={ "TT_amcatnlo_pythia8" } });
  bkg_ttbars.push_back({ .postfix="TTNLOHerwig",   .legend="t#bar{t} (MC@NLO, herwig)",              .color="620",/*DMagen*/ .dirs={ "TT_amcatnlo_herwig" } });
  bkg_ttbars.push_back({ .postfix="TTPowheg",      .legend="t#bar{t} (Powheg, pythia8)",             .color="803",/*DOran*/  .dirs={ "TT_powheg_pythia8" } });
  bkg_ttbars.push_back({ .postfix="TTPowhegmpiOFF",.legend="t#bar{t} (Powheg, pythia8, mpiOFF)",     .color="804",/*DOran*/  .dirs={ "TT_powheg_pythia8_mpiOFF" } });
  bkg_ttbars.push_back({ .postfix="TTPowhegnoCR",  .legend="t#bar{t} (Powheg, pythia8, noCR)",       .color="805",/*DOran*/  .dirs={ "TT_powheg_pythia8_noCR" } });
  bkg_ttbars.push_back({ .postfix="TTPowhegHerwig",.legend="t#bar{t} (Powheg, herwig)",              .color="803",/*DOran*/  .dirs={ "TT_powheg_herwig" } });

  std::vector<Sample> bkg_nonttbars;
  bkg_nonttbars.push_back({ .postfix="QCD",         .legend="QCD",                .color=  "4",/*Blue*/    
			      .dirs={ 
			      "QCD_HT100to200",  "QCD_HT200to300",   "QCD_HT300to500",   "QCD_HT500to700",
			      "QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"
			      //, "QCD_GenJets5_HT300to500",   "QCD_GenJets5_HT500to700",   "QCD_GenJets5_HT700to1000",
			      //  "QCD_GenJets5_HT1000to1500", "QCD_GenJets5_HT1500to2000", "QCD_GenJets5_HT2000toInf",
			    } });
  bkg_nonttbars.push_back({ .postfix="ZJets",       .legend="Z+jets", .color="401",/*Yellow*/  .dirs={ 
			      "ZJetsToNuNu_HT-100To200", "ZJetsToNuNu_HT-200To400", "ZJetsToNuNu_HT-400To600", "ZJetsToNuNu_HT-600ToInf",
			      "ZJetsToQQ"
			    } });
  bkg_nonttbars.push_back({ .postfix="DYJetsToLL",  .legend="DY#rightarrowll",    .color="403",/*DYellow*/ .dirs={ 
			      "DYJetsToLL_M-5to50_HT-100to200", "DYJetsToLL_M-5to50_HT-200to400", "DYJetsToLL_M-5to50_HT-400to600", "DYJetsToLL_M-5to50_HT-600toInf",
			      "DYJetsToLL_M-50_HT-100to200",    "DYJetsToLL_M-50_HT-200to400",    "DYJetsToLL_M-50_HT-400to600",    "DYJetsToLL_M-50_HT-600toInf"
			    } });
  bkg_nonttbars.push_back({ .postfix="WJets",       .legend="W+jets",   .color="418",/*Green*/   .dirs={ 
			      "WJetsToLNu_HT-100To200",  "WJetsToLNu_HT-200To400",   "WJetsToLNu_HT-400To600",
			      "WJetsToLNu_HT-600To800",  "WJetsToLNu_HT-800To1200", "WJetsToLNu_HT-1200To2500", "WJetsToLNu_HT-2500ToInf",
			      "WJetsToQQ"
			    } });
  bkg_nonttbars.push_back({ .postfix="TTX",         .legend="t#bar{t}+X",       .color="803",/*Brown*/   .dirs={ 
			      "TTWJetsToLNu", "TTWJetsToQQ",
			      "TTZToLLNuNu", "TTZToQQ",
			      "TTGJets",
			      "TTTT" 
			    } });
  bkg_nonttbars.push_back({ .postfix="Diboson",     .legend="Diboson",            .color="804",/*DOrange*/ .dirs={ 
			      "WWTo1L1Nu2Q", "WWTo2L2Nu", "WWTo4Q",
			      "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2L2Q", "WZTo3LNu",
			      "ZZTo2L2Nu", "ZZTo2L2Q", "ZZTo2Q2Nu", "ZZTo4L", "ZZTo4Q",
			      "ZHToTauTau", "ZH_HToBB_ZToLL", "ZH_HToBB_ZToNuNu"
			    } });
  bkg_nonttbars.push_back({ .postfix="Top",         .legend="Top",                .color="403",/*DYellow*/ .dirs={ 
			      "ST_s-channel_4f_leptonDecays",
			      "ST_t-channel_top_4f_inclusiveDecays", "ST_t-channel_antitop_4f_inclusiveDecays",
			      "ST_tW_top_5f_inclusiveDecays",        "ST_tW_antitop_5f_inclusiveDecays"
			    } });

  std::vector<Sample> bkg_all, bkg_selected;
  bkg_all.insert(bkg_all.end(), bkg_ttbars.begin(), bkg_ttbars.end());
  bkg_all.insert(bkg_all.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());
  bkg_selected.push_back(bkg_ttbars[5]); // powheg
  bkg_selected.insert(bkg_selected.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());

  std::vector<Sample> data_all, data_selected, single_ele, single_mu, met;
  //if ( runOnSkim ) {
  //  data_all.push_back({ .postfix="Data",      .legend="Data",             .color="1", .dirs={ "JetHT_25ns" } });
  //  data_all.push_back({ .postfix="SingleEle", .legend="Data (SingleEle)", .color="1", .dirs={ "SingleMuon_25ns" } });
  //  data_all.push_back({ .postfix="SingleMu",  .legend="Data (SingleMu)",  .color="1", .dirs={ "SingleElectron_25ns" } });
  //  data_all.push_back({ .postfix="MET",       .legend="Data (MET)",       .color="1", .dirs={ "MET_25ns" } });
  //} else {
  data_all.push_back({ .postfix="Data",      .legend="Data",             .color="1", .dirs={ "JetHT_25ns_2015C", "JetHT_25ns_2015D" } });
  data_all.push_back({ .postfix="SingleEle", .legend="Data (SingleEle)", .color="1", .dirs={ "SingleMuon_25ns_2015C", "SingleMuon_25ns_2015D" } });
  data_all.push_back({ .postfix="SingleMu",  .legend="Data (SingleMu)",  .color="1", .dirs={ "SingleElectron_25ns_2015C", "SingleElectron_25ns_2015D" } });
  data_all.push_back({ .postfix="MET",       .legend="Data (MET)",       .color="1", .dirs={ "MET_25ns_2015C", "MET_25ns_2015D" } });
  //}
  data_selected.push_back(data_all[0]);
  single_ele.push_back(data_all[1]);
  single_mu.push_back(data_all[2]);
  met.push_back(data_all[3]);

  std::vector<Sample> signal_all, signal_selected;
  signal_all.push_back({ .postfix="T1tttt",    .legend="T1tttt (M_{#tilde{g}}=1.5TeV, M_{#tilde{#chi}^{0}}=100GeV)",  .color="862",/*Azure*/ .dirs={ "SMS-T1tttt_mGluino-1500_mLSP-100_FullSim" } });
  signal_selected.push_back(signal_all[0]);

  //"T5ttttDeg (M_{#tilde{g}}=1TeV)","1",/*Black*/
  //"T1tttt (M_{#tilde{g}}=1.5TeV)", "862",/*Azure*/
  //"T1tttt (M_{#tilde{g}}=1.5TeV, M_{#tilde{#chi}^{0}}=100GeV)", "841",/*Teal*/     
  //"T1tttt (M_{#tilde{g}}=1.2TeV, M_{#tilde{#chi}^{0}}=800GeV)", "843",/*DarkTeal*/ 
  //"T5ttttDeg (M_{#tilde{g}}=1TeV, 2,3-body)",                   "12", /*DarkGrey*/ 
  //"T2ttDeg (M_{#tilde{t}}=350GeV)",                             "434",/*Cyan*/     

  // Sample postfixes
  static const PostfixOptions all_samples_opt=get_pf_opts_({data_all, bkg_all, signal_all}, dirname);
  sh.AddNewPostfix("AllSamples", [](){ return all_samples_opt.index; }, all_samples_opt.postfixes, all_samples_opt.legends, all_samples_opt.colors);

  static const PostfixOptions plot_samples_opt=get_pf_opts_({data_selected, signal_selected, bkg_selected}, dirname);
  sh.AddNewPostfix("PlotSamples", [](){ return plot_samples_opt.index; }, plot_samples_opt.postfixes, plot_samples_opt.legends, plot_samples_opt.colors);

  std::vector<Sample> background;
  std::vector<std::string> background_dirs;
  for (auto bkg : bkg_selected) for (auto dir : bkg.dirs) background_dirs.push_back(dir);
  background.push_back({ .postfix="Background", .legend="Background", .color="1", .dirs=background_dirs });
  static const PostfixOptions background_opt = get_pf_opts_({background}, dirname);
  sh.AddNewPostfix("Background",  [](){ return background_opt.index; }, background_opt.postfixes, background_opt.legends, background_opt.colors);

  static const PostfixOptions signals_background_opt = get_pf_opts_({signal_all, background}, dirname);
  sh.AddNewPostfix("Signals,Background",  [](){ return signals_background_opt.index; }, signals_background_opt.postfixes, signals_background_opt.legends, signals_background_opt.colors);

  static const PostfixOptions background_signal_opt = get_pf_opts_({background, signal_selected}, dirname);
  sh.AddNewPostfix("Background,Signal", [](){ return background_signal_opt.index; }, background_signal_opt.postfixes, background_signal_opt.legends, "633,601");

  static const PostfixOptions signals_ttbar_opt = get_pf_opts_({signal_all, bkg_ttbars}, dirname);
  sh.AddNewPostfix("Signals,TT",  [](){ return signals_ttbar_opt.index; }, signals_ttbar_opt.postfixes, signals_ttbar_opt.legends, signals_ttbar_opt.colors);

  static const PostfixOptions data_mc_opt = get_pf_opts_({data_selected, background}, dirname);
  sh.AddNewPostfix("Data,MC",  [](){ return data_mc_opt.index; }, data_mc_opt.postfixes, data_mc_opt.legends, "1,633");

  static const PostfixOptions single_lep_opt = get_pf_opts_({single_ele, single_mu}, dirname);
  sh.AddNewPostfix("SingleEle,SingleMu", [](){ return single_lep_opt.index; }, single_lep_opt.postfixes, single_lep_opt.legends, "1,633");

  //static const PostfixOptions triggers_opt = get_pf_opts_({data_selected, single_ele, single_mu, background}, dirname);
  //sh.AddNewPostfix("Triggers", [&d]()
  //                  {
  //                    bool Pass_any_PFHT = //(d.hlt.AK8PFJet360_TrimMass30==1) ||
  //                      (d.hlt.PFHT200==1) || (d.hlt.PFHT250==1) || (d.hlt.PFHT300==1) || (d.hlt.PFHT350==1) ||
  //                      (d.hlt.PFHT400==1) || (d.hlt.PFHT475==1) || (d.hlt.PFHT550_4Jet==1) ||
  //                      (d.hlt.PFHT600==1) || (d.hlt.PFHT650==1) || (d.hlt.PFHT650_4Jet==1) ||
  //                      (d.hlt.PFHT750_4Jet==1) ||(d.hlt.PFHT750_4JetPt50==1) ||
  //                      (d.hlt.PFHT800==1) || (d.hlt.PFHT900==1);
  //                    if (triggers_opt.index==0) {
  //                      // JetHT Dataset: Pass any low threshold HT Trigger
  //                      if (Pass_any_PFHT) return (size_t)0;
  //                    } else if (triggers_opt.index==1) {
  //                      // SingleElectron Dataset: Pass Ele22_eta2p1_WPLoose_Gsf
  //                      if (d.hlt.Ele22_eta2p1_WPLoose_Gsf==1) return (size_t)1;
  //                    } else if (triggers_opt.index==2) {
  //                      // SingleMuon Dataset: Pass Mu45_eta2p1
  //                      if (d.hlt.Mu45_eta2p1==1) return (size_t)2;
  //                    } else if (triggers_opt.index==3) {
  //                      // MC Datasets: Pass any low threshold HT Trigger
  //                      if (Pass_any_PFHT) return (size_t)3;
  //                    }
  //                    return (size_t)-1; 
  //                  }, "LowHT;Ele22;Mu45;MC", "Data: OR (HT triggers);Data: Ele22 (|#eta|<2.1, loose WP);Data: Mu45 (|#eta|<2.1);MC: OR (HT triggers)", "1,402,601,633");
  static const PostfixOptions triggers_opt = get_pf_opts_({data_selected, single_ele, single_mu, met, background}, dirname);
  sh.AddNewPostfix("Triggers", [&d]()
                    {
                      if (triggers_opt.index==0) {
			bool Pass_any_PFHT = //(d.hlt.AK8PFJet360_TrimMass30==1) ||
			  (d.hlt.PFHT400==1) || (d.hlt.PFHT475==1) || (d.hlt.PFHT600==1) || 
			  (d.hlt.PFHT650==1) || (d.hlt.PFHT800==1) ||
			  (d.hlt.PFHT550_4JetPt50==1) ||(d.hlt.PFHT650_4JetPt50==1) || (d.hlt.PFHT750_4JetPt50==1);
			// JetHT Dataset: Pass any low threshold HT Trigger
                        if (Pass_any_PFHT) return (size_t)0;
			/*
                      } else if (triggers_opt.index==1) {
                        // SingleElectron Dataset: Pass Ele22_eta2p1_WPLoose_Gsf
                        if (d.hlt.Ele22_eta2p1_WPLoose_Gsf==1) return (size_t)1;
                      } else if (triggers_opt.index==1) {
                        // SingleMuon Dataset: Pass Mu45_eta2p1
                        if (d.hlt.Mu45_eta2p1==1) return (size_t)1;
			*/
                      } else if (triggers_opt.index==1) {
                        // SingleElectron Dataset: Pass Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50
                        if (d.hlt.Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50==1) return (size_t)1;
                        //if (d.hlt.Ele22_eta2p1_WPLoose_Gsf==1) return (size_t)1;
                        //if (d.hlt.Ele15_IsoVVVL_PFHT350==1) return (size_t)1;
                      } else if (triggers_opt.index==2) {
                        // SingleMuon Dataset: Pass Mu30_eta2p1_PFJet150_PFJet5
			if (d.hlt.Mu30_eta2p1_PFJet150_PFJet50==1) return (size_t)1;
			//if (d.hlt.HLT_Mu45_eta2p1==1) return (size_t)1;
                        //if (d.hlt.Mu15_IsoVVVL_PFHT350==1) return (size_t)1;
                      } else if (triggers_opt.index==3) {
                        // MET Dataset: PFMET triggers not saved, just cut on events with MET>200
			if (d.met.Pt[0]>200) return (size_t)2;
                      } else if (triggers_opt.index==4) {
                        // MC Datasets: Pass any low threshold HT Trigger
                        return (size_t)3;
                      }
                      return (size_t)-1; 
                    }, "LowHT;SingleLep;MET;MC", "Data: JetHT (PFHT*);Data: SingleLep (lep + PFHT350);Data: MET (MET>200);Simulation", "1,601,618,633");
  
  //static const PostfixOptions triggers_opt = get_pf_opts_({data_selected, single_ele, single_mu, met, background}, dirname);
  //sh.AddNewPostfix("Triggers", [&d]()
  //                  {
  //      	      // Run on datasets without any trigger selection
  //      	      return triggers_opt.index;
  //                  }, "JetHT;SingleEle;SingleMu;MET;MC", "Data: JetHT;Data: SingleElectron;Data: SingeMuon;Data: MET;MC", "1,402,601,418,633");

  static const PostfixOptions abcd_signals_opt = get_pf_opts_({background, signal_all}, dirname);
  sh.AddNewPostfix("ABCD,Signals", [&d](){
                      return 
			abcd_signals_opt.index == (size_t)-1 ? (size_t)-1 :
			d.jetsAK8Puppi.size<2 ? (size_t)-1 : // For tau32 cuts
			d.evt.AK8Puppi_R<R_CUT_LOW ? (size_t)-1 : // Exclude low R region
                        abcd_signals_opt.index == 0 ? (d.evt.AK8Puppi_R>=R_CUT) + (passTau32Cuts)*2 :
                        abcd_signals_opt.index +3; }, std::string("A;B;C;D;")+signals_ttbar_opt.postfixes, std::string("A: 0.2<R<0.4, Fail any tau_{32} cut;B: R >= 0.4, Fail any tau_{32} cut;C: 0.2<R<0.4, Pass tau_{32} cuts;D (SR): R >= 0.4, Pass tau_{32} cuts")+signals_ttbar_opt.legends, std::string("858,602,628,634,")+signals_ttbar_opt.colors);

  // --------------------------------------------------------------------
  //                         Fill Parameters
  // --------------------------------------------------------------------

  // Muons
  sh.AddNewFillParam("MuEnergy",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.mu.E[d.mu.it];              }, .axis_title="Muon Energy (GeV)"});
  sh.AddNewFillParam("MuPt",              { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.mu.Pt[d.mu.it];             }, .axis_title="Muon p_{T} (GeV)"});
  //IMP sh.AddNewFillParam("MuDRJet",           { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.MuDRJet[d.mu.it];    }, .axis_title="#DeltaR (#mu, jet)"});
  //IMP sh.AddNewFillParam("MuRelPtJet",        { .nbin=  50, .bins={   0,    500}, .fill=[&d](){ return d.evt.MuRelPtJet[d.mu.it]; }, .axis_title="p_{T}^{rel} (#mu, jet) (GeV)"});
  //IMP sh.AddNewFillParam("MuJetCombMass",     { .nbin=  80, .bins={   0,   2000}, .fill=[&d](){ return d.evt.MuJetCombMass[d.mu.it]; }, .axis_title="Mass_{#mu+jet comb.} (GeV)"});

  // Electrons
  sh.AddNewFillParam("EleEnergy",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.ele.E[d.ele.it];              }, .axis_title="Electron Energy (GeV)"});
  sh.AddNewFillParam("ElePt",             { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.ele.Pt[d.ele.it];             }, .axis_title="Electron p_{T} (GeV)"});
  //IMP sh.AddNewFillParam("EleDRJet",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.EleDRJet[d.ele.it];   }, .axis_title="#DeltaR (e, jet)"});
  //IMP sh.AddNewFillParam("EleRelPtJet",       { .nbin=  50, .bins={   0,    500}, .fill=[&d](){ return d.evt.EleRelPtJet[d.ele.it];}, .axis_title="p_{T}^{rel} (e, jet) (GeV)"});
  //IMP sh.AddNewFillParam("EleJetCombMass",    { .nbin=  80, .bins={   0,   2000}, .fill=[&d](){ return d.evt.EleJetCombMass[d.ele.it]; }, .axis_title="Mass_{e+jet comb.} (GeV)"});

  // MET
  sh.AddNewFillParam("MetPt",             { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return d.met.Pt[0];                    }, .axis_title="MET p_{T} (GeV)"});

  // AK4 Jets
  sh.AddNewFillParam("AK4JetEnergy",      { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.E[d.jetsAK8Puppi.it];          }, .axis_title="AK4-jet Energy (GeV)"});
  sh.AddNewFillParam("AK4JetPt",          { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.Pt[d.jetsAK8Puppi.it];         }, .axis_title="AK4-jet p_{T} (GeV)"});
  //sh.AddNewFillParam("AK4JetMass",        { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK4.Mass[d.jetsAK8Puppi.it];       }, .axis_title="AK4-jet Mass (GeV)"});

  // Jets (AK8)
  sh.AddNewFillParam("JetEnergy",          { .nbin= 500, .bins={   0,  10000}, .fill=[&d](){ return d.jetsAK8Puppi.E[d.jetsAK8Puppi.it];                    }, .axis_title="AK8 (Puppi) jet Energy (GeV)", .def_range={0,3000} });
  sh.AddNewFillParam("JetPt",              { .nbin= 500, .bins={   0,  10000}, .fill=[&d](){ return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it];                   }, .axis_title="AK8 (Puppi) jet p_{T} (GeV)", .def_range={0,3000} });
  sh.AddNewFillParam("JetPtBins",          { .nbin=  15, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000, 5000}, .fill=[&d](){ return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet p_{T} (GeV)", .def_range={0,2000} });
  sh.AddNewFillParam("JetPtFewBins",       { .nbin=   6, .bins={0, 300, 400, 600, 1000, 2000, 5000}, .fill=[&d](){ return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet p_{T} (GeV)", .def_range={0,3000} });
  sh.AddNewFillParam("JetPtOneBin",        { .nbin=   1, .bins={ 400,    5000}, .fill=[&d](){ return d.jetsAK8Puppi.Pt[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet p_{T} (GeV)"});
  sh.AddNewFillParam("JetEta",             { .nbin=  40, .bins={  -4,      4},                  .fill=[&d](){ return d.jetsAK8Puppi.Eta[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet #eta"});
  sh.AddNewFillParam("JetPhi",             { .nbin=  16, .bins={-3.1416, 3.1416}, .fill=[&d](){ return d.jetsAK8Puppi.Phi[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet #phi"});
  sh.AddNewFillParam("JetNeutralHadronMultiplicity", { .nbin=  20, .bins={0,  20}, .fill=[&d](){ return d.jetsAK8Puppi.neutralHadronMultiplicity[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet Neutral Hadron Multiplicity"});
  sh.AddNewFillParam("JetChargedHadronMultiplicity", { .nbin=  50, .bins={0, 100}, .fill=[&d](){ return d.jetsAK8Puppi.ChargedHadronMultiplicity[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet Charged Hadron Multiplicity"});

  std::vector<double> Ebins = { 0, 100, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 5000, 10000};
  sh.AddNewFillParam("JetPhotonEnergy",              { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.PhotonEnergy[d.jetsAK8Puppi.it];              }, .axis_title="AK8 (Puppi) jet Photon Energy (GeV)"});
  sh.AddNewFillParam("JetElectronEnergy",            { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.ElectronEnergy[d.jetsAK8Puppi.it];            }, .axis_title="AK8 (Puppi) jet Electron Energy (GeV)"});
  sh.AddNewFillParam("JetMuonEnergy",                { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.MuonEnergy[d.jetsAK8Puppi.it];                }, .axis_title="AK8 (Puppi) jet Muon Energy (GeV)"});
  sh.AddNewFillParam("JetChargedMuEnergy",           { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.ChargeMuEnergy[d.jetsAK8Puppi.it];            }, .axis_title="AK8 (Puppi) jet Charged Mu Energy (GeV)"});
  sh.AddNewFillParam("JetChargedEmEnergy",           { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.chargedEmEnergy[d.jetsAK8Puppi.it];           }, .axis_title="AK8 (Puppi) jet Charged Em Energy (GeV)"});
  sh.AddNewFillParam("JetChargedHadronEnergy",       { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.chargedHadronEnergy[d.jetsAK8Puppi.it];       }, .axis_title="AK8 (Puppi) jet Charged Hadron Energy (GeV)"});
  sh.AddNewFillParam("JetNeutralHadronEnergy",       { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.neutralHadronEnergy[d.jetsAK8Puppi.it];       }, .axis_title="AK8 (Puppi) jet Neutral Hadron Energy (GeV)"});
  sh.AddNewFillParam("JetNeutralEmEnergy",           { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.neutralEmEnergy[d.jetsAK8Puppi.it];           }, .axis_title="AK8 (Puppi) jet Neutral Em Energy (GeV)"});
  sh.AddNewFillParam("JetHFHadronEnergy",            { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.HFHadronEnergy[d.jetsAK8Puppi.it];            }, .axis_title="AK8 (Puppi) jet HF Hadron Energy (GeV)"});
  sh.AddNewFillParam("JetHFEMEnergy",                { .nbin= Ebins.size()-1, .bins=Ebins, .fill=[&d](){ return d.jetsAK8Puppi.HFEMEnergy[d.jetsAK8Puppi.it];                }, .axis_title="AK8 (Puppi) jet HF EM Energy (GeV)"});
  sh.AddNewFillParam("JetNumberOfDaughters",         { .nbin= 120, .bins={    0,   120}, .fill=[&d](){ return d.jetsAK8Puppi.numberOfDaughters[d.jetsAK8Puppi.it];         }, .axis_title="AK8 (Puppi) jet Number Of Daughters"});
  sh.AddNewFillParam("JetPhotonMultiplicity",        { .nbin= 120, .bins={    0,   120}, .fill=[&d](){ return d.jetsAK8Puppi.photonMultiplicity[d.jetsAK8Puppi.it];        }, .axis_title="AK8 (Puppi) jet Photon Multiplicity"});
  sh.AddNewFillParam("JetElectronMultiplicity",      { .nbin=  10, .bins={    0,    10}, .fill=[&d](){ return d.jetsAK8Puppi.electronMultiplicity[d.jetsAK8Puppi.it];      }, .axis_title="AK8 (Puppi) jet Electron Multiplicity"});
  sh.AddNewFillParam("JetMuonMultiplicity",          { .nbin=  10, .bins={    0,    10}, .fill=[&d](){ return d.jetsAK8Puppi.muonMultiplicity[d.jetsAK8Puppi.it];          }, .axis_title="AK8 (Puppi) jet Muon Multiplicity"});
  sh.AddNewFillParam("JetNeutralMultiplicity",       { .nbin= 120, .bins={    0,   120}, .fill=[&d](){ return d.jetsAK8Puppi.neutralMultiplicity[d.jetsAK8Puppi.it];       }, .axis_title="AK8 (Puppi) jet Neutral Multiplicity"});
  sh.AddNewFillParam("JetChargedMultiplicity",       { .nbin= 120, .bins={    0,   120}, .fill=[&d](){ return d.jetsAK8Puppi.chargedMultiplicity[d.jetsAK8Puppi.it];       }, .axis_title="AK8 (Puppi) jet Charged Multiplicity"});
  sh.AddNewFillParam("JetChargedHadronMultiplicity", { .nbin= 120, .bins={    0,   120}, .fill=[&d](){ return d.jetsAK8Puppi.ChargedHadronMultiplicity[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet Charged Hadron Multiplicity"});
  sh.AddNewFillParam("JetNeutralHadronMultiplicity", { .nbin=  40, .bins={    0,    40}, .fill=[&d](){ return d.jetsAK8Puppi.neutralHadronMultiplicity[d.jetsAK8Puppi.it]; }, .axis_title="AK8 (Puppi) jet Neutral Hadron Multiplicity"});
  sh.AddNewFillParam("JetHFHadronMultiplicity",      { .nbin=  80, .bins={    0,    80}, .fill=[&d](){ return d.jetsAK8Puppi.HFHadronMultiplicity[d.jetsAK8Puppi.it];      }, .axis_title="AK8 (Puppi) jet HF Hadron Multiplicity"});
  sh.AddNewFillParam("JetHFEMMultiplicity",          { .nbin=  50, .bins={    0,    50}, .fill=[&d](){ return d.jetsAK8Puppi.HFEMMultiplicity[d.jetsAK8Puppi.it];          }, .axis_title="AK8 (Puppi) jet HF EM Multiplicity"});
  
  //sh.AddNewFillParam("JetMass",            { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.Mass[d.jetsAK8Puppi.it];                 }, .axis_title="AK8 (Puppi) jet Mass (GeV)"});
  sh.AddNewFillParam("JetPrunedMass",      { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.prunedMass[d.jetsAK8Puppi.it];           }, .axis_title="AK8 (Puppi) jet Pruned Mass (GeV)"});
  sh.AddNewFillParam("JetPrunedMassCoarse",{ .nbin= 200, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.prunedMass[d.jetsAK8Puppi.it];           }, .axis_title="AK8 (Puppi) jet Pruned Mass (GeV)"});
  sh.AddNewFillParam("JetFilteredMass",    { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.filteredMass[d.jetsAK8Puppi.it];         }, .axis_title="AK8 (Puppi) jet Filtered Mass (GeV)"});
  sh.AddNewFillParam("JetTrimmedMass",     { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.trimmedMass[d.jetsAK8Puppi.it];          }, .axis_title="AK8 (Puppi) jet Trimmed Mass (GeV)"});
  sh.AddNewFillParam("JetSoftDropMass",    { .nbin= 200, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.softDropMass[d.jetsAK8Puppi.it];         }, .axis_title="AK8 (Puppi) jet Soft-drop Mass (GeV)"});
  //sh.AddNewFillParam("JetTopMass",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.topMass[d.jetsAK8Puppi.it];              }, .axis_title="AK8 (Puppi) jet Top Mass (GeV)"});
  //sh.AddNewFillParam("JetMinMass",         { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return d.jetsAK8Puppi.minmass[d.jetsAK8Puppi.it];              }, .axis_title="AK8 (Puppi) jet Min. Subjet-pair Mass (GeV)"});
  //sh.AddNewFillParam("JetNSubJets",        { .nbin=  11, .bins={-0.5,   10.5}, .fill=[&d](){ return d.jetsAK8Puppi.nSubJets[d.jetsAK8Puppi.it];             }, .axis_title="AK8 (Puppi) jet N_{subjet}"});
  sh.AddNewFillParam("JetTau1",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8Puppi.tau1[d.jetsAK8Puppi.it];                 }, .axis_title="#tau_{1}"});
  sh.AddNewFillParam("JetTau2",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8Puppi.tau2[d.jetsAK8Puppi.it];                 }, .axis_title="#tau_{2}"});
  sh.AddNewFillParam("JetTau3",            { .nbin= 100, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8Puppi.tau3[d.jetsAK8Puppi.it];                 }, .axis_title="#tau_{3}"});
  sh.AddNewFillParam("JetTau21",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8Puppi.tau2[d.jetsAK8Puppi.it]/d.jetsAK8Puppi.tau1[d.jetsAK8Puppi.it];   }, .axis_title="#tau_{2}/#tau_{1}"});
  sh.AddNewFillParam("JetTau31",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8Puppi.tau3[d.jetsAK8Puppi.it]/d.jetsAK8Puppi.tau1[d.jetsAK8Puppi.it];   }, .axis_title="#tau_{3}/#tau_{1}"});
  sh.AddNewFillParam("JetTau32",           { .nbin=  50, .bins={   0,      1}, .fill=[&d](){ return d.jetsAK8Puppi.tau3[d.jetsAK8Puppi.it]/d.jetsAK8Puppi.tau2[d.jetsAK8Puppi.it];   }, .axis_title="#tau_{3}/#tau_{2}"});
  //IMP sh.AddNewFillParam("JetDRLep",           { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return d.evt.DRJetLep[d.jetsAK8Puppi.it];         }, .axis_title="#DeltaR (lepton, jet)"});
  //IMP sh.AddNewFillParam("JetRelPtLep",        { .nbin= 100, .bins={   0,    500}, .fill=[&d](){ return d.evt.RelPtJetLep[d.jetsAK8Puppi.it];      }, .axis_title="p_{T}^{rel} (lepton, jet) (GeV)"});
  sh.AddNewFillParam("SubJetBDisc",        { .nbin= 101, .bins={   0,   1.01}, .fill=[&d](){ return subjetBTagDiscr[d.jetsAK8Puppi.it];         }, .axis_title="AK8 (Puppi) jet - Max subjet B discr."});

  sh.AddNewFillParam("JetMatchedGenTopPt",        { .nbin= 400, .bins={   0,   2000}, .fill=[&d](){ return jetsAK8Puppi_PtNearGenTop[d.jetsAK8Puppi.it];         }, .axis_title="Gen. top p_{T} (GeV)"});
  sh.AddNewFillParam("JetMatchedGenTopPtCoarse",  { .nbin= 100, .bins={   0,   2000}, .fill=[&d](){ return jetsAK8Puppi_PtNearGenTop[d.jetsAK8Puppi.it];         }, .axis_title="Gen. top p_{T} (GeV)"});
  sh.AddNewFillParam("JetMatchedGenTopPtBins",    { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return jetsAK8Puppi_PtNearGenTop[d.jetsAK8Puppi.it]; }, .axis_title="Gen. top p_{T} (GeV)"});
  sh.AddNewFillParam("JetMatchedGenTopJetDR",     { .nbin=  50, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenTop[d.jetsAK8Puppi.it];      }, .axis_title="#DeltaR (Gen. top, jet)"});
  sh.AddNewFillParam("JetMatchedGenTopJetDRFine", { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenTop[d.jetsAK8Puppi.it];      }, .axis_title="#DeltaR (Gen. top, jet)"});
  sh.AddNewFillParam("GenBJetDR",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenBFromTop[d.jetsAK8Puppi.it];        }, .axis_title="#DeltaR (Gen. b, jet)"});
  sh.AddNewFillParam("GenBJetDRFine",      { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenBFromTop[d.jetsAK8Puppi.it];        }, .axis_title="#DeltaR (Gen. b, jet)"});
  sh.AddNewFillParam("GenWJetDR",          { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenWFromTop[d.jetsAK8Puppi.it];        }, .axis_title="#DeltaR (Gen. W, jet)"});
  sh.AddNewFillParam("GenWJetDRFine",      { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenWFromTop[d.jetsAK8Puppi.it];        }, .axis_title="#DeltaR (Gen. W, jet)"});
  sh.AddNewFillParam("GenLepJetDR",        { .nbin=  60, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenLepFromSLTop[d.jetsAK8Puppi.it];    }, .axis_title="#DeltaR (Gen. lep, jet)"});
  sh.AddNewFillParam("GenLepJetDRFine",    { .nbin= 600, .bins={   0,      6}, .fill=[&d](){ return jetsAK8Puppi_DRNearGenLepFromSLTop[d.jetsAK8Puppi.it];    }, .axis_title="#DeltaR (Gen. lep, jet)"});
  
  //IMP sh.AddNewFillParam("MaxSubJetCSV",       { .nbin=  9, .bins={ 0, 0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 }, .fill=[&d](){ return d.evt.maxSubjetCSV[d.jetsAK8Puppi.it]; }, .axis_title="Maximum Subjet CSV"});
  
  // Event variables
  sh.AddNewFillParam("NVertices",           { .nbin= 101, .bins={-0.5,   100.5}, .fill=[&d](){ return d.evt.NGoodVtx;                  }, .axis_title="N_{Vertices}"});
  //sh.AddNewFillParam("NJet",                { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.jetsAK8Puppi.size;                  }, .axis_title="N_{AK8 (Puppi) jet}"});
  sh.AddNewFillParam("NJet",                { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return nLooseJet;                            }, .axis_title="N_{AK8 (Puppi) jet, loose ID}"});
  sh.AddNewFillParam("NHadTopTag",          { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return nHadTopTag;                           }, .axis_title="N_{had. top tag}", .def_range={0,5} });
  //sh.AddNewFillParam("NJetSelected",        { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return nHadTopTag+d.evt.NTopLep;     }, .axis_title="N_{Hadronic AK8 (Puppi) jet}"});
  //sh.AddNewFillParam("NJetHadronic",        { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return nHadTopTag;                   }, .axis_title="N_{Hadronic AK8 (Puppi) jet}"});
  //sh.AddNewFillParam("NJetLeptonic",        { .nbin=  21, .bins={-0.5,    20.5}, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{Leptonic AK8 (Puppi) jet}"});
  //sh.AddNewFillParam("NLep",                { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NLep;                      }, .axis_title="N_{lepton}"});
  //IMP sh.AddNewFillParam("NLepTight",           { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neletight+d.evt.nmu;       }, .axis_title="N_{lepton}"});
  //IMP sh.AddNewFillParam("NMu",                 { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nmu;                       }, .axis_title="N_{muon}"});
  //IMP sh.AddNewFillParam("NEle",                { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nele;                      }, .axis_title="N_{electron}"});
  //IMP sh.AddNewFillParam("NEleTight",           { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neletight;                 }, .axis_title="N_{electron}"});
  //IMP sh.AddNewFillParam("NLepVeto",            { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neleveto+d.evt.nmuveto;    }, .axis_title="N_{lepton}"});
  //IMP sh.AddNewFillParam("NMuVeto",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.nmuveto;                   }, .axis_title="N_{muon}"});
  //IMP sh.AddNewFillParam("NEleVeto",            { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.neleveto;                  }, .axis_title="N_{electron}"});
  sh.AddNewFillParam("NHadTopPreTag",       { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return nHadTopPreTag;               }, .axis_title="N_{jet} (p_{T}>=400, pass Mass cut)"});
  sh.AddNewFillParam("NTopLep",             { .nbin=   6, .bins={-0.5,     5.5}, .fill=[&d](){ return d.evt.NTopLep;                   }, .axis_title="N_{top, leptonic}"});
  sh.AddNewFillParam("TTHadMR",             { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV)"});
  sh.AddNewFillParam("TTHadMRCoarse",       { .nbin=  20, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMR;                   }, .axis_title="M_{R,t#bar{t}} (GeV)"});
  sh.AddNewFillParam("TTHadMTR",            { .nbin=  50, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMTR;                  }, .axis_title="M_{T,t#bar{t}}^{R} (GeV)"});
  sh.AddNewFillParam("DPhiFine",            { .nbin=  64, .bins={   0,     3.2}, .fill=[](){ return dPhi;           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DPhiBins",            { .nbin=   9, .bins={ 0, 0.5, 1.0, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.15 }, .fill=[](){ return dPhi; }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DPhi",                { .nbin=  16, .bins={   0,     3.2}, .fill=[](){ return dPhi;           }, .axis_title="|#Delta#phi_{t#bar{t}}|"});
  sh.AddNewFillParam("DEta",                { .nbin=  20, .bins={   0,     5.0}, .fill=[](){ return dEta;           }, .axis_title="|#Delta#eta_{t#bar{t}}|"});
  sh.AddNewFillParam("DR",                  { .nbin=  24, .bins={   0,     6.0}, .fill=[](){ return dR;             }, .axis_title="|#DeltaR_{t#bar{t}}|"});
  sh.AddNewFillParam("Jet1Pt",              { .nbin=  15, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000, 5000}, .fill=[&d](){ return (d.jetsAK8Puppi.size<1) ? -9999. : d.jetsAK8Puppi.Pt[0]; }, .axis_title="Leading AK8 (Puppi) jet p_{T} (GeV)"});
  sh.AddNewFillParam("Jet2Pt",              { .nbin=  15, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000, 5000}, .fill=[&d](){ return (d.jetsAK8Puppi.size<2) ? -9999. : d.jetsAK8Puppi.Pt[1]; }, .axis_title="Subleading AK8 (Puppi) jet p_{T} (GeV)"});
  sh.AddNewFillParam("Jet1Mass",            { .nbin=  17, .bins={0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300, 500, 1000}, .fill=[&d](){ return (d.jetsAK8Puppi.size<1) ? -9999. : d.jetsAK8Puppi.softDropMass[0]; }, .axis_title="Leading AK8 (Puppi) jet M_{Soft-Drop} (GeV)"});
  sh.AddNewFillParam("Jet2Mass",            { .nbin=  17, .bins={0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300, 500, 1000}, .fill=[&d](){ return (d.jetsAK8Puppi.size<2) ? -9999. : d.jetsAK8Puppi.softDropMass[1]; }, .axis_title="Subleading AK8 (Puppi) jet M_{Soft-Drop} (GeV)"});
  sh.AddNewFillParam("SumPt",               { .nbin=  10, .bins={ 400, 600, 700, 750, 800, 850, 900, 1000, 1200, 2000, 10000}, .fill=[&d](){ if (d.jetsAK8Puppi.size<2) return -9999.; return double(d.jetsAK8Puppi.Pt[0]+d.jetsAK8Puppi.Pt[1]); }, .axis_title="Sum of leading + subleding jet p_{T} (GeV)", .def_range={400,2000} });
  sh.AddNewFillParam("SumPtOneBin",         { .nbin=   1, .bins={ 800,   10000}, .fill=[&d](){ if (d.jetsAK8Puppi.size<2) return -9999.; return double(d.jetsAK8Puppi.Pt[0]+d.jetsAK8Puppi.Pt[1]); }, .axis_title="Sum of leading + subleding jet p_{T} (GeV)"});
  sh.AddNewFillParam("TTHadSumPt",          { .nbin=  10, .bins={ 400, 600, 700, 750, 800, 850, 900, 1000, 1200, 2000, 10000}, .fill=[&d](){ return d.evt.TTHadSumPt; }, .axis_title="p_{T}^{Top1}+p_{T}^{Top2} (GeV)", .def_range={400,2000} });
  sh.AddNewFillParam("TTHadSumPtOneBin",    { .nbin=   1, .bins={ 800,   10000}, .fill=[&d](){ return d.evt.TTHadSumPt;                }, .axis_title="p_{T}^{Top1}+p_{T}^{Top2} (GeV)"});
  sh.AddNewFillParam("TTHadDEta",           { .nbin=  50, .bins={   0,       5}, .fill=[&d](){ return d.evt.TTHadDEta;                 }, .axis_title="#Delta#eta_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadDR",             { .nbin=  60, .bins={   0,       6}, .fill=[&d](){ return d.evt.TTHadDR;                   }, .axis_title="#DeltaR_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadPz",             { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadPz;                   }, .axis_title="p_{z,t#bar{t}} (GeV)"});
  sh.AddNewFillParam("TTHadDPz",            { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadDPz;                  }, .axis_title="#Deltap_{z,t#bar{t}} (GeV)"});
  sh.AddNewFillParam("TTHadHz",             { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadHz;                   }, .axis_title="H_{z,t#bar{t}} (GeV)"});
  sh.AddNewFillParam("TTHadMass",           { .nbin= 100, .bins={   0,    5000}, .fill=[&d](){ return d.evt.TTHadMass;                 }, .axis_title="M_{t#bar{t}} (GeV)"});
  sh.AddNewFillParam("TTHadR",              { .nbin=  24, .bins={   0,     1.2}, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadRFine",          { .nbin= 120, .bins={   0,     1.2}, .fill=[&d](){ return d.evt.TTHadR;                    }, .axis_title="R_{t#bar{t}}"});
  sh.AddNewFillParam("TTHadR2",             { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.TTHadR2;                   }, .axis_title="R_{t#bar{t}}^{2}"});
  sh.AddNewFillParam("AK4ChsR",                { .nbin=  24, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R (AK4)"});
  sh.AddNewFillParam("AK4ChsRFine",            { .nbin= 120, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.R;                         }, .axis_title="R (AK4)"});
  sh.AddNewFillParam("AK4ChsRBins",            { .nbin=  14, .bins={ 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2 }, .fill=[&d](){ return d.evt.R; }, .axis_title="R (AK4)"});
  sh.AddNewFillParam("AK8ChsR",             { .nbin=  24, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.AK4_R;                         }, .axis_title="R (AK4)"});
  sh.AddNewFillParam("AK8ChsRFine",         { .nbin= 120, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.AK4_R;                         }, .axis_title="R (AK4)"});
  sh.AddNewFillParam("AK8ChsRBins",         { .nbin=  14, .bins={ 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2 }, .fill=[&d](){ return d.evt.AK4_R; }, .axis_title="R (AK4)"});
  sh.AddNewFillParam("R",                   { .nbin=  24, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.AK8Puppi_R;                         }, .axis_title="R (AK8 Puppi)"});
  sh.AddNewFillParam("RFine",               { .nbin= 120, .bins={   0,    1.20}, .fill=[&d](){ return d.evt.AK8Puppi_R;                         }, .axis_title="R (AK8 Puppi)"});
  sh.AddNewFillParam("RBins",               { .nbin=  14, .bins={ 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2 }, .fill=[&d](){ return d.evt.AK8Puppi_R; }, .axis_title="R (AK8 Puppi)"});
  sh.AddNewFillParam("R2",                  { .nbin=  32, .bins={   0,     1.6}, .fill=[&d](){ return d.evt.AK8Puppi_R2;                        }, .axis_title="R^{2} (AK8 Puppi)"});
  sh.AddNewFillParam("MR",                  { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.AK8Puppi_MR;                        }, .axis_title="M_{R} (AK8 Puppi) (GeV)"});
  sh.AddNewFillParam("MTR",                 { .nbin= 100, .bins={   0,    2000}, .fill=[&d](){ return d.evt.AK8Puppi_MTR;                       }, .axis_title="M_{T}^{R} (AK8 Puppi) (GeV)"});
  sh.AddNewFillParam("HtTopFraction",       { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.HtTopFr;                   }, .axis_title="H_{T,tops}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV)"});
  sh.AddNewFillParam("HtExFraction",        { .nbin=  20, .bins={   0,       1}, .fill=[&d](){ return d.evt.HtExFr;                    }, .axis_title="H_{T,extra}/(H_{T}+H_{T,leptonic}+#slash{p}_{T}) (GeV)"});
  sh.AddNewFillParam("GenHt",               { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.Gen_Ht;                    }, .axis_title="H_{T}^{Gen} (GeV)",              .def_range={0, 6000} });
  sh.AddNewFillParam("AK8Ht",               { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return d.evt.Ht;                        }, .axis_title="H_{T}^{AK8 (CHS) jets} (GeV)",   .def_range={0, 6000} });
  sh.AddNewFillParam("AK4Ht",               { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return AK4_Ht;                          }, .axis_title="H_{T}^{AK4 (CHS) jets} (GeV)",   .def_range={0, 6000} });
  sh.AddNewFillParam("AK8PuppiHt",          { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return AK8Puppi_Ht;                     }, .axis_title="H_{T}^{AK8 (Puppi) jets} (GeV)", .def_range={0, 6000} });
  sh.AddNewFillParam("AK4PuppiHt",          { .nbin= 100, .bins={   0,   10000}, .fill=[&d](){ return AK4Puppi_Ht;                     }, .axis_title="H_{T}^{AK4 (Puppi) jets} (GeV)", .def_range={0, 6000} });
  sh.AddNewFillParam("AK8HtBins",           { .nbin=20,   .bins={ 0, 200, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1500, 2000, 3000, 6000, 10000}, .fill=[&d](){ return d.evt.Ht; }, .axis_title="H_{T}^{AK8 (Puppi) jets} (GeV)", .def_range={0, 6000} });
  sh.AddNewFillParam("AK8HtOneBin",            { .nbin=   1, .bins={800, 10000}, .fill=[&d](){ return d.evt.Ht;                        }, .axis_title="H_{T}^{AK8 (Puppi) jets} (GeV)"});
  sh.AddNewFillParam("HtAllCoarse",         { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HtAll;                     }, .axis_title="H_{T}^{AK8 (Puppi) jets}+H_{T}^{leptons}+#slash{p}_{T} (GeV)"});
  sh.AddNewFillParam("HtAll",               { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HtAll;                     }, .axis_title="H_{T}^{AK8 (Puppi) jets}+H_{T}^{leptons}+#slash{p}_{T} (GeV)"});
  sh.AddNewFillParam("HtTop",               { .nbin=  25, .bins={   0,    5000}, .fill=[&d](){ return d.evt.HtTop;                     }, .axis_title="H_{T}^{tops} (GeV)"});
  sh.AddNewFillParam("HtEx",                { .nbin=  50, .bins={   0,   10000}, .fill=[&d](){ return d.evt.HtEx;                      }, .axis_title="H_{T}^{extra} (GeV)"});
  
  sh.AddNewFillParam("FlaggoodVertices",                       { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.goodVertices;                       }, .axis_title="goodVertices"});
  sh.AddNewFillParam("FlageeBadScFilter",                      { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.eeBadScFilter;                      }, .axis_title="eeBadScFilter"});
  sh.AddNewFillParam("FlagEcalDeadCellTriggerPrimitiveFilter", { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.EcalDeadCellTriggerPrimitiveFilter; }, .axis_title="EcalDeadCellTriggerPrimitiveFilter"});
  sh.AddNewFillParam("FlagHBHENoiseFilter",                    { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.HBHENoiseFilter;                    }, .axis_title="HBHENoiseFilter"});
  sh.AddNewFillParam("FlagHBHENoiseIsoFilter",                 { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.HBHENoiseIsoFilter;                 }, .axis_title="HBHENoiseIsoFilter"});
  sh.AddNewFillParam("FlagCSCTightHalo2015Filter",             { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.CSCTightHalo2015Filter;             }, .axis_title="CSCTightHalo2015Filter"});
  sh.AddNewFillParam("FlagmuonBadTrackFilter",                 { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.muonBadTrackFilter;                 }, .axis_title="muonBadTrackFilter"});
  sh.AddNewFillParam("FlagchargedHadronTrackResolutionFilter", { .nbin=  2, .bins={-0.5,      1.5}, .fill=[&d](){ return d.filter.chargedHadronTrackResolutionFilter; }, .axis_title="chargedHadronTrackResolutionFilter"});
  
  // Gen Particles
  sh.AddNewFillParam("GenTopPtBins",      { .nbin=  14, .bins={0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 800, 1000, 1400, 2000}, .fill=[&d](){ return is_Gen_Top[d.gen.it] ? d.gen.Pt[d.gen.it] : -9999; }, .axis_title="Gen. top p_{T} (GeV)"});
  sh.AddNewFillParam("GenTopPtFewBins",   { .nbin=  4, .bins={0, 300, 400, 800, 2000}, .fill=[&d](){ return is_Gen_Top[d.gen.it] ? d.gen.Pt[d.gen.it] : -9999; }, .axis_title="Gen. top p_{T} (GeV)"});

  // Special Y/Z axis parameters:
  sh.AddSpecial({ .name="MergedTopFraction",         .name_plus_1d="IsGenTopMerged",          .axis="Fraction of Merged Tops",             .axis_plus_1d="Gen. b and W Merged in R<0.8 (bool)"});
  sh.AddSpecial({ .name="JetFindingEfficiency",      .name_plus_1d="HasJet",                  .axis="Jet finding Efficiency",              .axis_plus_1d="Found AK8 (Puppi) jet (bool)"});
  sh.AddSpecial({ .name="TopFindingEfficiency",      .name_plus_1d="HasHadTopTaggedJet",      .axis="Top finding Efficiency",              .axis_plus_1d="Has hadronic top-tagged jet (bool)"});
  sh.AddSpecial({ .name="TopTagEfficiency",          .name_plus_1d="JetIsHadTopTagged",       .axis="Top-tagging Efficiency",              .axis_plus_1d="Jet is hadronic top-tagged (bool)"});
  sh.AddSpecial({ .name="MisTagRate",                .name_plus_1d="JetHasNoGenTop",          .axis="Mis-tag Rate",                        .axis_plus_1d="Jet is mis-matched (bool)"});
  sh.AddSpecial({ .name="HLTEff_AK8PFJet360",        .name_plus_1d="HLTAK8PFJet360",          .axis="HLT_AK8PFJet360_TrimMass30",          .axis_plus_1d="#epsilon_{HLT_AK8PFJet360_TrimMass30}"});
  sh.AddSpecial({ .name="HLTEff_PFJet450",           .name_plus_1d="HLTPFJet450",             .axis="HLT_PFJet450",                        .axis_plus_1d="#epsilon_{HLT_PFJet450}"});
  sh.AddSpecial({ .name="HLTEff_DiPFJetAve500",      .name_plus_1d="HLTDiPFJetAve500",        .axis="HLT_DiPFJetAve500",                   .axis_plus_1d="#epsilon_{HLT_DiPFJetAve500}"});
  sh.AddSpecial({ .name="HLTEff_AK8PFHT650",         .name_plus_1d="HLTAK8PFHT650",           .axis="HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", .axis_plus_1d="#epsilon_{HLT_AK8PFHT650_TrimR0p1PT0p03Mass50}"});
  sh.AddSpecial({ .name="HLTEff_AK8PFHT700",         .name_plus_1d="HLTAK8PFHT700",           .axis="HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", .axis_plus_1d="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}"});
  sh.AddSpecial({ .name="HLTEff_PFHT750_4JetPt50",   .name_plus_1d="HLTPFHT750_4JetPt50",     .axis="HLT_PFHT750_4JetPt50",                .axis_plus_1d="#epsilon_{HLT_PFHT750_4JetPt50}"});
  sh.AddSpecial({ .name="HLTEff_PFHT800",            .name_plus_1d="HLTPFHT800",              .axis="HLT_PFHT800",                         .axis_plus_1d="#epsilon_{HLT_PFHT800}"});
  //sh.AddSpecial({ .name="HLTEff_PFHT350",            .name_plus_1d="HLTPFHT350",              .axis="HLT_PFHT350",                         .axis_plus_1d="#epsilon_{HLT_PFHT350}"});
  //sh.AddSpecial({ .name="HLTEff_PFHT900",            .name_plus_1d="HLTPFHT900",              .axis="HLT_PFHT900",                         .axis_plus_1d="#epsilon_{HLT_PFHT900}"});
  sh.AddSpecial({ .name="Counts",                    .name_plus_1d="Syst",                    .axis="Counts (Incl Syst Unc)",              .axis_plus_1d="Systematics variation index"});

  sh.AddNewFillParam("MergedTopFraction",         { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return top_Children_Within_Cone[d.jetsAK8Puppi.it]; }, .axis_title="Fraction of Merged Tops" });
  sh.AddNewFillParam("JetFindingEfficiency",      { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return has_Matched_Jet[d.gen.it]; }, .axis_title="Jet finding Efficiency" });
  sh.AddNewFillParam("TopFindingEfficiency",      { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return has_Matched_Tagged_Jet[d.gen.it]; }, .axis_title="Top finding Efficiency" });
  sh.AddNewFillParam("TopTagEfficiency",          { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return passHadTopTag[d.jetsAK8Puppi.it]; }, .axis_title="Top-tagging Efficiency" });
  sh.AddNewFillParam("HLTEff_AK8PFJet360",        { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.AK8PFJet360_TrimMass30;          }, .axis_title="#epsilon_{HLT_AK8PFJet360_TrimMass30}" });
  sh.AddNewFillParam("HLTEff_PFJet450",           { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.PFJet450;                        }, .axis_title="#epsilon_{HLT_PFJet450}" });
  sh.AddNewFillParam("HLTEff_DiPFJetAve500",      { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.DiPFJetAve500;                   }, .axis_title="#epsilon_{HLT_DiPFJetAve500}" });
  sh.AddNewFillParam("HLTEff_AK8PFHT650",         { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.AK8PFHT650_TrimR0p1PT0p03Mass50; }, .axis_title="#epsilon_{HLT_AK8PFHT650_TrimR0p1PT0p03Mass50}" });
  sh.AddNewFillParam("HLTEff_AK8PFHT700",         { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.AK8PFHT700_TrimR0p1PT0p03Mass50; }, .axis_title="#epsilon_{HLT_AK8PFHT700_TrimR0p1PT0p03Mass50}" });
  sh.AddNewFillParam("HLTEff_PFHT750_4JetPt50",   { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.PFHT750_4JetPt50;                }, .axis_title="#epsilon_{HLT_PFHT750_4JetPt50}" });
  sh.AddNewFillParam("HLTEff_PFHT800",            { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.PFHT800;                         }, .axis_title="#epsilon_{HLT_PFHT800}" });
  //sh.AddNewFillParam("HLTEff_PFHT350",            { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.PFHT350;                         }, .axis_title="#epsilon_{HLT_PFHT350}" });
  //sh.AddNewFillParam("HLTEff_PFHT900",            { .nbin=  2, .bins={ -0.5, 1.5}, .fill=[&d](){ return d.hlt.PFHT900;                         }, .axis_title="#epsilon_{HLT_PFHT900}" });
  sh.AddNewFillParam("Counts",                    { .nbin= 1+syst_nSyst, .bins={-0.5, syst_nSyst+0.5}, .fill=[&syst_index](){ return syst_index; }, .axis_title="Counts (Incl Syst Unc)"});

}


//_______________________________________________________
//                 List of Histograms
//TH1D* h_nvtx;
TH2D* h_abcd;

//_______________________________________________________
//              Define Histograms here
void
Analysis::init_analysis_histos(const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
  //h_nvtx = new TH1D("nvtx",";N_{Vertices}", 100,0,100);
  Double_t R_bins[4] = { 0, R_CUT_LOW, R_CUT, 2 };
  Double_t Ntop_bins[3] = { 0, 1, 2 };
  h_abcd = new TH2D("abcd",";R;Pass both tau32 cuts", 3, R_bins, 2, Ntop_bins );

  //__________________________________
  //        Define Smarthistos
  //__________________________________

  // Define histo types (for different object to loop on, and different cuts to apply)
  sh.AddHistoType("mu");
  sh.AddHistoType("ele");
  sh.AddHistoType("jetsAK4");
  sh.AddHistoType("jetsAK8Puppi");
  sh.AddHistoType("gen");
  sh.AddHistoType("all events");
  sh.AddHistoType("baseline notrigger");
  sh.AddHistoType("baseline events");
  // Plots including systematics variations
  sh.AddHistoType("jetsAK8Puppi syst");
  sh.AddHistoType("all events syst");
  sh.AddHistoType("baseline events syst");

  std::string stack_plot_opt = "LogSumw2Stack2AddRatio";
  //std::string stack_plot_opt = "LogSumw2Stack1";

  // --------------------------------------------------------------------------
  //                                   AK8 Jets

  std::vector<std::string> allcuts;
  for (int i=0; i<=8; ++i) { std::stringstream ss; ss<<i<<"Cut"; allcuts.push_back(ss.str()); }

  // Plots:
  // NJet - Cutflow, PlotSamples + Cuts

  // cuts:
  // 2 AK8 jets
  // id   1,2
  // eta  1,2
  // pt   1,2
  // mass 1,2
  // HLT

  // Things to check:
  // HLT->NJet
  // HLT->ID/Eta cutflow
  // HLT turnons

  // --------------------------------------
  //            Jet kinematics
  //         (pt, eta, phi, mass)

  bool Systematics_Only = false;

  if (!Systematics_Only) {

    // Cutflow plots
    //sh.AddHistos("all events", { .fill="NJet",              .pfs={"CutFlow","PassHLTCut","Data,MC"},                 .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="NHadTopTag",        .pfs={"CutFlow","PassHLTCut","Data,MC"},                 .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",           .pfs={"CutFlow","PassHLTCut","PtOrderedJets","Data,MC"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetPhi",          .pfs={"CutFlow","PassHLTCut","PtOrderedJets","Data,MC"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetEta",          .pfs={"CutFlow","PassHLTCut","PtOrderedJets","Data,MC"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",        .pfs={"CutFlow","PassHLTCut","PtOrderedJets","Data,MC"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass", .pfs={"CutFlow","PassHLTCut","PtOrderedJets","Data,MC"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0,      1.01e-3,1e7, 0.55,0.9} });

    // Data-MC Comparison - background composition
    for (size_t i=0; i<=13; ++i) { 
      std::stringstream cut; cut<<"Pass"<<i<<"Cuts";
      // Jet quantities
      sh.AddHistos("all events",        { .fill="NJet",            .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events",        { .fill="NHadTopTag",      .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPt",           .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPt",           .pfs={"PlotSamples",cut.str(),"PassHLTCut","PtOrderedJets"},     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPhi",          .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPhi",          .pfs={"PlotSamples",cut.str(),"PassHLTCut","PtOrderedJets"},     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetEta",          .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetEta",          .pfs={"PlotSamples",cut.str(),"PassHLTCut","PtOrderedJets"},     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetTau32",        .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetTau32",        .pfs={"PlotSamples",cut.str(),"PassHLTCut","PtOrderedJets"},     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetSoftDropMass", .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetSoftDropMass", .pfs={"PlotSamples",cut.str(),"PassHLTCut","PtOrderedJets"},     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      // Subjet BTags
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPt",           .pfs={"PlotSamples","NSubjetBTag",cut.str(),"PassHLTCut"},                       .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPt",           .pfs={"PlotSamples","NSubjetBTag",cut.str(),"PassHLTCut","PtOrderedJets"},       .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPt",           .pfs={"PlotSamples",cut.str(),"PassHLTCut","JetPassSubjetBTag"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="JetPt",           .pfs={"PlotSamples",cut.str(),"PassHLTCut","JetPassSubjetBTag","PtOrderedJets"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="SubJetBDisc",     .pfs={"PlotSamples",cut.str(),"PassHLTCut"},                                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("jetsAK8Puppi",      { .fill="SubJetBDisc",     .pfs={"PlotSamples",cut.str(),"PassHLTCut","PtOrderedJets"},                     .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    }
  
    // --------------------------------------
    //             Other variables
    //       (energies, multiplicities)

    //sh.AddHistos("jetsAK8Puppi", { .fill="JetPhotonEnergy",              .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetElectronEnergy",            .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetMuonEnergy",                .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetChargedMuEnergy",           .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetChargedEmEnergy",           .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetChargedHadronEnergy",       .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetNeutralHadronEnergy",       .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetNeutralEmEnergy",           .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetHFHadronEnergy",            .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetHFEMEnergy",                .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetNumberOfDaughters",         .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetPhotonMultiplicity",        .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetElectronMultiplicity",      .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetMuonMultiplicity",          .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetNeutralMultiplicity",       .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetChargedMultiplicity",       .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetChargedHadronMultiplicity", .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetNeutralHadronMultiplicity", .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetHFHadronMultiplicity",      .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetHFEMMultiplicity",          .pfs={"PlotSamples","PtOrderedJets","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={200,2000, 1.01e-3,1e7, 0.55,0.9} });

    // --------------------------------------------------------------------------
    //                                 Trigger

    sh.AddHistos("all events", { .fill="GenHt",      .pfs={"AllSamples"},                           .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK4Ht",      .pfs={"AllSamples"},                           .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK8Ht",      .pfs={"AllSamples"},                           .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK4PuppiHt", .pfs={"AllSamples"},                           .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK8PuppiHt", .pfs={"AllSamples"},                           .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="GenHt",      .pfs={"AllSamples","PassHLTCut"},              .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK4Ht",      .pfs={"AllSamples","PassHLTCut"},              .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK8Ht",      .pfs={"AllSamples","PassHLTCut"},              .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK4PuppiHt", .pfs={"AllSamples","PassHLTCut"},              .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="AK8PuppiHt", .pfs={"AllSamples","PassHLTCut"},              .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });
    sh.AddHistos("all events", { .fill="GenHt",      .pfs={"AllSamples","Pass9Cuts","PassHLTCut"},  .cuts={},  .draw="HISTE1", .opt="Sumw2Norm", .ranges={0,0, 0,0, 0.55,0.9} });

    // sumpt/ht distribution changes
    //sh.AddHistos("all events", { .fill="SumPt",           .pfs={"CutFlow","Triggers"},             .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="SumPt",           .pfs={"CutFlow","Data,MC","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK8Ht",           .pfs={"CutFlow","Triggers"},             .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK8Ht",           .pfs={"CutFlow","Data,MC","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK4Ht",           .pfs={"CutFlow","Triggers"},             .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK4Ht",           .pfs={"CutFlow","Data,MC","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK8PuppiHt",      .pfs={"CutFlow","Triggers"},             .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK8PuppiHt",      .pfs={"CutFlow","Data,MC","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK4PuppiHt",      .pfs={"CutFlow","Triggers"},             .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    //sh.AddHistos("all events", { .fill="AK4PuppiHt",      .pfs={"CutFlow","Data,MC","PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt="Sumw2Stack0KeepOrder", .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });

    // Event composition - under different cuts
    for (size_t i=1; i<=13; ++i) {
      std::stringstream cut; cut<<"Pass"<<i<<"Cuts";
      sh.AddHistos("all events", { .fill="SumPt",      .pfs={"PlotSamples",cut.str()},              .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="SumPt",      .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK8Ht",      .pfs={"PlotSamples",cut.str()},              .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK8Ht",      .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK4Ht",      .pfs={"PlotSamples",cut.str()},              .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK4Ht",      .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK8PuppiHt", .pfs={"PlotSamples",cut.str()},              .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK8PuppiHt", .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK4PuppiHt", .pfs={"PlotSamples",cut.str()},              .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="AK4PuppiHt", .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });

      sh.AddHistos("all events", { .fill="MetPt",   .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
      sh.AddHistos("all events", { .fill="MR",      .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
      sh.AddHistos("all events", { .fill="MTR",     .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
      sh.AddHistos("all events", { .fill="R",       .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("all events", { .fill="AK4ChsR", .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("all events", { .fill="AK8ChsR", .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("all events", { .fill="DPhi",    .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("all events", { .fill="DEta",    .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
      sh.AddHistos("all events", { .fill="DR",      .pfs={"PlotSamples",cut.str(),"PassHLTCut"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    }
    sh.AddHistos("all events",   { .fill="SumPt",      .pfs={"PlotSamples","Pass9Cuts","PassHLTCut","PassDPhiCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("all events",   { .fill="AK8Ht",      .pfs={"PlotSamples","Pass9Cuts","PassHLTCut","PassDPhiCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("all events",   { .fill="AK4Ht",      .pfs={"PlotSamples","Pass9Cuts","PassHLTCut","PassDPhiCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("all events",   { .fill="AK8PuppiHt", .pfs={"PlotSamples","Pass9Cuts","PassHLTCut","PassDPhiCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("all events",   { .fill="AK4PuppiHt", .pfs={"PlotSamples","Pass9Cuts","PassHLTCut","PassDPhiCut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });

    // --------------------------------------
    //               Efficiencies

    // 1D plots - sumpt - cutflow
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt", .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht", .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass1Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass3Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass5Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass7Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht", .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });

    // Check other triggers
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT650_vs_SumPt",       .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_PFHT800_vs_SumPt",          .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
  
    // 1 Bin
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT650_vs_SumPtOneBin",       .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0.95,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0.95,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPtOneBin",       .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0.95,1} });
    sh.AddHistos("all events", { .fill="HLTEff_PFHT800_vs_SumPtOneBin",          .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"}, .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_PFHT800_vs_SumPtOneBin",          .pfs={"Triggers","Pass9Cuts"},                   .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });
    sh.AddHistos("all events", { .fill="HLTEff_PFHT800_vs_SumPtOneBin",          .pfs={"Triggers","Pass9Cuts","PassDPhiCut"},     .cuts={}, .draw="PE1", .opt="Sumw2", .ranges={0,0, 0,1, 0.45,0.45} });

    // 2D plots
    // pt1 vs pt2
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Pt_vs_Jet2Pt",     .pfs={"Triggers","Pass1Cuts"},                                 .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,800, 0,800, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Pt_vs_Jet2Pt",     .pfs={"Triggers","Pass3Cuts"},                                 .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,800, 0,800, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Pt_vs_Jet2Pt",     .pfs={"Triggers","Pass5Cuts"},                                 .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,800, 0,800, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Pt_vs_Jet2Pt",     .pfs={"Triggers","Pass5Cuts","PassJetMassCuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,800, 0,800, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Pt_vs_Jet2Pt",     .pfs={"Triggers","Pass5Cuts","PassJetMassCuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,800, 0,800, 0,1} });
    // mass1 vs mass2
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Mass_vs_Jet2Mass", .pfs={"Triggers","Pass1Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 0,200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Mass_vs_Jet2Mass", .pfs={"Triggers","Pass3Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 0,200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Mass_vs_Jet2Mass", .pfs={"Triggers","Pass5Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 0,200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Mass_vs_Jet2Mass", .pfs={"Triggers","Pass7Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 0,200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_Jet1Mass_vs_Jet2Mass", .pfs={"Triggers","Pass7Cuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 0,200, 0,1} });
    // sumpt vs mass1
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass1Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass3Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_SumPt_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    // ak8 ht vs mass
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass1Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass3Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK8Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    // ak4 ht vs mass
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass1Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass3Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts"},               .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });
    sh.AddHistos("all events", { .fill="HLTEff_AK8PFHT700_vs_AK4Ht_vs_Jet1Mass",    .pfs={"Triggers","Pass5Cuts","PassDPhiCut"}, .cuts={}, .draw="COLZ", .opt="Sumw2", .ranges={0,200, 400,1200, 0,1} });

    // HLTjetPt vs AK8jetPt
  
    // Main variables (Shape, area, Data-MC agreement)
    // Hadronic top selection
    //   // No Cut
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",           .pfs={"PlotSamples"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetTau2",         .pfs={"PlotSamples"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetTau3",         .pfs={"PlotSamples"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",        .pfs={"PlotSamples"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} }); // Note
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass", .pfs={"PlotSamples"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note
    //   // Apply 1 Cut
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",           .pfs={"PlotSamples","JetMassCut"},  .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",           .pfs={"PlotSamples","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",        .pfs={"PlotSamples","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",        .pfs={"PlotSamples","JetMassCut"},  .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass", .pfs={"PlotSamples","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass", .pfs={"PlotSamples","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} });
    //   // Apply 2 Cuts (N-1)
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",           .pfs={"PlotSamples","JetMassCut","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",        .pfs={"PlotSamples","JetPtCut","JetMassCut"},    .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });  // Note
    //   sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass", .pfs={"PlotSamples","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,300, 1.01e-2,1e6, 0.55,0.9} }); // Note

    // 2D
    // No Cut
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass_vs_JetTau32",     .pfs={"AllSamples"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass_vs_JetPt",        .pfs={"AllSamples"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetPt_vs_JetTau3",                .pfs={"AllSamples"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });
    // Apply 1 Cut
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass_vs_JetTau32",     .pfs={"AllSamples","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass_vs_JetPt",        .pfs={"AllSamples","JetTau32Cut"}, .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetPt_vs_JetTau32",               .pfs={"AllSamples","JetMassCut"},  .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,0} });  

    // Same plots, but use Gen Particle Truth
    sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // No Cut
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau1",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau2",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau3",                               .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau21",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau31",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    //sh.AddHistos("jetsAK8Puppi", { .fill="JetNSubJets",                           .pfs={"JetGenTruth","Signals,Background"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // 1 Cut
    sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetMassCut"},  .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetPt",                                 .pfs={"JetGenTruth","Signals,Background","JetMassCut","JetTau32Cut"}, .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} }); // 2 Cuts
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32",                              .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau31",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau32_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetTau31_vs_JetTau21",                  .pfs={"JetGenTruth","Signals,Background","JetPtCut","JetMassCut"},    .cuts={},  .draw="COLZ",   .opt="",     .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass",                       .pfs={"JetGenTruth","Signals,Background","JetTau32Cut","JetPtCut"},   .cuts={},  .draw="HISTE1", .opt="Norm", .ranges={0,300, 0,0} });
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass_vs_JetTau32",           .pfs={"JetGenTruth","Signals,Background"},             .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 2D No Cut
    sh.AddHistos("jetsAK8Puppi", { .fill="JetSoftDropMass_vs_JetTau32",           .pfs={"JetGenTruth","Signals,Background","JetPtCut"},    .cuts={},  .draw="COLZ", .opt="Log",     .ranges={0,0, 0,500} }); // 1 Cut

    // Fraction of merged sub-jets
    sh.AddHistos("jetsAK8Puppi", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT"},                                   .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("jetsAK8Puppi", { .fill="MergedTopFraction_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT","JetMatchedGenTopType"},            .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("jetsAK8Puppi", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"Signals,TT"},                                   .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("jetsAK8Puppi", { .fill="MergedTopFraction_vs_JetPtBins",              .pfs={"Signals,TT","JetMatchedGenTopType"},            .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
  
    // Top Tag/Finding Efficiency
    sh.AddHistos("jetsAK8Puppi", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT"},              .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
    sh.AddHistos("jetsAK8Puppi", { .fill="TopTagEfficiency_vs_JetMatchedGenTopPtBins", .pfs={"Signals,TT","TopSizeCut"}, .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
    sh.AddHistos("jetsAK8Puppi", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"Signals,TT"},              .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });
    sh.AddHistos("jetsAK8Puppi", { .fill="TopTagEfficiency_vs_JetPtBins",              .pfs={"Signals,TT","TopSizeCut"}, .cuts={"JetHasMatchedGenTop"}, .draw="HISTE1", .opt="", .ranges={0,0, 0,1} });

    // --------------------------------------------------------------------------
    //                              Gen particles

    // Jet Finding Efficiency
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });  // Note
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="JetFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT"},              .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtBins",    .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });
    sh.AddHistos("gen",     { .fill="TopFindingEfficiency_vs_GenTopPtFewBins", .pfs={"Signals,TT","GenTopType"}, .cuts={"IsGenTop"}, .draw="HISTE1",     .opt="",        .ranges={0,0, 0,1} });

    // --------------------------------------------------------------------------
    //                                   HLT

    // // Trigger Efficiencies
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPt",       .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={800,2000, 0.8,1.05, 0.15,0.95} });
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPtOneBin", .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtBins",        .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,2000,   0.8,1.05, 0.15,0.95} });
    //sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtOneBin",      .pfs={"Data,MC"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });

    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPt",       .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={800,2000, 0.8,1.05, 0.15,0.95} });
    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_TTHadSumPtOneBin", .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });
    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtBins",        .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,2000,   0.8,1.05, 0.15,0.95} });
    // sh.AddHistos("evt",   { .fill="HLTEff_AK8PFHT700_vs_AK8HtOneBin",      .pfs={"Background,Signal"}, .cuts={"AllFilters","LowHTTrig","NHadTopPreTag>=2"}, .draw="PE1", .opt="Sumw2", .ranges={0,0,      0.8,1.05, 0.15,0.95} });

    // --------------------------------------------------------------------------
    //                                MET/Razor

    // MET, Razor Variables
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={"PlotSamples"},                              .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} }); // Note
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={"PlotSamples","Tau32Cuts"},                    .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={"PlotSamples","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={"PlotSamples","Tau32Cuts","RBands"},           .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="MetPt",            .pfs={"PlotSamples","Tau32Cuts","RBands","DPhiBands"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.55,0.9} });
    sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={"PlotSamples"},                              .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={"PlotSamples","Tau32Cuts"},                    .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="TTHadMR",          .pfs={"PlotSamples","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MR",               .pfs={"PlotSamples"},                              .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MR",               .pfs={"PlotSamples","Tau32Cuts"},                    .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MR",               .pfs={"PlotSamples","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MTR",              .pfs={"PlotSamples"},                              .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MTR",              .pfs={"PlotSamples","Tau32Cuts"},                    .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
    sh.AddHistos("baseline events", { .fill="MTR",              .pfs={"PlotSamples","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1e-3,1e6} });
  
    // 2D Correlation plots
    sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples","DPhiBands"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"AllSamples","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="R_vs_DPhi",        .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadR_vs_DPhi",   .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MR_vs_DPhi",       .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="TTHadMR_vs_DPhi",  .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_DPhi",      .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background"},                                  .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background","Tau32Cuts"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background","DPhiBands"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="MTR_vs_MR",        .pfs={"Background","Tau32Cuts","DPhiBands"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgR_vs_DPhi",        .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgTTHadR_vs_DPhi",   .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMR_vs_DPhi",       .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgTTHadMR_vs_DPhi",  .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_DPhi",      .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","DPhiBands","AllSamples"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","AllSamples"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Tau32Cuts","AllSamples"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgR_vs_DPhi",        .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgTTHadR_vs_DPhi",   .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMR_vs_DPhi",       .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgTTHadMR_vs_DPhi",  .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_DPhi",      .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"Tau32Cuts","DPhiBands","Background"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Background"},                      .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events", { .fill="AvgMTR_vs_MR",        .pfs={"DPhiBands","Tau32Cuts","Background"},          .cuts={}, .draw="COLZ", .opt="Log", .ranges={0,0, 0,0, 0,0} });
  
    // Signal selection (Apply loose pretag selection)
    sh.AddHistos("baseline events",   { .fill="R",                 .pfs={"PlotSamples"},                     .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={"PlotSamples"},                     .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="R",                 .pfs={"PlotSamples","Tau32Cuts"},           .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={"PlotSamples","Tau32Cuts"},           .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="R",                 .pfs={"PlotSamples","Tau32Cuts","DPhiBands"}, .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
    sh.AddHistos("baseline events",   { .fill="DPhi",              .pfs={"PlotSamples","Tau32Cuts","RBands"},  .cuts={}, .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e5, 0.15,0.9} });
  
    // 3D Plots to get best signal cuts (Maximize Smin) --> input for B2GAnalyzer
    sh.AddHistos("baseline events",   { .fill="MR_vs_DPhiFine_vs_RFine",           .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events",   { .fill="MR_vs_DPhiFine_vs_TTHadRFine",      .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
    sh.AddHistos("baseline events",   { .fill="TTHadMR_vs_DPhiFine_vs_RFine",      .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} });
    //sh.AddHistos("baseline events",   { .fill="HtAll_vs_DPhiFine_vs_RFine",        .pfs={"AllSamples","Tau32Cuts"}, .cuts={}, .draw="", .opt="", .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  
    // R plots - Distributions for All samples, All cut combinations, Ratios
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    // ABCD regions - DPhi<DPHI_CUT, define regions by: R and Ntop
    sh.AddHistos("baseline events",   { .fill="R",         .pfs={"ABCD,Signals","DPhiBands"},      .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="R",         .pfs={"ABCD","DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="R",         .pfs={"ABCD","DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    //sh.AddHistos("baseline events",   { .fill="R",         .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
    //sh.AddHistos("baseline events",   { .fill="R",         .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
    sh.AddHistos("baseline events",   { .fill="RFine",     .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RFine",     .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","Tau32Cuts","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="RBins",     .pfs={"DPhiBands","Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  
    // DPhi plots - Distributions for All samples, All cut combinations, Ratios
    sh.AddHistos("baseline events",   { .fill="DPhi",      .pfs={"ABCD,Signals","DPhiBands"},      .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands", "AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands", "Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",  "AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",  "Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    // Alternative ABCD regions - Ntop==2, define regions by: R and DPhi
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands",   "Tau32Cuts","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBands",   "Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",    "Tau32Cuts","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    //sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"RBins",    "Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","RBins",    "AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
    sh.AddHistos("baseline events",   { .fill="DPhiBins",  .pfs={"Tau32Cuts","RBins",    "Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });

    // ABCD method + systematics
    //sh.AddHistos("baseline events syst",   { .fill="Counts_vs_DPhiBins",  .pfs={"RBins",    "Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
    //sh.AddHistos("baseline events syst",   { .fill="Counts_vs_R",         .pfs={"Tau32Cuts","DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
    //sh.AddHistos("baseline events syst",   { .fill="Counts_vs_R",         .pfs={"Tau32Cuts","DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
    
  } // End all plots

  // --------------------------------------------------------------------------
  //                             Systematics

  // Data-MC Comparison - background composition
  std::vector<std::string> cuts = { "Pass5Cuts", "Pass9Cuts", "Pass10Cuts", "Pass12Cuts", "Pass13Cuts" };
  for (auto cut : cuts) { 
    sh.AddHistos("all events syst",   { .fill="Counts_vs_NJet",            .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("all events syst",   { .fill="Counts_vs_NHadTopTag",      .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPt",           .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPt",           .pfs={"PlotSamples",cut,"PassHLTCut","PtOrderedJets"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPhi",          .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetPhi",          .pfs={"PlotSamples",cut,"PassHLTCut","PtOrderedJets"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetEta",          .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetEta",          .pfs={"PlotSamples",cut,"PassHLTCut","PtOrderedJets"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetTau32",        .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetTau32",        .pfs={"PlotSamples",cut,"PassHLTCut","PtOrderedJets"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetSoftDropMass", .pfs={"PlotSamples",cut,"PassHLTCut"},                 .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
    sh.AddHistos("jetsAK8Puppi syst", { .fill="Counts_vs_JetSoftDropMass", .pfs={"PlotSamples",cut,"PassHLTCut","PtOrderedJets"}, .cuts={},  .draw="HISTE1", .opt=stack_plot_opt, .ranges={0,0, 1.01e-3,1e7, 0.55,0.9} });
  }

  // Background estimation
  sh.AddHistos("baseline events",      { .fill="HtAll_vs_DPhiFine_vs_RFine",      .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  sh.AddHistos("baseline events",      { .fill="AK8PuppiHt_vs_DPhiFine_vs_RFine", .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} }); // B2GAnalyzer
  sh.AddHistos("baseline events",      { .fill="DPhiBins",           .pfs={"RBins",     "Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_DPhiBins", .pfs={"RBins",     "Tau32Cuts","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 0
  sh.AddHistos("baseline events",      { .fill="R",                  .pfs={"Tau32Cuts", "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_R",        .pfs={"Tau32Cuts", "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 1/3
  sh.AddHistos("baseline events",      { .fill="R",                  .pfs={"Tau32Cuts", "DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_R",        .pfs={"Tau32Cuts", "DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // Method 2
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_RFine",    .pfs={"Tau32Cuts", "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} }); // python ABCD
  // Alternative Razor variables - AK4/8 CHS jets
  sh.AddHistos("baseline events",      { .fill="AK8PuppiHt_vs_DPhiFine_vs_AK4ChsRFine", .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events",      { .fill="AK8PuppiHt_vs_DPhiFine_vs_AK8ChsRFine", .pfs={"AllSamples","Tau32Cuts"},         .cuts={}, .draw="",       .opt="",         .ranges={0,0, 0,0, 0,0} });
  sh.AddHistos("baseline events",      { .fill="AK4ChsR",            .pfs={"Tau32Cuts",   "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_AK4ChsR",  .pfs={"Tau32Cuts",   "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",      { .fill="AK8ChsR",            .pfs={"Tau32Cuts",   "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_AK8ChsR",  .pfs={"Tau32Cuts",   "DPhiBands","AllSamples"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",      { .fill="AK4ChsR",            .pfs={"Tau32Cuts",   "DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_AK4ChsR",  .pfs={"Tau32Cuts",   "DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events",      { .fill="AK8ChsR",            .pfs={"Tau32Cuts",   "DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });
  sh.AddHistos("baseline events syst", { .fill="Counts_vs_AK8ChsR",  .pfs={"Tau32Cuts",   "DPhiBands","Background"}, .cuts={}, .draw="HISTE1", .opt="Sumw2Log", .ranges={0,0, 1e-3,1e5} });

}

//_______________________________________________________
//               Fill Histograms here
double weight0=0, weight1=0, weight2=0;
double count0=0, count1=0, count2=0, count1corr=0, count2corr=0;
void
Analysis::fill_analysis_histos(DataStruct& data, const unsigned int& syst_index, const double& weight)
{
  //h_nvtx->Fill(data.evt.NGoodVtx);
  if (syst_index == 0) {
    if ( _apply_ncut( analysis_cuts.size()-3 ) ) h_abcd->Fill(data.evt.AK8Puppi_R, passTau32Cuts, weight);

    //__________________________________
    //         Fill Smarthistos
    //__________________________________
  
    //while(data.mu.Loop())      sh.Fill("mu");
    //while(data.ele.Loop())     sh.Fill("ele");
    //while(data.jetsAK4.Loop()) sh.Fill("jetsAK4");
    //while(data.gen.Loop()) sh.Fill("gen");
    while(data.jetsAK8Puppi.Loop()) sh.Fill("jetsAK8Puppi");
    sh.Fill("all events");
    if (_apply_ncut(9)) sh.Fill("baseline notrigger");
    if (_apply_ncut(10)) sh.Fill("baseline events");
  }
  sh.Fill("all events syst");
  while(data.jetsAK8Puppi.Loop()) sh.Fill("jetsAK8Puppi syst");
  if (_apply_ncut(10)) sh.Fill("baseline events syst");
}

void
Analysis::load_analysis_histos(std::string inputfile)
{
  sh.Add(inputfile.c_str());
}

void
Analysis::save_analysis_histos(bool draw=0)
{
  if (draw) sh.DrawPlots();
  sh.Write();
}
