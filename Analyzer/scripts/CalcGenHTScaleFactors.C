void CalcGenHTScaleFactors() {

  // TTJets Madgraph-pythia8, binned/unbinned correction
  TChain *unbinned = new TChain("TTJets_unbinned");
  unbinned->Add("ntuple/grid18/Apr13_edm_Apr01/TTJets_madgraph/*.root/B2GTTreeMaker/B2GTree");
  TChain *binned_1 = new TChain("TTJets_binned_1");
  TChain *binned_2 = new TChain("TTJets_binned_2");
  TChain *binned_3 = new TChain("TTJets_binned_3");
  TChain *binned_4 = new TChain("TTJets_binned_4");
  TChain *binned_5 = new TChain("TTJets_binned_5");
  binned_1->Add("ntuple/grid18/Apr13_edm_Apr01/TTJets_HT-0to600/*.root/B2GTTreeMaker/B2GTree");
  binned_2->Add("ntuple/grid18/Apr13_edm_Apr01/TTJets_HT-600to800/*.root/B2GTTreeMaker/B2GTree");
  binned_3->Add("ntuple/grid18/Apr13_edm_Apr01/TTJets_HT-800to1200/*.root/B2GTTreeMaker/B2GTree");
  binned_4->Add("ntuple/grid18/Apr13_edm_Apr01/TTJets_HT-1200to2500/*.root/B2GTTreeMaker/B2GTree");
  binned_5->Add("ntuple/grid18/Apr13_edm_Apr01/TTJets_HT-2500toInf/*.root/B2GTTreeMaker/B2GTree");


  TH1::SetDefaultSumw2();
  Double_t bins[19] = { 0, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000 };
  TH1D* h  = new TH1D("ttjets_unbinned","HT^{Gen} - TTJets madgraphMLM-pythia8 unbinned (GeV)", 18,bins);
  TH1D* h1 = new TH1D("ttjets_binned_1","HT^{Gen} - TTJets madgraphMLM-pythia8 binned (GeV)",   18,bins);
  TH1D* h2 = new TH1D("ttjets_binned_2","HT^{Gen} - TTJets madgraphMLM-pythia8 binned (GeV)",   18,bins);
  TH1D* h3 = new TH1D("ttjets_binned_3","HT^{Gen} - TTJets madgraphMLM-pythia8 binned (GeV)",   18,bins);
  TH1D* h4 = new TH1D("ttjets_binned_4","HT^{Gen} - TTJets madgraphMLM-pythia8 binned (GeV)",   18,bins);
  TH1D* h5 = new TH1D("ttjets_binned_5","HT^{Gen} - TTJets madgraphMLM-pythia8 binned (GeV)",   18,bins);
  h1->SetLineColor(2);
  h2->SetLineColor(3);
  h3->SetLineColor(4);
  h4->SetLineColor(5);
  h5->SetLineColor(6);

  TCanvas *c = new TCanvas("ttjets","", 1200, 600);
  c->Divide(2);
  c->cd(1);
  unbinned->Draw("evt_Gen_Ht>>ttjets_unbinned",  "evt_Gen_Weight*831.76/1.02151e+07",  "HIST");
  binned_1->Draw("evt_Gen_Ht>>ttjets_binned_1",  "evt_Gen_Weight*827.793/9.34423e+06", "HIST SAME");
  binned_2->Draw("evt_Gen_Ht>>ttjets_binned_2",  "evt_Gen_Weight*2.66653/5.16352e+06", "HIST SAME");
  binned_3->Draw("evt_Gen_Ht>>ttjets_binned_3",  "evt_Gen_Weight*1.09808/3.55392e+06", "HIST SAME");
  binned_4->Draw("evt_Gen_Ht>>ttjets_binned_4",  "evt_Gen_Weight*0.198748/992490",     "HIST SAME");
  binned_5->Draw("evt_Gen_Ht>>ttjets_binned_5",  "evt_Gen_Weight*0.002368/517274",     "HIST SAME");
  unbinned->Draw("HIST");
  h1->Draw("HIST SAME");
  h2->Draw("HIST SAME");
  h3->Draw("HIST SAME");
  h4->Draw("HIST SAME");
  h5->Draw("HIST SAME");
  c->cd(2);
  TH1D* r = (TH1D*)h->Clone("ttjets_ratio");
  TH1D* den = (TH1D*)h1->Clone("ttjets_den");
  den->Add(h2);
  den->Add(h3);
  den->Add(h4);
  den->Add(h5);
  r->Divide(den);
  r->SetMarkerStyle(20);
  r->Draw("HIST");



  //    // WJets Madgraph-pythia8, binned/unbinned correction
  //    unbinned = new TChain("WJets_unbinned");
  //    unbinned->Add("ntuple/grid18/Jun08/WJetsToLNu/*.root/B2GTTreeMaker/B2GTree");
  //    binned_1 = new TChain("WJets_binned_1");
  //    binned_2 = new TChain("WJets_binned_2");
  //    binned_3 = new TChain("WJets_binned_3");
  //    binned_4 = new TChain("WJets_binned_4");
  //    binned_5 = new TChain("WJets_binned_5");
  //    TChain *binned_6 = new TChain("WJets_binned_5");
  //    TChain *binned_7 = new TChain("WJets_binned_5");
  //    binned_1->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-100To200/*.root/B2GTTreeMaker/B2GTree");
  //    binned_2->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-200To400/*.root/B2GTTreeMaker/B2GTree");
  //    binned_3->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-400To600/*.root/B2GTTreeMaker/B2GTree");
  //    binned_4->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-600To800/*.root/B2GTTreeMaker/B2GTree");
  //    binned_5->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-800To1200/*.root/B2GTTreeMaker/B2GTree");
  //    binned_6->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-1200To2500/*.root/B2GTTreeMaker/B2GTree");
  //    binned_7->Add("ntuple/grid18/Apr13_edm_Apr01/WJetsToLNu_HT-2500ToInf/*.root/B2GTTreeMaker/B2GTree");
  //    
  //    TH1::SetDefaultSumw2();
  //    Double_t bins2[19] = { 0, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000 };
  //    h  = new TH1D("wjets_unbinned","HT^{Gen} - WJets madgraphMLM-pythia8 unbinned (GeV)", 18,bins2);
  //    h1 = new TH1D("wjets_binned_1","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    h2 = new TH1D("wjets_binned_2","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    h3 = new TH1D("wjets_binned_3","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    h4 = new TH1D("wjets_binned_4","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    h5 = new TH1D("wjets_binned_5","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    TH1D* h6 = new TH1D("wjets_binned_6","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    TH1D* h7 = new TH1D("wjets_binned_7","HT^{Gen} - WJets madgraphMLM-pythia8 binned (GeV)",   18,bins2);
  //    h1->SetLineColor(2);
  //    h2->SetLineColor(3);
  //    h3->SetLineColor(4);
  //    h4->SetLineColor(5);
  //    h5->SetLineColor(6);
  //    h6->SetLineColor(7);
  //    h7->SetLineColor(8);
  //    
  //    TCanvas *c2 = new TCanvas("wjets","", 1200, 600);
  //    c2->Divide(2);
  //    c2->cd(1);
  //    unbinned->Draw("evt_Gen_Ht>>wjets_unbinned",  "evt_Gen_Weight*61526.5/4.71613e+07","HIST");
  //    binned_1->Draw("evt_Gen_Ht>>wjets_binned_1",  "evt_Gen_Weight*1632.53/1.02054e+07","HIST SAME");
  //    binned_2->Draw("evt_Gen_Ht>>wjets_binned_2",  "evt_Gen_Weight*436.597/4.94957e+06","HIST SAME");
  //    binned_3->Draw("evt_Gen_Ht>>wjets_binned_3",  "evt_Gen_Weight*59.366/1.94366e+06", "HIST SAME");
  //    binned_4->Draw("evt_Gen_Ht>>wjets_binned_4",  "evt_Gen_Weight*14.626/3.76777e+06", "HIST SAME");
  //    binned_5->Draw("evt_Gen_Ht>>wjets_binned_5",  "evt_Gen_Weight*6.677/1.56828e+06",  "HIST SAME");
  //    binned_6->Draw("evt_Gen_Ht>>wjets_binned_6",  "evt_Gen_Weight*1.61311/246239",     "HIST SAME");
  //    binned_7->Draw("evt_Gen_Ht>>wjets_binned_7",  "evt_Gen_Weight*0.039035/251982",    "HIST SAME");
  //    unbinned->Draw("HIST");
  //    h1->Draw("HIST SAME");
  //    h2->Draw("HIST SAME");
  //    h3->Draw("HIST SAME");
  //    h4->Draw("HIST SAME");
  //    h5->Draw("HIST SAME");
  //    h6->Draw("HIST SAME");
  //    h7->Draw("HIST SAME");
  //    c2->cd(2);
  //    TH1D* r2 = (TH1D*)h->Clone("wjets_ratio");
  //    TH1D* den2 = (TH1D*)h1->Clone("wjets_den");
  //    den2->Add(h2);
  //    den2->Add(h3);
  //    den2->Add(h4);
  //    den2->Add(h5);
  //    den2->Add(h6);
  //    den2->Add(h7);
  //    r2->Divide(den);
  //    r2->SetMarkerStyle(20);
  //    r2->Draw("HIST");

}
