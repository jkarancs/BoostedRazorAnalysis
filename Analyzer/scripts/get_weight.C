void get_weight(std::string files) {
  TChain c("c");
  c.Add(files.c_str());
  TH1D *h = new TH1D("h","",2,-1,1);
  float NEvent_Corr, XSec, Lumi_Weight, Gen_Weight, prev_Gen_Weight;
  bool ok = 1;
  for (size_t i=0; i<c.GetListOfFiles()->GetEntries(); ++i) {
    TFile *f = TFile::Open(c.GetListOfFiles()->At(i)->GetTitle());
    h->Add((TH1D*)f->Get("EventCounter/NEventNoFilter"));
    TTree *tree = (TTree*)f->Get("B2GTTreeMaker/B2GTree");
    tree->GetBranch("evt_NEvent_Corr")->SetAddress(&NEvent_Corr);
    tree->GetBranch("evt_XSec")->SetAddress(&XSec);
    tree->GetBranch("evt_Lumi_Weight")->SetAddress(&Lumi_Weight);
    tree->GetBranch("evt_Gen_Weight")->SetAddress(&Gen_Weight);
    tree->GetEntry(0);
    // Check that abs(evt_Gen_Weight) is constant
    prev_Gen_Weight = fabs(Gen_Weight);
    for (int i=0; i<10; i++) {
      tree->GetEntry(i);
      if (fabs(Gen_Weight)!=prev_Gen_Weight) { ok=0; break; }
    }
  }
  double NEvents = h->GetEntries();
  double NegFrac = h->GetBinContent(1)/h->GetEntries();
  double SumWeights = h->GetBinContent(2) - h->GetBinContent(1);
  double Lumi_Weight_calc = 1000 * XSec / SumWeights;
  printf("%-50s   NEvents: %10.0lf   NegFrac: %1.4lf   Corr: %1.6lf   SumWeights: %10.0lf   NEventCorr(file): %10.0lf   XSec(file): %lf   LumiWeight(file): %lf   LumiWeight(calc): %lf\n", files.c_str(), NEvents, NegFrac, NEvents/SumWeights, SumWeights, NEvent_Corr, XSec, Lumi_Weight, Lumi_Weight_calc);
  if (!ok) std::cout<<"!!! Warning fabs(evt_gen_Weight) is not a constant !!!"<<std::endl;
  gApplication->Terminate();
}
