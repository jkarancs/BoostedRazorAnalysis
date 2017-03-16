void get_totweight(std::string files) {
  TChain *c = new TChain("c");
  c->Add(files.c_str());
  double total = 0;
  for (size_t i=0, n=c->GetListOfFiles()->GetEntries(); i<n; ++i) {
    TFile *f = TFile::Open(c->GetListOfFiles()->At(i)->GetTitle());
    TH1D* h = (TH1D*)f->Get("EventCounter/totweight");
    total += h->GetBinContent(1);
    f->Close();
  }
  std::cout<<"totweight: "<<std::setprecision(10)<<total<<std::endl;
  gApplication->Terminate();
}
