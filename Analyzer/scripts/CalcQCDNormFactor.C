void CalcQCDNormFactor() {
  //TFile *f = TFile::Open("results/Plotter_out_2016_05_29_22h19m32.root"); // 76X Silver JSON
  TFile *f = TFile::Open("results/Plotter_out_2016_06_21_15h27m59.root"); // 76X Golden JSON
  TCanvas *c = (TCanvas*)f->Get("NJet/PlotSamples_Pass5Cuts_PassHLT");
  THStack *s = (THStack*)c->GetListOfPrimitives()->At(1);
  TH1D *data = (TH1D*)c->GetListOfPrimitives()->At(3);
  double MC_integral = 0;
  double QCD_count = 0;
  for (int i=s->GetHists()->GetEntries()-1; i>=0; --i) {
    TH1D* h = (TH1D*)s->GetHists()->At(i);
    if (i==s->GetHists()->GetEntries()-1) QCD_count = h->Integral();
    std::cout<<h->GetName()<<" "<<h->Integral()<<std::endl;
    MC_integral += h->Integral();
  }
  double data_QCD_estimate = data->Integral()- (MC_integral-QCD_count);
  double QCD_scale = data_QCD_estimate/QCD_count;
  std::cout<<"MC:                  "<<MC_integral<<std::endl;
  std::cout<<"Data:                "<<data->Integral()<<std::endl;
  std::cout<<"MC   (QCD only):     "<<QCD_count<<std::endl;
  std::cout<<"Data (QCD only est): "<<data_QCD_estimate<<std::endl;
  std::cout<<"QCD Scale: "<<QCD_scale<<std::endl;

  TH1D* qcd = (TH1D*)s->GetHists()->At(s->GetHists()->GetEntries()-1);
  qcd->Scale(QCD_scale);
  c->Draw();
  
  //TCanvas *c = (TCanvas*)f->Get("NJet/PlotSamples_Pass5Cuts_PassHLT_Ratio");
  //
  //TH1D* ratio = (TH1D*)((TVirtualPad*)c->cd(2))->GetListOfPrimitives()->At(0);
  //TF1* fit = new TF1("fit","pol1", 400, 2000);
  //ratio->Fit("fit","RBQ0");
  //fit->SetLineColor(2);
  //fit->SetLineWidth(1);
  //TF1* fit_up   = (TF1*)fit->Clone("fit_up");
  //TF1* fit_down = (TF1*)fit->Clone("fit_down");
  //fit_up  ->SetParameter(0,fit->GetParameter(0)+fit->GetParError(0)*1);
  //fit_down->SetParameter(0,fit->GetParameter(0)-fit->GetParError(0)*1);
  //fit_up  ->SetParameter(1,fit->GetParameter(1)+fit->GetParError(1)*1);
  //fit_down->SetParameter(1,fit->GetParameter(1)-fit->GetParError(1)*1);
  //fit_up  ->SetLineColor(4); fit_up  ->Draw("SAME");
  //fit_down->SetLineColor(4); fit_down->Draw("SAME");
  //fit->Draw("SAME");
  //
  //std::cout<<"Fit result:"<<std::endl;
  //std::cout<<"p0: "<<fit->GetParameter(0)<<" +- "<<fit->GetParError(0)*1<<std::endl;
  //std::cout<<"p1: "<<fit->GetParameter(1)<<" +- "<<fit->GetParError(1)*1<<std::endl;
  //f->Close();  
}
