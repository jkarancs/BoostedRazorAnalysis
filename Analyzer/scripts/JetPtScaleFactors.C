void JetPtScaleFactors() {
  //TFile *f = TFile::Open("results/Plotter_out_Mar_29_10h20m53_CEST_2016.root"); // No top scale factor
  TFile *f = TFile::Open("results/Plotter_out_Mar_31_14h17m43_CEST_2016.root"); // After using top scale factors
  TCanvas *c = (TCanvas*)f->Get("JetPt/PlotSamples_Pass5Cuts_PassHLT_Ratio");
  c->Draw();
  TH1D* ratio = (TH1D*)((TVirtualPad*)c->cd(2))->GetListOfPrimitives()->At(0);
  TF1* fit = new TF1("fit","pol1", 400, 2000);
  ratio->Fit("fit","RBQ0");
  fit->SetLineColor(2);
  fit->SetLineWidth(1);
  TF1* fit_up   = (TF1*)fit->Clone("fit_up");
  TF1* fit_down = (TF1*)fit->Clone("fit_down");
  fit_up  ->SetParameter(0,fit->GetParameter(0)+fit->GetParError(0));
  fit_down->SetParameter(0,fit->GetParameter(0)-fit->GetParError(0));
  fit_up  ->SetParameter(1,fit->GetParameter(1)+fit->GetParError(1));
  fit_down->SetParameter(1,fit->GetParameter(1)-fit->GetParError(1));
  fit_up  ->SetLineColor(4); fit_up  ->Draw("SAME");
  fit_down->SetLineColor(4); fit_down->Draw("SAME");
  fit->Draw("SAME");

  std::cout<<"Fit result:"<<std::endl;
  std::cout<<"p0: "<<fit->GetParameter(0)<<" +- "<<fit->GetParError(0)<<std::endl;
  std::cout<<"p1: "<<fit->GetParameter(1)<<" +- "<<fit->GetParError(1)<<std::endl;
  f->Close();  
}
