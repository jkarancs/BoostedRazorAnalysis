void CalcHTScaleFactors() {
  //TFile *f = TFile::Open("results/Plotter_out_Mar_29_10h20m53_CEST_2016.root"); // No top scale factor
  //TFile *f = TFile::Open("results/Plotter_out_Mar_31_14h17m43_CEST_2016.root"); // After using top scale factors
  //TFile *f = TFile::Open("results/Plotter_out_2016_05_24_15h23m41.root"); // Test run in 76X
  //TFile *f = TFile::Open("results/Plotter_out_2016_05_31_08h48m57_replot.root"); // 76X, after rescaling QCD to match data
  //TFile *f = TFile::Open("results/Plotter_out_2016_06_07_15h57m16_replot.root"); // Same, added Puppi HT
  //TFile *f = TFile::Open("results/Plotter_out_2016_06_13_12h27m53.root");  // 76X Silver JSON, After QCD scale
  TFile *f = TFile::Open("results/Plotter_out_2016_06_21_15h27m59.root"); // 76X Golden JSON, Before QCD scale

  bool pt = false;
  bool ranges = true;
  bool scale_qcd = true;
  float qcd_scale_factor = 0.776458;
  if (pt) {
    // Pt Reweighting

    TCanvas *c = (TCanvas*)f->Get("JetPt/PlotSamples_Pass5Cuts_PassHLT_Ratio");
    c->Draw();
    TH1D* ratio = (TH1D*)((TVirtualPad*)c->cd(2))->GetListOfPrimitives()->At(0);
    if (ranges) {
      TF1* fit1 = new TF1("fit1","pol1", 400, 1000);
      ratio->Fit("fit1","RBQ0");
      fit1->SetLineColor(2);
      fit1->SetLineWidth(1);
      TF1* fit1_up   = (TF1*)fit1->Clone("fit1_up");
      TF1* fit1_down = (TF1*)fit1->Clone("fit1_down");
      fit1_up  ->SetParameter(0,fit1->GetParameter(0)+fit1->GetParError(0));
      fit1_down->SetParameter(0,fit1->GetParameter(0)-fit1->GetParError(0));
      fit1_up  ->SetParameter(1,fit1->GetParameter(1)+fit1->GetParError(1));
      fit1_down->SetParameter(1,fit1->GetParameter(1)-fit1->GetParError(1));
      fit1_up  ->SetLineColor(4); fit1_up  ->Draw("SAME");
      fit1_down->SetLineColor(4); fit1_down->Draw("SAME");
      fit1->Draw("SAME");
      
      TF1* fit2 = new TF1("fit2","pol1", 1000, 2000);
      ratio->Fit("fit2","RBQ0");
      fit2->SetLineColor(2);
      fit2->SetLineWidth(1);
      TF1* fit2_up   = (TF1*)fit2->Clone("fit2_up");
      TF1* fit2_down = (TF1*)fit2->Clone("fit2_down");
      fit2_up  ->SetParameter(0,fit2->GetParameter(0)+fit2->GetParError(0));
      fit2_down->SetParameter(0,fit2->GetParameter(0)-fit2->GetParError(0));
      fit2_up  ->SetParameter(1,fit2->GetParameter(1)+fit2->GetParError(1));
      fit2_down->SetParameter(1,fit2->GetParameter(1)-fit2->GetParError(1));
      fit2_up  ->SetLineColor(4); fit2_up  ->Draw("SAME");
      fit2_down->SetLineColor(4); fit2_down->Draw("SAME");
      fit2->Draw("SAME");
      
      TF1* fit3 = new TF1("fit3","pol0", 2000, 3000);
      ratio->Fit("fit3","RBQ0");
      fit3->SetLineColor(2);
      fit3->SetLineWidth(1);
      TF1* fit3_up   = (TF1*)fit3->Clone("fit3_up");
      TF1* fit3_down = (TF1*)fit3->Clone("fit3_down");
      fit3_up  ->SetParameter(0,fit3->GetParameter(0)+fit3->GetParError(0));
      fit3_down->SetParameter(0,fit3->GetParameter(0)-fit3->GetParError(0));
      fit3_up  ->SetParameter(1,fit3->GetParameter(1)+fit3->GetParError(1));
      fit3_down->SetParameter(1,fit3->GetParameter(1)-fit3->GetParError(1));
      fit3_up  ->SetLineColor(4); fit3_up  ->Draw("SAME");
      fit3_down->SetLineColor(4); fit3_down->Draw("SAME");
      fit3->Draw("SAME");
      
      std::cout<<"Fit result:"<<std::endl;
      std::cout<<" fit [400,1000]:"<<std::endl;
      std::cout<<"  p0: "<<fit1->GetParameter(0)<<" +- "<<fit1->GetParError(0)<<std::endl;
      std::cout<<"  p1: "<<fit1->GetParameter(1)<<" +- "<<fit1->GetParError(1)<<std::endl;
      std::cout<<" fit [1000,2000]:"<<std::endl;
      std::cout<<"  p0: "<<fit2->GetParameter(0)<<" +- "<<fit2->GetParError(0)<<std::endl;
      std::cout<<"  p1: "<<fit2->GetParameter(1)<<" +- "<<fit2->GetParError(1)<<std::endl;
      std::cout<<" fit [2000,3000]:"<<std::endl;
      std::cout<<"  p0: "<<fit3->GetParameter(0)<<" +- "<<fit3->GetParError(0)<<std::endl;
      //std::cout<<"  p1: "<<fit3->GetParameter(1)<<" +- "<<fit3->GetParError(1)<<std::endl;
    } else {
      TF1* fit = new TF1("fit","pol1", 400, 2500);
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
      std::cout<<" fit [400,2500]:"<<std::endl;
      std::cout<<"  p0: "<<fit->GetParameter(0)<<" +- "<<fit->GetParError(0)<<std::endl;
      std::cout<<"  p1: "<<fit->GetParameter(1)<<" +- "<<fit->GetParError(1)<<std::endl;
    }
  } else {
    
    // HT Reweighting
    TCanvas *c = (TCanvas*)f->Get("AK8PuppiHt/PlotSamples_Pass5Cuts_PassHLT_Ratio");
    c->Draw();
    TH1D* ratio = (TH1D*)((TVirtualPad*)c->GetListOfPrimitives()->At(1))->GetListOfPrimitives()->At(0);
    if (scale_qcd) {
      //THStack* ratio = (TH1D*)((TVirtualPad*)c->GetListOfPrimitives()->At(1))->GetListOfPrimitives()->At(0);
      TH1D* data = (TH1D*)((TVirtualPad*)c->GetListOfPrimitives()->At(0))->GetListOfPrimitives()->At(1);
      TList* vh = ((THStack*)((TVirtualPad*)c->GetListOfPrimitives()->At(0))->GetListOfPrimitives()->At(2))->GetHists();
      TH1D* mc;
      for (int i=0; i<vh->GetEntries(); ++i) {
	TH1D* h = (TH1D*)((TH1D*)vh->At(i))->Clone();
	if (TString(h->GetName()).Contains("QCD")) h->Scale(qcd_scale_factor);
	if (i==0) mc = (TH1D*)h->Clone("mc");
	else mc->Add(h);
      }
      TH1D *ratio_scaled = (TH1D*)data->Clone();
      ratio_scaled->Divide(mc);
      ratio_scaled->GetYaxis()->SetTitle("Data/MC");
      TVirtualPad *pad = (TVirtualPad*)c->GetListOfPrimitives()->At(1);
      pad->cd();
      for (int bin=1; bin<=ratio->GetNbinsX(); ++bin) {
	ratio->SetBinContent(bin, ratio_scaled->GetBinContent(bin));
	ratio->SetBinError(bin, ratio_scaled->GetBinError(bin));
      }
      gPad->Update();
    }
    
    // QCD dominant region
    if (ranges) {
      TF1* fit1 = new TF1("fit1","pol1", 800, 2000);
      ratio->Fit("fit1","RBQ0");
      fit1->SetLineColor(2);
      fit1->SetLineWidth(1);
      TF1* fit1_up   = (TF1*)fit1->Clone("fit1_up");
      TF1* fit1_down = (TF1*)fit1->Clone("fit1_down");
      fit1_up  ->SetParameter(0,fit1->GetParameter(0)+fit1->GetParError(0));
      fit1_down->SetParameter(0,fit1->GetParameter(0)-fit1->GetParError(0));
      fit1_up  ->SetParameter(1,fit1->GetParameter(1)+fit1->GetParError(1));
      fit1_down->SetParameter(1,fit1->GetParameter(1)-fit1->GetParError(1));
      fit1_up  ->SetLineColor(4); fit1_up  ->Draw("SAME");
      fit1_down->SetLineColor(4); fit1_down->Draw("SAME");
      fit1->Draw("SAME");
      
      TF1* fit2 = new TF1("fit2","pol1", 2000, 6000);
      ratio->Fit("fit2","RBQ0");
      fit2->SetLineColor(2);
      fit2->SetLineWidth(1);
      TF1* fit2_up   = (TF1*)fit2->Clone("fit2_up");
      TF1* fit2_down = (TF1*)fit2->Clone("fit2_down");
      fit2_up  ->SetParameter(0,fit2->GetParameter(0)+fit2->GetParError(0));
      fit2_down->SetParameter(0,fit2->GetParameter(0)-fit2->GetParError(0));
      fit2_up  ->SetParameter(1,fit2->GetParameter(1)+fit2->GetParError(1));
      fit2_down->SetParameter(1,fit2->GetParameter(1)-fit2->GetParError(1));
      fit2_up  ->SetLineColor(4); fit2_up  ->Draw("SAME");
      fit2_down->SetLineColor(4); fit2_down->Draw("SAME");
      fit2->Draw("SAME");
      
      std::cout<<"Fit result:"<<std::endl;
      std::cout<<" fit [800,2000] Chi2/Ndof: "<<fit1->GetChisquare()/fit1->GetNDF()<<std::endl;
      std::cout<<"  p0: "<<fit1->GetParameter(0)<<" +- "<<fit1->GetParError(0)<<std::endl;
      std::cout<<"  p1: "<<fit1->GetParameter(1)<<" +- "<<fit1->GetParError(1)<<std::endl;
      std::cout<<" fit [2000,6000] Chi2/Ndof: "<<fit2->GetChisquare()/fit2->GetNDF()<<std::endl;
      std::cout<<"  p0: "<<fit2->GetParameter(0)<<" +- "<<fit2->GetParError(0)<<std::endl;
      std::cout<<"  p1: "<<fit2->GetParameter(1)<<" +- "<<fit2->GetParError(1)<<std::endl;
      std::cout<<std::endl;
      std::cout<<"const double p0[2]     = { "<<fit1->GetParameter(0)<<", "<<fit2->GetParameter(0)<<" };"<<std::endl;
      std::cout<<"const double p0_err[2] = { "<<fit1->GetParError(0) <<", "<<fit2->GetParError(0) <<" };"<<std::endl;
      std::cout<<"const double p1[2]     = { "<<fit1->GetParameter(1)<<", "<<fit2->GetParameter(1)<<" };"<<std::endl;
      std::cout<<"const double p1_err[2] = { "<<fit1->GetParError(1) <<", "<<fit2->GetParError(1) <<" };"<<std::endl;
	
    } else {

      TF1* fit = new TF1("fit","pol1", 800, 6000);
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
      std::cout<<" fit [800,6000]:"<<std::endl;
      std::cout<<"  p0: "<<fit->GetParameter(0)<<" +- "<<fit->GetParError(0)<<std::endl;
      std::cout<<"  p1: "<<fit->GetParameter(1)<<" +- "<<fit->GetParError(1)<<std::endl;
      std::cout<<std::endl;
      std::cout<<"const double p0     = "<<fit->GetParameter(0)<<";"<<std::endl;
      std::cout<<"const double p0_err = "<<fit->GetParError(0) <<";"<<std::endl;
      std::cout<<"const double p1     = "<<fit->GetParameter(1)<<";"<<std::endl;
      std::cout<<"const double p1_err = "<<fit->GetParError(1) <<";"<<std::endl;
    }
  }
  gPad->Update();
  f->Close();
}
