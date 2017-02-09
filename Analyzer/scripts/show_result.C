void show_result(const char* file="Latest_plots.root") {
  TFile *f = TFile::Open(file);
  gStyle->SetOptStat(0);
  new TBrowser("Browser", "Plotter output", 1200, 800);
}
