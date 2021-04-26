#include<PostAnalysis/Plotting.h>

void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style)
{
  std::vector<Color_t> Colors = {kRed, kGreen, kBlue, kCyan, kViolet, kOrange, kCoffee, kAurora}; 
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = Hists[i];
    H -> SetLineColor(Colors[i]);
    H -> SetLineStyle(Style);
    H -> SetLineWidth(1);
    H -> Draw("SAMEHIST");
    len -> AddEntry(H, H -> GetTitle());
    len -> Draw("SAME");
    can -> Update(); 
  }
}

void PlotHists(TH1F* Hist, TCanvas* can)
{
  can -> SetLogy(); 
  Hist -> Draw("SAMEHIST"); 
  can -> Update(); 
}

void PlotHists(std::vector<TH1F*> Hists, std::vector<TString> Legend_Titles, TCanvas* can)
{
  gStyle -> SetOptStat(0); 
  for (int i(0); i < Hists.size(); i++){Hists[i] -> SetTitle(Legend_Titles[i]); }
  
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  can -> cd(1); 
  Populate(Hists, can, len, kSolid);
}

void PlotHists(std::vector<TH1F*> Hists, TCanvas* can)
{
  float m; 
  for (TH1F* H : Hists)
  {
    int b = H -> GetMaximumBin(); 
    float m_n = H -> GetBinContent(b+1); 
    if (m_n > m){ m = m_n; }
  }
  
  can -> Clear();
  Hists[0] -> Draw("HIST");
  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  Hists[0] -> GetYaxis() -> SetRangeUser(1, m*2); 
  Populate(Hists, can, len, kSolid); 
  can -> Update();
}

void PlotHists(TH1F* Data, std::vector<TH1F*> Hists, TCanvas* can)
{
  int bin = Data -> GetMaximumBin(); 
  float m = Data -> GetBinContent(bin+1); 
  
  can -> Clear();
  gStyle -> SetOptStat(0); 
  Data -> GetYaxis() -> SetRangeUser(1, m);  
  Data -> GetXaxis() -> SetRangeUser(0, 10);
  Data -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(Hists, can, len, kSolid); 
}

void PlotHists(TH1F* Data, std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TCanvas* can)
{
  int bin = Data -> GetMaximumBin(); 
  float m = Data -> GetBinContent(bin+1); 
  
  can -> Clear();
  gStyle -> SetOptStat(0); 
  Data -> GetYaxis() -> SetRangeUser(1, m*2);
  Data -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(truth, can, len, kSolid); 
  Populate(prediction, can, len, kDashed); 
}

void PlotHists(std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TCanvas* can)
{
  float sum = 0; 
  for (int i(0); i < prediction.size(); i++)
  {
    int bin = prediction[i] -> GetMaximumBin(); 
    float m = prediction[i] -> GetBinContent(bin+1); 
    sum += m;  
  }

  can -> Clear();
  gStyle -> SetOptStat(0); 
  truth[0] -> GetYaxis() -> SetRangeUser(1, sum*2);
  truth[0] -> GetXaxis() -> SetRangeUser(0, 10);
  truth[0] -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(truth, can, len, kDashed); 
  Populate(prediction, can, len, kSolid); 

}

void PlotHists(std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TString title, TCanvas* can)
{
  float sum = 0; 
  for (int i(0); i < prediction.size(); i++)
  {
    int bin = prediction[i] -> GetMaximumBin(); 
    float m = prediction[i] -> GetBinContent(bin+1); 
    sum += m;  
  }

  can -> Clear();
  gStyle -> SetOptStat(0); 
  TH1F* empty = (TH1F*)truth[0] -> Clone(title); 
  empty -> SetTitle(title); 
  empty -> Reset();

  empty -> GetYaxis() -> SetRangeUser(1, sum*2);
  empty -> GetXaxis() -> SetRangeUser(0, 10);
  empty -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(truth, can, len, kDashed); 
  Populate(prediction, can, len, kSolid); 
}

void PlotHists(TH1F* Data, std::vector<TH1F*> Prediction, std::vector<TH1F*> Truth, TString Title, float FLost_P, float FLost_T, TCanvas* can)
{
  can -> SetLogy(); 
  Data -> SetTitle(Title);
 	Data -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^2]");
  if (Truth.size() == 0){PlotHists(Data, Prediction, can);}
  else{PlotHists(Data, Prediction, Truth, can);}
  
  TPaveText* pt = new TPaveText(.7, .5, 0.9, .7, "NDCNB"); 
  pt -> AddText("#int TrackN:");
  pt -> SetFillStyle(4000);
  for (int i(0); i < Prediction.size(); i++)
  {
    float x = Prediction[i] -> Integral(); 
    std::stringstream stream; 
    stream << std::fixed << std::setprecision(2) << x; 
    TString Line = "Track-"; Line += (i+1); Line += (": "); Line += (stream.str()); 
    pt -> AddText(Line); 
  }
  std::stringstream stream; 
  stream << std::fixed << std::setprecision(4) << FLost_P; 
  TString l = "FLost: "; l += (stream.str()); 
  pt -> AddText(l); 
  pt -> Draw();

  if(FLost_T != 0)
  {
    TPaveText* t = new TPaveText(.5, .5, 0.7, .7, "NDCNB"); 
    t -> AddText("#int TrackN:");
    t -> SetFillStyle(4000);
    for (int i(0); i < Truth.size(); i++)
    {
      float x = Truth[i] -> Integral(); 
      std::stringstream stream; 
      stream << std::fixed << std::setprecision(2) << x; 
      TString Line = "Track-"; Line += (i+1); Line += (": "); Line += (stream.str()); 
      t -> AddText(Line); 
    }
    std::stringstream stream; 
    stream << std::fixed << std::setprecision(4) << FLost_T; 
    TString l = "FLost T: "; l += (stream.str()); 
    t -> AddText(l); 
    t -> Draw();
  }
  can -> Update(); 

}

void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can)
{
  can -> Clear(); 
  
  TString name = "Ratio Plot: "; name += (H1 -> GetTitle()); name += ("/"); name += H2 -> GetTitle(); 
  TH1F* Ratio = (TH1F*)H1 -> Clone(name);
  Ratio -> SetTitle(name);
  Ratio -> Divide(H2); 
  
  for (int i(0); i < Ratio -> GetNbinsX(); i++)
  {
    if (std::isnan(Ratio -> GetBinContent(i+1)) || std::isnan(H1 -> GetBinContent(i+1)))
    {
      Ratio -> SetBinContent(i+1, 1e-10);
      H1 -> SetBinContent(i+1, 1e-10); 
    }
  }
  
  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
	H1 -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^2]");

  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1.0);
  P1 -> Draw(); 
  P1 -> cd(); 
  H1 -> SetLineColor(kBlack); 
  H1 -> GetYaxis() -> SetRangeUser(1, 2*GetMaxValue(H1));
  H1 -> Draw("HIST"); 
  len -> AddEntry(H1, H1 -> GetTitle() ); 
  len -> Draw("SAME");
  H1 -> SetStats(0);  
  H2 -> SetLineColor(kRed);
  H2 -> Draw("SAMEHIST"); 
  len -> AddEntry(H2, H2 -> GetTitle() );
  len -> Draw("SAME");  
  P1 -> SetLogy();  
  can -> cd();

  TPad *P2 = new TPad("P2", "P2", 0.0, 0.05, 1, 0.3);  
  P2 -> Draw();
  P2 -> cd(); 
	Ratio -> SetStats(0); 
  Ratio -> SetLineColor(kWhite);
  Ratio -> GetYaxis() -> SetRangeUser(0.5, 1.5);
  Ratio -> SetMarkerSize(1); 
  Ratio -> Draw("epl"); 
}

void RooFitPullPlot(RooAddPdf model, RooRealVar* Domain, std::vector<RooFFTConvPdf*> PDFs, RooDataHist* Data, TString Name)
{
  RooPlot* xframe = Domain -> frame(RooFit::Title("Figure")); 
  Data -> plotOn(xframe); 
  model.plotOn(xframe);  
  xframe -> chiSquare();

  RooHist* resid = xframe -> residHist(); 
  RooHist* pull = xframe -> pullHist(); 
  

  RooPlot* frame2 = Domain -> frame(RooFit::Title("Residual")); 
  frame2 -> addPlotable(resid, "P"); 

  RooPlot* frame3 = Domain -> frame(RooFit::Title("Pull")); 
  frame3 -> addPlotable(pull, "P"); 
  
  TCanvas* can = new TCanvas(); 
  
  can -> Divide(3); 
  can -> cd(1); xframe -> Draw(); 
  can -> cd(2); frame2 -> Draw(); 
  can -> cd(3); frame3 -> Draw(); 
  can -> Update();
  can -> Print(Name); 
  delete can; 
}

void PlotLikelyHood(RooAbsReal* nll, RooRealVar* var, TString name)
{
  // Define the frame 
  RooPlot* frame1 = var -> frame(RooFit::Bins(1000), RooFit::Range(-1, 1), RooFit::Title(name)); 
  nll -> plotOn(frame1, RooFit::ShiftToZero());

  RooAbsReal* pll_var = nll -> createProfile(*var); 
  pll_var -> plotOn(frame1, RooFit::LineColor(kRed)); 
  
  frame1 -> SetMinimum(-5); 
  frame1 -> SetMaximum(5); 
  TCanvas* can = new TCanvas();   
  frame1 -> Draw(); 
  can -> Update(); 
  can -> Print(name); 
  delete can;
}


void RooFitPullPlot(RooAddPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data, TString Name)
{
  RooPlot* xframe = Domain -> frame(RooFit::Title("Figure")); 
  Data -> plotOn(xframe); 
  model.plotOn(xframe);  
  xframe -> chiSquare();

  RooHist* resid = xframe -> residHist(); 
  RooHist* pull = xframe -> pullHist(); 

  RooPlot* frame2 = Domain -> frame(RooFit::Title("Residual")); 
  frame2 -> addPlotable(resid, "P"); 

  RooPlot* frame3 = Domain -> frame(RooFit::Title("Pull")); 
  frame3 -> addPlotable(pull, "P"); 
  
  TCanvas* can = new TCanvas(); 
  
  can -> Divide(3); 
  can -> cd(1); xframe -> Draw(); 
  can -> cd(2); frame2 -> Draw(); 
  can -> cd(3); frame3 -> Draw(); 
  can -> Update();
  can -> Print(Name); 
  delete can; 
  delete xframe; 
  delete frame2; 
  delete frame3; 
}

void PlotRooFit(RooAddPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data)
{
  RooPlot* xframe = Domain -> frame(RooFit::Title("Figure")); 
  Data -> plotOn(xframe); 
  
  Data -> plotOn(xframe, RooFit::Name("Data")); 
  std::vector<Color_t> Colors = {kRed, kGreen, kOrange, kBlue}; 
  for (int i(0); i < PDFs.size(); i++)
  {
    TString name = "trk-"; name += (i+1); 
    model.plotOn(xframe, RooFit::Name(name), RooFit::Components(*PDFs[i]), RooFit::LineStyle(kDotted), RooFit::LineColor(Colors[i])); 
  }
  TCanvas* can = new TCanvas(); 
  gPad -> SetLogy(); 
  xframe -> SetMinimum(1e-6); 
  xframe -> Draw(); 

  can -> Update();
  can -> Print("debug.pdf"); 
  delete can; 
}

// ====================== Proper histogram plotting for the presentation ======================== //
void GeneratePlot(TH1F* H, TString Title, TCanvas* can, Color_t color, ELineStyle style, TString DrawOption, float Intensity)
{
  gStyle -> SetOptStat(0); 
  if (Title != "")
  {
    H -> SetTitle(Title); 
  }
  H -> SetLineColor(color); 
  H -> SetLineColorAlpha(color, Intensity); 
  H -> SetLineStyle(style); 
  H -> SetLineWidth(1); 
  H -> Draw(DrawOption);
  can -> Update();  
}

TLegend* GenerateLegend(std::vector<TH1F*> Hist_V, TCanvas* can)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.5, 0.8); 
  for (TH1F* H : Hist_V)
  {
    len -> AddEntry(H); 
  }
  len -> Draw("SAME");
  can -> Update(); 
  return len; 
}

void GenerateRatioPlot(TH1F* H1, TH1F* H2, TCanvas* can, TString Title, TString Xaxis)
{
  TString name = "Ratio Plot: "; name += (H1 -> GetTitle()); name += (" / "); name += H2 -> GetTitle(); 
  TH1F* Ratio = (TH1F*)H1 -> Clone(name);
  Ratio -> SetTitle(name);
  Ratio -> Divide(H2); 

  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.5, 0.8); 

  TH1F* Empty = (TH1F*)H1 -> Clone("empty"); 
  Empty -> Reset();
  Empty -> SetTitle(Title);
  
  len -> AddEntry(H1); 
  len -> AddEntry(H2);  
  
  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1.0);

  if (Xaxis == "LOG")
  {
    P1 -> SetLogy(); 
    Empty -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
  }
  else { Empty -> GetXaxis() -> SetTitle(Xaxis);}

  P1 -> Draw(); 
  P1 -> cd();
  Empty -> Draw("HIST");   
   
  H1 -> Draw("SAMEHIST"); 
  H2 -> Draw("SAMEHIST"); 
  len -> Draw("SAME");  

  can -> cd();

  TPad *P2 = new TPad("P2", "P2", 0.0, 0.05, 1, 0.3);  
  P2 -> Draw();
  P2 -> cd(); 
	Ratio -> SetStats(0); 
  Ratio -> SetLineColor(kWhite);
  Ratio -> GetYaxis() -> SetRangeUser(0.25, 1.75);
  Ratio -> SetMarkerSize(1); 
  Ratio -> Draw("epl"); 
}
