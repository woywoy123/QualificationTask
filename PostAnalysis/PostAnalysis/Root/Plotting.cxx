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
  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  Hists[0] -> GetYaxis() -> SetRangeUser(1e-6, 1e6); 
  Populate(Hists, can, len, kSolid); 
}

void PlotHists(TH1F* Data, std::vector<TH1F*> Hists, TCanvas* can)
{
  gStyle -> SetOptStat(0); 
  Data -> SetLineColor(kBlack); 
  Data -> GetYaxis() -> SetRangeUser(0.1, Data -> Integral());
  Data -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(Hists, can, len, kSolid); 
}

void PlotHists(TH1F* Data, std::vector<TH1F*> truth, std::vector<TH1F*> prediction, TCanvas* can)
{
  can -> Clear();
  gStyle -> SetOptStat(0); 
  Data -> GetYaxis() -> SetRangeUser(1e-9, Data -> Integral());
  Data -> GetXaxis() -> SetRangeUser(0, 12);
  Data -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(truth, can, len, kSolid); 
  Populate(prediction, can, len, kDashed); 
}

void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can)
{
  TString name = "Ratio Plot: "; name += (H1 -> GetTitle()); name += ("/"); name += H2 -> GetTitle(); 
  TH1F* Ratio = (TH1F*)H1 -> Clone(name);
  Ratio -> SetTitle(name);
  Ratio -> Divide(H2); 
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

void PlotRooFit(RooAddPdf model, RooRealVar* Domain, RooDataHist* Data)
{
  RooPlot* xframe = Domain -> frame(RooFit::Title("Figure")); 
  Data -> plotOn(xframe, RooFit::Name("Data")); 
  model.plotOn(xframe); 
  TCanvas* can = new TCanvas(); 
  gPad -> SetLogy(); 
  xframe -> SetMinimum(1); 
  xframe -> Draw(); 
  can -> Update();
  can -> Print("debug.pdf"); 
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
