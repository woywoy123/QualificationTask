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


void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can)
{

  TH1F* Ratio = (TH1F*)H1 -> Clone("Ratio");
  Ratio -> Clear();  
  for (int i(0); i < H1 -> GetNbinsX(); i++)
  {
    float e1 = H1 -> GetBinContent(i+1); 
    float e2 = H2 -> GetBinContent(i+1); 
    float r; 
    if (e2 != 0){r = e1/e2;}
    else {r = 1; }
    Ratio -> SetBinContent(i+1, r); 
  }
	H1 -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^2]");

  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1.0);
  P1 -> Draw(); 
  P1 -> cd(); 
  P1 -> SetLogy(); 
  H1 -> SetLineColor(kBlack); 
  H1 -> Draw("SAMEHIST"); 
  H1 -> SetStats(0);  
  H2 -> SetLineColor(kRed);
  H2 -> Draw("SAMEHIST"); 
  
  can -> cd();

  TPad *P2 = new TPad("P2", "P2", 0.0, 0.05, 1, 0.3);  
  P2 -> Draw();
  P2 -> cd(); 
	Ratio -> SetStats(0); 
  Ratio -> Draw("SAMEHIST"); 
}
