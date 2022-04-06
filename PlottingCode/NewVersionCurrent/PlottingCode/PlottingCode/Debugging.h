#ifndef DEBUGGING_H
#define DEBUGGING_H

#include "../PlottingCode/IO.h"
#include "../PlottingCode/SharedFunctions.h"

typedef std::map<_TS, std::map<_TS, std::map<_TS, std::map<_TS, TH1F*>>>> maps; 
typedef std::map<_TS, std::map<_TS, std::map<_TS, std::map<_TS, TH1F*>>>>::iterator maps_i; 

void Residual(std::vector<TH1F*> pred, std::vector<TH1F*> truth)
{
  
  TH1F* R1 = (TH1F*)pred[0] -> Clone("Residual"); 
  R1 -> Reset();
  R1 -> GetYaxis() -> SetRangeUser(0.5, 1.5);
  R1 -> SetFillColor(kWhite);
  R1 -> SetTitle("");
  R1 -> GetYaxis() -> SetTitle("Fit/Truth");
  R1 -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
  R1 -> GetXaxis() -> CenterTitle();
  R1 -> GetYaxis() -> CenterTitle();
  R1 -> GetXaxis() -> SetLabelSize(0.075);
  R1 -> GetYaxis() -> SetLabelSize(0.075);
  R1 -> SetTitleOffset(0.4, "Y");
  R1 -> GetXaxis() -> SetTitleSize(0.075);
  R1 -> GetYaxis() -> SetTitleSize(0.075);
  R1 -> SetMarkerStyle(kCircle); 
  R1 -> SetMarkerSize(0.3); 
  R1 -> SetLineColor(kBlack);

  for (int j(0); j < R1 -> GetNbinsX(); j++)
  {
    float sum_t = 0; 
    float sum_f = 0; 
    for (int i(0); i < pred.size(); i++)
    {
      sum_t += truth[i] -> GetBinContent(j+1); 
      sum_f += pred[i] -> GetBinContent(j+1);    
    }

    R1 -> SetBinError(j+1, 0.0001);
    if (sum_t == 0 || sum_f == 0){continue;}
    R1 -> SetBinContent(j+1, sum_f/sum_t);

  }
  R1 -> Draw("PE3");
}



void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style)
{
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = Hists[i];
    
    TString name = H -> GetTitle(); 
    if (name.Contains("ntrk1-tru1")){ H -> SetTitle("Truth-1 (Monte Carlo)"); }
    if (name.Contains("ntrk1-tru2")){ H -> SetTitle("Truth-2 (Monte Carlo)"); }

    if (name.Contains("Case") && name.Contains("ntru_1")){ H -> SetTitle("Truth-1 (Fit)"); }
    if (name.Contains("Case") && name.Contains("ntru_2")){ H -> SetTitle("Truth-2 (Fit)"); }

    H -> SetLineColor(Colors_H[i]);
    H -> SetFillColorAlpha(Colors_H[i], 0.2); 
    H -> SetLineStyle(Style);
    H -> SetLineWidth(2);
    H -> Draw("SAMEHIST");
    len -> AddEntry(H, H -> GetTitle());
    len -> Draw("SAME");
    can -> Update(); 
  }
}

void Printer(TCanvas* can, _TS dir, std::vector<TH1F*> trks, std::vector<TH1F*> truth, _TS Title)
{
  can -> Clear();
  gStyle -> SetOptTitle(1); 
  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1);
  P1 -> SetBottomMargin(0.0001); 
  P1 -> SetBorderMode(0); 
  
  TPad *P2 = new TPad("P2", "P2", 0, 0, 1, 0.3);
  P2 -> SetTopMargin(0.0001); 
  P2 -> SetBorderMode(0);
  P2 -> SetBottomMargin(0.25); 

  P1 -> Draw(); 
  P2 -> Draw(); 
  P1 -> cd(); 
  
  P1 -> SetLogy();
  TLegend* len = new TLegend(0.8, 0.9, 0.65, 0.75); 
  TH1F* sum = (TH1F*)truth[0] -> Clone("Summed");
  sum -> Add(truth[1]); 
  sum -> SetLineColor(kBlack);
  sum -> SetLineStyle(kSolid);
  sum -> SetLineWidth(2); 
  len -> AddEntry(sum, "Summed Distribution");
  
  TH1F* empty = (TH1F*)truth[0] -> Clone("empty");
  empty -> Reset();
  empty -> SetTitle(Title); 
  empty -> SetLineWidth(2); 
  empty -> GetYaxis() -> SetTitle("Clusters");
  empty -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^{2}]"); 
  empty -> GetXaxis() -> CenterTitle();
  empty -> GetYaxis() -> CenterTitle();
  empty -> GetYaxis() -> SetRangeUser(1, sum -> Integral()); 
  empty -> Draw("HIST");
  
  Populate(trks, can, len, kSolid); 
  Populate(truth, can, len, kDashed); 
  sum -> Draw("SAME");

  P2 -> cd();
  P2 -> SetGrid();
  Residual(trks, truth); 
  can -> cd();

  can -> Print(dir);

};

TGraph* MakeGraphSteps(std::vector<float> Values, _TS Title, float step_s, float start)
{
  
  TGraph* gr = new TGraph(); 
 
  for (int i(0); i < Values.size(); i++)
  {
    gr -> SetPoint(i, i+0.5, Values[i]); 
  }

  TH1F* H = new TH1F(Title, Title, Values.size(), 0, Values.size());
  gr -> SetHistogram(H); 
  gr -> SetTitle(Title); 
  H -> SetTitle(Title); 

  H -> GetYaxis() -> SetTitle("Prediction/Truth");
  H -> GetYaxis() -> CenterTitle("Prediction/Truth");

  H -> GetXaxis() -> SetTitle("Scaling of Original 2-Truth Distribution (%)");
  H -> GetXaxis() -> CenterTitle("Scaling of Original 2-Truth Distribution (%)");
 
  for (int i(0); i < Values.size(); i++)
  {
    float cur = start + step_s*(i+1); 
    TString name = ""; name += (int(cur));  
    gr -> GetXaxis() -> SetBinLabel(i+1, name); 
  }
  return gr;
};


TLegend* GenerateLegend(std::vector<TGraph*> Hist_V, TCanvas* can, float x_size, float x_position, float box_size, float size)
{
  TLegend* len = new TLegend(x_size, box_size, x_position, size); 
  len -> SetBorderSize(0); 
  for (TGraph* H : Hist_V)
  {
    len -> AddEntry(H);
  }
  return len; 
}

void CombineStepGraphs(std::vector<TGraph*> GR, _TS Title, TCanvas* can, _TS Dir, _TS XTitle, _TS YTitle, float min, float max)
{
  auto Conform = [&](TGraph* g, ELineStyle Style, int index)
  {
    g -> SetFillColorAlpha(Colors_H[index], 0.2); 
    g -> SetLineColor(Colors_H[index]);
    g -> SetLineWidth(2.);
    g -> SetLineStyle(Style);
    g -> Draw("SAME");
  };

  TLegend* len = GenerateLegend(GR, can, 0.8, 0.5, 0.89, 0.79);
  can -> SetLogy(false);
  TGraph* H = (TGraph*)GR[0] -> Clone("EMPTY"); 
  H -> SetTitle(Title);
  H -> SetLineColor(kWhite); 
  H -> SetLineWidth(0.0001);
  H -> GetYaxis() -> SetTitle(YTitle);
  H -> GetYaxis() -> CenterTitle(YTitle);

  H -> GetXaxis() -> SetTitle(XTitle);
  H -> GetXaxis() -> CenterTitle(XTitle);
  H -> GetYaxis() -> SetRangeUser(min, max); 
  H -> Draw("ALP"); 
  
  Conform(GR[0], kDashed, 0); 
  Conform(GR[1], kDashed, 1);
  Conform(GR[2], kDashed, 2);

  Conform(GR[3], kSolid, 0); 
  Conform(GR[4], kSolid, 1);
  Conform(GR[5], kSolid, 2);
  
  can -> SetLogy(); 
  can -> Update(); 
  can -> SetLogx(false);
  can -> Update(); 
  len -> Draw("SAME");
  
  can -> Print(Dir);

  delete H;
}

TGraph* MakeGraphJets(std::vector<float> Entries, _TS Title)
{
  
  TGraph* gr = new TGraph(); 
 
  for (int i(0); i < Entries.size(); i++)
  {
    gr -> SetPoint(i, i+0.5, float(Entries[i])); 
  }
  
  TH1F* H = new TH1F(Title, Title, Energy_H.size(), 0, Energy_H.size());
  gr -> SetHistogram(H); 
  gr -> SetTitle(Title); 
  H -> SetTitle(Title); 

  for (int i(0); i < Energy_H.size(); i++)
  {
    _TS name = Energy_H[i]; name.ReplaceAll("_GeV", "").ReplaceAll("_", "-");
    gr -> GetXaxis() -> SetBinLabel(i+1, name); 
  }
  return gr;
}



#endif
