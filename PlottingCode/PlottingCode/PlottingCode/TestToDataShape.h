#ifndef TEST_TO_DATA_SHAPE_H
#define TEST_TO_DATA_SHAPE_H

#include "../PlottingCode/SharedFunctions.h"
#include "../PlottingCode/IO.h"

void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style)
{
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = Hists[i];
    
    TString name = H -> GetTitle(); 
    if (name.Contains("005_ntru_1")){ H -> SetTitle("Truth-1 (Monte Carlo)"); }
    if (name.Contains("005_ntru_2")){ H -> SetTitle("Truth-2 (Monte Carlo)"); }
    if (name.Contains("005_ntru_3")){ H -> SetTitle("Truth-3 (Monte Carlo)"); }
    if (name.Contains("005_ntru_4")){ H -> SetTitle("Truth-4 (Monte Carlo)"); }

    if (!name.Contains("rless") && name.Contains("ntru_1")){ H -> SetTitle("Truth-1 (Fit)"); }
    if (!name.Contains("rless") && name.Contains("ntru_2")){ H -> SetTitle("Truth-2 (Fit)"); }
    if (!name.Contains("rless") && name.Contains("ntru_3")){ H -> SetTitle("Truth-3 (Fit)"); }
    if (!name.Contains("rless") && name.Contains("ntru_4")){ H -> SetTitle("Truth-4 (Fit)"); }


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

void Residual(std::vector<TH1F*> pred, TH1F* Data)
{
  
  TH1F* R1 = (TH1F*)pred[0] -> Clone("Residual"); 
  R1 -> Reset();
  R1 -> GetYaxis() -> SetRangeUser(-1.2, 1.2);
  R1 -> SetFillColor(kWhite);
  R1 -> SetTitle("");
  R1 -> GetYaxis() -> SetTitle("log(Fit/Truth)");
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
    float sum_f = 0; 
    float sum_t = Data -> GetBinContent(j+1); 
    for (int i(0); i < pred.size(); i++)
    {
      sum_f += pred[i] -> GetBinContent(j+1);    
    }

    R1 -> SetBinError(j+1, 0.0001);
    if (sum_t == 0 || sum_f == 0){continue;}
    R1 -> SetBinContent(j+1, std::log(sum_f/sum_t));

  }
  R1 -> Draw("PE3");
}


void PlotHists(TH1F* Data, std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TString title, TCanvas* can)
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
  TH1F* empty = (TH1F*)truth[0] -> Clone(title); 
  empty -> SetTitle(title); 
  empty -> Reset();

  empty -> GetYaxis() -> SetRangeUser(1, sum*1.2);
  empty -> GetYaxis() -> SetTitle("Clusters");
  empty -> GetYaxis() -> CenterTitle();

  empty -> GetXaxis() -> SetRangeUser(0, 10);
 	empty -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
  empty -> GetXaxis() -> CenterTitle();

  empty -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.89, 0.75, 0.65); 
  len -> SetBorderSize(0); 

  Data -> SetTitle("Cluster Measurements"); 
  Data -> SetLineColor(kBlack); 
  Data -> SetLineWidth(2); 
  len -> AddEntry(Data, Data -> GetTitle());
  Data -> Draw("SAMEHIST"); 

  Populate(truth, can, len, kDashed); 
  Populate(prediction, can, len, kSolid); 

  P2 -> cd();
  P2 -> SetGrid();
  Residual(prediction, Data); 
  can -> cd();
}

TGraph* MakeGraph(std::map<_TS, float> Entries, _TS Title, _TS Layer, _TS ntrk)
{
  
  TGraph* gr = new TGraph(); 
 
  for (int i(0); i < Energy_H.size(); i++)
  {
    gr -> SetPoint(i, i+0.5, Entries[ntrk +"/" + Layer + "/" + Energy_H[i] + "/" + Title]); 
  }

  TH1F* H = new TH1F(Title, Title, Energy_H.size(), 0, Energy_H.size());
  gr -> SetHistogram(H); 
  gr -> SetTitle(Title); 
  H -> SetTitle(Title); 

  H -> GetYaxis() -> SetTitle("Weighted Shape Error");
  H -> GetYaxis() -> CenterTitle("Weighted Shape Error");

  H -> GetXaxis() -> SetTitle("Jet Energy (GeV)");
  H -> GetXaxis() -> CenterTitle("Jet Energy (GeV)");
 
  for (int i(0); i < Energy_H.size(); i++)
  {
    _TS name = Energy_H[i]; name.ReplaceAll("_GeV", "").ReplaceAll("_", "-");
    gr -> GetXaxis() -> SetBinLabel(i+1, name); 
  }
  return gr;
}

void PlotMultiGraph(std::vector<TGraph*>* GR, _TS Title, float min, float max, _TS Output)
{
  TLegend* len = new TLegend(0.8, 0.9, 0.65, 0.75); 
  TGraph* n; 
  TCanvas* can = new TCanvas(); 
  ConformCanvas(can); 

  std::vector<TGraph*> out; 
  n = (TGraph*)(*GR)[0] -> Clone("TMP");
  n -> SetTitle(Title); 
  n -> SetLineColor(kWhite);
  n -> GetYaxis() -> SetRangeUser(min, max);  
  n -> Draw("ALP");   
  out.push_back(n);  

  for (int i(0); i < (*GR).size(); i++)
  {
    TGraph* gr = (*GR)[i]; 
    gr -> SetLineColor(Colors_H[i]);
    gr -> SetLineWidth(5);
    gr -> GetYaxis() -> SetRangeUser(min, max);  
    gr -> Draw("SAME");
    len -> AddEntry(gr);
    out.push_back(gr); 
  }
  can -> SetLogy(true); 
  can -> SetLogx(false);
  can -> Update(); 
  len -> Draw("SAME");
 
  can -> Print(Output);
  BulkDelete(out); 
  delete can; 
  (*GR)={};
}




















#endif
