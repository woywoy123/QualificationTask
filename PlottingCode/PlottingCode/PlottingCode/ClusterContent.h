#ifndef CLUSTER_CONTENT_H
#define CLUSTER_CONTENT_H

#include "../PlottingCode/IO.h"
#include "../PlottingCode/SharedFunctions.h"

TGraph* MakeGraph(std::map<_TS, std::vector<float>> Entries, _TS Title, int ntrk)
{
  
  TGraph* gr = new TGraph(); 
 
  for (int i(0); i < Energy_H.size(); i++)
  {
    gr -> SetPoint(i, i+0.5, float(Entries[Energy_H[i]][ntrk])); 
  }

  TH1F* H = new TH1F(Title, Title, Energy_H.size(), 0, Energy_H.size());
  gr -> SetHistogram(H); 
  gr -> SetTitle(Title); 
  H -> SetTitle(Title); 

  H -> GetYaxis() -> SetTitle("Entries");
  H -> GetYaxis() -> CenterTitle("Entries");

  H -> GetXaxis() -> SetTitle("Jet Energy (GeV)");
  H -> GetXaxis() -> CenterTitle("Jet Energy (GeV)");
 
  for (int i(0); i < Energy_H.size(); i++)
  {
    _TS name = Energy_H[i]; name.ReplaceAll("_GeV", "").ReplaceAll("_", "-");
    gr -> GetXaxis() -> SetBinLabel(i+1, name); 
  }
  return gr;
}

std::vector<TGraph*> PlotMultiGraph(std::map<_TS, std::vector<float>> GR, _TS Title, TCanvas* can, int ntrks, float min, float max)
{
  TLegend* len = new TLegend(0.8, 0.9, 0.65, 0.75); 
  TGraph* n; 
  can -> SetLogy(false);
  
  std::vector<TGraph*> out; 
  for (int i(0); i < ntrks; i++)
  {
    _TS name = "Track-"; name +=(i+1); 
    TGraph* gr = MakeGraph(GR, name, i); 
    gr -> SetLineColor(Colors_H[i]);
    gr -> SetLineWidth(5);
    gr -> GetYaxis() -> SetRangeUser(min, max);  
    out.push_back(gr); 
    if ( i == 0)
    {
      n = (TGraph*)gr -> Clone("TMP");
      n -> SetTitle(Title); 
      n -> Draw("ALP");   
      out.push_back(n);  
    }
    gr -> Draw("SAME");
    len -> AddEntry(gr);
  }
  
  can -> SetLogy(); 
  can -> Update(); 
  can -> SetLogx(false);
  can -> Update(); 
  len -> Draw("SAME");
  
  return out;
}




#endif 
