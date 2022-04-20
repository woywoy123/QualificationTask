#ifndef BEST_FLOST_PREDICTION
#define BEST_FLOST_PREDICTION

#include "../PlottingCode/IO.h"
#include "../PlottingCode/SharedFunctions.h"
#include<thread>

_TS Flost; 
bool ErrorMode;
float min = -1; 
float max = -1; 
TGraph* MakeGraph(std::vector<float> Entries, _TS Title)
{
  
  TGraph* gr = new TGraph(); 
 
  for (int i(0); i < Energy_H.size(); i++)
  {
    gr -> SetPoint(i, i+0.5, Entries[i]); 
  }

  TH1F* H = new TH1F(Title, Title, Energy_H.size(), 0, Energy_H.size());
  gr -> SetHistogram(H); 
  gr -> SetTitle(Title); 
  H -> SetTitle(Title); 
  
  _TS x_Title; 
  if (ErrorMode){ x_Title = "#frac{|#Delta F_{Lost_{" + Flost +"}}|}{F_{Lost_{" + Flost + "}(Truth)}} (%)"; }
  else { x_Title = "F_{Lost_{" + Flost +"}}"; }

  H -> GetYaxis() -> SetTitle(x_Title);
  H -> GetYaxis() -> CenterTitle();
  H -> GetYaxis() -> SetTitleSize(0.025);
  H -> GetYaxis() -> SetTitleOffset(1.5);

  H -> GetXaxis() -> SetTitle("Jet Energy (GeV)");
  H -> GetXaxis() -> CenterTitle("Jet Energy (GeV)");
  H -> GetXaxis() -> SetTitleOffset(1); 
 
  for (int i(0); i < Energy_H.size(); i++)
  {
    _TS name = Energy_H[i]; name.ReplaceAll("_GeV", "").ReplaceAll("_", "-");
    gr -> GetXaxis() -> SetBinLabel(i+1, name); 
  }

  if (Title.Contains("Minimizer")){ gr -> SetLineStyle(kSolid); }
  if (Title.Contains("FitTo")){ gr -> SetLineStyle(kDashed); }
  
  return gr;
};


void PlotMultiGraph(std::vector<TGraph*>* GR, _TS Title, _TS Output)
{

  TLegend* len = new TLegend(0.8, 0.9, 0.65, 0.75); 
  TGraph* n; 
  TCanvas* can = new TCanvas(); 
  ConformCanvas(can); 
  can -> SetLeftMargin(0.105); 

  std::vector<TGraph*> out; 
  n = (TGraph*)(*GR)[0] -> Clone("TMP");
  n -> SetTitle(Title); 
  n -> SetLineColor(kWhite);
  n -> GetYaxis() -> SetRangeUser(min*0.5, max*2);  
  n -> Draw("ALP");   
  out.push_back(n);  

  int f, m = 0; 
  for (int i(0); i < (*GR).size(); i++)
  {
    TGraph* gr = (*GR)[i]; 
    _TS T = gr -> GetTitle(); 
    if (T.Contains("FitTo"))
    {
      gr -> SetLineColor(Colors_H[f]);
      gr -> SetFillColorAlpha(Colors_H[f], 0.2);
      f++; 
    }
    else if (T.Contains("Minimizer"))
    {
      gr -> SetLineColor(Colors_H[m]);
      gr -> SetFillColorAlpha(Colors_H[m], 0.2);
      m++; 
    }
    else { gr -> SetLineColor(Colors_H[i]); }

    gr -> SetLineWidth(2);
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
};




#endif
