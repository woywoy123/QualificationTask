#include<TStyle.h>
#include<TH1F.h>
#include<TLegend.h>
#include<TCanvas.h>

void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style);
void PlotHists(TH1F* Hist, TCanvas* can); 
void PlotHists(std::vector<TH1F*> Hists, TCanvas* can); 
void PlotHists(std::vector<TH1F*> Hists, std::vector<TString> Legend_Titles, TCanvas* can); 
void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can); 



