#include<TStyle.h>
#include<TH1F.h>
#include<TLegend.h>
#include<TCanvas.h>
#include<TPaveText.h>
#include<sstream>
#include<iomanip>
#include<RooHist.h>
#include<TMultiGraph.h>
#include<TGraph.h>
#include<THStack.h>

#ifndef PLOTTING_H 
#define PLOTTING_H

#include "BaseFunctions.h"

TLegend* GenerateLegend(std::vector<TH1F*> Hist_V, TCanvas* can, float x_size = 0.9, float x_position = 0.5, float box_size = 0.89, float size = 0.8); 
TLegend* GenerateLegend(std::vector<TGraph*> Hist_V, TCanvas* can, float x_size = 0.9, float x_position = 0.5, float box_size = 0.89, float size = 0.8);

void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style);
void PlotHists(TH1F* Hist, TCanvas* can); 
void PlotHists(std::vector<TH1F*> truth, std::vector<TH1F*> prediction, TCanvas* can); 
void PlotHists(std::vector<TH1F*> Hists, TCanvas* can); 
void PlotHists(std::vector<TH1F*> Hists, std::vector<TString> Legend_Titles, TCanvas* can); 
void PlotHists(TH1F* Data, std::vector<TH1F*> truth, std::vector<TH1F*> prediction, TCanvas* can); 
void PlotHists(TH1F* Data, std::vector<TH1F*> truth, TCanvas* can); 
void PlotHists(TH1F* Data, std::vector<TH1F*> Prediction, std::vector<TH1F*> Truth, TString Title, float FLost_P, float FLost_T, TCanvas* can);
void PlotHists(std::vector<TH1F*> truth, std::vector<TH1F*> prediction, TString title, TCanvas* can); 

void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can); 

void GeneratePlot(TH1F* H, TString Title, TCanvas* can, Color_t color, ELineStyle style, TString DrawOption, float Intensity); 
void GenerateRatioPlot(TH1F* H1, TH1F* H2, TCanvas* can, TString Title, TString Xaxis);

void PlotGraphs(std::vector<TH1F*>, TString Title, TCanvas* can);
void GenerateNiceStacks(std::vector<TH1F*> vec, TString Title, TCanvas* can, TString x_axis, TString y_axis, TString options); 
TGraph* GenerateGraph(std::vector<float> Input, TString name); 
TGraph* GenerateGraph(std::map<TString, float> Input, TString name, TString yname = "Percent Error %"); 
std::vector<TGraph*> GenerateMultiGraph(MMF ntru, TString Title, TCanvas* can, float min = 0.01, float max = 100, TString yname = "Percent Error %"); 

#endif
