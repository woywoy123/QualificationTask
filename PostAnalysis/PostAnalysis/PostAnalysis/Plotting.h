#include<TStyle.h>
#include<TH1F.h>
#include<TLegend.h>
#include<TCanvas.h>
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooAddPdf.h>
#include<RooPlot.h>
#include<PostAnalysis/BaseFunctions.h>
#include<RooFFTConvPdf.h>
#include<TPaveText.h>
#include<sstream>
#include<iomanip>
#include<RooHist.h>
#include<RooMinimizer.h>
#include<RooRealVar.h>
#include<RooFit.h>
#include<RooProfileLL.h>
#include<TMultiGraph.h>
#include<TGraph.h>
#include<RooRealSumPdf.h>
#include<THStack.h>

#ifndef PLOTTING_H 
#define PLOTTING_H


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
void PlotHists(TH1F* Data, TH1F* Hist, TString Title, TCanvas* can);

void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can); 
void RooFitPullPlot(RooAddPdf model, RooRealVar* Domain, std::vector<RooFFTConvPdf*> PDFs, RooDataHist* Data, TString Name);
void RooFitPullPlot(RooAddPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data, TString Name);

void RooFitPullPlot(RooRealSumPdf model, RooRealVar* Domain, std::vector<RooFFTConvPdf*> PDFs, RooDataHist* Data, TString Name);
void RooFitPullPlot(RooRealSumPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data, TString Name);

void PlotRooFit(RooAddPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data);

void GeneratePlot(TH1F* H, TString Title, TCanvas* can, Color_t color, ELineStyle style, TString DrawOption, float Intensity); 
void GenerateRatioPlot(TH1F* H1, TH1F* H2, TCanvas* can, TString Title, TString Xaxis);

void PlotLikelyHood(RooAbsReal* nll, RooRealVar* var, TString name); 
void PlotGraphs(std::vector<TH1F*>, TString Title, TCanvas* can);
void GenerateNiceStacks(std::vector<TH1F*> vec, TString Title, TCanvas* can, TString x_axis, TString y_axis, TString options); 
void GeneratePerformanceGraphs(std::vector<TGraph*> Graphs, TString title, TString x_axis, TString y_axis, float min, float max, TCanvas* can); 
TGraph* GenerateGraph(std::vector<float> Input, TString name); 
TGraph* GenerateGraph(std::map<TString, float> Input, TString name); 

static std::vector<Color_t> Colors_F = {kRed, kGreen, kBlue, kCyan, kViolet, kOrange, kCoffee, kAurora}; 


#endif
