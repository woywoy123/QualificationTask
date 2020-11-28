#ifndef PRESENTATIONFIGURES_H
#define PRESENTATIONFIGURES_H
#include<TFile.h>
#include<PostAnalysis/Plotting.h>

void FigureCompiler(TFile* F); 
void PlotLandau(TCanvas* can, std::vector<TH1F*> Hist_Vi, TString filename); 
void PlotGaussian(TCanvas* can, std::vector<TH1F*> Hist_Vi, TString filename); 
void PlotGaussianXGaussian(TCanvas* can, std::vector<TH1F*> Hist_Vi, TString filename); 
void PlotLandauXLandau(TCanvas* can, std::vector<TH1F*> Hist_Vi, TString filename); 
void PlotLandauXGaussian(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename); 
void PlotDeconvGausXGaus(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename); 
void PlotDeconvLandauXLandau(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename); 
void PlotDeconvLandauXGaussian(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename); 
void PlotGaussianDeconvolutionFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename); 
void PlotLandauXGausFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotNLandauXNGausFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotDeconvolutionFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotComparisonBinCenteringLandauXLandau(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotOscillationLucyRichardson(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotAlgorithm(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotTestReadFile(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);
void PlotMonteCarloMatchConvolution(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename);

#endif
