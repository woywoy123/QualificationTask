#include<TH1F.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<TStyle.h>
#include<TRandom.h>
#include<TFile.h>
#include<RooAddPdf.h>
#include<RooHistPdf.h>
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooPlot.h>
#include<RooFFTConvPdf.h>
#include<TStyle.h>
#include<PostAnalysis/Constants.h>

#ifndef PLOTTING_H
#define PLOTTING_H

class Plotting
{
  public:
    TCanvas* SimplePlot(TH1F* Hist);
    void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style = kDashed);
    TCanvas* PlotHists(std::vector<TH1F*> Hists);
    void PlotHists(std::vector<TH1F*> Hists, TCanvas* can);
    TCanvas* PlotHists(std::vector<TH1F*> Hists, TH1F* Data);
    void PlotHists(std::vector<TH1F*> Hists, TH1F* Data, TCanvas* can);
    TCanvas* PlotHists(TH1F* Hists, TH1F* Data);
    TCanvas* PlotHists(std::vector<std::vector<TH1F*>> Hists, std::vector<TH1F*> Data);
    void PlotHists(std::vector<std::vector<TH1F*>> Hists, std::vector<TH1F*> Data, TCanvas* can); 
    void PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Closure, TH1F* Data, TCanvas* can);
    void PlotHists(std::vector<std::vector<TH1F*>> Hists, std::vector<std::vector<TH1F*>> Closure, std::vector<TH1F*> Data, TCanvas* can);
    void PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Closure, std::vector<TH1F*> Data, TCanvas* can);

    TCanvas* PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Data);
    TCanvas* PlotHists(std::vector<std::vector<TH1F*>> Hists);
    TCanvas* PlotHists(RooAddPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data);
    TCanvas* PlotHists(RooAddPdf model, RooRealVar* Domain, std::vector<RooFFTConvPdf*> PDFs, RooDataHist* Data);
    TCanvas* PlotHists(RooHistPdf model, RooRealVar Domain, RooDataHist Data);

};

class DistributionGenerators
{
  public:
    // Landau Generator 
    void Landau(std::vector<TH1F*> Hists, 
                std::vector<float> COMP, 
                std::vector<float> Parameters, 
                int Number, float min, float max);
    // Gaussian Generator 
    void Gaussian(float mean, float stdev, int Number, TH1F* Hist);

    // Monte Carlo
    std::vector<TH1F*> FillTH1F(std::vector<TString>, TString dir);
    TH1F* FillTH1F(TString name, std::vector<TString> SubDir, TString dir, std::vector<TString> Detector = Constants::Detector);

  private:
    TH1F* TH1FFromFile(TString Name, TString Layer, TFile* file);  

};

#endif
