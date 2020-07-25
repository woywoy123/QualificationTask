// Include all the required ROOT libraries
#include<TString.h>
#include<TFile.h>
#include<TH1F.h>
#include<TImage.h>
#include<TApplication.h>
#include<TSystem.h>
#include<TCanvas.h>

// Include the standard C++ library 
#include<iostream>

// Include RooFit classes 
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooHistPdf.h>
#include<RooFit.h>
#include<RooAddPdf.h>
#include<TROOT.h>
#include<TStyle.h>
#include<TLegend.h>
#include<RooPlot.h>
#include<TVirtualFFT.h>
#include<TComplex.h>

using namespace RooFit;

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

class Functions
{
  public:
    std::vector<TH1F*> MakeTH1F(std::vector<TString> Names,  int bins, int min, int max, TString Extension = "");

    void FillTH1F_From_File(std::vector<TH1F*> Histograms, TFile* File, TString DetectorLayer, TString Extension = "");
    void FillTH1F_From_File(TH1F* Histogram, TFile* File, TString DetectorLayer, std::vector<TString> List, TString Extension = "");

    std::vector<TString> SplitString(TString token, TString Split);
    TString RemoveExtension(TString Name, TString Extension);

    template <typename T> std::vector<T> AppendVectors(std::vector<std::vector<T>> Input)
    {
      std::vector<T> Output; 
      for (std::vector<T> i : Input)
      {
        for (T x : i){ Output.push_back(x); }
      }
      return Output;
    }

    TH1F* VectorToTH1F(std::vector<float> Vec, TString name, int bins, int min, int max);
    void VectorToTH1F(std::vector<float> Vec, TH1F* hist);
};

class Fit_Functions
{
  public:
    std::vector<RooDataHist*> ConvertTH1toDataHist(std::vector<TH1F*> Histograms, RooRealVar* domain);
    std::vector<RooHistPdf*> ConvertTH1FtoPDF(std::vector<TH1F*> Histograms, RooRealVar* domain);
    std::vector<RooRealVar*> GenerateVariables(std::vector<TString>, std::vector<double> Begin, std::vector<double> End); 
    
    RooHistPdf* ConvertTH1FtoPDF(RooDataHist* Histogram, TString Name, RooRealVar* domain);
    RooDataHist* ConvertTH1toDataHist(TH1F* Hist, RooRealVar* domain);
    RooDataHist* ConvertTH1toDataHist(TH1* Hist, RooRealVar* domain);
    RooRealVar* GenerateVariable(TString name, double begin, double end);   
  
    RooArgList VectorToArgList(std::vector<RooRealVar*> Vector);
    RooArgList VectorToArgList(std::vector<RooHistPdf*> Vector);

    std::vector<float> LRDeconvolution(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y);
    void ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv, int offset);  

  private:
    std::vector<float> ConvolveHists(std::vector<float> Hist1, std::vector<float> Hist2);
   
};

class Plot_Functions
{
  public: 
    TCanvas* GeneratePlot(TString Title, RooRealVar* range, RooDataHist* Data, RooAddPdf Model, std::vector<RooHistPdf*> PDFs, std::vector<TString> pdf_titles); 
};

class Benchmark
{
  public:
    float WeightedEuclidean(std::vector<float> v1, std::vector<float> v2);
};

namespace Constants
{
  const std::vector<TString> Detector = {"IBL", "Blayer", "layer1", "layer2"};
  const std::vector<TString> Pure_Names = {"dEdx_ntrk_1_ntru_1", 
                                           "dEdx_ntrk_2_ntru_2", 
                                           "dEdx_ntrk_3_ntru_3", 
                                           "dEdx_ntrk_4_ntru_4"};

  const std::vector<TString> trk_1 = {"dEdx_ntrk_1_ntru_1", 
                                      "dEdx_ntrk_1_ntru_2", 
                                      "dEdx_ntrk_1_ntru_3", 
                                      "dEdx_ntrk_1_ntru_4"};

  const std::vector<TString> trk_2 = {"dEdx_ntrk_2_ntru_1", 
                                      "dEdx_ntrk_2_ntru_2", 
                                      "dEdx_ntrk_2_ntru_3", 
                                      "dEdx_ntrk_2_ntru_4"};

  const std::vector<TString> trk_3 = {"dEdx_ntrk_3_ntru_1", 
                                      "dEdx_ntrk_3_ntru_2", 
                                      "dEdx_ntrk_3_ntru_3", 
                                      "dEdx_ntrk_3_ntru_4"};

  const std::vector<TString> trk_4 = {"dEdx_ntrk_4_ntru_1", 
                                      "dEdx_ntrk_4_ntru_2", 
                                      "dEdx_ntrk_4_ntru_3", 
                                      "dEdx_ntrk_4_ntru_4"};

  const std::vector<TString> Variable_Names = {"ntrk_1", "ntrk_2", "ntrk_3", "ntrk_4"};
  const std::vector<double> Begin = {0, 0, 0, 0};
  const std::vector<double> End = {1e8, 1e8, 1e8, 1e8};
  const std::vector<Color_t> Colors = {kRed, kBlue, kOrange, kCyan, kGreen, kYellow, kViolet, kAzure};
}


#endif
