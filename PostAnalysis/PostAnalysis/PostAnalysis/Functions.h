// Include all the required ROOT libraries
#include<TString.h>
#include<TFile.h>
#include<TH1F.h>
#include<TImage.h>
#include<TApplication.h>
#include<TSystem.h>
#include<TCanvas.h>
#include<TRandom.h>

// Include the standard C++ library 
#include<iostream>
#include<numeric>

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
#include<RooFormulaVar.h>
#include<RooGaussian.h>
#include<RooFFTConvPdf.h>

using namespace RooFit;

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

class Functions
{
  public:
    std::vector<TH1F*> MakeTH1F(std::vector<TString> Names,  int bins, float min, float max, TString Extension = "");

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
    void VectorToTH1F(std::vector<float> Vec, TH1F* hist, int StartBin = 0);
    std::vector<float> TH1FToVector(TH1F* hist, int CustomLength = 0, int start = 0);
    void CutTH1F(TH1F* Input, TH1F* Output);
    void ExpandTH1F(TH1F* Input, TH1F* Output, int Start = 0);
};

class Fit_Functions
{
  public: 
    // Bulk variable generation.
    std::vector<RooDataHist*> ConvertTH1toDataHist(std::vector<TH1F*> Histograms, RooRealVar* domain);
    std::vector<RooHistPdf*> ConvertTH1FtoPDF(std::vector<TH1F*> Histograms, RooRealVar* domain);
    std::vector<RooRealVar*> GenerateVariables(std::vector<TString>, std::vector<float> Begin, std::vector<float> End); 
    std::vector<RooGaussian*> GaussianVariables(std::vector<TString> Names, std::vector<RooRealVar*> Mean, std::vector<RooRealVar*> Stdev, RooRealVar* Domain);
    std::vector<RooFFTConvPdf*> ConvolveVariables(std::vector<TString> Names, std::vector<RooHistPdf*> PDFs, std::vector<RooGaussian*> Gaus, RooRealVar* domain);

    // Base variables used for vector derivation 
    RooHistPdf* ConvertTH1FtoPDF(RooDataHist* Histogram, TString Name, RooRealVar* domain);
    RooDataHist* ConvertTH1toDataHist(TH1F* Hist, RooRealVar* domain);
    RooDataHist* ConvertTH1toDataHist(TH1* Hist, RooRealVar* domain);
    RooRealVar* GenerateVariable(TString name, float begin, float end);   
    RooGaussian* GenerateGaussian(TString name, RooRealVar* mean, RooRealVar* stdev, RooRealVar* domain);
    RooFFTConvPdf* GenerateConvolve(TString name, RooHistPdf* PDF, RooGaussian* Gaus, RooRealVar* domain);

    // Generate the RooArgList varibles for a model 
    RooArgList VectorToArgList(std::vector<RooRealVar*> Vector);
    RooArgList VectorToArgList(std::vector<RooHistPdf*> Vector);

    // Deconvolution algorithm Lucy Richardson 
    std::vector<float> LRDeconvolution(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y);
    
    // Fast Fourier Transformation of two TH1F histograms  
    void ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv, int StartBin = 0); 
   
    // Replace the tail of a histogram with another after fitting  
    std::vector<float> TailReplaceClosure(TH1F* hist, std::vector<float> deconv);
    std::vector<float> TailReplace(TH1F* hist, std::vector<float> deconv);

    // Normalize the histograms
    void Normalizer(TH1F* Hist);

    // Fitting PDFs to data
    std::vector<RooRealVar*> FitPDFtoData(std::vector<TH1F*> PDFs, TH1F* Data, float min, float max);
   
    // Cleaning up the cross contamination via subtraction 
    void Subtraction(std::vector<TH1F*> nTrk, TH1F* Target, int ntrk, std::vector<RooRealVar*> var);
    
    // Fractionalizer
    std::vector<float> Fractionalizer(std::vector<RooRealVar*> vars, TH1F* Data);

    // Removes the weird artifacts before the peak from the hists
    void ArtifactRemove(TH1F* Hist, TString mode = "b");

    // Creates the datavector with the mirror tail 
    std::vector<float> TH1FDataVector(TH1F* Data, float Offset);

    // Produce a gaussian distribution
    void GaussianGenerator(float mean, float std, int N, TH1F* Hist);

    // Gaussian convolution fit 
    std::vector<RooRealVar*> GaussianConvolutionFit(std::vector<TH1F*> PDFs, TH1F* trk2, float min, float max, float offset, float stdev_s, float stdev_e, float mean_s, float mean_e); 

  private:
    std::vector<float> ConvolveHists(std::vector<float> Hist1, std::vector<float> Hist2);
    int ArtifactRemoveForward(TH1F* Hist);
    int ArtifactRemoveBackward(TH1F* Hist);   
};

class Algorithms
{
  public:
    void MinimalAlgorithm(TH1F* trk1_D, TH1F* trk2_D, std::vector<TH1F*> O_PDFs, float min, float max, float offset, int iter, int trkn = 2);
    std::vector<float> GaussianAlgorithm(TH1F* trk1_D, TH1F* trk2_D, std::vector<TH1F*> O_PDFs, float min, float max, float offset, float mean_s, float mean_e, float stdev_s, float stdev_e, int iter, int trkn = 2);
};


class Benchmark
{
  public:
    float WeightedEuclidean(std::vector<float> v1, std::vector<float> v2);
    float PythagoreanDistance(std::vector<float> v1, std::vector<float> v2);
    TCanvas* ClosurePlot(TString Name, std::vector<TH1F*> Data, std::vector<TH1F*> PDFs, std::vector<std::vector<float>> Closure);
};

namespace Constants
{
  const std::vector<TString> Detector = {"Blayer"};
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
  const std::vector<Color_t> Colors = {kOrange, kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kGray, kBlack, kCoffee, kAurora};
}

#endif
