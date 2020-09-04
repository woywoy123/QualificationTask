#include<TString.h>
#include<TH1F.h>
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooHistPdf.h>
#include<RooAddPdf.h>
#include<TVirtualFFT.h>
#include<TComplex.h>
#include<RooGaussian.h>
#include<RooFFTConvPdf.h>
#include<RooFitResult.h>
#include<RooConstVar.h>

// Functions used for testing.
#include<RooPolynomial.h>
#include<RooDataSet.h>
#include<RooProdPdf.h>

#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

class BaseFunctions
{
  public: 
    // ==== Histogram Generators
    std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");
    std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, TH1F* Input);
    std::vector<TH1F*> CopyTH1F(std::vector<TH1F*> Hists, TString Extension);

    // ==== Data Generators and Properties 
    std::vector<float> Ratio(std::vector<TH1F*> Hists, TH1F* Data); 
    std::vector<float> Ratio(std::vector<RooRealVar*> Vars, TH1F* Data);
    std::vector<float> ClosureAndData(std::vector<TH1F*> Hists, TH1F* Data);  
    std::vector<float> TH1FDataVector(TH1F* Data, float offset = 0);
    void ToTH1F(std::vector<float> Vector, TH1F* Hist);

    // ==== Operational Functions  
    void Normalize(TH1F* Hist);
    void Normalize(std::vector<TH1F*> Hist);
    void Subtraction(std::vector<TH1F*> ntrk, TH1F* Data, int Exclude, std::vector<float> Ratios);  
    void ShiftExpandTH1F(TH1F* In, TH1F* Out, int start = 0);
    void ShiftExpandTH1F(std::vector<TH1F*> In, std::vector<TH1F*> Out, int start = 0); 
    void Scale(std::vector<TH1F*> PDFs, std::vector<RooRealVar*> Vars);

    // ==== RooFit functions 
    // Simple Scale Fits
    std::vector<RooRealVar*> RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End);
    std::vector<RooRealVar*> RooVariables(std::vector<TString> Names, std::vector<float> Var1, std::vector<float> Var2, std::vector<float> Var3); 
    std::vector<RooDataHist*> RooData(std::vector<TH1F*> Hist, RooRealVar* Domain); 
    RooDataHist* RooData(TH1F* Hist, RooRealVar* Domain);
    std::vector<RooHistPdf*> RooPDF(std::vector<TH1F*> Hist, RooRealVar* Domain); 
    RooArgList RooList(std::vector<RooRealVar*> Vector);
    RooArgList RooList(std::vector<RooHistPdf*> Vector);
    
    // Gaussian Convolution Fit
    std::vector<RooGaussian*> RooVariables(std::vector<TString> Names, std::vector<RooRealVar*> Mean, std::vector<RooRealVar*> Stdev, RooRealVar* Domain); 
    std::vector<RooGaussian*> RooVariables(std::vector<TString> Names, std::vector<RooRealVar*> V1, std::vector<Double_t> V2, std::vector<Double_t> V3); 
    std::vector<RooFFTConvPdf*> RooVariables(std::vector<TString> Names, std::vector<RooHistPdf*> PDFs, std::vector<RooGaussian*> Gaus, RooRealVar* Domain);
   
    // ==== Benchmarks  
    float ChiSquare(std::vector<float> V1, std::vector<float> V2);
    void PredictionTruthPrint(std::vector<float> Truth, std::vector<float> Prediction);

    // ==== Core Algorithm Functions 
    std::vector<float> LucyRichardson(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y = 0.75); 
    void ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv);
    std::vector<float> ConvolveHists(std::vector<float> Hist1, std::vector<float> Hist2);
    void ResidualRemove(TH1F* Hist);
};

#endif
