#include<TString.h>
#include<TH1F.h>
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooHistPdf.h>
#include<RooAddPdf.h>

#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

class BaseFunctions
{
  public: 
    std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");
    std::vector<float> Ratio(std::vector<TH1F*> Hists, TH1F* Data); 
    std::vector<float> Ratio(std::vector<RooRealVar*>, TH1F* Data);
    std::vector<float> ClosureAndData(std::vector<TH1F*> Hists, TH1F* Data); 
    void Normalize(TH1F* Hist);
    void Normalize(std::vector<TH1F*> Hist);
    void Subtraction(std::vector<TH1F*> ntrk, TH1F* Data, int Exclude, std::vector<float> Ratios);  
    std::vector<RooRealVar*> RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End);
    std::vector<RooDataHist*> RooData(std::vector<TH1F*> Hist, RooRealVar* Domain); 
    RooDataHist* RooData(TH1F* Hist, RooRealVar* Domain);
    std::vector<RooHistPdf*> RooPDF(std::vector<TH1F*> Hist, RooRealVar* Domain); 
    RooArgList RooList(std::vector<RooRealVar*> Vector);
    RooArgList RooList(std::vector<RooHistPdf*> Vector);
  
    float ChiSquare(std::vector<float> V1, std::vector<float> V2);
    void PredictionTruthPrint(std::vector<float> Truth, std::vector<float> Prediction);

};

#endif
