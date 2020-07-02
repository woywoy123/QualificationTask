// Include all the required ROOT libraries
#include<TString.h>
#include<TFile.h>
#include<TH1F.h>
#include<TImage.h>
#include<TApplication.h>
#include<TSystem.h>

// Include the standard C++ library 
#include<iostream>

// Include RooFit classes 
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooHistPdf.h>
#include<RooFit.h>
#include<RooAddPdf.h>

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

class Functions
{
  public:
    std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, int min, int max);
    void FillTH1F_From_File(std::vector<TH1F*> Histograms, TFile* File, TString DetectorLayer);
};

class Fitting
{
  public:
    std::vector<RooDataHist*> ConvertTH1FtoDataHist(std::vector<TH1F*> Histograms, RooRealVar* domain);
    std::vector<RooHistPdf*> ConvertTH1FtoPDF(std::vector<TH1F*> Histograms, RooRealVar* domain);
    std::vector<RooRealVar*> GenerateVariables(std::vector<TString>, std::vector<double> Begin, std::vector<double> End); 
    
    RooHistPdf* ConvertTH1FtoPDF(RooDataHist* Histogram, TString Name, RooRealVar* domain);
    RooDataHist* ConvertTH1FtoDataHist(TH1F* Hist, RooRealVar* domain);
    RooRealVar* GenerateVariable(TString name, double begin, double end);   
  
    RooArgList VectorToArgList(std::vector<RooRealVar*> Vector);
    RooArgList VectorToArgList(std::vector<RooHistPdf*> Vector);
   
};

namespace Constants
{
  std::vector<TString> Detector = {"IBL", "Blayer", "layer1", "layer2"};
  std::vector<TString> Pure_Names = {"dEdx_ntrk_1_ntru_1", 
                                     "dEdx_ntrk_2_ntru_2", 
                                     "dEdx_ntrk_3_ntru_3", 
                                     "dEdx_ntrk_4_ntru_4"};

  std::vector<TString> Variable_Names = {"ntrk_1", "ntrk_2", "ntrk_3", "ntrk_4"};
  std::vector<double> Begin = {1e6, 0., 0., 0.};
  std::vector<double> End = {1e8, 1e1, 1e1, 1e1};
}

#endif
