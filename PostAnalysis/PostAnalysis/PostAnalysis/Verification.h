
#ifndef VERIFCATION_H
#define VERIFCATION_H

#include<RooAddPdf.h>
#include<RooHistPdf.h>
#include<RooFit.h>
#include<iostream>

using namespace RooFit;

class Verification
{
  public: 
    void RecoverScaling(RooAddPdf model, 
                                  std::vector<TH1F*> Histograms, 
                                  std::vector<RooHistPdf*> PDF, 
                                  RooRealVar* range, 
                                  std::vector<RooRealVar*> Variables,
                                  float S1, float S2);

    void Subtraction(std::vector<TH1F*> Histograms, 
                                   RooRealVar* range, 
                                   RooAddPdf model, 
                                   std::vector<RooHistPdf*> PDF, 
                                   std::vector<RooRealVar*> Variables);

    void Reconstruction(std::vector<TH1F*> trk1, 
                                  std::vector<TH1F*> trk2,
                                  std::vector<TH1F*> trk3,
                                  std::vector<TH1F*> trk4,
                                  std::vector<RooRealVar*> Variables, 
                                  std::vector<RooHistPdf*> PDFs,
                                  RooAddPdf model, 
                                  RooRealVar* range);
    
    void FLostLayer(std::map<TString, std::vector<TH1F*>> Layer_Track, 
                                  float lower, float upper);

    void NormalizedSubtraction(float lower, float upper, std::vector<TH1F*> Hist1, std::vector<TH1F*> Hist2, RooAddPdf model, std::vector<RooRealVar*> Variables, RooRealVar* Range);

    void Debug(std::vector<TH1F*> Hist, std::vector<std::vector<float>> Params);


};
#endif 
