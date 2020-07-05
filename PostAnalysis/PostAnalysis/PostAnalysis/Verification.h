
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

};
#endif 
