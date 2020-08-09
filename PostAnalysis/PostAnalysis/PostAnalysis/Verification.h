
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
    void UnitTesting();
    void MainAlgorithm(std::vector<TH1F*> Data, TH1F* Target, std::vector<TH1F*> Closure);
    void MainGaussianUnfolding(std::vector<TH1F*> Data, TH1F* Target, std::vector<TH1F*> Closure);
    void Debug(TH1F* trk1, TH1F* trk2); 
    void CalibrationDataConvolution();
    void NewLRTesting(TH1F* trk);
};
#endif 
