#ifndef UNITCLOSURES_H
#define UNITCLOSURES_H

#include<RooAddPdf.h>
#include<RooHistPdf.h>
#include<RooFit.h>
#include<iostream>

using namespace RooFit;

class UnitClosures
{
  public:
    void TestFit(std::vector<TH1F*> PDF, std::vector<TH1F*> Data, float min, float max, std::vector<std::vector<float>> Closure);
    void TestTailAndDeconv(TH1F* trk1, TH1F* trk2, int iter, float min, float max);
    void TestDeconvolution(TH1F* h1, TH1F* PSF, int iter);
    
    // Need to write test units for:
    // - Subtraction 
    // - MainAlgorithm
    // - GaussianUnfold (with data being convolved with gaussian)
    // - Getting MC data for closure

};

class Debug
{
  public: 
    void NewLRTesting(TH1F* trk);
};


#endif
