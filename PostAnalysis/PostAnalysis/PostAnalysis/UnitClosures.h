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
    void TestSubtraction(TH1F* Data, int trk, std::vector<TH1F*> PDFs, float min, float max, std::vector<float> Closure); 

    // Need to write test units for:
    // - GaussianUnfold (with data being convolved with gaussian)
    // - Getting MC data for closure

};

class Presentation
{
  public:
    void Threshold(TString DataDir);
    void TestMinimalAlgorithm(std::vector<TH1F*> Data, float min, float max, float offset, std::vector<TH1F*> Pure, std::vector<std::vector<float>> Closure);
    void TestGaussianAlgorithm(std::vector<TH1F*> Data, float min, float max, float offset, std::vector<TH1F*> Pure, std::vector<std::vector<float>> Closure);

};

class Debug
{
  public: 

};


#endif
