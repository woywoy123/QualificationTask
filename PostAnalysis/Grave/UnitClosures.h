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
    void TestFit(std::vector<std::vector<TH1F*>> PDF, std::vector<TH1F*> Data, float min, float max, std::vector<std::vector<float>> Closure);
    void TestTailAndDeconv(TH1F* trk1, TH1F* trk2, int iter, float min, float max);
    void TestDeconvolution(TH1F* h1, TH1F* PSF, int iter);
    void TestSubtraction(TH1F* Data, int trk, std::vector<TH1F*> PDFs, float min, float max, std::vector<float> Closure); 

  private:
    void ClosureBaseFit(std::vector<TH1F*> PDFs, TH1F* d_trk, std::vector<float> Clos, float min, float max);

};

class Presentation
{
  public:
    void Threshold(TString DataDir);
};

class DataGeneration
{
  public:
    void MonteCarlo(std::vector<TH1F*> Hists, TString dir, TString Extension = "");
    void IdealLandau(std::vector<TH1F*> Hists, std::vector<float> COMP, std::vector<float> LandauParams);
    std::vector<float> MergeToys(std::vector<TH1F*> Hists, TH1F* trk);    
};

class Debug
{
  public: 

};


#endif
