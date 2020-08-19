#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DerivedFunctions.h>

#ifndef UNITTEST_H
#define UNITTEST_H

class BaseFunctionTest
{
  public:
    void Subtraction();
    void NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max);  
    void Convolve(TH1F* Hist1, TH1F* Hist2, TH1F* Expection);
    void Deconvolve(TH1F* Trk2, TH1F* Trk1, float offset, int iter); 
};

class DerivedFunctionTest
{
  public:
    void NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max);


};

#endif
