#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DerivedFunctions.h>

#ifndef UNITTEST_H
#define UNITTEST_H

class BaseFunctionTest
{
  public:
    void Subtraction();
    void NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max);  
};

class DerivedFunctionTest
{
  public:
    void NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max);


};

#endif
