
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
    void Debug(std::vector<TH1F*> Hist, std::vector<float> Params);
    void UnitTesting();


};
#endif 
