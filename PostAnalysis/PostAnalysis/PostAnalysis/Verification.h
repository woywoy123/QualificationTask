
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
    void MainAlgorithm(std::vector<TH1F*> Data, TH1F* Target);
};
#endif 
