#ifndef DISTRIBUTIONGENERATOR_H
#define DISTRIBUTIONGENERATOR_H

#include<TF1.h>
#include<iostream>
#include<TH1F.h>

void Landau(std::vector<TH1F*> Hists, std::vector<float> COMP, std::vector<float> Parameters, int Number, float min, float max);
TH1F* Gaussian(float mean, float stdev, int bins, float min, float max, TString Extension = "");

#endif
