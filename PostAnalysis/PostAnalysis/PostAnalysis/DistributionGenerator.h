#ifndef DISTRIBUTIONGENERATOR_H
#define DISTRIBUTIONGENERATOR_H

#include<TF1.h>
#include<iostream>
#include<TH1F.h>

std::vector<TH1F*> Landau(std::vector<TString> Hists, std::vector<float> COMP, std::vector<float> Parameters, int Number, int bins, float min, float max);
std::vector<TH1F*> WrongLandau(std::vector<TString> Hists, std::vector<float> COMP, std::vector<float> Parameters, int Number, int bins, float min, float max);
TH1F* Gaussian(float mean, float stdev, int bins, float min, float max, TString Extension = "");

#endif
