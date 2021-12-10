#ifndef DISTRIBUTIONGENERATOR_H
#define DISTRIBUTIONGENERATOR_H

#include<TF1.h>
#include<iostream>
#include<TH1F.h>

TH1F* Landau(float mean, float stdev, int bins, float min, float max, TString Extension = "");
TH1F* Gaussian(float mean, float stdev, int bins, float min, float max, TString Extension = "");

#endif
