#ifndef DISTRIBUTIONGENERATOR_H
#define DISTRIBUTIONGENERATOR_H

#include<PostAnalysis/BaseFunctions.h>
#include<TF1.h>

void Landau(std::vector<TH1F*> Hists, std::vector<float> COMP, std::vector<float> Parameters, int Number, float min, float max);
void Gaussian(float mean, float stdev, TH1F* Hist); 



#endif
