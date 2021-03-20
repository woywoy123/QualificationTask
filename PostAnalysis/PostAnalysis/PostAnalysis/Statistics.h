#ifndef STATISTICS_H
#define STATISTICS_H

#include<iostream>
#include<TH1F.h>
#include<TString.h>

// Benchmarking 
float Pythagoras(std::vector<float> v1, std::vector<float> v2); 
float SquareError(TH1F* Hist1, TH1F* Hist2, float x_min = 0, float x_max = 0);  
float ErrorByIntegral(TH1F* Hist1, TH1F* Hist2, float x_min = 0, float x_max = 0); 
void Statistics(TH1F* H1, TH1F* H2, float x_min = 0, float x_max = 0); 


#endif
