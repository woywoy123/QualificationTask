#ifndef STATISTICS_H
#define STATISTICS_H

#include<iostream>
#include<TH1F.h>
#include<TString.h>
#include<PostAnalysis/BaseFunctions.h>

// Benchmarking 
float ErrorByIntegral(TH1F* Hist1, TH1F* Hist2, float x_min = 0, float x_max = 0); 
std::vector<float> Flost2(std::vector<std::vector<TH1F*>> ntrk, std::vector<TH1F*> Err); 
std::vector<float> Flost3(std::vector<std::vector<TH1F*>> ntrk, std::vector<TH1F*> Err); 

#endif
