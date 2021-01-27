#ifndef EXPERIMENTAL_H
#define EXPERIMENTAL_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/RooFitFunctions.h>
#include<PostAnalysis/Plotting.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth, int trk_Data); 
void AlgorithmMonteCarlo(); 
void Shifting(TH1F* H); 

#endif
