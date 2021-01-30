#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

#include<TString.h>
#include<TH1F.h>
#include<iostream> 
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/RooFitFunctions.h>

// Histogram Generators
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");
std::vector<TH1F*> CloneTH1F(TH1F* Hist, std::vector<TString> Names); 
TH1F* SumHists(std::vector<TH1F*> Hists_V, TString name);

// Histogram manipulation 
void Normalize(TH1F* Hist); 
void Normalize(std::vector<TH1F*> Hists); 
void Shift(TH1F* Hist, int shift);

// Variable Name Generator 
std::vector<TString> NameGenerator(int number, TString shorty); 

// vector normalization 
std::vector<float> Normalize(std::vector<float> V1); 

// Convert TH1F to vector
std::vector<float> ToVector(TH1F* Hist); 
void ToTH1F(std::vector<float> Input, TH1F* Hist); 

// Benchmarking 
float Pythagoras(std::vector<float> v1, std::vector<float> v2); 
float SquareError(TH1F* Hist1, TH1F* Hist2, float x_min = 0, float x_max = 0);  
float ErrorByIntegral(TH1F* Hist1, TH1F* Hist2, float x_min = 0, float x_max = 0); 
void Statistics(TH1F* H1, TH1F* H2, float x_min = 0, float x_max = 0); 

// Simple repetitive functions 
float GetMaxValue(TH1F* H); 
void BulkWrite(std::vector<TH1F*> Hist_V); 
void BulkDelete(std::vector<TH1F*> Hist_V); 

#endif
