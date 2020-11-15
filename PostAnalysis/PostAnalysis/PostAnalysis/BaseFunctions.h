#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

#include<TString.h>
#include<TH1F.h>
#include<iostream> 

// Histogram Generators
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");

// Histogram manipulation 
void Normalize(TH1F* Hist); 
void Normalize(std::vector<TH1F*> Hists); 
void Shift(TH1F* Hist, int shift); 

// vector normalization 
std::vector<float> Normalize(std::vector<float> V1); 

// Translate bin to domain range
float BinToDomain(int Bin_Number, int bins, float min, float max);
#endif
