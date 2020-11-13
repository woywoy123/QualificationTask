#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

#include<TString.h>
#include<TH1F.h>
#include<iostream> 
#include<TVirtualFFT.h>
#include<TComplex.h>

// Histogram Generators
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");


// Histogram manipulation 
void Normalize(TH1F* Hist); 
void Normalize(std::vector<TH1F*> Hists); 
void Shift(TH1F* Hist, int shift); 

// Convolution - Deconvolution 
void ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* H1xH2);
void Convolution(TH1F* H1, TH1F* H2, TH1F* Out); 
void ArtifactRemove(TH1F* Hist); 
void ArtifactRemove(std::vector<TH1F*> Hists); 






















#endif
