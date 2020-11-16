#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include<TVirtualFFT.h>
#include<TComplex.h>
#include<TH1F.h>
#include<iostream>
#include<PostAnalysis/BaseFunctions.h>

// Delete after
#include<TCanvas.h>


// Convolution 
std::vector<float> ConvolutionFFT(std::vector<float> V1, std::vector<float> V2, int ZeroPointBin); 
void Convolution(TH1F* H1, TH1F* H2, TH1F* Out); 
void ResidualRemove(TH1F* Hist); 

// Deconvolution 
std::vector<float> LucyRichardson(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y); 
void Deconvolution(TH1F* PDF, TH1F* PSF, TH1F* Output, int Max_Iter); 

#endif
