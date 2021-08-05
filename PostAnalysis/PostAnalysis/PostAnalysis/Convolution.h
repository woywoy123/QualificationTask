#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include<TVirtualFFT.h>
#include<TComplex.h>
#include<TH1F.h>
#include<iostream>
#include<PostAnalysis/BaseFunctions.h>
#include<thread>



// Convolution 
void Convolution(TH1F* H1, TH1F* H2, TH1F* Out); 
std::vector<float> ConvolutionSum(std::vector<float> V1, std::vector<float> V2, int ZeroPointBin); 
std::vector<TH1F*> ConvolveNTimes(TH1F* Start, int n, std::vector<TString> base, TString extension); 

// Deconvolution 
std::vector<float> LucyRichardson(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y); 
void Deconvolution(TH1F* PDF, TH1F* PSF, TH1F* Output, int Max_Iter); 
void MultiThreadingDeconvolution(std::vector<TH1F*> Data, std::vector<TH1F*> PSF, std::vector<TH1F*> Result, int Iter); 

void Smooth1Trk(TH1F* trk1_start, int iter); 

#endif
