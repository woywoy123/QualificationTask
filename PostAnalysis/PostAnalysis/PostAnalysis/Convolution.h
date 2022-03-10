#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include<TH1F.h>
#include<iostream>
#include<PostAnalysis/BaseFunctions.h>


// Convolution 
void Convolution(TH1F* H1, TH1F* H2, TH1F* Out); 
std::vector<float> ConvolutionSum(std::vector<float> V1, std::vector<float> V2, int ZeroPointBin); 
std::vector<TH1F*> ConvolveNTimes(TH1F* Start, int n, std::vector<TString> base, TString extension); 

#endif
