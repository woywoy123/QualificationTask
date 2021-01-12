#ifndef ROOFITFUNCTIONS_H
#define ROOFITFUNCTIONS_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/Plotting.h>
#include<RooRealVar.h>
#include<RooGaussian.h>
#include<RooFFTConvPdf.h>
#include<RooDataHist.h>
#include<RooAddPdf.h>
#include<PostAnalysis/DistributionGenerator.h>

// Constraint classes 
#include<RooConstVar.h>




std::vector<TH1F*> FitDeconvolution(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache = 0, int cache = 0); 

std::vector<std::pair<TH1F*, std::vector<float>>> FitDeconvolutionPerformance(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache = 0, int cache = 0); 

#endif
