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
#include<TF1.h>
#include<RooArgList.h>
#include<RooNumConvPdf.h>
#include<RooConstVar.h>
#include<RooProdPdf.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<RooFormulaVar.h>
#include<RooFitResult.h> 

std::vector<TH1F*> FitDeconvolution(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache = 0, int cache = 0); 
std::vector<std::pair<TH1F*, std::vector<float>>> FitDeconvolutionPerformance(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache = 0, int cache = 0); 
std::vector<std::pair<TH1F*, std::vector<float>>> FitWithConstraint(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache, int cache);
TH1F* ExplicitConstraining(TH1F* Data, TH1F* PDF, std::map<TString, std::vector<float>> Params); 
TH1F* ExplicitConstrainingExternal(TH1F* Data, TH1F* PDF, std::map<TString, std::vector<float>> Params); 
std::vector<TH1F*> IterativeFitting(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache, int cache); 
std::map<TString, std::vector<float>> ScalingFit(TH1F* Data, std::vector<TH1F*> PDF_H); 
std::map<TString, std::vector<float>> ScalingShift(TH1F* Data, std::vector<TH1F*> ntrk); 
#endif
