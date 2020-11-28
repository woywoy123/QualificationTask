#ifndef BASEFUNCTIONTEST_H
#define BASEFUNCTIONTEST_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/RooFitFunctions.h>
#include<TFile.h>

void TestLandau(TFile* F); 
void TestGaussian(TFile* F); 
void TestGaussianXGaussian(TFile* F); 
void TestLandauXLandau(TFile* F); 
void TestLandauXGaussian(TFile* F);
void TestDeconvGausXGaus(TFile* F); 
void TestDeconvLandauXLandau(TFile* F); 
void TestDeconvLandauXGaussian(TFile* F); 
void TestLandauXGausFit(TFile* F); 
void TestNLandauXNGausFit(TFile* F); 
void TestGaussianDeconvolutionFit(TFile* F); 
void TestDeconvolutionFit(TFile* F); 
void TestComparisonBinCenteringLandauXLandau(TFile* F); 
void TestOscillationLucyRichardson(TFile* F); 
void TestAlgorithm(TFile* F); 
void TestReadFile(TFile* F); 
void TestMonteCarloMatchConvolution(TFile* F); 
#endif
