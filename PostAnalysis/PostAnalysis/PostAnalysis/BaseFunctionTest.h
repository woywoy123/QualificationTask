#ifndef BASEFUNCTIONTEST_H
#define BASEFUNCTIONTEST_H


#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>

void PlotLandau(); 
void PlotGaussian(); 
void PlotGaussianXGaussian(); 
void PlotLandauXLandau(); 
void PlotLandauXGaussian();
void PlotDeconvGausXGaus(); 
void PlotDeconvLandauXLandau(); 
void PlotDeconvLandauXGaussian(); 

#endif
