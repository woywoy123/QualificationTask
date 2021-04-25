#ifndef ALGORITHM_FUNCTIONS_H
#define ALGORITHM_FUNCTIONS_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/IO.h>

std::vector<TH1F*> Normalization_Fit(std::vector<TH1F*> Data, TH1F* trk1_Start, std::map<TString, std::vector<float>> Params, TString JE); 
std::vector<TH1F*> NormalizationShift_Fit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE);
std::vector<TH1F*> NormalizationShiftFFT_Fit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE);
std::vector<TH1F*> NormalizationShiftWidthFFT_Fit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE);

#endif 
