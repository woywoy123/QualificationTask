#ifndef ALGORITHM_FUNCTIONS_H
#define ALGORITHM_FUNCTIONS_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/IO.h>

std::vector<std::vector<TH1F*>> Normalization_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_Start, std::map<TString, std::vector<float>> Params, TString JE); 
std::vector<std::vector<TH1F*>> NormalizationShift_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE);
std::vector<std::vector<TH1F*>> NormalizationShiftFFT_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE);
std::vector<std::vector<TH1F*>> NormalizationShiftWidthFFT_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE);
std::vector<std::vector<TH1F*>> Simultaneous_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE); 
std::vector<std::vector<TH1F*>> IncrementalFit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE); 

void Experimental(std::vector<TH1F*> Data, std::vector<std::vector<TH1F*>> ntrk_mtru, std::map<TString, std::vector<float>> Params);
std::vector<std::vector<TH1F*>> BuildNtrkMtru(int n, TH1F* trk1_start, TString extension, int tru = 4);



#endif 
