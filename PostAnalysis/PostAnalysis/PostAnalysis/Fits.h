#ifndef FITS_H
#define FITS_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/RooFitFunctions.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/IO.h>
#include<iostream>
#include<fstream>

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> OnlyNormal(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name); 
std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShift(std::vector<TH1F*> Inside_Using, TH1F* Outside, std::map<TString, std::vector<float>> Params, TString Name);
std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> ShiftNormalFFT(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name);
std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalWidthDeconvShiftFFT(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name);
void AnalysisLoop(); 
void Evaluation(); 

#endif
