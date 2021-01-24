#ifndef ALGORITHM_FUNCTIONS_H
#define ALGORITHM_FUNCTIONS_H

#include<PostAnalysis/BaseFunctions.h>

void ShapeSigmoid(TH1F* trk_Fit, TH1F* ntrk_Conv); 
void Scale(TH1F* Data, std::vector<TH1F*> ntrk); 
void ScaleShape(TH1F* Data, std::vector<TH1F*> ntrk); 
void Flush(std::vector<TH1F*> F_C, std::vector<TH1F*> ntrk_Conv, bool sig = true); 
void Average(TH1F* Data); 
std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, std::vector<TH1F*> PSF, TH1F* Data, std::map<TString, std::vector<float>> Params); 
float FitError(TH1F* Data, std::vector<TH1F*> ntrk);




std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, int trk_Data);

#endif 
