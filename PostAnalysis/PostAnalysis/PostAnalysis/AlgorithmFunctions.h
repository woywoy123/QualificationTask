#ifndef ALGORITHM_FUNCTIONS_H
#define ALGORITHM_FUNCTIONS_H

#include<PostAnalysis/BaseFunctions.h>
typedef std::map<TString, std::vector<float>>::iterator it; 

std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params); 
std::vector<std::pair<TH1F*, std::vector<float>>> LoopGenAll(std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params); 

#endif 
