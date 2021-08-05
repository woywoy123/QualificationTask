#ifndef SIMPLEEXAMPLEFITS_H
#define SIMPLEEXAMPLEFITS_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/Evaluate.h>

void FitTemplateToTruth( std::vector<std::vector<TH1F*>> Truth, TH1F* trk1_start, std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, TString Mode, TString JE); 



#endif
