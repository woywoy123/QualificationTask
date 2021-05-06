#ifndef EVALUATE_H
#define EVALUATE_H

#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/Statistics.h>

void CompareToTruth(TString dir);
void CompileCout(std::vector<TString> JetEnergy, std::vector<TString> Algo_Strings, std::vector<std::map<TString, std::map<TString, std::vector<float>>>> Collector); 
void MultiTrackTruthComparison(TString dir); 



#endif
