#ifndef EVALUATE_H
#define EVALUATE_H
#include "BaseFunctions.h"
#include "IO.h"
#include "Plotting.h"

void ShapePerformance(TString dir, TString Mode = "Test");
MMVTH1F FLostPerformance(TString dir, TString Mode = "Test", bool Plotting = true);
void WriteFLostFile(TString Input, TString Output); 
void ErrorFlostAlgo(TString dir); 
void BestMinimizerAlgo(TString dir);
void Clusters(TString dir); 
void TemplateVariants(TString dir); 
void Debug_Cases(TString dir);

static std::vector<TString> Ag = {"Normalization", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT", "Incremental", "Experimental"};
#endif
