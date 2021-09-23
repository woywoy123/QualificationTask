#ifndef EVALUATE_H
#define EVALUATE_H
#include "BaseFunctions.h"
#include "IO.h"
#include "Plotting.h"

void ShapePerformance(TString dir, TString Mode = "Test");
MMVTH1F FLostPerformance(TString dir, TString Mode = "Test", bool Plotting = true);
void WriteFLostFile(TString Input, TString Output); 
void ErrorFlostAlgo(TString dir); 

#endif
