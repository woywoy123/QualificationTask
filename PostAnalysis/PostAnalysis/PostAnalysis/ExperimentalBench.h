#ifndef EXPERIMENTALBENCH_H
#define EXPERIMENTALBENCH_H
#include<TString.h>

void IOTest(); 
void RooFitBaseFunctionTest(); 
void TestFits_NTruth_NTrack(); 
void TestRead(); 
void TestFits_AllTruth_ToTrack(TString JE = "", TString Mode = ""); 

#endif
