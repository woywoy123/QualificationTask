#ifndef EVALUATE_FITS_H
#define EVALUATE_FITS_H

#include<iostream>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/Statistics.h>

typedef std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>>::iterator it; 
typedef std::map<TString, std::map<TString, std::vector<TH1F*>>>::iterator ip; 
typedef std::map<TString, std::vector<TH1F*>>::iterator iz; 
typedef std::map<TString, std::vector<float>>::iterator ix; 
typedef std::map<TString, std::vector<std::pair<TString, std::map<TString, std::vector<float>>>>>::iterator xi; 
typedef std::map<TString, std::vector<int>>::iterator ui; 
void Evaluate_Fits(TString Filename); 
void EvaluateErrorImpact(TString Filename); 
void Evaluate_nTrackFits(TString Filename); 

#endif
