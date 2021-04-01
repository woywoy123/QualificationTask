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

void Evaluate_Fits(TString Filename); 



#endif
