#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

#include<TString.h>
#include<TH1F.h>
#include<iostream> 
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/DistributionGenerator.h>

// Histogram Generators
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");
std::vector<TH1F*> CloneTH1F(TH1F* Hist, std::vector<TString> Names); 
TH1F* SumHists(std::vector<TH1F*> Hists_V, TString name);

// Histogram manipulation 
void Normalize(TH1F* Hist); 
void Normalize(std::vector<TH1F*> Hists); 
void Shift(TH1F* Hist, int shift);

// Variable Name Generator 
std::vector<TString> NameGenerator(int number, TString shorty); 
std::vector<TString> NameGenerator(std::vector<TH1F*> Hists, TString append); 

// vector normalization 
std::vector<float> Normalize(std::vector<float> V1); 

// Convert TH1F to vector
std::vector<float> ToVector(TH1F* Hist); 
void ToTH1F(std::vector<float> Input, TH1F* Hist); 

// Simple repetitive functions 
float GetMaxValue(TH1F* H); 
void BulkWrite(std::vector<TH1F*> Hist_V); 
void BulkDelete(std::vector<TH1F*> Hist_V); 
void BulkDelete(std::map<TString, std::vector<TH1F*>> Hists); 
std::vector<TH1F*> BulkClone(std::vector<TH1F*> Hists, std::vector<TString> Name);


// Others
void Flush(std::vector<TH1F*> F_C, std::vector<TH1F*> ntrk_Conv, bool DeletePointer = false);
void Average(TH1F* Data); 
void Average(std::vector<TH1F*> Data); 

// Cout text based stuff 
void CoutText(TString *Input, int v, TString Text = "");
TString PrecisionString(float number, int precision, bool sci = false);

// Iterators that are used throughout the package 
typedef std::map<TString, std::map<TString, std::vector<TH1F*>>>::iterator MMVi; 
typedef std::map<TString, std::vector<TH1F*>>::iterator MVi; 
typedef std::map<TString, std::vector<float>>::iterator MVFi; 
typedef std::map<TString, std::map<TString, std::vector<float>>>::iterator MMVFi; 



#endif
