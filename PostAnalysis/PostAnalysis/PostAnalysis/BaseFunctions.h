#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

#include<TRandom2.h>
#include<TF1.h>
#include<TMatrixD.h>
#include<TMatrixDSym.h>
#include<TVectorD.h>

#include<TString.h>
#include<TH1F.h>
#include<iostream> 
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<TFile.h>
#include<TGraph.h>
#include<TGraphSmooth.h>

// Histogram Generators
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension = "");
std::vector<TH1F*> CloneTH1F(TH1F* Hist, std::vector<TString> Names); 
TH1F* SumHists(std::vector<TH1F*> Hists_V, TString name);

// Histogram manipulation 
void Normalize(TH1F* Hist); 
void Normalize(std::vector<TH1F*> Hists); 
void Shift(TH1F* Hist, int shift);
void MatchBins(std::vector<TH1F*> In, TH1F* Data); 
void SubtractData(std::vector<TH1F*> In, TH1F* Data, int trk, bool trutrk = false); 
void Smooth(TH1F* Hist, float kern); 
void SmoothHist(TH1F* Hist, int order, float sigma); 
TH1F* ExpandTH1F(TH1F* Hist, float min, float max); 
std::vector<TH1F*> ExpandTH1F(std::vector<TH1F*> Hists, float min, float max); 

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

// Similarity metric
float ChiSquareError(TH1F* Truth, TH1F* Pred); 



// Iterators that are used throughout the package 
typedef std::map<TString, std::map<TString, std::map<TString, std::vector<float>>>>::iterator MMMVFi;
typedef std::map<TString, std::map<TString, std::map<TString, int>>>::iterator MMMSIi; 
typedef std::map<TString, std::map<TString, std::vector<TH1F*>>>::iterator MMVi; 
typedef std::map<TString, std::map<TString, std::vector<float>>>::iterator MMVFi; 
typedef std::map<TString, std::map<TString, std::map<int, int>>>::iterator MMIIi; 
typedef std::map<TString, std::vector<std::vector<TH1F*>>> MVVi; 
typedef std::map<TString, std::vector<float>>::iterator MVFi; 
typedef std::map<TString, std::vector<TH1F*>>::iterator MVi; 
typedef std::map<TString, std::map<int, int>>::iterator MIIi; 
typedef std::map<TString, std::map<TString, int>>::iterator MMSIi; 
typedef std::map<TString, std::map<TString, float>>::iterator MMFi; 
typedef std::map<TString, float>::iterator MFi;
typedef std::map<TString, int>::iterator MIi; 
typedef std::map<int, int>::iterator IIi; 

// Non iterators used in this package
typedef std::map<TString, std::map<TString, std::map<TString, std::vector<float>>>> MMMVF;
typedef std::map<TString, std::map<TString, std::map<TString, int>>> MMMSI; 
typedef std::map<TString, std::map<TString, std::vector<TH1F*>>> MMVT; 
typedef std::map<TString, std::map<TString, std::vector<float>>> MMVF;

typedef std::map<TString, std::map<TString, std::map<int, int>>> MMII; 

typedef std::map<TString, std::map<TString, int>> MMSI; 
typedef std::map<TString, std::map<int, int>> MII; 
typedef std::map<TString, std::vector<float>> MVF; 
typedef std::map<TString, std::vector<TH1F*>> MVT; 

typedef std::vector<std::vector<TH1F*>> VVT; 
typedef std::vector<std::vector<float>> VVF; 

typedef std::map<TString, int> MI;
typedef std::map<int, int> II; 
typedef std::map<TString, float> MF; 
typedef std::vector<TH1F*> VT; 
typedef std::vector<float> VF;
typedef std::vector<TString> VS; 
#endif
