#ifndef BASEFUNCTIONS_H
#define BASEFUNCTIONS_H

#include<TFile.h>
#include<TGraph.h>
#include<TKey.h>
#include<TTree.h>
#include<TH1F.h>
#include<TString.h>
#include<iostream>
#include<sstream>
#include<fstream>

const std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
const std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 
const std::vector<TString> Algos = {"Normalization_TRUTH", "ShiftNormal_TRUTH", "ShiftNormalFFT_TRUTH", "ShiftNormalWidthFFT_TRUTH", "Incremental_TRUTH", 
                                    "Normalization", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT", "Incremental"};



typedef std::map<TString, std::map<TString, std::vector<TH1F*>>> MMVTH1F; 
typedef MMVTH1F::iterator MMVTH1Fi; 

typedef std::map<TString, std::vector<TH1F*>> MVTH1F; 
typedef MVTH1F::iterator MVTH1Fi;

typedef std::vector<TH1F*> VTH1F; 
typedef std::vector<float> VF;

float ChiSquareDistanceSym(TH1F* H1, TH1F* H2);
VTH1F GetSubVector(std::vector<TH1F*> Input, int start, int end);
VF ShapeComparison(VTH1F Hists, VTH1F Truth); 
float WeightedComparisonToTruth(VTH1F Hists, VTH1F Truth); 

typedef std::map<TString, std::vector<std::vector<float>>> MVVF; 
typedef std::map<TString, std::vector<std::vector<float>>>::iterator MVVFi; 

typedef std::map<TString, float> MF; 
typedef std::map<TString, float>::iterator MFi; 

typedef std::map<TString, std::vector<float>> MVF; 


typedef std::map<TString, std::map<TString, std::map<TString, float>>> MMMF; 
typedef std::map<TString, std::map<TString, std::map<TString, float>>>::iterator MMMFi; 

typedef std::map<TString, std::map<TString, float>> MMF; 
typedef std::map<TString, std::map<TString, float>>::iterator MMFi; 

typedef std::map<TString, std::map<TString, std::map<TString, std::map<TString, float>>>> MMMMF; 
typedef std::map<TString, std::map<TString, std::map<TString, std::map<TString, float>>>>::iterator MMMMFi; 


typedef std::map<TString, std::map<TString, std::map<TString, std::map<TString, std::map<TString, float>>>>> MMMMMF; 
typedef std::map<TString, std::map<TString, std::map<TString, std::map<TString, std::map<TString, float>>>>>::iterator MMMMMFi; 

typedef std::map<TString, std::map<TString, std::vector<TH1F*>>> MMTF; 
typedef std::map<TString, std::vector<TH1F*>> MTF;

typedef std::map<TString, std::map<TString, std::vector<TH1F*>>>::iterator MMTFi; 
typedef std::map<TString, std::vector<TH1F*>>::iterator MTFi;

typedef std::vector<TH1F*> VTF;
typedef std::vector<float> VF;
typedef std::vector<VTF> VVTF; 
typedef std::vector<std::vector<float>> VVF; 

typedef std::map<TString, TString> MT; 


static void BulkDelete(std::vector<TGraph*> graph)
{
  for (int i(0); i < graph.size(); i++){ delete graph[i]; }
}

float Flost2(std::vector<std::vector<TH1F*>> ntrk); 
float Flost3(std::vector<std::vector<TH1F*>> ntrk); 
MVTH1F PrintFlost(MMMF FLostMap, TString Layer, TString dir); 
TString PrecisionString(float number, int precision, bool sci); 
void CoutText(TString *Input, int v, TString Text); 
void BulkWrite(std::vector<TH1F*> H); 
MF ReadTH1F(TH1F* H); 


#endif 
