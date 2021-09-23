#ifndef IO_H
#define IO_H

#include<TFile.h>
#include<TKey.h>
#include<TTree.h>
#include<TH1F.h>
#include<TString.h>
#include<iostream>

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

typedef std::map<TString, std::map<TString, std::map<TString, std::map<TString, std::map<TString, float>>>>> MMMMMF; 


typedef std::map<TString, std::map<TString, std::vector<TH1F*>>> MMTF; 
typedef std::map<TString, std::vector<TH1F*>> MTF;

typedef std::map<TString, std::map<TString, std::vector<TH1F*>>>::iterator MMTFi; 
typedef std::map<TString, std::vector<TH1F*>>::iterator MTFi;

typedef std::vector<TH1F*> VTF;
typedef std::vector<float> VF;
typedef std::vector<VTF> VVTF; 
typedef std::vector<std::vector<float>> VVF; 


const std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
const std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 
const std::vector<TString> Algos = {"Normalization", "ShiftNormalWidthFFT", "ShiftNormalFFT", "Incremental", "ShiftNormal"};
std::map<TString, std::map<TString, std::vector<TH1F*>>> ReadAlgorithmResults(TString dir); 

VF ShapeComparison(VTF Hists, VTF Truth); 
VF IntegralComparison(VTF Hists, VTF Truth); 
VF AdjustedShape(VTF Hists, VTF Truth); 
VF TruthIntegrals(VTF Truth); 


#endif 

