#ifndef IO_H
#define IO_H
#include<PostAnalysis/BaseFunctions.h>
#include<TFile.h>
#include<TKey.h>
#include<TTree.h>

const std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
const std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 

std::map<TString, std::map<TString, std::vector<TH1F*>>> ReadCTIDE(TString dir);
void WriteHistsToFile(std::vector<TH1F*> ntrk_ntru, TString dir); 
void WriteOutputMapToFile(std::map<TString, std::vector<float>> Map, TString dir, TString name);


#endif 
