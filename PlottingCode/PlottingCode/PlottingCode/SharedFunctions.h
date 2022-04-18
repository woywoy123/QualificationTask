#ifndef SHARED_FUNCTIONS_H
#define SHARED_FUNCTIONS_H

#include<iostream>
#include<sstream>
#include<fstream>

#include<TString.h>
#include<TFile.h>
#include<TH1F.h>
#include<TGraph.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TPad.h>


typedef TString _TS;
typedef std::map<_TS, std::map<_TS, std::vector<TH1F*>>>::iterator MMVTFi; 
typedef std::map<_TS, std::vector<TH1F*>>::iterator MVTFi; 
typedef std::map<_TS, std::vector<float>>::iterator MVFi; 

const std::vector<_TS> Layer_H = {"IBL", "Blayer", "layer1", "layer2"}; 
const std::vector<_TS> Energy_H = {"200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                        "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV",
                                        "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV"}; 
const std::vector<_TS> Algos_H = {"Normalization", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT", "Incremental", "Experimental"}; //, "Simultaneous"};

static std::vector<_TS> Layer_Energy_H()
{
  std::vector<_TS> Output; 
  for (_TS L : Layer_H)
  {
    for (_TS E : Energy_H)
    {
      Output.push_back(L + "-" + E);  
    }
  }
  return Output;
};

const std::vector<Color_t> Colors_H = {kRed, kGreen, kBlue, kCyan, kViolet, kOrange, kCoffee, kAurora}; 


std::vector<_TS> Split(_TS Input, _TS Sub); 
void BulkDelete(std::vector<TGraph*> gr); 
static void ConformCanvas(TCanvas* can)
{
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(3); 
  can -> SetTopMargin(0.1);
};

float WeightedShapeError(std::vector<TH1F*> Pred, std::vector<TH1F*> Truth);
void Normalize(std::vector<TH1F*> Hists); 
void Normalize(TH1F* Hist);

class Table
{
  public:
    int x_width;
    int y_width; 

    int x_i = -1; 
    int y_i = -1;

    int Round = 3; 
    bool Scientific = false; 
    bool Integer = false; 

    _TS Name = "Table.txt";
    _TS Outdir = ""; 

    std::vector<std::vector<float>> values; 
    std::vector<_TS> x_cellNames; 
    std::vector<_TS> y_cellNames; 

    std::vector<_TS> Output; 
    std::vector<std::vector<float>> TMP_Score; 


    _TS Spacer = " | "; 
 
    int LargestElement(std::vector<_TS> c); 
    _TS Rounded(float l); 
    void CompileTable();
    void AppendLine(std::vector<std::vector<float>> val);
    void BestValues(); 
    void MakeHeader();
    void ResultList(); 
    void OtherStats();



    void Print(){ for (_TS O : Output){ std::cout << O << std::endl; } }; 
    void Print_Debug() 
    { 
      for (std::vector<float> i : values)
      { 
        for (float c : i) {std::cout << c << std::endl; }
      }
    };

    void WriteToFile()
    {
      std::ofstream print; 
      print.open(Outdir + "/" + Name); 
      for (_TS a : Output){ print << a << "\n"; }
      print.close(); 
    };
    
    void RepeatString(_TS* inp, int ntimes, _TS string)
    { 
      for (int i(0); i < ntimes; i++){*inp += string;}
    }; 

    void Newline()
    { 
      values.push_back({}); 
      y_i++; 
      x_i = -1; 
    };

    void AddValue(float x)
    {
      values[y_i].push_back(x); 
      x_i++; 
    };

};





#endif
