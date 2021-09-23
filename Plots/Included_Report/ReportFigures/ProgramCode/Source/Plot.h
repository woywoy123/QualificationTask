#include<TCanvas.h>
#include "Plotting.h"


#ifndef PLOT_H
#define PLOT_H

MMMMMF ResultsMap; 
MMMMMF ResultsMapAdjustedShape; 
MMMMMF ResultsMapIntegralShape; 

void GetResults(VTF Hists, VTF Hists_T, TString input, TString trk, MVVF* Output, TString metric, TString name)
{
  if (input.Contains(trk))
  {
    VTF Use; 
    for (int i(0); i < Hists.size(); i++)
    { 
      TString Title = (Hists[i] -> GetTitle()); 
      if (Title.Contains(name)){ Use.push_back(Hists[i]); }
    }

    TString al = input.ReplaceAll(trk, "");
    VF r; 
    if (metric == "Shape"){ r = ShapeComparison(Use, Hists_T); }
    if (metric == "Integral"){ r = IntegralComparison(Use, Hists_T);}
    if (metric == "AdjustedShape"){ r = AdjustedShape(Use, Hists_T); }
    if (metric == "TruthIntegral"){ r = TruthIntegrals(Hists_T); }
    
    (*Output)[al].push_back(r); 
  }
}

void FillShapeErrorMaps(TString LJE, MVVF Collect, MMMMMF* ResultsMap)
{
  
  TString Layer; 
  if (LJE.Contains("IBL")){Layer = "IBL";}
  if (LJE.Contains("Blayer")){Layer = "Blayer";}
  if (LJE.Contains("layer1")){Layer = "layer1";}
  if (LJE.Contains("layer2")){Layer = "layer2";}

  TString JE = LJE.ReplaceAll(Layer + "_", ""); 
  for (MVVFi x = Collect.begin(); x != Collect.end(); x++)
  {
    TString alg = x -> first; 
    VVF ntrk_mtru = Collect[alg]; 
    
    for (int i(0); i < ntrk_mtru.size(); i++)
    {
      TString ntrk = "ntrk_"; ntrk += (i+1); 
      
      for (int k(0); k < ntrk_mtru[i].size(); k++)
      {
        TString ntru = "ntru_"; ntru += (k+1); 
        (*ResultsMap)[Layer][alg][ntrk][ntru][JE] = 100*ntrk_mtru[i][k]; 
      }
    }
  }
}

void CompileGraphs(TString Layer, TString Algo, MMMMMF ResultsMap, TString dir, TString name)
{
  MMMF Results = ResultsMap[Layer][Algo]; 
  
  for (MMMFi x = Results.begin(); x != Results.end(); x++)
  {
    TCanvas* can = new TCanvas(); 
    
    std::cout << x -> first << std::endl;
    gStyle -> SetOptStat(0); 
    gStyle -> SetOptTitle(0);
    gStyle -> SetImageScaling(4);
    can -> SetRightMargin(0.01); 
    can -> SetTopMargin(0.1); 

    TString ntrk = x -> first; 
    MMF ntru = Results[ntrk];
    
    TString Normalized = "Template "; 
    Normalized += (name); 
    Normalized += (" Matching (" + ntrk + "): " + Layer + " Using: " + Algo); 
    std::vector<TGraph*> gr = GenerateMultiGraph(ntru, Normalized, can); 
    TPaveText* T = new TPaveText(0.1, 0.9, 0.9, 1, "NDCNB"); 
    T -> AddText(Normalized);
    T -> Draw();
    can -> Update();
    can -> Print(dir + "/" + Layer + "_" + Algo + "_" + ntrk + ".png"); 
    can -> Clear();

    for (TGraph* g : gr){ delete g; }
    delete can;
  }
}


void CompileGraphs(MMMF ResultsMap, TString dir, TString Layer, TString name = "", float min = 1, float max = 200, TString yname = "")
{
  for (MMMFi x = ResultsMap.begin(); x != ResultsMap.end(); x++)
  {
    TCanvas* can = new TCanvas(); 
 
    gStyle -> SetOptStat(0); 
    gStyle -> SetOptTitle(0);
    gStyle -> SetImageScaling(4);
    can -> SetRightMargin(0.01); 
    can -> SetTopMargin(0.1); 

    TString trk = x -> first; 
    MMF L = ResultsMap[trk];
    
    TString Normalized = name + Layer + " ("+trk+")"; 
    std::vector<TGraph*> gr; 
    if (yname == ""){gr = GenerateMultiGraph(L, Normalized, can, min, max);}
    else {gr = GenerateMultiGraph(L, Normalized, can, min, max, yname);}

    TPaveText* T = new TPaveText(0.1, 0.9, 0.9, 1, "NDCNB"); 
    T -> AddText(Normalized);
    T -> Draw();
    can -> Update();
    can -> Print(dir + "/" + Layer + "_" + trk + ".png"); 
    can -> Clear();

    for (TGraph* g : gr){ delete g; }
    delete can;
  }
}








#endif
