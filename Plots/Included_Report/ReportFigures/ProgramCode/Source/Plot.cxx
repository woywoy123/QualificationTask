#include "IO.h"
#include "Plot.h"
#include "Plotting.h"

void ShapePerformance()
{

  TString d = "ShapePerformance"; 
  TString dir = "MultiTrackFit.root"; 
  MMTF Hists = ReadAlgorithmResults(dir);
  
  for (MMTFi x = Hists.begin(); x != Hists.end(); x++)
  {
    TString LJE = x -> first; 
    MTF Algs = x -> second; 
    
    VTF ntrk1_T = Algs["ntrk_1_Truth"]; 
    VTF ntrk2_T = Algs["ntrk_2_Truth"]; 
    VTF ntrk3_T = Algs["ntrk_3_Truth"];     
    VTF ntrk4_T = Algs["ntrk_4_Truth"]; 
    
    std::cout << LJE << std::endl;
    MVVF Collect;  
    for (MTFi y = Algs.begin(); y != Algs.end(); y++)
    {
      TString alg = (y -> first);
      if (alg.Contains("Truth")){continue;}

      GetResults(y -> second, ntrk1_T, alg, "_ntrk_1", &Collect, "Shape", "Test"); 
      GetResults(y -> second, ntrk2_T, alg, "_ntrk_2", &Collect, "Shape", "Test"); 
      GetResults(y -> second, ntrk3_T, alg, "_ntrk_3", &Collect, "Shape", "Test"); 
      GetResults(y -> second, ntrk4_T, alg, "_ntrk_4", &Collect, "Shape", "Test"); 
    }
  
    FillShapeErrorMaps(LJE, Collect, &ResultsMap); 
  }
  
  for (TString L : Layer)
  {
    for (TString a : Algos)
    {
      CompileGraphs(L, a, ResultsMap, d, "Shape"); 
    }
  }
}


void IntegralPerformance()
{

  TString d = "IntegralPerformance"; 
  TString dir = "MultiTrackFit.root"; 
  MMTF Hists = ReadAlgorithmResults(dir);
  
  for (MMTFi x = Hists.begin(); x != Hists.end(); x++)
  {
    TString LJE = x -> first; 
    MTF Algs = x -> second; 
    
    VTF ntrk1_T = Algs["ntrk_1_Truth"]; 
    VTF ntrk2_T = Algs["ntrk_2_Truth"]; 
    VTF ntrk3_T = Algs["ntrk_3_Truth"];     
    VTF ntrk4_T = Algs["ntrk_4_Truth"]; 
    
    MVVF Collect;  
    for (MTFi y = Algs.begin(); y != Algs.end(); y++)
    {
      TString alg = (y -> first);
      if (alg.Contains("Truth")){continue;}

      GetResults(y -> second, ntrk1_T, alg, "_ntrk_1", &Collect, "Integral", "Template");
      GetResults(y -> second, ntrk2_T, alg, "_ntrk_2", &Collect, "Integral", "Template");
      GetResults(y -> second, ntrk3_T, alg, "_ntrk_3", &Collect, "Integral", "Template");
      GetResults(y -> second, ntrk4_T, alg, "_ntrk_4", &Collect, "Integral", "Template");
    }
  
    FillShapeErrorMaps(LJE, Collect, &ResultsMap); 
  }

  for (TString L : Layer)
  {
    for (TString a : Algos)
    {
      CompileGraphs(L, a, ResultsMap, d, "Integral"); 
    }
  }
}

void Histograms()
{

  TString d = "Histograms"; 
  TString dir = "MultiTrackFit.root"; 
  MMTF Hists = ReadAlgorithmResults(dir);
  TCanvas* can = new TCanvas(); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(4);
  can -> SetRightMargin(0.01); 
  can -> SetTopMargin(0.1); 
  can -> SetLogy(); 
  
  for (MMTFi x = Hists.begin(); x != Hists.end(); x++)
  {
    TString LJE = x -> first; 
    MTF Algs = x -> second; 
    
    VTF ntrk1_T = Algs["ntrk_1_Truth"]; 
    VTF ntrk2_T = Algs["ntrk_2_Truth"]; 
    VTF ntrk3_T = Algs["ntrk_3_Truth"];     
    VTF ntrk4_T = Algs["ntrk_4_Truth"]; 
    
    MVVF Collect;  
    for (MTFi y = Algs.begin(); y != Algs.end(); y++)
    {
      TString alg = (y -> first);
      if (alg.Contains("Truth")){continue;}

      VTF Alg = y -> second;
      
      VTF ntrk1; 
      VTF ntrk2; 
      VTF ntrk3; 
      VTF ntrk4; 
      float sum = 0; 
      for (int i(0); i < Alg.size(); i++)
      {
        TString name = Alg[i] -> GetTitle(); 
        
        if (!name.Contains("Test")){continue;}
        if (name.Contains("ntrk_1_ntru")){ ntrk1.push_back(Alg[i]);}
        if (name.Contains("ntrk_2_ntru")){ ntrk2.push_back(Alg[i]);}
        if (name.Contains("ntrk_3_ntru")){ ntrk3.push_back(Alg[i]);}
        if (name.Contains("ntrk_4_ntru")){ ntrk4.push_back(Alg[i]);}
        
        sum += Alg[i] -> Integral(); 
      }
     
      TString Layer; 
      if (LJE.Contains("IBL")){Layer = "IBL";}
      if (LJE.Contains("Blayer")){Layer = "Blayer";}
      if (LJE.Contains("layer1")){Layer = "layer1";}
      if (LJE.Contains("layer2")){Layer = "layer2";}

      TString LJE_T = LJE; 
      TString alg_T = alg; 
      if (Layer != ""){LJE_T.ReplaceAll(Layer+ "_", "").ReplaceAll("_GeV", "").ReplaceAll("_", "-");}
      else {LJE_T.ReplaceAll("_GeV", "").ReplaceAll("_", "-");}
      
      std::cout << LJE_T << std::endl;
      TString name; 
      if (LJE_T != Layer){ name = Layer + " at " + LJE_T;}
      else {name = Layer;}

      if (sum <= 5){ continue; }
      if (ntrk1.size() != 0)
      {
        name = alg_T.ReplaceAll("_ntrk_1", ": ") + name + " Track-1";
        PlotHists(ntrk1, ntrk1_T, name, can); 
      }

      if (ntrk2.size() != 0)
      {
        name = alg_T.ReplaceAll("_ntrk_2", ": ") + name + " Track-2";
        PlotHists(ntrk2, ntrk2_T, name, can); 
      }

      if (ntrk3.size() != 0)
      {
        name = alg_T.ReplaceAll("_ntrk_3", ": ") + name + " Track-3";
        PlotHists(ntrk3, ntrk3_T, name, can); 
      }

      if (ntrk4.size() != 0)
      {
        name = alg_T.ReplaceAll("_ntrk_4", ": ") + name + " Track-4";
        PlotHists(ntrk4, ntrk4_T, name, can); 
      }

      can -> Print(d + "/" + name + ".png"); 
    }
  }
 
}

void AdjustedShapePerformance()
{
  TString d = "ShapePerformance"; 
  TString dir = "MultiTrackFit.root"; 
  MMTF Hists = ReadAlgorithmResults(dir);
  
  for (MMTFi x = Hists.begin(); x != Hists.end(); x++)
  {
    TString LJE = x -> first; 
    MTF Algs = x -> second; 
    
    VTF ntrk1_T = Algs["ntrk_1_Truth"]; 
    VTF ntrk2_T = Algs["ntrk_2_Truth"]; 
    VTF ntrk3_T = Algs["ntrk_3_Truth"];     
    VTF ntrk4_T = Algs["ntrk_4_Truth"]; 
    
    MVVF Collect;  
    MVVF Integral; 
    for (MTFi y = Algs.begin(); y != Algs.end(); y++)
    {
      TString alg = (y -> first);
      if (alg.Contains("Truth")){continue;}

      GetResults(y -> second, ntrk1_T, alg, "_ntrk_1", &Collect, "AdjustedShape", "Test"); 
      GetResults(y -> second, ntrk2_T, alg, "_ntrk_2", &Collect, "AdjustedShape", "Test"); 
      GetResults(y -> second, ntrk3_T, alg, "_ntrk_3", &Collect, "AdjustedShape", "Test"); 
      GetResults(y -> second, ntrk4_T, alg, "_ntrk_4", &Collect, "AdjustedShape", "Test"); 


      GetResults(y -> second, ntrk1_T, alg, "_ntrk_1", &Integral, "TruthIntegral", "Test"); 
      GetResults(y -> second, ntrk2_T, alg, "_ntrk_2", &Integral, "TruthIntegral", "Test"); 
      GetResults(y -> second, ntrk3_T, alg, "_ntrk_3", &Integral, "TruthIntegral", "Test"); 
      GetResults(y -> second, ntrk4_T, alg, "_ntrk_4", &Integral, "TruthIntegral", "Test"); 
    }
  
    FillShapeErrorMaps(LJE, Collect, &ResultsMapAdjustedShape); 
    FillShapeErrorMaps(LJE, Integral, &ResultsMapIntegralShape); 
  }
  
  MMMMF AdjustedShape; // Layer / Trk /Algo /JE - Error 
  MMMF Integral; // Layer / trk / JE - Integral
  for (TString L : Layer)
  {
    for (TString a : Algos)
    {
      MMMF ntrk_mtru= ResultsMapAdjustedShape[L][a]; 
      
      for (MMMFi ntrk = ntrk_mtru.begin(); ntrk != ntrk_mtru.end(); ntrk++)
      {
        TString trk = ntrk -> first; 
        for (MMFi ntru = (ntrk -> second).begin(); ntru != (ntrk -> second).end(); ntru++)
        {
          TString tru = ntru -> first; 
          for (MFi JE = (ntru -> second).begin(); JE != (ntru -> second).end(); JE++)
          {
            TString je = JE -> first; 
            if (!je.Contains("_GeV")){continue;}
            AdjustedShape[L][trk][a][je] += ResultsMapAdjustedShape[L][a][trk][tru][je]; 
            
            if (a != "Normalization"){continue;}
            Integral[L][trk][je] += 0.01*ResultsMapIntegralShape[L][a][trk][tru][je];
          }
        }
      }
    }
  }
  for (TString L : Layer)
  {
    CompileGraphs(AdjustedShape[L], "AdjustedShapePerformance", L, "Algorithm Performance Over All Jet Energies in Layer: ");
  }
  
  CompileGraphs(Integral, "AdjustedShapePerformance", "", "Number of Clusters Sampled for each Jet Energy: ", 1, 1e9, "Number of Clusters");
  
}



int main()
{
  ShapePerformance();  
  //IntegralPerformance();  
  //Histograms(); 
  //AdjustedShapePerformance();
  return 0;
}
