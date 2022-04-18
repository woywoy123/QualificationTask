#include "../PlottingCode/TestToDataShape.h"

void Analysis_TestFits(_TS Layer, std::map<_TS, std::vector<TH1F*>> CTIDE_C, 
                                  std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST_C, 
                                  bool Plot, 
                                  std::map<_TS, float>* Output
)
{
  
  TCanvas* can = new TCanvas(); 
  for (_TS E : Energy_H)
  {
    std::vector<TH1F*> ntrk1_T = CTIDE_C[Layer + "/" + E + "/ntrk1/InsideTruth"]; 
    std::vector<TH1F*> ntrk2_T = CTIDE_C[Layer + "/" + E + "/ntrk2/InsideTruth"]; 
    std::vector<TH1F*> ntrk3_T = CTIDE_C[Layer + "/" + E + "/ntrk3/InsideTruth"]; 
    std::vector<TH1F*> ntrk4_T = CTIDE_C[Layer + "/" + E + "/ntrk4/InsideTruth"]; 

    TH1F* ntrk1 = CTIDE_C[Layer + "/" + E + "/ntrk1/InsideData"][0];   
    TH1F* ntrk2 = CTIDE_C[Layer + "/" + E + "/ntrk2/InsideData"][0];   
    TH1F* ntrk3 = CTIDE_C[Layer + "/" + E + "/ntrk3/InsideData"][0];   
    TH1F* ntrk4 = CTIDE_C[Layer + "/" + E + "/ntrk4/InsideData"][0];   

    for (_TS Alg : Algos_H)
    {
      _TS Mode = "Template"; 

      std::vector<TH1F*> ntrk1_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk1"]; 
      std::vector<TH1F*> ntrk2_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk2"]; 
      std::vector<TH1F*> ntrk3_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk3"]; 
      std::vector<TH1F*> ntrk4_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk4"]; 
      if (ntrk1_Fit.size() == 0){std::cout << "-----> SKIPPING @ " + Layer + " " + E + " " + Alg << std::endl; continue;}

      if (Plot)
      {
        ConformCanvas(can);
        _TS tmp = E; 
        tmp.ReplaceAll("_GeV", " GeV").ReplaceAll("_", "-");
        _TS name1 = "n-Track Templates Fitted to 1-Track Data (Cluster Measurements) Distributions: " + Alg + " " + Layer + " " + tmp; 
        _TS name2 = "n-Track Templates Fitted to 2-Track Data (Cluster Measurements) Distributions: " + Alg + " " + Layer + " " + tmp; 
        _TS name3 = "n-Track Templates Fitted to 3-Track Data (Cluster Measurements) Distributions: " + Alg + " " + Layer + " " + tmp; 
        _TS name4 = "n-Track Templates Fitted to 4-Track Data (Cluster Measurements) Distributions: " + Alg + " " + Layer + " " + tmp; 

        tmp.ReplaceAll(" GeV", "");
        _TS dir1 = "Distributions/"+Alg+"/Track-1/"+Layer+"_"+ tmp +".png"; 
        _TS dir2 = "Distributions/"+Alg+"/Track-2/"+Layer+"_"+ tmp +".png"; 
        _TS dir3 = "Distributions/"+Alg+"/Track-3/"+Layer+"_"+ tmp +".png"; 
        _TS dir4 = "Distributions/"+Alg+"/Track-4/"+Layer+"_"+ tmp +".png"; 

        PlotHists(ntrk1, ntrk1_Fit, ntrk1_T, name1, can);
        can -> Print(dir1);

        PlotHists(ntrk2, ntrk2_Fit, ntrk2_T, name2, can);
        can -> Print(dir2);

        PlotHists(ntrk3, ntrk3_Fit, ntrk3_T, name3, can);
        can -> Print(dir3);

        PlotHists(ntrk4, ntrk4_Fit, ntrk4_T, name4, can);
        can -> Print(dir4);
      }
      else
      {
        (*Output)["trk1/" + Layer + "/" + E + "/" + Alg] = WeightedShapeError(ntrk1_Fit, ntrk1_T);
        (*Output)["trk2/" + Layer + "/" + E + "/" + Alg] = WeightedShapeError(ntrk2_Fit, ntrk2_T);
        (*Output)["trk3/" + Layer + "/" + E + "/" + Alg] = WeightedShapeError(ntrk3_Fit, ntrk3_T);
        (*Output)["trk4/" + Layer + "/" + E + "/" + Alg] = WeightedShapeError(ntrk4_Fit, ntrk4_T);
      }
    }
  }
}


int main(int argc, char* argv[])
{
  TString Fit_ROOT = argv[1]; 
  TString Merged_ROOT = argv[2]; 
  
  std::map<_TS, std::vector<TH1F*>> CTIDE = ReadCTIDE(Merged_ROOT);
  std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST = ReadPostAnalysis(Fit_ROOT); 
  std::map<_TS, float> Map; 
  float min = 0.0001; 
  float max = 10; 
  bool Plotting = true;
  
  for (_TS L : Layer_H)
  { 
    Analysis_TestFits(L, CTIDE, POST, Plotting, &Map); 
    if (!Plotting){continue;}
    Analysis_TestFits(L, CTIDE, POST, false, &Map); 
  }
  

  for (_TS L : Layer_H)
  {
    _TS trk1 = "ShapePerformance/" + L + "/1_Track.png"; 
    _TS trk2 = "ShapePerformance/" + L + "/2_Track.png"; 
    _TS trk3 = "ShapePerformance/" + L + "/3_Track.png"; 
    _TS trk4 = "ShapePerformance/" + L + "/4_Track.png"; 

    _TS name1 = "1-Track Shape Performance of Algorithms Fitted to Data (Cluster Measurements) " + L; 
    _TS name2 = "2-Track Shape Performance of Algorithms Fitted to Data (Cluster Measurements) " + L; 
    _TS name3 = "3-Track Shape Performance of Algorithms Fitted to Data (Cluster Measurements) " + L; 
    _TS name4 = "4-Track Shape Performance of Algorithms Fitted to Data (Cluster Measurements) " + L; 

    std::vector<TGraph*> grs; 
    for (_TS Alg : Algos_H)
    {
      TGraph* gr1 = MakeGraph(Map, Alg, L, "trk1");
      grs.push_back(gr1); 
    }
    PlotMultiGraph(&grs, name1, min, max, trk1); 
    
    for (_TS Alg : Algos_H)
    {
      TGraph* gr1 = MakeGraph(Map, Alg, L, "trk2");
      grs.push_back(gr1); 
    }
    PlotMultiGraph(&grs, name2, min, max, trk2); 

    for (_TS Alg : Algos_H)
    {
      TGraph* gr1 = MakeGraph(Map, Alg, L, "trk3");
      grs.push_back(gr1); 
    }
    PlotMultiGraph(&grs, name3, min, max, trk3); 

    for (_TS Alg : Algos_H)
    {
      TGraph* gr1 = MakeGraph(Map, Alg, L, "trk4");
      grs.push_back(gr1); 
    }
    PlotMultiGraph(&grs, name4, min, max, trk4); 
  }  
  return 0; 
}
