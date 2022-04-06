#include "../PlottingCode/TestToTruthShape.h"

void Analysis_TestFits(_TS Layer, std::map<_TS, std::vector<TH1F*>> CTIDE_C, 
                                  std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST_C, 
                                  bool Plot
)
{

  TCanvas* can = new TCanvas(); 
  for (_TS E : Energy_H)
  {
    std::vector<TH1F*> ntrk1_T = CTIDE_C[Layer + "/" + E + "/ntrk1/InsideTruth"]; 
    std::vector<TH1F*> ntrk2_T = CTIDE_C[Layer + "/" + E + "/ntrk2/InsideTruth"]; 
    std::vector<TH1F*> ntrk3_T = CTIDE_C[Layer + "/" + E + "/ntrk3/InsideTruth"]; 
    std::vector<TH1F*> ntrk4_T = CTIDE_C[Layer + "/" + E + "/ntrk4/InsideTruth"]; 

    for (_TS Alg : Algos_H)
    {
      _TS Mode = "Test"; 
      if (Alg == "Experimental"){ Mode = "Template"; }

      std::vector<TH1F*> ntrk1_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk1"]; 
      std::vector<TH1F*> ntrk2_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk2"]; 
      std::vector<TH1F*> ntrk3_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk3"]; 
      std::vector<TH1F*> ntrk4_Fit = POST_C[Layer + "/" + E + "/" + Alg][Mode + "_ntrk4"]; 
      if (Plot)
      {
        ConformCanvas(can);
        _TS tmp = E; 
        tmp.ReplaceAll("_GeV", " GeV").ReplaceAll("_", "-");
        _TS name1 = "n-Track Templates Fitted to 1-Track Truth Distributions: " + Alg + " " + Layer + " " + tmp; 
        _TS name2 = "n-Track Templates Fitted to 2-Track Truth Distributions: " + Alg + " " + Layer + " " + tmp; 
        _TS name3 = "n-Track Templates Fitted to 3-Track Truth Distributions: " + Alg + " " + Layer + " " + tmp; 
        _TS name4 = "n-Track Templates Fitted to 4-Track Truth Distributions: " + Alg + " " + Layer + " " + tmp; 

        tmp.ReplaceAll(" GeV", "");
        _TS dir1 = "Distributions/"+Alg+"/Track-1/"+Layer+"_"+ tmp +".png"; 
        _TS dir2 = "Distributions/"+Alg+"/Track-2/"+Layer+"_"+ tmp +".png"; 
        _TS dir3 = "Distributions/"+Alg+"/Track-3/"+Layer+"_"+ tmp +".png"; 
        _TS dir4 = "Distributions/"+Alg+"/Track-4/"+Layer+"_"+ tmp +".png"; 

        PlotHists(ntrk1_Fit, ntrk1_T, name1, can);
        can -> Print(dir1);

        PlotHists(ntrk2_Fit, ntrk2_T, name2, can);
        can -> Print(dir2);

        PlotHists(ntrk3_Fit, ntrk3_T, name3, can);
        can -> Print(dir3);

        PlotHists(ntrk4_Fit, ntrk4_T, name4, can);
        can -> Print(dir4);
      }

      //std::cout << WeightedShapeError(ntrk1_Fit, ntrk1_T) << " " << Alg<< std::endl;
    }
  }
}




int main(int argc, char* argv[])
{
  TString Fit_ROOT = argv[1]; 
  TString Merged_ROOT = argv[2]; 
  
  std::map<_TS, std::vector<TH1F*>> CTIDE = ReadCTIDE(Merged_ROOT);
  std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST = ReadPostAnalysis(Fit_ROOT); 
  
  bool Plotting = true;
  
  for (_TS L : Layer_H){ Analysis_TestFits(L, CTIDE, POST, Plotting); }

  return 0; 
}
