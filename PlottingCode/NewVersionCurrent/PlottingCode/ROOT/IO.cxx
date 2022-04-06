#include "../PlottingCode/IO.h"

std::vector<_TS> ReturnCurrentDirs()
{
  std::vector<_TS> Output; 
  TDirectory* dir = gDirectory; 
  for (TObject* key : *dir -> GetListOfKeys())
  {
    auto k = dynamic_cast<TObject*>(key); 
    TString dir = (_TS)k->GetName(); 
    Output.push_back(dir);
  }
  return Output; 
}


std::map<_TS, std::vector<TH1F*>> ReadCTIDE(_TS dir)
{
  TFile* f = new TFile(dir, "READ");
  std::vector<_TS> Detector_Layer = ReturnCurrentDirs();
  std::map<_TS, std::vector<TH1F*>> Output; 

  for (_TS L : Detector_Layer)
  {
    f -> cd(L);
    std::vector<_TS> Energy = ReturnCurrentDirs();
    for (_TS E : Energy)
    {
      if (!E.Contains("radius")){continue;}

      f -> cd(L + "/" +E);
      std::vector<_TS> ntrks = ReturnCurrentDirs();
      for (_TS n : ntrks)
      {
        TH1F* H = (TH1F*)gDirectory -> Get(n);

        TString d = L + "/" + E + "/"; 
        d.ReplaceAll("_radius", "");  
        
        if (n.Contains("rless_") && n.Contains("ntru")){ Output[d + n + "/InsideTruth"].push_back(H); }
        if (n.Contains("rless_") && !n.Contains("ntru")){ Output[d + n + "/InsideData"].push_back(H); }
        if (n.Contains("rgreater_") && n.Contains("ntru")){ Output[d + n + "/OutsideTruth"].push_back(H); }
        if (n.Contains("rgreater_") && !n.Contains("ntru")){ Output[d + n + "/OutsideData"].push_back(H); }
      }
    }
  }
  
  return Output;
}

std::map<_TS, std::map<_TS, std::map<_TS, std::map<_TS, TH1F*>>>> ReadDebugging(_TS dir)
{
  TFile* f = new TFile(dir); 
  std::vector<_TS> Layer_Energy = ReturnCurrentDirs(); 
  
  // Layer_Energy - Step - Algo_Case - Result/truth
  std::map<_TS, std::map<_TS, std::map<_TS, std::map<_TS, TH1F*>>>> Output; 
  for (_TS le : Layer_Energy)
  {
    f -> cd(le); 
    std::vector<_TS> AlgosTestCases = ReturnCurrentDirs(); 
    for (_TS agc : AlgosTestCases)
    {
      if (agc.Contains("Results")){continue;}
      f -> cd(le + "/" + agc); 
      std::vector<_TS> Steps = ReturnCurrentDirs(); 
      for (_TS i : Steps)
      {
        f -> cd(le + "/" + agc + "/" + i); 
        _TS S = Split(i, "_")[1];
        std::vector<_TS> Hists = ReturnCurrentDirs(); 
        
        for (_TS H : Hists)
        {
          _TS T = ""; 
          if (H.Contains("dEdx_ntrk_1_ntru_1")){ T = "ntrk1-ntru1-D"; }
          if (H.Contains("dEdx_ntrk_1_ntru_2")){ T = "ntrk1-ntru2-D"; }
          if (H.Contains("ntrk1-tru1")){ T = "ntrk1-ntru1-T"; }
          if (H.Contains("ntrk1-tru2")){ T = "ntrk1-ntru2-T"; }
          
          Output[le][S][agc][T] = (TH1F*)gDirectory -> Get(H); 
        }
      }
    } 
  }
  return Output;
}

void ReadPostAnalysis(_TS dir)
{
  TFile* f = new TFile(dir);
  std::vector<_TS> Layer_Energy = ReturnCurrentDirs();
  for (_TS LE : Layer_Energy)
  {
    f -> cd(LE); 
    std::vector<_TS> Algos = ReturnCurrentDirs();
    for (_TS A : Algos)
    {
      f -> cd(LE + "/" + A); 
      std::vector<_TS> Hists = ReturnCurrentDirs(); 

      std::cout << LE + "/" + A << std::endl;
      for (_TS H : Hists)
      {
        std::cout << "-> " << H << std::endl;
        

      }


    }
  }
}
