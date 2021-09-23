#include "IO.h"

std::vector<TString> ReturnCurrentDirs(bool FolderOnly = true)
{
  std::vector<TString> Output; 
  TDirectory* dir = gDirectory; 
  for (TObject* key : *dir -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key); 
    TString dir = (TString)k -> GetName();
    
    if (dir.Contains("_t")){ continue;}
    Output.push_back(dir); 
  }
  return Output; 
}

MMVTH1F ReadAlgorithmResults(TString dir)
{
  MMVTH1F Output; 
  TFile* F = new TFile(dir, "READ");
  std::vector<TString> JetEnergy_Layer_Folder = ReturnCurrentDirs(); 
  
  for (TString Folder : JetEnergy_Layer_Folder)
  {
    F -> cd(Folder);
    std::map<TString, std::vector<TH1F*>> Algorithms_Map;  
    
    std::vector<TString> Algorithm_Folder = ReturnCurrentDirs(); 
    for (TString Alg_Folder : Algorithm_Folder)
    {
      F -> cd(Folder + "/" + Alg_Folder); 
      
      std::vector<TString> Alg_V = ReturnCurrentDirs(); 
      for (TString H_TS : Alg_V)
      {
        if (H_TS.Contains("_error")){continue;}
        TH1F* H = (TH1F*)gDirectory -> Get(H_TS);

        // Truth inside the jet core. 
        if (!Alg_Folder.Contains("ntrk_1_T") && H_TS.Contains("ntrk_1")){ Algorithms_Map[Alg_Folder + "_ntrk_1"].push_back(H); }
        if (!Alg_Folder.Contains("ntrk_2_T") && H_TS.Contains("ntrk_2")){ Algorithms_Map[Alg_Folder + "_ntrk_2"].push_back(H); }
        if (!Alg_Folder.Contains("ntrk_3_T") && H_TS.Contains("ntrk_3")){ Algorithms_Map[Alg_Folder + "_ntrk_3"].push_back(H); }
        if (!Alg_Folder.Contains("ntrk_4_T") && H_TS.Contains("ntrk_4")){ Algorithms_Map[Alg_Folder + "_ntrk_4"].push_back(H); }
        
        // Truth for multitrack fit:
        if (Alg_Folder.Contains("ntrk_")){ Algorithms_Map[Alg_Folder + "ruth"].push_back(H); }
      }
    }
    Output[Folder] = Algorithms_Map;
    F -> cd(); 
  }

  return Output;  
}

MMMMMF ReadMinimizers(TString dir)
{
  // Minimizer - Layer - *_FLostX - Algo - JE
  MMMMMF Out; 

  TFile* F = new TFile(dir);
  std::vector<TString> Mini = ReturnCurrentDirs();
  for (TString m : Mini)
  {
    F -> cd(m); 
    std::vector<TString> Lay = ReturnCurrentDirs();
    
    for (TString L : Lay)
    {
      F -> cd(m + "/" + L); 
      std::vector<TString> FL = ReturnCurrentDirs();

      for (TString fl : FL)
      {
        F -> cd(m + "/" + L + "/" + fl); 
        std::vector<TString> alg = ReturnCurrentDirs();
        
        for (TString a : alg)
        {
          TH1F* H = (TH1F*)gDirectory -> Get(a); 
          TString t = H -> GetTitle();
          
          if (t.Contains("Truth")){continue;}
          
          TString a_k; 
          for (TString ag : Algos){ if (t.Contains(ag)){ a_k = ag; } }
          
          TString d = m;
          d.ReplaceAll("_", "");
          Out[fl][L][a_k][d] = ReadTH1F(H); 
        }
      }
    }
    F -> cd(); 
  }
  return Out;
}
