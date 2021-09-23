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

std::map<TString, std::map<TString, std::vector<TH1F*>>> ReadAlgorithmResults(TString dir)
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Output; 
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

float ChiSquareDistanceSym(TH1F* H1, TH1F* H2)
{
  float chi2 = 0; 
  for (int i = 0; i < H1 -> GetNbinsX(); i++)
  {
    float x = H1 -> GetBinContent(i+1); 
    float y = H2 -> GetBinContent(i+1); 
    float sum = x + y; 
    float dif = x-y; 
    if (sum == 0){continue;}
    chi2 += std::pow(dif, 2)/sum; 
  }
  chi2 = 0.5*chi2; 
  return chi2;
}

VF ShapeComparison(VTF Hists, VTF Truth)
{
  VF Output; 
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = (TH1F*)Hists[i] -> Clone("H"); 
    TH1F* H_T = (TH1F*)Truth[i] -> Clone("T"); 
    
    float l_H = H -> Integral();
    H -> Scale(1./l_H); 

    float l_T = H_T -> Integral();
    H_T -> Scale(1./l_T); 

    float chi2 = ChiSquareDistanceSym(H, H_T);  
    Output.push_back(chi2); 
    
    delete H, H_T; 
  }
  return Output;
}

VF AdjustedShape(VTF Hists, VTF Truth)
{
  VF Output; 
  VF Compare = ShapeComparison(Hists, Truth); 
  for (int i(0) ; i < Hists.size(); i++)
  {
    float r = Compare[i]; 
    float l_p = Hists[i] -> Integral(); 
    float l_t = Truth[i] -> Integral();
    float adj = r* ((l_p)/(l_t)); 
    if (std::isnan(adj)){ adj = 0; }
    Output.push_back(adj); 
  }
  return Output;
}

VF IntegralComparison(VTF Hists, VTF Truth)
{
  VF Output; 
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = (TH1F*)Hists[i] -> Clone("H"); 
    TH1F* H_T = (TH1F*)Truth[i] -> Clone("T"); 
   
    float e = std::abs((H -> Integral()) - (H_T -> Integral())) / (H_T -> Integral()); 
    if (H_T -> Integral() <= 0) { e = 1; }
    if (H -> Integral() <= 2) { e = 1;}

    Output.push_back( e); 
    
    delete H, H_T; 
  }
  return Output;
}

VF TruthIntegrals(VTF Truth)
{
  VF Output; 
  for (int i(0); i < Truth.size(); i++)
  {
    float e = Truth[i] -> GetEntries(); 
    Output.push_back(e);  
  }
  return Output; 
}





