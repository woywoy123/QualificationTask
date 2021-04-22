#include<PostAnalysis/IO.h>

std::vector<TString> ReturnCurrentDirs()
{
  std::vector<TString> Output; 
  TDirectory* dir = gDirectory; 
  for (TObject* key : *dir -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key); 
    TString dir = (TString)k -> GetName(); 
    Output.push_back(dir); 
  }
  return Output; 
}

std::map<TString, std::map<TString, std::vector<TH1F*>>> ReadCTIDE(TString dir)
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Output; 
  TFile* F = new TFile(dir, "READ");
  for (TString L1 : Layer)
  {
    for (TString E1 : JetEnergy)
    {
      std::map<TString, std::vector<TH1F*>> Trk_Tru; 

      F -> cd(L1 + "/" + E1 + "_radius");  
      std::vector<TString> Folders = ReturnCurrentDirs(); 
      std::vector<TString> CollectHists; 
      for (TString x : Folders)
      {
        TH1F* H = (TH1F*)gDirectory -> Get(x); 
        
        // Truth inside the jet core. 
        if (x.Contains("ntrk_1") && x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_1_T_I"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_2_T_I"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_3_T_I"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_4_T_I"].push_back(H); }

        // Non Truth measurement. 
        if (x.Contains("ntrk_1") && !x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_1_M_I"].push_back(H); }
        if (x.Contains("ntrk_2") && !x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_2_M_I"].push_back(H); }
        if (x.Contains("ntrk_3") && !x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_3_M_I"].push_back(H); }
        if (x.Contains("ntrk_4") && !x.Contains("ntru") && x.Contains("rless") ){ Trk_Tru["ntrk_4_M_I"].push_back(H); }

        // Outside jetcore measurement 
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && !x.Contains("ntru")){ Trk_Tru["ntrk_1_M_O"].push_back(H); }
      }
      
      Output[L1 + "_" + E1] = Trk_Tru; 
      F -> cd(); 
    }
  }

  return Output; 
}

void TestReadCTIDE(TString dir)
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(dir); 
  
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    std::cout << x -> first << std::endl;
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 
    
    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    std::vector<TH1F*> ntrk_2_T = M["ntrk_2_T_I"]; 
    std::vector<TH1F*> ntrk_3_T = M["ntrk_3_T_I"]; 
    std::vector<TH1F*> ntrk_4_T = M["ntrk_4_T_I"]; 

    std::vector<TH1F*> ntrk_1_M = M["ntrk_1_M_I"]; 
    std::vector<TH1F*> ntrk_2_M = M["ntrk_2_M_I"]; 
    std::vector<TH1F*> ntrk_3_M = M["ntrk_3_M_I"]; 
    std::vector<TH1F*> ntrk_4_M = M["ntrk_4_M_I"]; 
    
    std::vector<TH1F*> ntrk_1 = M["ntrk_1_M_O"]; 
    
    for (int i(0); i < ntrk_1_M.size(); i++)
    {
      std::cout << ntrk_1_M[i] -> GetTitle() << " " << ntrk_2_M[i] -> GetTitle() << " " << ntrk_3_M[i] -> GetTitle() << " " << ntrk_4_M[i] -> GetTitle() << std::endl; 
    }

    for (int i(0); i < ntrk_1_T.size(); i++)
    {
      std::cout << ntrk_1_T[i] -> GetTitle() << " " << ntrk_2_T[i] -> GetTitle() << " " << ntrk_3_T[i] -> GetTitle() << " " << ntrk_4_T[i] -> GetTitle() << std::endl; 
    }
  }
}























