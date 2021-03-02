#include<PostAnalysis/IO.h>

std::map<TString, std::vector<TH1F*>> ReadEntries(TFile* F)
{

  std::map<TString, std::vector<TH1F*>> Map; 
  for (TObject* key1 : *F -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key1);
    Map[(TString)k -> GetName()] = {}; 
    F -> cd(k -> GetName()); 
    TDirectory* subdir = gDirectory; 
    for (TObject* key2 : *subdir -> GetListOfKeys())
    {
      auto k2 = dynamic_cast<TKey*>(key2); 
      TString Histname = (TString)k2 -> GetName(); 
      TH1F* H = (TH1F*)gDirectory-> Get(Histname); 
      Map[(TString)k -> GetName()].push_back(H); 
    } 
    F -> cd(); 
  }
  return Map; 
}

std::map<TString, std::vector<TH1F*>> GetHist(std::map<TString, std::vector<TString>> Map, TString dir, TString scnd = "")
{ 
  TFile* F = new TFile(dir, "READ");
  std::map<TString, std::vector<TH1F*>> Out;  
  
  std::map<TString, std::vector<TString>>::iterator M; 
  for (M = Map.begin(); M != Map.end(); M++)
  {
    TString layer = M -> first; 
    std::vector<TString> H = M -> second; 
    std::vector<TH1F*> Hist_V; 
    
    F -> cd(layer + scnd); 
    for (TString H_N : H)   
    {
      TH1F* Hist = (TH1F*)gDirectory -> Get(H_N);
      Hist_V.push_back(Hist);  
    }
    F -> cd(); 
    Out[layer] = Hist_V; 
  }
  return Out; 

}

std::map<TString, std::vector<TH1F*>> MC_Reader(TString Dir)
{
  auto FillingMap =[] (TString Key, std::vector<TString> Names, std::map<TString, std::vector<TH1F*>> Output)
  {
    if (Output[Key].size() == 0)
    {
      for (TString name : Names)
      {
        TH1F* H = (TH1F*)gDirectory -> Get(name); 
        TH1F* H_N = (TH1F*)H -> Clone(name+"_"+Key); 
        H_N -> SetTitle(name+"_"+Key); 
        Output[Key].push_back(H_N);  
        delete H;
      }
    }
    else
    {
      for(int n(0); n < Names.size(); n++)
      {
        TH1F* H = (TH1F*)gDirectory -> Get(Names[n]); 
        Output[Key][n] -> Add(H, 1); 
        delete H; 
      }
    }
    return Output; 
  };


  TFile* F = new TFile(Dir, "READ"); 
  std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 

  int nsdo = 4; 
  int ntrk = 4; 
  std::vector<TString> Names; 
  for (int i(0); i < ntrk; i++)
  {
    for (int j(0); j < nsdo; j++)
    {
      TString n = "dEdx_ntrk_"; n+=(i+1); n+= ("_ntru_"); n+= (j+1);  
      Names.push_back(n); 
    }
  }
  
  std::vector<TString> Data_Names; 
  for (int i(0); i < ntrk; i++)
  {
    TString n = "dEdx_ntrk_"; n+=(i+1); 
    Data_Names.push_back(n);  
  }


  std::vector<TString> R_Gre_Names; 
  std::vector<TString> R_Les_Names; 
  for (int i(0); i < 4; i++)
  {
    TString base = "dEdx_ntrk_1_"; 
    TString rgre = base + "rgreater_1_ntru_"; rgre += (i+1); 
    TString rles = base + "rless_005_ntru_"; rles += (i+1); 
    R_Gre_Names.push_back(rgre); 
    R_Les_Names.push_back(rles); 
  }
  
  std::map<TString, std::vector<TH1F*>> Output; 

  for (int i(0); i < JetEnergy.size(); i++)
  {
    TString JE = JetEnergy[i]; 
    TString RE = JE + "_radius";
    for (int j(0); j < Layer.size(); j++)
    {
      TString L = Layer[j]; 
      TString D = L + "/" + JE + "/"; 
      F -> cd(D); 

      Output = FillingMap(JE, Names, Output); 
      Output = FillingMap(L, Names, Output); 
      Output = FillingMap("All", Names, Output); 
      Output = FillingMap(JE + "_Data", Data_Names, Output); 
      Output = FillingMap(L + "_Data", Data_Names, Output); 
      Output = FillingMap("All_Data", Data_Names, Output); 

      F -> cd();
      
      F -> cd(L + "/" + RE + "/");  
      Output = FillingMap(RE + "_Less", R_Les_Names, Output); 
      Output = FillingMap(RE + "_Greater", R_Gre_Names, Output); 
  
      Output = FillingMap(L + "_radius_Less", R_Les_Names, Output); 
      Output = FillingMap(L + "_radius_Greater", R_Gre_Names, Output); 
 
      Output = FillingMap("All_radius_Less", R_Les_Names, Output); 
      Output = FillingMap("All_radius_Greater", R_Gre_Names, Output); 
    
      F -> cd();
    
    }
  }
  return Output;
}

std::map<TString, std::vector<TH1F*>> Result_Reader(TString Dir)
{
  auto Names =[] (int trk, int tru)
  {
    std::vector<TString> Output; 
    for (int i(0); i < trk; i++)
    {
      for (int x(0); x < tru; x++)
      {
        TString title = "dEdx_ntrk_"; title += (i+1); title += ("_ntru_"); title +=(x+1); 
        Output.push_back(title);  
      }
    }
    return Output; 
  };

  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV", "IBL", "Blayer", "layer1", "layer2", "All"}; 
  std::vector<TString> Tracks = {"dEdx_ntrk_1", "dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4"}; 


  TFile* F = new TFile(Dir, "READ"); 
  F -> ReOpen("UPDATE"); 
  std::map<TString, std::vector<TH1F*>> Files = ReadEntries(F);  
  
  int trk = 4; 
  int tru = 4; 

  typedef std::map<TString, std::vector<TH1F*>>::iterator  it; 
  for (it p = Files.begin(); p != Files.end(); p++)
  {
    TH1F* H = (TH1F*)gDirectory -> Get(p -> first);
    Files[p -> first][0] = (TH1F*)H -> Clone(p -> first);
  }
  
  std::vector<TString> trkTitles = Names(trk, tru);
  std::map<TString, std::vector<TH1F*>> Map; 
  for (TString JE : JetEnergy)
  {
    for (TString trk : trkTitles)
    {
      for (it p = Files.begin(); p != Files.end(); p++)
      {
        TString key = p -> first; 
        if (key.Contains(JE) && key.Contains(trk)){Map[JE].push_back(Files[p -> first][0]);}
      }
    }

    for (it p = Files.begin(); p != Files.end(); p++)
    {
      TString key = p -> first; 
      if (key.Contains(JE) && key.Contains("_R") == false && key.Contains("ntrk_1")){Map[JE +"_Data"].push_back(Files[p -> first][0]);}
      if (key.Contains(JE) && key.Contains("_R") == false && key.Contains("ntrk_2")){Map[JE +"_Data"].push_back(Files[p -> first][0]);}
      if (key.Contains(JE) && key.Contains("_R") == false && key.Contains("ntrk_3")){Map[JE +"_Data"].push_back(Files[p -> first][0]);}
      if (key.Contains(JE) && key.Contains("_R") == false && key.Contains("ntrk_4")){Map[JE +"_Data"].push_back(Files[p -> first][0]);}
    }
  }
  
  return Map; 
}
