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

std::map<TString, std::vector<TH1F*>> FillingMap(TString Key, std::vector<TString> Names, std::map<TString, std::vector<TH1F*>> Output)
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
}

std::vector<TString> NameGenerator(int ntrk, int nsdo, bool tru)
{
  std::vector<TString> n_trk_tru;
  std::vector<TString> n_trk; 
  for (int i(0); i < ntrk; i++)
  {
    TString n = "dEdx_ntrk_"; n+=(i+1); 
    for (int j(0); j < nsdo; j++)
    {
      TString t = n; 
      t+= ("_ntru_"); t+= (j+1);  
      n_trk_tru.push_back(t); 
    }
    n_trk.push_back(n); 
  }

  if (tru == true){return n_trk_tru;}
  else {return n_trk;}
}


std::map<TString, std::vector<TH1F*>> MC_Reader(TString Dir)
{
  TFile* F = new TFile(Dir, "READ"); 
  std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 
  int nsdo = 4; 
  int ntrk = 4; 
  std::vector<TString> Names = NameGenerator(ntrk, nsdo, true); 
  std::vector<TString> Data_Names = NameGenerator(ntrk, nsdo, false); 
  
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

std::map<TString, std::vector<TH1F*>> MC_Reader_All(TString Dir)
{
  TFile* F = new TFile(Dir, "READ"); 
  std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"};  
   
  int nsdo = 4; 
  int ntrk = 4; 
  std::vector<TString> Names = NameGenerator(ntrk, nsdo, true); 
  std::vector<TString> Data_Names = NameGenerator(ntrk, nsdo, false); 
  
  std::vector<TString> R_Gre_Names_T; 
  std::vector<TString> R_Les_Names_T; 
  std::vector<TString> R_Gre_Names; 
  std::vector<TString> R_Les_Names; 
  for (int i(0); i < ntrk; i++)
  {
    TString base = "dEdx_ntrk_"; base += (i+1); base += ("_"); 
    for (int j(0); j < nsdo; j++)
    {
      TString rgre = base + "rgreater_1_ntru_"; rgre += (j+1); 
      TString rles = base + "rless_005_ntru_"; rles += (j+1); 
      R_Gre_Names_T.push_back(rgre); 
      R_Les_Names_T.push_back(rles); 
    }
    TString rg = base + "rgreater_1"; 
    TString rl = base + "rless_005";
    R_Gre_Names.push_back(rg); 
    R_Les_Names.push_back(rl); 
  }
  
  std::map<TString, std::vector<TH1F*>> Output; 
  for (int i(0); i < Layer.size(); i++)
  {
    TString L = Layer[i]; 
    for (int j(0); j < JetEnergy.size(); j++)
    {
      TString J = JetEnergy[j]; 
      TString LJ = L + "/" + J + "/"; 
      F -> cd(LJ); 
      Output = FillingMap(L + "_" + J + "_Truth", Names, Output); 
      Output = FillingMap(L + "_" + J, Data_Names, Output); 
      
      F -> cd(); 
      TString RE = L + "/" + J + "_radius/"; 
      F -> cd(RE);
      Output = FillingMap(L + "_" + J + "_" + "radius_Less", R_Les_Names, Output); 
      Output = FillingMap(L + "_" + J + "_" + "radius_Greater", R_Gre_Names, Output); 
      Output = FillingMap(L + "_" + J + "_" + "radius_Less_Truth", R_Les_Names_T, Output); 
      Output = FillingMap(L + "_" + J + "_" + "radius_Greater_Truth", R_Gre_Names_T, Output); 
     
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

std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> ReadResults(TString File_Name)
{
  TFile* F = new TFile(File_Name, "READ");
  
  std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> Maps; 
  for (TObject* key : *F -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key); 
    TString dir = (TString)k -> GetName();
    F -> cd(dir); 
    TDirectory* subdir = gDirectory;  
    
    std::vector<std::pair<TString, std::vector<TH1F*>>> M; 
    for (TObject* subkey : *subdir -> GetListOfKeys())
    {
      auto k_sub = dynamic_cast<TKey*>(subkey); 
      TString sub_dir = (TString)k_sub -> GetName(); 
      
      F -> cd(dir + "/" + sub_dir); 
      TDirectory* hist_dir = gDirectory;
      
      std::vector<TH1F*> Hi; 
      for (TObject* hists : *hist_dir -> GetListOfKeys())
      {
        auto h = dynamic_cast<TKey*>(hists); 
        TString H_name = (TString)h -> GetName(); 
        TH1F* h_temp = (TH1F*)gDirectory -> Get(H_name); 
        TString n = h_temp -> GetTitle(); 
        
        TH1F* H_Copy; 
        if (n.Contains("Error"))
        { 
          H_Copy = (TH1F*)h_temp -> Clone(H_name + "_" + dir); 
          H_Copy -> SetTitle(H_name + "_" + dir);
        } 
        else
        {
          H_Copy = (TH1F*)h_temp -> Clone(H_name); 
        }
        Hi.push_back(H_Copy); 
      }
      M.push_back(std::pair<TString, std::vector<TH1F*>>(sub_dir, Hi)); 
    }
    Maps[dir] = M;
    F -> cd(); 
  }
  return Maps; 
}
