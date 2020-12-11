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

std::map<TString, std::vector<TH1F*>> GetHist(std::map<TString, std::vector<TString>> Map, TString dir)
{ 
  TFile* F = new TFile(dir, "READ");
  std::map<TString, std::vector<TH1F*>> Out;  
  
  std::map<TString, std::vector<TString>>::iterator M; 
  for (M = Map.begin(); M != Map.end(); M++)
  {
    TString layer = M -> first; 
    std::vector<TString> H = M -> second; 
    std::vector<TH1F*> Hist_V; 
    
    F -> cd(layer); 
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

std::map<TString, std::vector<TH1F*>> MonteCarlo(TString dir)
{
  auto Name_Generator = [](int ntrk, int i)
  {
    std::vector<TString> Names; 
    for (int x(0); x < i; x++)
    {
      TString name = "DEDX_NEW_TRUTH_DEF/dEdx_ntrk_"; name +=(x+1); name += ("_nsdo_"); name += (x+1); name +=("_1keV"); 
      Names.push_back(name); 
    }
    return Names; 
  };

  auto Name_Layer_Track_Generator =[] (std::vector<TString> Layer, TString trk)
  {
    std::vector<TString> Out; 
    for (TString l : Layer)
    {
      TString n = trk + "_" + l; 
      Out.push_back(n);  
    }
    return Out; 
  };
  auto Merge_All = [&](std::map<TString, std::vector<TH1F*>> m, TString trk)
  {
    std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
    std::vector<TString> n = Name_Layer_Track_Generator(Layer, trk); 
   
    std::vector<TH1F*> ntrk = m[n[0]];  
    std::vector<TH1F*> trk_n; 
    for (int i(0); i < ntrk.size(); i++)
    {
      TH1F* h = (TH1F*)ntrk[i] -> Clone(trk); 
      h -> Reset(); 
      h -> SetTitle(trk); 
      trk_n.push_back(h);   
    }
    
    for (int i(0); i < n.size(); i++)
    {
      std::vector<TH1F*> Hists = m[n[i]]; 
      for (int x(0); x < Hists.size(); x++)
      {
        TH1F* H = trk_n[x];  
        TH1F* H2 = Hists[x]; 
        H -> Add(H2);  
      }
    }
    return trk_n; 
  };
 
  std::map<TString, std::vector<TString>> trk1_m; 
  trk1_m["IBL"] =  Name_Generator(1, 4); 
  trk1_m["Blayer"] =  Name_Generator(1, 4); 
  trk1_m["layer1"] =  Name_Generator(1, 4); 
  trk1_m["layer2"] =  Name_Generator(1, 4); 
  std::map<TString, std::vector<TH1F*>> trk1_H = GetHist(trk1_m, dir);  

  std::map<TString, std::vector<TString>> trk2_m; 
  trk2_m["IBL"] =  Name_Generator(2, 4); 
  trk2_m["Blayer"] =  Name_Generator(2, 4); 
  trk2_m["layer1"] =  Name_Generator(2, 4); 
  trk2_m["layer2"] =  Name_Generator(2, 4); 
  std::map<TString, std::vector<TH1F*>> trk2_H = GetHist(trk2_m, dir);  

  std::map<TString, std::vector<TString>> trk3_m; 
  trk3_m["IBL"] =  Name_Generator(3, 4); 
  trk3_m["Blayer"] =  Name_Generator(3, 4); 
  trk3_m["layer1"] =  Name_Generator(3, 4); 
  trk3_m["layer2"] =  Name_Generator(3, 4); 
  std::map<TString, std::vector<TH1F*>> trk3_H = GetHist(trk3_m, dir);  

  std::map<TString, std::vector<TString>> trk4_m; 
  trk4_m["IBL"] =  Name_Generator(4, 4); 
  trk4_m["Blayer"] =  Name_Generator(4, 4); 
  trk4_m["layer1"] =  Name_Generator(4, 4); 
  trk4_m["layer2"] =  Name_Generator(4, 4); 
  std::map<TString, std::vector<TH1F*>> trk4_H = GetHist(trk4_m, dir);  

  std::map<TString, std::vector<TString>> trk_dEdx; 
  trk_dEdx["IBL"] = {"dEdx_ntrk_1", "dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4"}; 
  trk_dEdx["Blayer"] = {"dEdx_ntrk_1", "dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4"}; 
  trk_dEdx["layer1"] = {"dEdx_ntrk_1", "dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4"}; 
  trk_dEdx["layer2"] = {"dEdx_ntrk_1", "dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4"}; 
  std::map<TString, std::vector<TH1F*>> dEdx_H = GetHist(trk_dEdx, dir);  

  std::map<TString, std::vector<TH1F*>> Map; 
  Map["trk1_IBL"] = trk1_H["IBL"]; 
  Map["trk2_IBL"] = trk2_H["IBL"]; 
  Map["trk3_IBL"] = trk3_H["IBL"]; 
  Map["trk4_IBL"] = trk4_H["IBL"]; 
  Map["trk1_Blayer"] = trk1_H["Blayer"]; 
  Map["trk2_Blayer"] = trk2_H["Blayer"]; 
  Map["trk3_Blayer"] = trk3_H["Blayer"]; 
  Map["trk4_Blayer"] = trk4_H["Blayer"]; 
  Map["trk1_layer1"] = trk1_H["layer1"]; 
  Map["trk2_layer1"] = trk2_H["layer1"]; 
  Map["trk3_layer1"] = trk3_H["layer1"]; 
  Map["trk4_layer1"] = trk4_H["layer1"]; 
  Map["trk1_layer2"] = trk1_H["layer2"]; 
  Map["trk2_layer2"] = trk2_H["layer2"]; 
  Map["trk3_layer2"] = trk3_H["layer2"]; 
  Map["trk4_layer2"] = trk4_H["layer2"]; 
  Map["trk1_All"] = Merge_All(Map, "trk1"); 
  Map["trk2_All"] = Merge_All(Map, "trk2"); 
  Map["trk3_All"] = Merge_All(Map, "trk3"); 
  Map["trk4_All"] = Merge_All(Map, "trk4"); 
  Map["dEdx_IBL"] = dEdx_H["IBL"]; 
  Map["dEdx_Blayer"] = dEdx_H["Blayer"]; 
  Map["dEdx_layer1"] = dEdx_H["layer1"]; 
  Map["dEdx_layer2"] = dEdx_H["layer2"]; 
  Map["dEdx_ntrk_All"] = Merge_All(Map, "dEdx"); 

  return Map;

}
