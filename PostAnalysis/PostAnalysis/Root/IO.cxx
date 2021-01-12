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

std::map<TString, std::vector<TH1F*>> MonteCarlo(TString dir)
{
  auto Name_Generator = [](int ntrk, int i)
  {
    std::vector<TString> Names; 
    for (int x(0); x < i; x++)
    {
      TString name = "dEdx_ntrk_"; name +=(ntrk); name += ("_nsdo_"); name += (x+1); name +=("_1keV"); 
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
      TString name = ntrk[i] -> GetTitle(); name += ("_All"); 
      TH1F* h = (TH1F*)ntrk[i] -> Clone(name); 
      h -> Reset(); 
      h -> SetTitle(name); 
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
  std::map<TString, std::vector<TH1F*>> trk1_H = GetHist(trk1_m, dir, "/DEDX_NEW_TRUTH_DEF");  

  std::map<TString, std::vector<TString>> trk2_m; 
  trk2_m["IBL"] =  Name_Generator(2, 4); 
  trk2_m["Blayer"] =  Name_Generator(2, 4); 
  trk2_m["layer1"] =  Name_Generator(2, 4); 
  trk2_m["layer2"] =  Name_Generator(2, 4); 
  std::map<TString, std::vector<TH1F*>> trk2_H = GetHist(trk2_m, dir, "/DEDX_NEW_TRUTH_DEF");  

  std::map<TString, std::vector<TString>> trk3_m; 
  trk3_m["IBL"] =  Name_Generator(3, 4); 
  trk3_m["Blayer"] =  Name_Generator(3, 4); 
  trk3_m["layer1"] =  Name_Generator(3, 4); 
  trk3_m["layer2"] =  Name_Generator(3, 4); 
  std::map<TString, std::vector<TH1F*>> trk3_H = GetHist(trk3_m, dir, "/DEDX_NEW_TRUTH_DEF");  

  std::map<TString, std::vector<TString>> trk4_m; 
  trk4_m["IBL"] =  Name_Generator(4, 4); 
  trk4_m["Blayer"] =  Name_Generator(4, 4); 
  trk4_m["layer1"] =  Name_Generator(4, 4); 
  trk4_m["layer2"] =  Name_Generator(4, 4); 
  std::map<TString, std::vector<TH1F*>> trk4_H = GetHist(trk4_m, dir, "/DEDX_NEW_TRUTH_DEF");  

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

std::map<TString, std::vector<TH1F*>> MonteCarloLayerEnergy(TString dir)
{
  std::vector<TString> Layers = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 
  auto GenerateNames = [](int trks, int tru, TString Energy = "", int term = -1)
  {
    std::vector<TString> Output; 
    for (int i(0); i < trks; i++)
    {
      if ( term != -1 )
      {
        i = term; 
      }
      for (int y(0); y < tru; y++)
      {
        TString name = "dEdx_ntrk_"; name += (i+1); name += ("_ntru_"); name += (y+1); name += (Energy); 
        Output.push_back(name); 
      }
      if (term != -1 )
      {
        break; 
      }
    }
    return Output; 
  };

  auto Sort = [&](int trk, int tru, std::vector<TString> Dirs, std::map<TString, std::vector<TH1F*>> map, TString Energy)
  {
    int s_index = (trk-1)*tru; 
    int e_index = trk*tru;    
   
    std::vector<TString> Names = GenerateNames(trk, tru, Energy); 
    std::vector<TH1F*> Hists_V; 

    for (int p(0); p < Dirs.size(); p++)
    {
      TString n = Dirs[p]; 
      std::vector<TH1F*> Hist = map[n];  
      if (p == 0)
      {
        for (int i(s_index); i < e_index; i++)
        {
          TH1F* Temp = (TH1F*)Hist[i] -> Clone(Names[i]); 
          Temp -> SetTitle(Names[i]);  
          Hists_V.push_back(Temp);  
        }
      }
      else
      {
        int x(0); 
        for (int i(s_index); i < e_index; i++)
        {
          Hists_V[x] -> Add(Hist[i]); 
          x++; 
        }
      }
    }
    return Hists_V; 
  };

  auto Merge =[] (std::vector<TH1F*> Hout, std::vector<TH1F*> Hin)
  {
    for (int i(0); i < Hout.size(); i++)
    {
      Hout[i] -> Add(Hin[i]); 
    }
  }; 


  int ntru = 4; 
  int ntrk = 4; 
  
  // ================= Getting the eventWeight based histograms ========================= // 
  std::vector<TString> Names = GenerateNames(ntrk, ntru); 
  std::map<TString, std::vector<TString>> ToGet; 
  for (int i(0); i < Layers.size(); i++)
  {
    for (int y(0); y < JetEnergy.size(); y++)
    {
      TString name = Layers[i]; name += ("/"); name += (JetEnergy[y]); 
      ToGet[name] = Names; 
    }
  }

  std::map<TString, std::vector<TH1F*>> Dump = GetHist(ToGet, dir); 
  std::map<TString, std::vector<TH1F*>> Map; 
  for (int i(0); i < JetEnergy.size(); i++)
  {
    std::vector<TString> Dirs; 
    for (int y(0); y < Layers.size(); y++)
    {
      TString name = Layers[y]; name += ("/"); name += (JetEnergy[i]); 
      Dirs.push_back(name); 
    }
    
    Map["Track_1_" + JetEnergy[i]] = Sort(1, ntru, Dirs, Dump, "_" + JetEnergy[i]); 
    Map["Track_2_" + JetEnergy[i]] = Sort(2, ntru, Dirs, Dump, "_" + JetEnergy[i]); 
    Map["Track_3_" + JetEnergy[i]] = Sort(3, ntru, Dirs, Dump, "_" + JetEnergy[i]); 
    Map["Track_4_" + JetEnergy[i]] = Sort(4, ntru, Dirs, Dump, "_" + JetEnergy[i]); 
  }
  
  std::vector<TH1F*> trk1_all; 
  std::vector<TH1F*> trk2_all; 
  std::vector<TH1F*> trk3_all; 
  std::vector<TH1F*> trk4_all; 
  for (int i(0); i < JetEnergy.size(); i++)
  {
    std::vector<TH1F*> t1 = Map["Track_1_" + JetEnergy[i]]; 
    std::vector<TH1F*> t2 = Map["Track_2_" + JetEnergy[i]]; 
    std::vector<TH1F*> t3 = Map["Track_3_" + JetEnergy[i]]; 
    std::vector<TH1F*> t4 = Map["Track_4_" + JetEnergy[i]]; 
    
    if (i == 0)
    {
      std::vector<TString> trk1_n = GenerateNames(ntrk, ntru, "", 0); 
      std::vector<TString> trk2_n = GenerateNames(ntrk, ntru, "", 1); 
      std::vector<TString> trk3_n = GenerateNames(ntrk, ntru, "", 2); 
      std::vector<TString> trk4_n = GenerateNames(ntrk, ntru, "", 3); 
      
      trk1_all = CloneTH1F(t1[0], trk1_n); 
      trk2_all = CloneTH1F(t2[0], trk2_n);   
      trk3_all = CloneTH1F(t3[0], trk3_n);   
      trk4_all = CloneTH1F(t4[0], trk4_n);   
    }
    else 
    {
      Merge(trk1_all, t1); 
      Merge(trk2_all, t2); 
      Merge(trk3_all, t3); 
      Merge(trk4_all, t4); 
    }
  }
  Map["Track_1_All"] = trk1_all; 
  Map["Track_2_All"] = trk2_all; 
  Map["Track_3_All"] = trk3_all; 
  Map["Track_4_All"] = trk4_all; 


  // ==================== Get the eventWeight = 1 based histograms ======================== //
  std::vector<TString> Names_eW1 = GenerateNames(ntrk, ntru, "_EW_1"); 
  std::map<TString, std::vector<TString>> ToGet_eW1; 
  for (int i(0); i < Layers.size(); i++)
  {
    for (int y(0); y < JetEnergy.size(); y++)
    {
      TString name = Layers[i]; name += ("/"); name += (JetEnergy[y]); 
      ToGet[name] = Names; 
    }
  }

  std::map<TString, std::vector<TH1F*>> Dump_eW1 = GetHist(ToGet, dir); 
  for (int i(0); i < JetEnergy.size(); i++)
  {
    std::vector<TString> Dirs; 
    for (int y(0); y < Layers.size(); y++)
    {
      TString name = Layers[y]; name += ("/"); name += (JetEnergy[i]); 
      Dirs.push_back(name); 
    }
    
    Map["Track_1_" + JetEnergy[i] + "_eW1"] = Sort(1, ntru, Dirs, Dump_eW1, "_" + JetEnergy[i] + "_EW_1"); 
    Map["Track_2_" + JetEnergy[i] + "_eW1"] = Sort(2, ntru, Dirs, Dump_eW1, "_" + JetEnergy[i] + "_EW_1"); 
    Map["Track_3_" + JetEnergy[i] + "_eW1"] = Sort(3, ntru, Dirs, Dump_eW1, "_" + JetEnergy[i] + "_EW_1"); 
    Map["Track_4_" + JetEnergy[i] + "_eW1"] = Sort(4, ntru, Dirs, Dump_eW1, "_" + JetEnergy[i] + "_EW_1"); 
  }
  
  std::vector<TH1F*> trk1_all_eW1; 
  std::vector<TH1F*> trk2_all_eW1; 
  std::vector<TH1F*> trk3_all_eW1; 
  std::vector<TH1F*> trk4_all_eW1; 
  for (int i(0); i < JetEnergy.size(); i++)
  {
    std::vector<TH1F*> t1 = Map["Track_1_" + JetEnergy[i] + "_eW1"]; 
    std::vector<TH1F*> t2 = Map["Track_2_" + JetEnergy[i] + "_eW1"]; 
    std::vector<TH1F*> t3 = Map["Track_3_" + JetEnergy[i] + "_eW1"]; 
    std::vector<TH1F*> t4 = Map["Track_4_" + JetEnergy[i] + "_eW1"]; 
    
    if (i == 0)
    {
      std::vector<TString> trk1_n = GenerateNames(ntrk, ntru, "_EW_1", 0); 
      std::vector<TString> trk2_n = GenerateNames(ntrk, ntru, "_EW_1", 1); 
      std::vector<TString> trk3_n = GenerateNames(ntrk, ntru, "_EW_1", 2); 
      std::vector<TString> trk4_n = GenerateNames(ntrk, ntru, "_EW_1", 3); 
      
      trk1_all = CloneTH1F(t1[0], trk1_n); 
      trk2_all = CloneTH1F(t2[0], trk2_n);   
      trk3_all = CloneTH1F(t3[0], trk3_n);   
      trk4_all = CloneTH1F(t4[0], trk4_n);   
    }
    else 
    {
      Merge(trk1_all, t1); 
      Merge(trk2_all, t2); 
      Merge(trk3_all, t3); 
      Merge(trk4_all, t4); 
    }
  }
  Map["Track_1_All_eW1"] = trk1_all; 
  Map["Track_2_All_eW1"] = trk2_all; 
  Map["Track_3_All_eW1"] = trk3_all; 
  Map["Track_4_All_eW1"] = trk4_all; 

  return Map; 
}




