#include<PostAnalysis/IO.h>

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

std::map<TString, std::map<TString, std::vector<TH1F*>>> ReadCTIDE(TString dir)
{
 
  auto Merging =[] (std::vector<TH1F*> In, std::vector<TH1F*>* out, std::vector<TString> Names, bool Contains = false)
  {
    if (out -> size() == 0)
    {
      for (int j(0); j < Names.size(); j++)
      {
        TH1F* H = (TH1F*)In[0] -> Clone(Names[j]); 
        H -> Reset();
        H -> SetTitle(Names[j]); 
        out -> push_back(H); 
      }
    }
    else 
    {
      for (int i(0); i < In.size(); i++)
      {
        for (int k(0); k < out -> size(); k++)
        {
          TString temp = (*out)[k] -> GetTitle();
          if (temp == (In[i] -> GetTitle()) && !Contains){(*out)[k] -> Add(In[i], 1);}
          if (temp.Contains(In[i] -> GetTitle()) && Contains && !temp.Contains("Split")){(*out)[k] -> Add(In[i], 1);}
        }
      }
    }
  };  

  auto GenerateNames =[&] ( TString L, TString radius, TString Ext = "")
  {
    std::vector<std::vector<TString>> Names; 
    for (int i(0); i < 4; i++)
    {
      std::vector<TString> Temp; 
      TString nam = "dEdx_ntrk_"; nam +=(i+1); nam +=(radius); nam+=("_"+L); nam += (Ext);
      if (!radius.Contains("ntru")){Temp.push_back(nam); Names.push_back(Temp); Temp.clear(); continue;}

      for (int v(0); v < 4; v++)
      {
        TString na = "dEdx_ntrk_"; na += (i+1); na += (radius); na += (v+1); na += ("_" + L); na += (Ext);
        Temp.push_back(na); 
      }   
      Names.push_back(Temp); 
    }
    return Names;
  }; 

  std::map<TString, std::map<TString, std::vector<TH1F*>>> Output; 
  TFile* F = new TFile(dir, "READ");
  for (TString L1 : Layer)
  {
    for (TString E1 : JetEnergy)
    {
      std::map<TString, std::vector<TH1F*>> Trk_Tru; 

      F -> cd(L1 + "/" + E1 + "_radius");  
      std::vector<TString> Folders = ReturnCurrentDirs(false); 
      for (TString x : Folders)
      {
        TH1F* H = (TH1F*)gDirectory -> Get(x); 

        // Non Truth measurement. 
        if (x.Contains("ntrk_1") && !x.Contains("ntru") && x.Contains("rless") && !x.Contains("Split")){ Trk_Tru["ntrk_1_M_I"].push_back(H); }
        if (x.Contains("ntrk_2") && !x.Contains("ntru") && x.Contains("rless") && !x.Contains("Split")){ Trk_Tru["ntrk_2_M_I"].push_back(H); }
        if (x.Contains("ntrk_3") && !x.Contains("ntru") && x.Contains("rless") && !x.Contains("Split")){ Trk_Tru["ntrk_3_M_I"].push_back(H); }
        if (x.Contains("ntrk_4") && !x.Contains("ntru") && x.Contains("rless") && !x.Contains("Split")){ Trk_Tru["ntrk_4_M_I"].push_back(H); }

        // Outside jetcore measurement 
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && !x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_1_M_O"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && !x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_2_M_O"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && !x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_3_M_O"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && !x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_4_M_O"].push_back(H); }


        // Truth inside the jet core. 
        if (x.Contains("ntrk_1") && x.Contains("ntru") && x.Contains("rless")  && !x.Contains("Split")){ Trk_Tru["ntrk_1_T_I"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("ntru") && x.Contains("rless")  && !x.Contains("Split")){ Trk_Tru["ntrk_2_T_I"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("ntru") && x.Contains("rless")  && !x.Contains("Split")){ Trk_Tru["ntrk_3_T_I"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("ntru") && x.Contains("rless")  && !x.Contains("Split")){ Trk_Tru["ntrk_4_T_I"].push_back(H); }

        // Outside jetcore truth
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_1_T_O"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_2_T_O"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_3_T_O"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("Split")){ Trk_Tru["ntrk_4_T_O"].push_back(H); }

        // Split
        // Non Truth measurement. 
        if (x.Contains("ntrk_1") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("IsSplit")){ Trk_Tru["ntrk_1_M_I_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("IsSplit")){ Trk_Tru["ntrk_2_M_I_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("IsSplit")){ Trk_Tru["ntrk_3_M_I_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("IsSplit")){ Trk_Tru["ntrk_4_M_I_IsSplit"].push_back(H); }

        // Outside jetcore measurement 
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("IsSplit")){ Trk_Tru["ntrk_1_M_O_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("IsSplit")){ Trk_Tru["ntrk_2_M_O_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("IsSplit")){ Trk_Tru["ntrk_3_M_O_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("IsSplit")){ Trk_Tru["ntrk_4_M_O_IsSplit"].push_back(H); }

        // Truth inside the jet core. 
        if (x.Contains("ntrk_1") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("IsSplit")){ Trk_Tru["ntrk_1_T_I_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("IsSplit")){ Trk_Tru["ntrk_2_T_I_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("IsSplit")){ Trk_Tru["ntrk_3_T_I_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("IsSplit")){ Trk_Tru["ntrk_4_T_I_IsSplit"].push_back(H); }

        // Outside jetcore truth
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("IsSplit")){ Trk_Tru["ntrk_1_T_O_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("IsSplit")){ Trk_Tru["ntrk_2_T_O_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("IsSplit")){ Trk_Tru["ntrk_3_T_O_IsSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("IsSplit")){ Trk_Tru["ntrk_4_T_O_IsSplit"].push_back(H); }


        // No Split
        // Non Truth measurement. 
        if (x.Contains("ntrk_1") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("NotSplit")){ Trk_Tru["ntrk_1_M_I_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("NotSplit")){ Trk_Tru["ntrk_2_M_I_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("NotSplit")){ Trk_Tru["ntrk_3_M_I_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && !x.Contains("ntru") && x.Contains("rless") && x.Contains("NotSplit")){ Trk_Tru["ntrk_4_M_I_NotSplit"].push_back(H); }

        // Outside jetcore measurement 
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("NotSplit")){ Trk_Tru["ntrk_1_M_O_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("NotSplit")){ Trk_Tru["ntrk_2_M_O_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("NotSplit")){ Trk_Tru["ntrk_3_M_O_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && !x.Contains("ntru") && x.Contains("NotSplit")){ Trk_Tru["ntrk_4_M_O_NotSplit"].push_back(H); }

        // Truth inside the jet core. 
        if (x.Contains("ntrk_1") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("NotSplit")){ Trk_Tru["ntrk_1_T_I_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("NotSplit")){ Trk_Tru["ntrk_2_T_I_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("NotSplit")){ Trk_Tru["ntrk_3_T_I_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("ntru") && x.Contains("rless")  && x.Contains("NotSplit")){ Trk_Tru["ntrk_4_T_I_NotSplit"].push_back(H); }

        // Outside jetcore truth
        if (x.Contains("ntrk_1") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("NotSplit")){ Trk_Tru["ntrk_1_T_O_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("NotSplit")){ Trk_Tru["ntrk_2_T_O_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("NotSplit")){ Trk_Tru["ntrk_3_T_O_NotSplit"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && x.Contains("ntru") && !x.Contains("NotSplit")){ Trk_Tru["ntrk_4_T_O_NotSplit"].push_back(H); }       
        
      }
      Output[L1 + "_" + E1] = Trk_Tru; 
      F -> cd(); 
    }
  }
  
  for (MMVi n = Output.begin(); n != Output.end(); n++)
  {
    for (MVi h = Output[n -> first].begin(); h != Output[n -> first].end(); h++)
    {
      std::vector<TH1F*> H = Output[n -> first][h -> first]; 
      
      if ((h -> first).Contains("NotSplit") || (h -> first).Contains("IsSplit")){continue;}
      for (TH1F* k : H){k -> Rebin(5);}
    }
  }
  
  return Output; // <<<< --------------------------- Premature close

  std::map<TString, std::vector<TH1F*>> All_Hists; 
  for (MMVi x = Output.begin(); x != Output.end(); x++)
  {
    TString title = x -> first; 
    std::map<TString, std::vector<TH1F*>> Trk_Tru = x -> second; 

    for (TString L : Layer)
    {
      if ( title.Contains(L) )
      {
        // Nominal 
        std::vector<std::vector<TString>> Names_T_I = GenerateNames( L, "_rless_005_ntru_" ); 
        std::vector<std::vector<TString>> Names_M_I = GenerateNames( L, "_rless_005" ); 
        std::vector<std::vector<TString>> Names_T_O = GenerateNames( L, "_rgreater_1_ntru_" );  
        std::vector<std::vector<TString>> Names_M_O = GenerateNames( L, "_rgreater_1" );  
        
        // Split clusters
        std::vector<std::vector<TString>> Names_T_I_IsSplit = GenerateNames( L, "_rless_005_ntru_", "_IsSplit" ); 
        std::vector<std::vector<TString>> Names_M_I_IsSplit = GenerateNames( L, "_rless_005" , "_IsSplit"); 
        std::vector<std::vector<TString>> Names_T_O_IsSplit = GenerateNames( L, "_rgreater_1_ntru_" , "_IsSplit");  
        std::vector<std::vector<TString>> Names_M_O_IsSplit = GenerateNames( L, "_rgreater_1" , "_IsSplit");  
      
        // Non Split clusters
        std::vector<std::vector<TString>> Names_T_I_NotSplit = GenerateNames( L, "_rless_005_ntru_" , "_NotSplit"); 
        std::vector<std::vector<TString>> Names_M_I_NotSplit = GenerateNames( L, "_rless_005" , "_NotSplit"); 
        std::vector<std::vector<TString>> Names_T_O_NotSplit = GenerateNames( L, "_rgreater_1_ntru_" , "_NotSplit");  
        std::vector<std::vector<TString>> Names_M_O_NotSplit = GenerateNames( L, "_rgreater_1" , "_NotSplit");  
        
        // Nominal
        Merging(Trk_Tru["ntrk_1_T_I"], &Output[L]["ntrk_1_T_I"], Names_T_I[0], true); 
        Merging(Trk_Tru["ntrk_2_T_I"], &Output[L]["ntrk_2_T_I"], Names_T_I[1], true); 
        Merging(Trk_Tru["ntrk_3_T_I"], &Output[L]["ntrk_3_T_I"], Names_T_I[2], true); 
        Merging(Trk_Tru["ntrk_4_T_I"], &Output[L]["ntrk_4_T_I"], Names_T_I[3], true); 

        Merging(Trk_Tru["ntrk_1_M_I"], &Output[L]["ntrk_1_M_I"], Names_M_I[0], true); 
        Merging(Trk_Tru["ntrk_2_M_I"], &Output[L]["ntrk_2_M_I"], Names_M_I[1], true); 
        Merging(Trk_Tru["ntrk_3_M_I"], &Output[L]["ntrk_3_M_I"], Names_M_I[2], true); 
        Merging(Trk_Tru["ntrk_4_M_I"], &Output[L]["ntrk_4_M_I"], Names_M_I[3], true); 

        Merging(Trk_Tru["ntrk_1_T_O"], &Output[L]["ntrk_1_T_O"], Names_T_O[0], true);      
        Merging(Trk_Tru["ntrk_2_T_O"], &Output[L]["ntrk_2_T_O"], Names_T_O[1], true);
        Merging(Trk_Tru["ntrk_3_T_O"], &Output[L]["ntrk_3_T_O"], Names_T_O[2], true);     
        Merging(Trk_Tru["ntrk_4_T_O"], &Output[L]["ntrk_4_T_O"], Names_T_O[3], true);     

        Merging(Trk_Tru["ntrk_1_M_O"], &Output[L]["ntrk_1_M_O"], Names_M_O[0], true);      
        Merging(Trk_Tru["ntrk_2_M_O"], &Output[L]["ntrk_2_M_O"], Names_M_O[1], true);
        Merging(Trk_Tru["ntrk_3_M_O"], &Output[L]["ntrk_3_M_O"], Names_M_O[2], true);     
        Merging(Trk_Tru["ntrk_4_M_O"], &Output[L]["ntrk_4_M_O"], Names_M_O[3], true);     

        // IsSplit
        Merging(Trk_Tru["ntrk_1_T_I_IsSplit"], &Output[L]["ntrk_1_T_I_IsSplit"], Names_T_I_IsSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_T_I_IsSplit"], &Output[L]["ntrk_2_T_I_IsSplit"], Names_T_I_IsSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_T_I_IsSplit"], &Output[L]["ntrk_3_T_I_IsSplit"], Names_T_I_IsSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_T_I_IsSplit"], &Output[L]["ntrk_4_T_I_IsSplit"], Names_T_I_IsSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_M_I_IsSplit"], &Output[L]["ntrk_1_M_I_IsSplit"], Names_M_I_IsSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_M_I_IsSplit"], &Output[L]["ntrk_2_M_I_IsSplit"], Names_M_I_IsSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_M_I_IsSplit"], &Output[L]["ntrk_3_M_I_IsSplit"], Names_M_I_IsSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_M_I_IsSplit"], &Output[L]["ntrk_4_M_I_IsSplit"], Names_M_I_IsSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_T_O_IsSplit"], &Output[L]["ntrk_1_T_O_IsSplit"], Names_T_O_IsSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_T_O_IsSplit"], &Output[L]["ntrk_2_T_O_IsSplit"], Names_T_O_IsSplit[1], true);
        Merging(Trk_Tru["ntrk_3_T_O_IsSplit"], &Output[L]["ntrk_3_T_O_IsSplit"], Names_T_O_IsSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_T_O_IsSplit"], &Output[L]["ntrk_4_T_O_IsSplit"], Names_T_O_IsSplit[3], true);     

        Merging(Trk_Tru["ntrk_1_M_O_IsSplit"], &Output[L]["ntrk_1_M_O_IsSplit"], Names_M_O_IsSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_M_O_IsSplit"], &Output[L]["ntrk_2_M_O_IsSplit"], Names_M_O_IsSplit[1], true);
        Merging(Trk_Tru["ntrk_3_M_O_IsSplit"], &Output[L]["ntrk_3_M_O_IsSplit"], Names_M_O_IsSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_M_O_IsSplit"], &Output[L]["ntrk_4_M_O_IsSplit"], Names_M_O_IsSplit[3], true);     

        // NotSplit
        Merging(Trk_Tru["ntrk_1_T_I_NotSplit"], &Output[L]["ntrk_1_T_I_NotSplit"], Names_T_I_NotSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_T_I_NotSplit"], &Output[L]["ntrk_2_T_I_NotSplit"], Names_T_I_NotSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_T_I_NotSplit"], &Output[L]["ntrk_3_T_I_NotSplit"], Names_T_I_NotSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_T_I_NotSplit"], &Output[L]["ntrk_4_T_I_NotSplit"], Names_T_I_NotSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_M_I_NotSplit"], &Output[L]["ntrk_1_M_I_NotSplit"], Names_M_I_NotSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_M_I_NotSplit"], &Output[L]["ntrk_2_M_I_NotSplit"], Names_M_I_NotSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_M_I_NotSplit"], &Output[L]["ntrk_3_M_I_NotSplit"], Names_M_I_NotSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_M_I_NotSplit"], &Output[L]["ntrk_4_M_I_NotSplit"], Names_M_I_NotSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_T_O_NotSplit"], &Output[L]["ntrk_1_T_O_NotSplit"], Names_T_O_NotSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_T_O_NotSplit"], &Output[L]["ntrk_2_T_O_NotSplit"], Names_T_O_NotSplit[1], true);
        Merging(Trk_Tru["ntrk_3_T_O_NotSplit"], &Output[L]["ntrk_3_T_O_NotSplit"], Names_T_O_NotSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_T_O_NotSplit"], &Output[L]["ntrk_4_T_O_NotSplit"], Names_T_O_NotSplit[3], true);     

        Merging(Trk_Tru["ntrk_1_M_O_NotSplit"], &Output[L]["ntrk_1_M_O_NotSplit"], Names_M_O_NotSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_M_O_NotSplit"], &Output[L]["ntrk_2_M_O_NotSplit"], Names_M_O_NotSplit[1], true);
        Merging(Trk_Tru["ntrk_3_M_O_NotSplit"], &Output[L]["ntrk_3_M_O_NotSplit"], Names_M_O_NotSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_M_O_NotSplit"], &Output[L]["ntrk_4_M_O_NotSplit"], Names_M_O_NotSplit[3], true);     

      }
    }

    for (TString L : JetEnergy)
    {
      if ( title.Contains(L) )
      {
        // Nominal 
        std::vector<std::vector<TString>> Names_T_I = GenerateNames( L, "_rless_005_ntru_" ); 
        std::vector<std::vector<TString>> Names_M_I = GenerateNames( L, "_rless_005" ); 
        std::vector<std::vector<TString>> Names_T_O = GenerateNames( L, "_rgreater_1_ntru_" );  
        std::vector<std::vector<TString>> Names_M_O = GenerateNames( L, "_rgreater_1" );  
        
        // Split clusters
        std::vector<std::vector<TString>> Names_T_I_IsSplit = GenerateNames( L, "_rless_005_ntru_", "_IsSplit" ); 
        std::vector<std::vector<TString>> Names_M_I_IsSplit = GenerateNames( L, "_rless_005" , "_IsSplit"); 
        std::vector<std::vector<TString>> Names_T_O_IsSplit = GenerateNames( L, "_rgreater_1_ntru_" , "_IsSplit");  
        std::vector<std::vector<TString>> Names_M_O_IsSplit = GenerateNames( L, "_rgreater_1" , "_IsSplit");  
      
        // Non Split clusters
        std::vector<std::vector<TString>> Names_T_I_NotSplit = GenerateNames( L, "_rless_005_ntru_" , "_NotSplit"); 
        std::vector<std::vector<TString>> Names_M_I_NotSplit = GenerateNames( L, "_rless_005" , "_NotSplit"); 
        std::vector<std::vector<TString>> Names_T_O_NotSplit = GenerateNames( L, "_rgreater_1_ntru_" , "_NotSplit");  
        std::vector<std::vector<TString>> Names_M_O_NotSplit = GenerateNames( L, "_rgreater_1" , "_NotSplit");  
        
        // Nominal
        Merging(Trk_Tru["ntrk_1_T_I"], &Output[L]["ntrk_1_T_I"], Names_T_I[0], true); 
        Merging(Trk_Tru["ntrk_2_T_I"], &Output[L]["ntrk_2_T_I"], Names_T_I[1], true); 
        Merging(Trk_Tru["ntrk_3_T_I"], &Output[L]["ntrk_3_T_I"], Names_T_I[2], true); 
        Merging(Trk_Tru["ntrk_4_T_I"], &Output[L]["ntrk_4_T_I"], Names_T_I[3], true); 

        Merging(Trk_Tru["ntrk_1_M_I"], &Output[L]["ntrk_1_M_I"], Names_M_I[0], true); 
        Merging(Trk_Tru["ntrk_2_M_I"], &Output[L]["ntrk_2_M_I"], Names_M_I[1], true); 
        Merging(Trk_Tru["ntrk_3_M_I"], &Output[L]["ntrk_3_M_I"], Names_M_I[2], true); 
        Merging(Trk_Tru["ntrk_4_M_I"], &Output[L]["ntrk_4_M_I"], Names_M_I[3], true); 

        Merging(Trk_Tru["ntrk_1_T_O"], &Output[L]["ntrk_1_T_O"], Names_T_O[0], true);      
        Merging(Trk_Tru["ntrk_2_T_O"], &Output[L]["ntrk_2_T_O"], Names_T_O[1], true);
        Merging(Trk_Tru["ntrk_3_T_O"], &Output[L]["ntrk_3_T_O"], Names_T_O[2], true);     
        Merging(Trk_Tru["ntrk_4_T_O"], &Output[L]["ntrk_4_T_O"], Names_T_O[3], true);     

        Merging(Trk_Tru["ntrk_1_M_O"], &Output[L]["ntrk_1_M_O"], Names_M_O[0], true);      
        Merging(Trk_Tru["ntrk_2_M_O"], &Output[L]["ntrk_2_M_O"], Names_M_O[1], true);
        Merging(Trk_Tru["ntrk_3_M_O"], &Output[L]["ntrk_3_M_O"], Names_M_O[2], true);     
        Merging(Trk_Tru["ntrk_4_M_O"], &Output[L]["ntrk_4_M_O"], Names_M_O[3], true);     

        // IsSplit
        Merging(Trk_Tru["ntrk_1_T_I_IsSplit"], &Output[L]["ntrk_1_T_I_IsSplit"], Names_T_I_IsSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_T_I_IsSplit"], &Output[L]["ntrk_2_T_I_IsSplit"], Names_T_I_IsSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_T_I_IsSplit"], &Output[L]["ntrk_3_T_I_IsSplit"], Names_T_I_IsSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_T_I_IsSplit"], &Output[L]["ntrk_4_T_I_IsSplit"], Names_T_I_IsSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_M_I_IsSplit"], &Output[L]["ntrk_1_M_I_IsSplit"], Names_M_I_IsSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_M_I_IsSplit"], &Output[L]["ntrk_2_M_I_IsSplit"], Names_M_I_IsSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_M_I_IsSplit"], &Output[L]["ntrk_3_M_I_IsSplit"], Names_M_I_IsSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_M_I_IsSplit"], &Output[L]["ntrk_4_M_I_IsSplit"], Names_M_I_IsSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_T_O_IsSplit"], &Output[L]["ntrk_1_T_O_IsSplit"], Names_T_O_IsSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_T_O_IsSplit"], &Output[L]["ntrk_2_T_O_IsSplit"], Names_T_O_IsSplit[1], true);
        Merging(Trk_Tru["ntrk_3_T_O_IsSplit"], &Output[L]["ntrk_3_T_O_IsSplit"], Names_T_O_IsSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_T_O_IsSplit"], &Output[L]["ntrk_4_T_O_IsSplit"], Names_T_O_IsSplit[3], true);     

        Merging(Trk_Tru["ntrk_1_M_O_IsSplit"], &Output[L]["ntrk_1_M_O_IsSplit"], Names_M_O_IsSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_M_O_IsSplit"], &Output[L]["ntrk_2_M_O_IsSplit"], Names_M_O_IsSplit[1], true);
        Merging(Trk_Tru["ntrk_3_M_O_IsSplit"], &Output[L]["ntrk_3_M_O_IsSplit"], Names_M_O_IsSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_M_O_IsSplit"], &Output[L]["ntrk_4_M_O_IsSplit"], Names_M_O_IsSplit[3], true);     

        // NotSplit
        Merging(Trk_Tru["ntrk_1_T_I_NotSplit"], &Output[L]["ntrk_1_T_I_NotSplit"], Names_T_I_NotSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_T_I_NotSplit"], &Output[L]["ntrk_2_T_I_NotSplit"], Names_T_I_NotSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_T_I_NotSplit"], &Output[L]["ntrk_3_T_I_NotSplit"], Names_T_I_NotSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_T_I_NotSplit"], &Output[L]["ntrk_4_T_I_NotSplit"], Names_T_I_NotSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_M_I_NotSplit"], &Output[L]["ntrk_1_M_I_NotSplit"], Names_M_I_NotSplit[0], true); 
        Merging(Trk_Tru["ntrk_2_M_I_NotSplit"], &Output[L]["ntrk_2_M_I_NotSplit"], Names_M_I_NotSplit[1], true); 
        Merging(Trk_Tru["ntrk_3_M_I_NotSplit"], &Output[L]["ntrk_3_M_I_NotSplit"], Names_M_I_NotSplit[2], true); 
        Merging(Trk_Tru["ntrk_4_M_I_NotSplit"], &Output[L]["ntrk_4_M_I_NotSplit"], Names_M_I_NotSplit[3], true); 

        Merging(Trk_Tru["ntrk_1_T_O_NotSplit"], &Output[L]["ntrk_1_T_O_NotSplit"], Names_T_O_NotSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_T_O_NotSplit"], &Output[L]["ntrk_2_T_O_NotSplit"], Names_T_O_NotSplit[1], true);
        Merging(Trk_Tru["ntrk_3_T_O_NotSplit"], &Output[L]["ntrk_3_T_O_NotSplit"], Names_T_O_NotSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_T_O_NotSplit"], &Output[L]["ntrk_4_T_O_NotSplit"], Names_T_O_NotSplit[3], true);     

        Merging(Trk_Tru["ntrk_1_M_O_NotSplit"], &Output[L]["ntrk_1_M_O_NotSplit"], Names_M_O_NotSplit[0], true);      
        Merging(Trk_Tru["ntrk_2_M_O_NotSplit"], &Output[L]["ntrk_2_M_O_NotSplit"], Names_M_O_NotSplit[1], true);
        Merging(Trk_Tru["ntrk_3_M_O_NotSplit"], &Output[L]["ntrk_3_M_O_NotSplit"], Names_M_O_NotSplit[2], true);     
        Merging(Trk_Tru["ntrk_4_M_O_NotSplit"], &Output[L]["ntrk_4_M_O_NotSplit"], Names_M_O_NotSplit[3], true);    

      }
    }  


    std::vector<std::vector<TString>> Names_T_I_A = GenerateNames("All", "_rless_005_ntru_" ); 
    std::vector<std::vector<TString>> Names_M_I_A = GenerateNames("All", "_rless_005" ); 
    std::vector<std::vector<TString>> Names_T_O_A = GenerateNames("All", "_rgreater_1_ntru_" );  
    std::vector<std::vector<TString>> Names_M_O_A = GenerateNames("All", "_rgreater_1" );  

    // Split clusters
    std::vector<std::vector<TString>> Names_T_I_A_IsSplit = GenerateNames("All", "_rless_005_ntru_", "_IsSplit" ); 
    std::vector<std::vector<TString>> Names_M_I_A_IsSplit = GenerateNames("All", "_rless_005" , "_IsSplit"); 
    std::vector<std::vector<TString>> Names_T_O_A_IsSplit = GenerateNames("All", "_rgreater_1_ntru_" , "_IsSplit");  
    std::vector<std::vector<TString>> Names_M_O_A_IsSplit = GenerateNames("All", "_rgreater_1" , "_IsSplit");  
    
    // Non Split clusters
    std::vector<std::vector<TString>> Names_T_I_A_NotSplit = GenerateNames("All", "_rless_005_ntru_" , "_NotSplit"); 
    std::vector<std::vector<TString>> Names_M_I_A_NotSplit = GenerateNames("All", "_rless_005" , "_NotSplit"); 
    std::vector<std::vector<TString>> Names_T_O_A_NotSplit = GenerateNames("All", "_rgreater_1_ntru_" , "_NotSplit");  
    std::vector<std::vector<TString>> Names_M_O_A_NotSplit = GenerateNames("All", "_rgreater_1" , "_NotSplit");  

    Merging(Trk_Tru["ntrk_1_T_I"], &Output["All"]["ntrk_1_T_I"], Names_T_I_A[0], true); 
    Merging(Trk_Tru["ntrk_2_T_I"], &Output["All"]["ntrk_2_T_I"], Names_T_I_A[1], true); 
    Merging(Trk_Tru["ntrk_3_T_I"], &Output["All"]["ntrk_3_T_I"], Names_T_I_A[2], true); 
    Merging(Trk_Tru["ntrk_4_T_I"], &Output["All"]["ntrk_4_T_I"], Names_T_I_A[3], true); 

    Merging(Trk_Tru["ntrk_1_M_I"], &Output["All"]["ntrk_1_M_I"], Names_M_I_A[0], true); 
    Merging(Trk_Tru["ntrk_2_M_I"], &Output["All"]["ntrk_2_M_I"], Names_M_I_A[1], true); 
    Merging(Trk_Tru["ntrk_3_M_I"], &Output["All"]["ntrk_3_M_I"], Names_M_I_A[2], true); 
    Merging(Trk_Tru["ntrk_4_M_I"], &Output["All"]["ntrk_4_M_I"], Names_M_I_A[3], true); 

    Merging(Trk_Tru["ntrk_1_T_O"], &Output["All"]["ntrk_1_T_O"], Names_T_O_A[0], true);      
    Merging(Trk_Tru["ntrk_2_T_O"], &Output["All"]["ntrk_2_T_O"], Names_T_O_A[1], true);
    Merging(Trk_Tru["ntrk_3_T_O"], &Output["All"]["ntrk_3_T_O"], Names_T_O_A[2], true);     
    Merging(Trk_Tru["ntrk_4_T_O"], &Output["All"]["ntrk_4_T_O"], Names_T_O_A[3], true);     

    Merging(Trk_Tru["ntrk_1_M_O"], &Output["All"]["ntrk_1_M_O"], Names_M_O_A[0], true);      
    Merging(Trk_Tru["ntrk_2_M_O"], &Output["All"]["ntrk_2_M_O"], Names_M_O_A[1], true);
    Merging(Trk_Tru["ntrk_3_M_O"], &Output["All"]["ntrk_3_M_O"], Names_M_O_A[2], true);     
    Merging(Trk_Tru["ntrk_4_M_O"], &Output["All"]["ntrk_4_M_O"], Names_M_O_A[3], true);     


    // IsSplit
    Merging(Trk_Tru["ntrk_1_T_I_IsSplit"], &Output["All"]["ntrk_1_T_I_IsSplit"], Names_T_I_A_IsSplit[0], true); 
    Merging(Trk_Tru["ntrk_2_T_I_IsSplit"], &Output["All"]["ntrk_2_T_I_IsSplit"], Names_T_I_A_IsSplit[1], true); 
    Merging(Trk_Tru["ntrk_3_T_I_IsSplit"], &Output["All"]["ntrk_3_T_I_IsSplit"], Names_T_I_A_IsSplit[2], true); 
    Merging(Trk_Tru["ntrk_4_T_I_IsSplit"], &Output["All"]["ntrk_4_T_I_IsSplit"], Names_T_I_A_IsSplit[3], true); 

    Merging(Trk_Tru["ntrk_1_M_I_IsSplit"], &Output["All"]["ntrk_1_M_I_IsSplit"], Names_M_I_A_IsSplit[0], true); 
    Merging(Trk_Tru["ntrk_2_M_I_IsSplit"], &Output["All"]["ntrk_2_M_I_IsSplit"], Names_M_I_A_IsSplit[1], true); 
    Merging(Trk_Tru["ntrk_3_M_I_IsSplit"], &Output["All"]["ntrk_3_M_I_IsSplit"], Names_M_I_A_IsSplit[2], true); 
    Merging(Trk_Tru["ntrk_4_M_I_IsSplit"], &Output["All"]["ntrk_4_M_I_IsSplit"], Names_M_I_A_IsSplit[3], true); 

    Merging(Trk_Tru["ntrk_1_T_O_IsSplit"], &Output["All"]["ntrk_1_T_O_IsSplit"], Names_T_O_A_IsSplit[0], true);      
    Merging(Trk_Tru["ntrk_2_T_O_IsSplit"], &Output["All"]["ntrk_2_T_O_IsSplit"], Names_T_O_A_IsSplit[1], true);
    Merging(Trk_Tru["ntrk_3_T_O_IsSplit"], &Output["All"]["ntrk_3_T_O_IsSplit"], Names_T_O_A_IsSplit[2], true);     
    Merging(Trk_Tru["ntrk_4_T_O_IsSplit"], &Output["All"]["ntrk_4_T_O_IsSplit"], Names_T_O_A_IsSplit[3], true);     

    Merging(Trk_Tru["ntrk_1_M_O_IsSplit"], &Output["All"]["ntrk_1_M_O_IsSplit"], Names_M_O_A_IsSplit[0], true);      
    Merging(Trk_Tru["ntrk_2_M_O_IsSplit"], &Output["All"]["ntrk_2_M_O_IsSplit"], Names_M_O_A_IsSplit[1], true);
    Merging(Trk_Tru["ntrk_3_M_O_IsSplit"], &Output["All"]["ntrk_3_M_O_IsSplit"], Names_M_O_A_IsSplit[2], true);     
    Merging(Trk_Tru["ntrk_4_M_O_IsSplit"], &Output["All"]["ntrk_4_M_O_IsSplit"], Names_M_O_A_IsSplit[3], true);     

    // NotSplit
    Merging(Trk_Tru["ntrk_1_T_I_NotSplit"], &Output["All"]["ntrk_1_T_I_NotSplit"], Names_T_I_A_NotSplit[0], true); 
    Merging(Trk_Tru["ntrk_2_T_I_NotSplit"], &Output["All"]["ntrk_2_T_I_NotSplit"], Names_T_I_A_NotSplit[1], true); 
    Merging(Trk_Tru["ntrk_3_T_I_NotSplit"], &Output["All"]["ntrk_3_T_I_NotSplit"], Names_T_I_A_NotSplit[2], true); 
    Merging(Trk_Tru["ntrk_4_T_I_NotSplit"], &Output["All"]["ntrk_4_T_I_NotSplit"], Names_T_I_A_NotSplit[3], true); 

    Merging(Trk_Tru["ntrk_1_M_I_NotSplit"], &Output["All"]["ntrk_1_M_I_NotSplit"], Names_M_I_A_NotSplit[0], true); 
    Merging(Trk_Tru["ntrk_2_M_I_NotSplit"], &Output["All"]["ntrk_2_M_I_NotSplit"], Names_M_I_A_NotSplit[1], true); 
    Merging(Trk_Tru["ntrk_3_M_I_NotSplit"], &Output["All"]["ntrk_3_M_I_NotSplit"], Names_M_I_A_NotSplit[2], true); 
    Merging(Trk_Tru["ntrk_4_M_I_NotSplit"], &Output["All"]["ntrk_4_M_I_NotSplit"], Names_M_I_A_NotSplit[3], true); 

    Merging(Trk_Tru["ntrk_1_T_O_NotSplit"], &Output["All"]["ntrk_1_T_O_NotSplit"], Names_T_O_A_NotSplit[0], true);      
    Merging(Trk_Tru["ntrk_2_T_O_NotSplit"], &Output["All"]["ntrk_2_T_O_NotSplit"], Names_T_O_A_NotSplit[1], true);
    Merging(Trk_Tru["ntrk_3_T_O_NotSplit"], &Output["All"]["ntrk_3_T_O_NotSplit"], Names_T_O_A_NotSplit[2], true);     
    Merging(Trk_Tru["ntrk_4_T_O_NotSplit"], &Output["All"]["ntrk_4_T_O_NotSplit"], Names_T_O_A_NotSplit[3], true);     

    Merging(Trk_Tru["ntrk_1_M_O_NotSplit"], &Output["All"]["ntrk_1_M_O_NotSplit"], Names_M_O_A_NotSplit[0], true);      
    Merging(Trk_Tru["ntrk_2_M_O_NotSplit"], &Output["All"]["ntrk_2_M_O_NotSplit"], Names_M_O_A_NotSplit[1], true);
    Merging(Trk_Tru["ntrk_3_M_O_NotSplit"], &Output["All"]["ntrk_3_M_O_NotSplit"], Names_M_O_A_NotSplit[2], true);     
    Merging(Trk_Tru["ntrk_4_M_O_NotSplit"], &Output["All"]["ntrk_4_M_O_NotSplit"], Names_M_O_A_NotSplit[3], true);    


  }

  return Output; 
}

void WriteHistsToFile(std::vector<TH1F*> ntrk_ntru, TString dir)
{
  gDirectory -> cd("/");
  gDirectory -> mkdir(dir); 
  gDirectory -> cd(dir); 
  BulkWrite(ntrk_ntru); 
  gDirectory -> cd("/"); 
}


void WriteOutputMapToFile(std::map<TString, std::vector<float>> Map, TString dir, TString name)
{
  gDirectory -> cd("/");
  gDirectory -> cd(dir); 

  TTree* tree = new TTree(name, name); 
  for(MVFi x = Map.begin(); x != Map.end(); x++)
  {
    TString key = name + "_" + (x -> first); 
    tree -> Branch(key, &Map[x -> first]);
  }
  tree -> Fill(); 
  tree -> Write();

  delete tree;
}
