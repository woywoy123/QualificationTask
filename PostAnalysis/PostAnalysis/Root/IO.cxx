#include<PostAnalysis/IO.h>

std::vector<TString> ReturnCurrentDirs(bool FolderOnly = true)
{
  std::vector<TString> Output; 
  TDirectory* dir = gDirectory; 
  for (TObject* key : *dir -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key); 
    TString dir = (TString)k -> GetName();
    
    if (dir.Contains("CX") || dir.Contains("temp")){ continue;}
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
      std::vector<TString> Folders = ReturnCurrentDirs(false); 
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
        if (x.Contains("ntrk_2") && x.Contains("rgreater") && !x.Contains("ntru")){ Trk_Tru["ntrk_2_M_O"].push_back(H); }
        if (x.Contains("ntrk_3") && x.Contains("rgreater") && !x.Contains("ntru")){ Trk_Tru["ntrk_3_M_O"].push_back(H); }
        if (x.Contains("ntrk_4") && x.Contains("rgreater") && !x.Contains("ntru")){ Trk_Tru["ntrk_4_M_O"].push_back(H); }
      }
      Output[L1 + "_" + E1] = Trk_Tru; 
      F -> cd(); 
    }
  }

  std::map<TString, std::vector<TH1F*>> All_Hists; 
  for (MMVi x = Output.begin(); x != Output.end(); x++)
  {
    auto Merging =[] (std::vector<TH1F*> In, std::vector<TH1F*>* out, std::vector<TString> Names)
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
            if (temp.Contains(In[i] -> GetTitle()))
            {
              (*out)[k] -> Add(In[i], 1); 
            }
          }
        }
      }
    }; 



    TString title = x -> first; 
    std::map<TString, std::vector<TH1F*>> Trk_Tru = x -> second; 

    for (TString L : Layer)
    {
      if ( title.Contains(L) )
      {
        std::vector<std::vector<TString>> Names_T; 
        std::vector<std::vector<TString>> Names_M; 
        for (int t(0); t < 4; t++)
        {
          std::vector<TString> Temp_T;  
          for (int v(0); v < 4; v++)
          {
            // Inside truth
            TString na = "dEdx_ntrk_"; na +=(t+1); na +=("_rless_005_ntru_"); na +=(v+1); na+=("_"+L); 
            Temp_T.push_back(na); 
          }
          Names_T.push_back(Temp_T);
          Temp_T.clear();
          
          TString nam = "dEdx_ntrk_"; nam +=(t+1); nam +=("_rless_005"); nam+=("_"+L); 
          Temp_T.push_back(nam);  
          Names_M.push_back(Temp_T);
        }

        Merging(Trk_Tru["ntrk_1_T_I"], &Output[L]["ntrk_1_T_I"], Names_T[0]); 
        Merging(Trk_Tru["ntrk_2_T_I"], &Output[L]["ntrk_2_T_I"], Names_T[1]); 
        Merging(Trk_Tru["ntrk_3_T_I"], &Output[L]["ntrk_3_T_I"], Names_T[2]); 
        Merging(Trk_Tru["ntrk_4_T_I"], &Output[L]["ntrk_4_T_I"], Names_T[3]); 

        Merging(Trk_Tru["ntrk_1_M_I"], &Output[L]["ntrk_1_M_I"], Names_M[0]); 
        Merging(Trk_Tru["ntrk_2_M_I"], &Output[L]["ntrk_2_M_I"], Names_M[1]); 
        Merging(Trk_Tru["ntrk_3_M_I"], &Output[L]["ntrk_3_M_I"], Names_M[2]); 
        Merging(Trk_Tru["ntrk_4_M_I"], &Output[L]["ntrk_4_M_I"], Names_M[3]); 
       
        TString na = "dEdx_ntrk_1_rgreater_1_" + L; 
        Merging(Trk_Tru["ntrk_1_M_O"], &Output[L]["ntrk_1_M_O"], {na});      
      }
    }

    for (TString L : JetEnergy)
    {
      if ( title.Contains(L) )
      {
        std::vector<std::vector<TString>> Names_T; 
        std::vector<std::vector<TString>> Names_M; 
        for (int t(0); t < 4; t++)
        {
          std::vector<TString> Temp_T;  
          for (int v(0); v < 4; v++)
          {
            // Inside truth
            TString na = "dEdx_ntrk_"; na +=(t+1); na +=("_rless_005_ntru_"); na +=(v+1); na+=("_"+L); 
            Temp_T.push_back(na); 
          }
          Names_T.push_back(Temp_T);
          Temp_T.clear(); 
          
          TString nam = "dEdx_ntrk_"; nam +=(t+1); nam +=("_rless_005"); nam+=("_"+L); 
          Temp_T.push_back(nam);  
          Names_M.push_back(Temp_T);
        }

        Merging(Trk_Tru["ntrk_1_T_I"], &Output[L]["ntrk_1_T_I"], Names_T[0]); 
        Merging(Trk_Tru["ntrk_2_T_I"], &Output[L]["ntrk_2_T_I"], Names_T[1]); 
        Merging(Trk_Tru["ntrk_3_T_I"], &Output[L]["ntrk_3_T_I"], Names_T[2]); 
        Merging(Trk_Tru["ntrk_4_T_I"], &Output[L]["ntrk_4_T_I"], Names_T[3]); 

        Merging(Trk_Tru["ntrk_1_M_I"], &Output[L]["ntrk_1_M_I"], Names_M[0]); 
        Merging(Trk_Tru["ntrk_2_M_I"], &Output[L]["ntrk_2_M_I"], Names_M[1]); 
        Merging(Trk_Tru["ntrk_3_M_I"], &Output[L]["ntrk_3_M_I"], Names_M[2]); 
        Merging(Trk_Tru["ntrk_4_M_I"], &Output[L]["ntrk_4_M_I"], Names_M[3]); 
        
        TString na = "dEdx_ntrk_1_rgreater_1_" + L; 
        Merging(Trk_Tru["ntrk_1_M_O"], &Output[L]["ntrk_1_M_O"], {na});      
      }
    }  

    std::vector<std::vector<TString>> Names_AT; 
    std::vector<std::vector<TString>> Names_AM; 
    for (int t(0); t < 4; t++)
    {
      std::vector<TString> Temp_T;  
      for (int v(0); v < 4; v++)
      {
        // Inside truth
        TString na = "dEdx_ntrk_"; na +=(t+1); na +=("_rless_005_ntru_"); na +=(v+1); na+=("_All"); 
        Temp_T.push_back(na); 
      }
      Names_AT.push_back(Temp_T);
      Temp_T.clear(); 
      
      TString nam = "dEdx_ntrk_"; nam +=(t+1); nam +=("_rless_005"); nam+=("_All"); 
      Temp_T.push_back(nam);  
      Names_AM.push_back(Temp_T);
    }

    Merging(Trk_Tru["ntrk_1_T_I"], &Output["All"]["ntrk_1_T_I"], Names_AT[0]); 
    Merging(Trk_Tru["ntrk_2_T_I"], &Output["All"]["ntrk_2_T_I"], Names_AT[1]); 
    Merging(Trk_Tru["ntrk_3_T_I"], &Output["All"]["ntrk_3_T_I"], Names_AT[2]); 
    Merging(Trk_Tru["ntrk_4_T_I"], &Output["All"]["ntrk_4_T_I"], Names_AT[3]); 
  
    Merging(Trk_Tru["ntrk_1_M_I"], &Output["All"]["ntrk_1_M_I"], Names_AM[0]); 
    Merging(Trk_Tru["ntrk_2_M_I"], &Output["All"]["ntrk_2_M_I"], Names_AM[1]); 
    Merging(Trk_Tru["ntrk_3_M_I"], &Output["All"]["ntrk_3_M_I"], Names_AM[2]); 
    Merging(Trk_Tru["ntrk_4_M_I"], &Output["All"]["ntrk_4_M_I"], Names_AM[3]); 
    
    TString na = "dEdx_ntrk_1_rgreater_1_All"; 
    Merging(Trk_Tru["ntrk_1_M_O"], &Output["All"]["ntrk_1_M_O"], {na});      
  }

  return Output; 
}

void TestReadCTIDE(TString dir)
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(dir); 
  
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
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
    
  }
}


void WriteHistsToFile(std::vector<TH1F*> ntrk_ntru, TString dir)
{
  gDirectory -> cd("/");
  gDirectory -> mkdir(dir); 
  gDirectory -> cd(dir); 
  BulkWrite(ntrk_ntru); 
  gDirectory -> cd("/"); 
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

void TestReadAlgorithm()
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Map = ReadAlgorithmResults();
  for (MMVi x = Map.begin(); x != Map.end(); x++)
  {
    TString current = x -> first; 
    std::map<TString, std::vector<TH1F*>> Algs = x -> second; 
    
    for ( MVi p = Algs.begin(); p != Algs.end(); p++)
    {
      std::vector<TH1F*> Hists = Algs[p -> first]; 
    }
  }
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

std::map<TString, std::map<TString, std::map<TString, std::vector<float>>>> ReadOutputFileToMap(TString dir) 
{
  std::map<TString, std::map<TString, std::map<TString, std::vector<float>>>> Output; 
  TFile* F = new TFile(dir, "READ");
  std::vector<TString> JetEnergy_Layer_Folder = ReturnCurrentDirs(); 
  
  for (TString Folder : JetEnergy_Layer_Folder)
  {
    F -> cd(Folder);
    std::map<TString, std::map<TString, std::vector<float>>> Algorithms_Map;  
    
    std::vector<TString> Algorithm_Folder = ReturnCurrentDirs(); 
    for (TString Alg_Folder : Algorithm_Folder)
    {
      F -> cd(Folder + "/" + Alg_Folder); 
      
      std::vector<TString> Alg_V = ReturnCurrentDirs(false); 
      std::map<TString, std::vector<float>> key_map; 
      for (TString H_TS : Alg_V)
      {
        if (H_TS.Contains("_error")) 
        {
          TTree* T = (TTree*)F -> Get(Folder + "/" + Alg_Folder + "/" + H_TS); 
         
          std::vector<TString> keys; 
          TObjArray* List = T -> GetListOfBranches(); 
          for (int i(0); i < List -> GetEntries(); i++){keys.push_back( (List -> At(i)) -> GetName() );}
          
          for (TString k : keys)
          {
            std::vector<float>* v = 0; 
            T -> SetBranchAddress(k, &v); 
            for (int t(0); t < (int)T -> GetEntries(); t++)
            {
              T -> GetEntry(t); 
              key_map[k] = *v; 
            }
            delete v;  
          }
        }
      }
      if (Alg_Folder.Contains("ntrk")){continue;}
      Algorithms_Map[Alg_Folder] = key_map; 
    }
    Output[Folder] = Algorithms_Map; 
    F -> cd(); 
  }
  delete F;

  return Output;  
}
  
















