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


void WriteHistsToFile(std::vector<TH1F*> ntrk_ntru, TString dir)
{
  gDirectory -> mkdir(dir); 
  gDirectory -> cd(dir); 

  BulkWrite(ntrk_ntru); 
  
  gDirectory -> cd("../"); 
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
        TH1F* H = (TH1F*)gDirectory -> Get(H_TS);

        // Truth inside the jet core. 
        if (H_TS.Contains("ntrk_1")){ Algorithms_Map[Alg_Folder + "_ntrk_1"].push_back(H); }
        if (H_TS.Contains("ntrk_2")){ Algorithms_Map[Alg_Folder + "_ntrk_2"].push_back(H); }
        if (H_TS.Contains("ntrk_3")){ Algorithms_Map[Alg_Folder + "_ntrk_3"].push_back(H); }
        if (H_TS.Contains("ntrk_4")){ Algorithms_Map[Alg_Folder + "_ntrk_4"].push_back(H); }
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
      std::cout << current << "/" << p -> first << std::endl; 
      std::vector<TH1F*> Hists = Algs[p -> first]; 
      for (TH1F* J : Hists)
      {
        std::cout << "----- " << J -> GetTitle() << std::endl;
      }
    }
  }
}

void WriteOutputMapToFile(std::map<TString, std::vector<float>> Map, TString dir, TString name)
{
  gDirectory -> mkdir(dir);
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

  gDirectory -> cd("../"); 
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
      
      std::vector<TString> Alg_V = ReturnCurrentDirs(); 
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

              std::cout << k << std::endl;
            }
            delete v;  
          }
        }
      }
      Algorithms_Map[Alg_Folder] = key_map; 
    }
    Output[Folder] = Algorithms_Map; 
    F -> cd(); 
  }

  return Output;  
}
  
















