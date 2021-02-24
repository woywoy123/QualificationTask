#include<PostAnalysis/AlgorithmTest.h>

void TestAlgorithmMonteCarlo()
{
  float nsdo = 4; 
  float m = 0.1; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params["iterations"] = {5}; 
  Params["LR_iterations"] = {25}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 
  Params["x_range"] = {0.1, 11.5}; 
  
  TString Dir = "Merged_MC.root"; 
  std::map<TString, std::vector<TH1F*>> Output = MC_Reader(Dir);
 
  TFile* F = new TFile("Results_MC.root", "RECREATE"); 
  F -> ReOpen("UPDATE"); 
 
  typedef std::map<TString, std::vector<TH1F*>>::iterator M;
  for (M x = Output.begin(); x != Output.end(); x++)
  {
    TString name = x -> first; 
    std::vector<TH1F*> Trks = x -> second; 
    
    if (name.Contains("_Less") || name.Contains("Greater")){continue;}
    std::vector<TH1F*> Outside_Core = Output[name + "_radius_Greater"]; 
    
    std::vector<TH1F*> T1;  
    std::vector<TH1F*> T2;  
    std::vector<TH1F*> T3;  
    std::vector<TH1F*> T4;  
    for (int i(0); i < nsdo; i++){T1.push_back(Trks[i]);}
    for (int i(4); i < nsdo+4; i++){T2.push_back(Trks[i]);}
    for (int i(8); i < nsdo+8; i++){T3.push_back(Trks[i]);}
    for (int i(12); i < nsdo+12; i++){T4.push_back(Trks[i]);} 
    
    TH1F* Outside_JetCore = SumHists(Outside_Core, name + "_Outside_JetCore");  
    TH1F* Trk1 = SumHists(T1, x -> first + "Track1");
    TH1F* Trk2 = SumHists(T2, x -> first + "Track2");
    TH1F* Trk3 = SumHists(T3, x -> first + "Track3");
    TH1F* Trk4 = SumHists(T4, x -> first + "Track4");
    std::vector<TH1F*> Data = {Trk1, Trk2, Trk3, Trk4}; 
  
    std::cout << " :::: " << name << std::endl;  
    bool kill = false; 
    for (TH1F* O : Data)
    {
      std::cout << O -> GetEntries() << std::endl;
      if (O -> GetEntries() < 1000){kill = true; }
    }
    if (kill == true)
    {
      BulkDelete(Trks); 
      BulkDelete(Outside_Core); 
      delete Outside_JetCore;
      continue;
    }

    std::map<TString, std::vector<TH1F*>> trk1_res = MainAlgorithm(Outside_JetCore, Data, Params, 0); 
    std::map<TString, std::vector<TH1F*>> trk2_res = MainAlgorithm(Outside_JetCore, Data, Params, 1);    
    std::map<TString, std::vector<TH1F*>> trk3_res = MainAlgorithm(Outside_JetCore, Data, Params, 2); 
    std::map<TString, std::vector<TH1F*>> trk4_res = MainAlgorithm(Outside_JetCore, Data, Params, 3); 

    M p = trk1_res.end(); 
    p--; 
    p--; 

    std::vector<TH1F*> trk1_r = trk1_res[p -> first]; 
    std::vector<TH1F*> trk2_r = trk2_res[p -> first]; 
    std::vector<TH1F*> trk3_r = trk3_res[p -> first]; 
    std::vector<TH1F*> trk4_r = trk4_res[p -> first]; 
      
    BulkWrite(trk1_r); 
    BulkWrite(trk2_r); 
    BulkWrite(trk3_r); 
    BulkWrite(trk4_r); 
    trk1_res["Data"][0] -> Write();
    trk2_res["Data"][0] -> Write();
    trk3_res["Data"][0] -> Write();
    trk4_res["Data"][0] -> Write();
    
    BulkDelete(trk1_res); 
    BulkDelete(trk2_res); 
    BulkDelete(trk3_res); 
    BulkDelete(trk4_res); 

    delete trk1_res["Data"][0];
    delete trk2_res["Data"][0];
    delete trk3_res["Data"][0];
    delete trk4_res["Data"][0];
   
    F -> cd(); 
  }
  F -> Close(); 
}

void DataAlgorithm()
{
  TFile* F = new TFile("Results_Data.root", "RECREATE"); 
  F -> ReOpen("UPDATE"); 










}

float FLost_Track_1(std::vector<TH1F*> Hists)
{
  float L1 = Hists[0] -> Integral(); 
  float L2 = Hists[1] -> Integral(); 
  float L3 = Hists[2] -> Integral(); 
  float L4 = Hists[3] -> Integral(); 

  float FLost = (L2 + 2*L3 + 3*L4) / (2*L2 + 3*L3 + 4*L4); 
  return FLost;
}

float FLost_Track_2(std::vector<TH1F*> Hists)
{
  float L2 = Hists[1] -> Integral(); 
  float L3 = Hists[2] -> Integral(); 
  float L4 = Hists[3] -> Integral(); 

  float FLost = (L3 + 2*L4) / (3*L3 + 4*L4);  
  return FLost; 
}

void ProcessDataResults()
{






}





void ProcessMonteCarloResults()
{





} 




