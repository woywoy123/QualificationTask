#include<PostAnalysis/AlgorithmTest.h>

void TestAlgorithmMonteCarlo()
{
  int nsdo = 4;
  float m = 0.0001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.04, 0.04, 0.04, 0.04};  
  Params["iterations"] = {2}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.02, 0.02, 0.02, 0.02}; 
  Params["cache"] = {10000}; 
  Params["x_range"] = {0., 10.};  

  TString Dir = "Merged_MC.root"; 
  std::map<TString, std::vector<TH1F*>> Output = MC_Reader(Dir);
 
  TFile* F = new TFile("Results_MC.root", "RECREATE"); 
  F -> ReOpen("UPDATE"); 
 
  typedef std::map<TString, std::vector<TH1F*>>::iterator M;
  for (M x = Output.begin(); x != Output.end(); x++)
  {
    TString name = x -> first; 
    std::vector<TH1F*> Trks = x -> second; 
   
    if (name.Contains("_Less") || name.Contains("Greater") || name.Contains("Data")){continue;}
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
    TH1F* Trk1 = Output[name + "_Data"][0]; 
    TH1F* Trk2 = Output[name + "_Data"][1]; 
    TH1F* Trk3 = Output[name + "_Data"][2]; 
    TH1F* Trk4 = Output[name + "_Data"][3]; 
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
    
    std::vector<TH1F*> trk1_r = trk1_res["Result"]; 
    std::vector<TH1F*> trk2_r = trk2_res["Result"]; 
    std::vector<TH1F*> trk3_r = trk3_res["Result"]; 
    std::vector<TH1F*> trk4_r = trk4_res["Result"];


    std::cout << trk1_r[0] -> GetTitle() << std::endl;
      
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

float FLost_Track_2(std::vector<std::vector<TH1F*>> Hists)
{
  float L1_2 = Hists[0][1] -> Integral(); // Track-1, Truth-2 
  float L2_2 = Hists[1][1] -> Integral(); // Track-2, Truth-2  
  float FLost = L1_2 / (2 * L1_2 + 2 * L2_2); 
  return FLost; 
}

void ProcessDataResults()
{






}

void ProcessMonteCarloResults()
{
  auto ChangeNames =[] (TString Track, std::vector<TH1F*> Hists, TString Layer)
  {
    for (int i(0); i < Hists.size(); i++)
    {
      TString name = Track + "_ntru_"; name += (i+1); name += ("_"); name+=(Layer); 
      Hists[i] -> SetTitle(name);
    }
  }; 


  std::map<TString, std::vector<TH1F*>> Map = Result_Reader("Results_MC.root");
  std::map<TString, std::vector<TH1F*>> MC = MC_Reader("Merged_MC.root"); 
  typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
 
  TCanvas* can = new TCanvas(); 
  can -> Print("Results_MC.pdf["); 
  can -> SetLogy();
  for (it p = Map.begin(); p != Map.end(); p++)
  {
    TString n = p -> first; 
    std::vector<TH1F*> Truth = MC[p -> first]; 
    if (Truth.size() == 0 || n.Contains("_Data")){continue;} 
    
    std::cout << " :::: " <<  p -> first << std::endl;
    std::vector<TH1F*> Track1_T = {Truth[0], Truth[1], Truth[2], Truth[3]}; 
    std::vector<TH1F*> Track2_T = {Truth[0+4], Truth[1+4], Truth[2+4], Truth[3+4]}; 
    std::vector<TH1F*> Track3_T = {Truth[0+8], Truth[1+8], Truth[2+8], Truth[3+8]}; 
    std::vector<TH1F*> Track4_T = {Truth[0+12], Truth[1+12], Truth[2+12], Truth[3+12]}; 
    std::vector<std::vector<TH1F*>> Trk_Truth = {Track1_T, Track2_T, Track3_T, Track4_T}; 
    float FLost_T = FLost_Track_2(Trk_Truth);

    std::vector<TH1F*> Pred = Map[p -> first]; 
    std::vector<TH1F*> Track1_P = {Pred[0], Pred[1], Pred[2], Pred[3]}; 
    std::vector<TH1F*> Track2_P = {Pred[0+4], Pred[1+4], Pred[2+4], Pred[3+4]}; 
    std::vector<TH1F*> Track3_P = {Pred[0+8], Pred[1+8], Pred[2+8], Pred[3+8]}; 
    std::vector<TH1F*> Track4_P = {Pred[0+12], Pred[1+12], Pred[2+12], Pred[3+12]}; 
    ChangeNames("P_trk_1", Track1_P, p -> first); 
    ChangeNames("P_trk_2", Track2_P, p -> first);    
    ChangeNames("P_trk_3", Track3_P, p -> first); 
    ChangeNames("P_trk_4", Track4_P, p -> first); 
    
    std::vector<std::vector<TH1F*>> Trk_Pred = {Track1_P, Track2_P, Track3_P, Track4_P}; 
    float FLost_P = FLost_Track_2(Trk_Pred);   
    
    TH1F* Data_trk1 = Map[p -> first + "_Data"][0]; 
    TH1F* Data_trk2 = Map[p -> first + "_Data"][1]; 
    TH1F* Data_trk3 = Map[p -> first + "_Data"][2]; 
    TH1F* Data_trk4 = Map[p -> first + "_Data"][3]; 
    Data_trk1 -> SetLineColor(kBlack); 
    Data_trk2 -> SetLineColor(kBlack); 
    Data_trk3 -> SetLineColor(kBlack); 
    Data_trk4 -> SetLineColor(kBlack); 
    
    PlotHists(Data_trk1, Track1_P, Track1_T, can); 
    can -> Print("Results_MC.pdf"); 

    PlotHists(Data_trk2, Track2_P, Track2_T, p -> first + " Track2", FLost_P, FLost_T, can); 
    can -> Print("Results_MC.pdf"); 
   
    PlotHists(Data_trk3, Track3_P, Track3_T, can); 
    can -> Print("Results_MC.pdf"); 
   
    PlotHists(Data_trk4, Track4_P, Track4_T, can); 
    can -> Print("Results_MC.pdf"); 
    
    std::cout << "__________________" << std::endl;
  
  }
  can -> Print("Results_MC.pdf]"); 
} 

