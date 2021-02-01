#include<PostAnalysis/AlgorithmTest.h>

void TestAlgorithmMonteCarlo()
{
  std::map<TString, std::vector<TH1F*>> MC_Layers = MonteCarlo("Merged_MC.root"); 
  std::vector<TString> layers = {"IBL", "Blayer", "layer1", "layer2", "All"}; 
  
  TFile* F = new TFile("Results_MC.root", "RECREATE"); 
  F -> ReOpen("UPDATE"); 
  
  float m = 0.001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.03, 0.03, 0.03, 0.03};  
  Params["x_range"] = {0.05, 11.5}; 
  Params["iterations"] = {5}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.01, 0.01, 0.01, 0.01}; 
  Params["cache"] = {10000};  
  
  for (int i(0); i < layers.size(); i++)
  {
    F -> mkdir(layers[i]); 
    F -> cd(layers[i]); 
    std::cout << ":::::::::::::: Current -> " << layers[i] << std::endl;
    std::vector<TH1F*> trk1 = MC_Layers["trk1_"+layers[i]]; 
    std::vector<TH1F*> trk2 = MC_Layers["trk2_"+layers[i]]; 
    std::vector<TH1F*> trk3 = MC_Layers["trk3_"+layers[i]]; 
    std::vector<TH1F*> trk4 = MC_Layers["trk4_"+layers[i]]; 
    TH1F* Data_trk1 = SumHists(trk1, "Data_trk1_"+layers[i]); 
    TH1F* Data_trk2 = SumHists(trk2, "Data_trk2_"+layers[i]); 
    TH1F* Data_trk3 = SumHists(trk3, "Data_trk3_"+layers[i]); 
    TH1F* Data_trk4 = SumHists(trk4, "Data_trk4_"+layers[i]); 
    std::vector<TH1F*> Data = {Data_trk1, Data_trk2, Data_trk3, Data_trk4}; 

    bool kill = false; 
    for (TH1F* H : Data)
    {
      if (H -> Integral() == 0)
      {
        kill = true;
      }
    }
    if (kill == true)
    {
      for (TH1F* H : Data){delete H;}
      continue;
    }

    std::map<TString, std::vector<TH1F*>> trk1_res = MainAlgorithm(Data, Params, 0); 
    std::map<TString, std::vector<TH1F*>> trk2_res = MainAlgorithm(Data, Params, 1);    
    std::map<TString, std::vector<TH1F*>> trk3_res = MainAlgorithm(Data, Params, 2); 
    std::map<TString, std::vector<TH1F*>> trk4_res = MainAlgorithm(Data, Params, 3); 

    typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
    it p = trk1_res.end(); 
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

  std::map<TString, std::vector<TH1F*>> MC_Energy = MonteCarloLayerEnergy("Merged_MC.root"); 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 

  for (int i(0); i < JetEnergy.size(); i++)
  {
    F -> mkdir(JetEnergy[i]); 
    F -> cd(JetEnergy[i]); 
    std::cout << ":::::::::::::: Current -> " << JetEnergy[i] << std::endl;
    std::vector<TH1F*> trk1 = MC_Energy["Track_1_"+JetEnergy[i]]; 
    std::vector<TH1F*> trk2 = MC_Energy["Track_2_"+JetEnergy[i]]; 
    std::vector<TH1F*> trk3 = MC_Energy["Track_3_"+JetEnergy[i]]; 
    std::vector<TH1F*> trk4 = MC_Energy["Track_4_"+JetEnergy[i]]; 
    TH1F* Data_trk1 = SumHists(trk1, "Data_trk1_"+JetEnergy[i]); 
    TH1F* Data_trk2 = SumHists(trk2, "Data_trk2_"+JetEnergy[i]); 
    TH1F* Data_trk3 = SumHists(trk3, "Data_trk3_"+JetEnergy[i]); 
    TH1F* Data_trk4 = SumHists(trk4, "Data_trk4_"+JetEnergy[i]); 
    std::vector<TH1F*> Data = {Data_trk1, Data_trk2, Data_trk3, Data_trk4}; 

    bool kill = false; 
    for (TH1F* H : Data)
    {
      if (H -> Integral() == 0)
      {
        kill = true;
      }
    }
    if (kill == true)
    {
      for (TH1F* H : Data){delete H;}
      continue;
    }

    std::map<TString, std::vector<TH1F*>> trk1_res = MainAlgorithm(Data, Params, 0); 
    std::map<TString, std::vector<TH1F*>> trk2_res = MainAlgorithm(Data, Params, 1);    
    std::map<TString, std::vector<TH1F*>> trk3_res = MainAlgorithm(Data, Params, 2); 
    std::map<TString, std::vector<TH1F*>> trk4_res = MainAlgorithm(Data, Params, 3); 

    typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
    it p = trk1_res.end(); 
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
  
  float m = 0.001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.03, 0.03, 0.03, 0.03};  
  Params["x_range"] = {0.05, 11.5}; 
  Params["iterations"] = {2}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.01, 0.01, 0.01, 0.01}; 
  Params["cache"] = {10000};  
 
  std::map<TString, TH1F*> MC_Energy = Data("Merged_Data.root"); 
  std::vector<TString> Names = {"Track_1_", "Track_2_", "Track_3_", "Track_4_"};  
  std::vector<TString> Layers = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 

  std::vector<TString> Ext = {"W", "U"}; 
  typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
  for (TString w : Ext)
  {
    for (int i(0); i < JetEnergy.size(); i++)
    {
      std::vector<TH1F*> Data;
      
      for (int t(0); t < Layers.size(); t++)
      { 
        for (int x(0); x < Names.size(); x++)
        {
          TString name = Layers[t] + "/" + JetEnergy[i] + "/" + Names[x] + w; 
          TH1F* H = MC_Energy[name]; 
          if (t == 0)
          {
            TH1F* N = (TH1F*)H -> Clone(JetEnergy[i] + "_" + Names[x] + "_" + w);
            N -> SetTitle(JetEnergy[i] + "_" + Names[x] + "_" + w); 
            Data.push_back(N); 
          }
          else{Data[x] -> Add(H,1);}
        }
      }

      F -> mkdir(JetEnergy[i] + w); 
      F -> cd(JetEnergy[i] + w); 

      bool kill = false;
      std::vector<TH1F*> Data_T;
      for (int x(0); x < Names.size(); x++)
      {
        TH1F* H = Data[x]; 
        if (H -> GetEntries() < 1000){kill = true;}
        Data_T.push_back(H);
      }
      if (kill == true)
      {
        BulkDelete(Data_T); 
        continue;
      }

      std::cout << ":::::::::::::: Current -> " << JetEnergy[i] << std::endl; 
      for (int p(0); p < Data_T.size(); p++)
      {
        std::map<TString, std::vector<TH1F*>> trk_res = MainAlgorithm(Data_T, Params, p); 
        it x = trk_res.end(); 
        x--; 
        std::vector<TH1F*> trk_r = trk_res[x -> first]; 
        std::vector<TH1F*> D = trk_res["Data"]; 
        BulkWrite(trk_r); 
        D[0] -> Write();
        BulkDelete(trk_res);
      }
      BulkDelete(Data_T);  

      F -> cd(); 

    }

    for (int i(0); i < Layers.size(); i++)
    {
      F -> mkdir(Layers[i] + w); 
      F -> cd(Layers[i] + w); 

      std::vector<TH1F*> Data;  
      bool kill = false;
      for (int x(0); x < Names.size(); x++)
      {
        TH1F* H = MC_Energy[Layers[i] + "/" + Names[x] + w];
        if (H -> GetEntries() < 1000){kill = true;}
        Data.push_back(H); 
      }
      if (kill == true)
      {
        BulkDelete(Data);
        continue;
      }

      std::cout << ":::::::::::::: Current Layer -> " << Layers[i] << std::endl; 
      for (int p(0); p < Data.size(); p++)
      {
        std::map<TString, std::vector<TH1F*>> trk_res = MainAlgorithm(Data, Params, p); 
        it k = trk_res.end(); 
        k--; 
        std::vector<TH1F*> trk_r = k -> second; 
        std::vector<TH1F*> D = trk_res["Data"]; 
        BulkWrite(trk_r); 
        D[0] -> Write();
        BulkDelete(trk_r);
      }
      BulkDelete(Data);  

      F -> cd(); 
    }
  }
  F -> Close(); 
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
  auto Rename =[](std::vector<TH1F*> Hists)
  {
    for (int i(0); i < Hists.size(); i++)
    {
      TString name = "Track-"; name +=(i+1); 
      Hists[i] -> SetTitle(name); 
    }
  };

  TFile* F = new TFile("Results_Data.root", "READ"); 
  std::map<TString, std::vector<TH1F*>> Out = ReadEntries(F); 
  typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print("FLost.pdf[");
  for (it p = Out.begin(); p != Out.end(); p++)
  {
    TString name = p -> first; 
    std::vector<TH1F*> Hists_V = Out[p -> first]; 
    std::vector<TH1F*> trk1; 
    std::vector<TH1F*> trk2; 
    std::vector<TH1F*> trk3; 
    for (TH1F* H : Hists_V)
    {
      TString histname = H -> GetTitle(); 
      if (histname.Contains("ntrk_1")){trk1.push_back(H);}
      if (histname.Contains("ntrk_2")){trk2.push_back(H);}
      if (histname.Contains("ntrk_3")){trk3.push_back(H);}
    }
    //Rename(trk1); 
    //Rename(trk2);
    //Rename(trk3); 
    
    if (trk1.size() != 3){continue;}
    float FLost_1 = FLost_Track_1(trk1); 
    float FLost_2 = FLost_Track_2(trk2);

    TH1F* Data1 = SumHists(trk1, "DATA");  // <<===== replace afterwards with real data
    Data1 -> SetLineColor(kBlack); 
    PlotHists(Data1, trk1, trk2, "Track-1 Distribution FLost", FLost_1, FLost_2, can); 
    can -> Print("FLost.pdf");
    can -> Clear();
  }
  can -> Print("FLost.pdf]");
  delete can;
}





void ProcessMonteCarloResults()
{
   auto Rename =[](std::vector<TH1F*> Hists)
  {
    for (int i(0); i < Hists.size(); i++)
    {
      TString name = "Track-"; name +=(i+1); 
      Hists[i] -> SetTitle(name); 
    }
  };

  TFile* F = new TFile("Results_MC.root", "READ"); 
  std::map<TString, std::vector<TH1F*>> Out = ReadEntries(F); 

  typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
  
  TCanvas* can = new TCanvas(); 
  can -> Print("FLost_MonteCarlo.pdf["); 
  for (it p = Out.begin(); p != Out.end(); p++)
  {
    TString name = p -> first; 
    std::vector<TH1F*> Hists_V = Out[name]; 
   
    std::vector<TH1F*> trk1; 
    std::vector<TH1F*> trk2; 
    std::vector<TH1F*> trk3; 
    for (int i(0); i < Hists_V.size(); i++)
    {
      TString histname = Hists_V[i] -> GetTitle(); 
      if (histname.Contains("trk_1")){trk1.push_back(Hists_V[i]);} 
      if (histname.Contains("trk_2")){trk2.push_back(Hists_V[i]);} 
      if (histname.Contains("trk_3")){trk3.push_back(Hists_V[i]);} 
    }
    
    std::vector<TH1F*> trk1_MC; 
    std::vector<TH1F*> trk2_MC; 
    std::vector<TH1F*> trk3_MC; 
    std::map<TString, std::vector<TH1F*>> MC; 
    TString ext; 
    if (name.Contains("_GeV"))
    {
      MC = MonteCarloLayerEnergy("Merged_MC.root");
      ext = "Track_"; 
    }
    else
    {
      MC = MonteCarlo("Merged_MC.root");
      ext = "trk";
    }

    for (it i = MC.begin(); i != MC.end(); i++)
    {
      TString histname = i -> first;
      if (histname.Contains("dEdx_") || histname.Contains("_eW1")){continue;}
      if (histname.Contains(name) && histname.Contains(ext+"1")){trk1_MC = MC[histname];}
      if (histname.Contains(name) && histname.Contains(ext+"2")){trk2_MC = MC[histname];}
      if (histname.Contains(name) && histname.Contains(ext+"3")){trk3_MC = MC[histname];}
    }
    if (trk1.size() != trk1_MC.size()){continue;}

    float FLost_1 = FLost_Track_1(trk1); 
    float FLost_2 = FLost_Track_2(trk2);

    float FLost_1_MC = FLost_Track_1(trk1_MC); 
    float FLost_2_MC = FLost_Track_1(trk2_MC); 

    for (int i(0); i < trk1.size(); i++)
    {
      Statistics(trk1[i], trk1_MC[i], 0.1, 10); 
    }
    
    TH1F* Data_trk1 = SumHists(trk1_MC, "Data Track-1 Monte Carlo"); 
    Data_trk1 -> SetLineColor(kBlack); 
    PlotHists(Data_trk1, trk1, trk1_MC, "Track-1 Distribution FLost: " + name , FLost_1, FLost_1_MC, can); 
    can -> Print("FLost_MonteCarlo.pdf");

    std::cout << "======================================" << std::endl;
    for (int i(0); i < trk1.size(); i++)
    {
      Statistics(trk1[i], trk1_MC[i], 0.1, 10); 
    }

    TH1F* Data_trk2 = SumHists(trk2_MC, "Data Track-2 Monte Carlo"); 
    Data_trk2 -> SetLineColor(kBlack); 
    PlotHists(Data_trk2, trk2, trk2_MC, "Track-2 Distribution FLost: " + name , FLost_2, FLost_2_MC, can); 
    can -> Print("FLost_MonteCarlo.pdf");
  }
  can -> Print("FLost_MonteCarlo.pdf]");
  delete can;
} 




