#include<PostAnalysis/Fits.h>

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> OnlyNormal(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  
  std::vector<std::vector<TH1F*>> Conv; 

  for (int i(0); i < Data.size(); i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, Data.size(), N + Name + "_Normal"); 
    Conv.push_back(ntrk_Conv);
  }

  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  for (int i(0); i < Data.size(); i++)
  {
    TH1F* trk = Data[i]; 
    std::vector<TH1F*> ntrk_Conv = Conv[i]; 
    std::map<TString, std::vector<float>> Result = ScalingFit(trk, ntrk_Conv); 

    typedef std::map<TString, std::vector<float>>::iterator it; 
    std::vector<float> Error;  
    for (it p = Result.begin(); p != Result.end(); p++)
    {
      std::vector<float> R = p -> second; 
      Error.push_back(R[1]);  
    }
    Output[trk -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Error, ntrk_Conv); 
  }
  return Output; 
}

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShift(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  std::vector<std::vector<TH1F*>> Conv; 

  for (int i(0); i < Data.size(); i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, Data.size(), N + Name + "_NormalShift"); 
    Conv.push_back(ntrk_Conv);
  }

  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  for (int i(0); i < Data.size(); i++)
  {
    TH1F* trk = Data[i]; 
    std::vector<TH1F*> ntrk_Conv = Conv[i]; 
    std::map<TString, std::vector<float>> Result = ScalingShift(trk, ntrk_Conv); 

    typedef std::map<TString, std::vector<float>>::iterator it; 
    std::vector<float> Error;  
    for (it p = Result.begin(); p != Result.end(); p++)
    {
      std::vector<float> R = p -> second; 
      Error.push_back(R[1]);  
    }
    

    Output[trk -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Error, ntrk_Conv); 
  }
  return Output; 
}

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> ShiftNormalFFT(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  std::vector<std::vector<TH1F*>> Conv; 
  for (int i(0); i < Data.size(); i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, Data.size(), N + Name + "_ShiftNormalFFT"); 
    Conv.push_back(ntrk_Conv);
  }

  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  for (int i(0); i < Data.size(); i++)
  {
    TH1F* trk = Data[i]; 
    std::vector<TH1F*> ntrk_Conv = Conv[i]; 
    std::vector<std::pair<TH1F*, std::vector<float>>> Result = FitDeconvolutionPerformance(trk, ntrk_Conv, Params, Params["cache"][0], Params["cache"][0]); 
    
    std::vector<float> Error; 
    std::vector<TH1F*> Temp_H; 
    for (std::pair<TH1F*, std::vector<float>> H : Result)
    {
      Temp_H.push_back(H.first); 
      Error.push_back(H.second[5]); 
    }
    Flush(Temp_H, ntrk_Conv); 
    Output[trk -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Error, ntrk_Conv); 
  }
  return Output; 
}

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalWidthDeconvShiftFFT(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  std::vector<std::vector<TH1F*>> Conv; 
  for (int i(0); i < Data.size(); i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, Data.size(), N + Name + "_NormalWidthDeconvShiftFFT"); 
    Conv.push_back(ntrk_Conv);
  }

  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  for (int i(0); i < Data.size(); i++)
  {
    TH1F* trk = Data[i]; 
    std::vector<TH1F*> ntrk_Conv = Conv[i]; 
    std::vector<std::pair<TH1F*, std::vector<float>>> Result = LoopGenAll(ntrk_Conv, trk, Params); 
    
    std::vector<float> Error; 
    std::vector<TH1F*> Temp_H; 
    for (std::pair<TH1F*, std::vector<float>> H : Result)
    {
      Temp_H.push_back(H.first); 
      Error.push_back(H.second[5]); 
    }
    Flush(Temp_H, ntrk_Conv); 
    Output[trk -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Error, ntrk_Conv); 
  }
  return Output; 
}

void AnalysisPlots(std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Input, std::vector<std::vector<TH1F*>> Truth, TString Name)
{
  typedef std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>>::iterator it; 
  
  int t = 0; 
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print(Name + ".pdf[");  
  for (it i = Input.begin(); i != Input.end(); i++)
  {
    TString Name_trk = i -> first; 
    std::vector<float> Error_V = i -> second.first; 
    std::vector<TH1F*> ntru_p = i -> second.second; 
    std::vector<TH1F*> ntru_t = Truth[t]; 
       
    PlotHists(ntru_p, ntru_t, can); 
    can -> Print(Name + ".pdf"); 
    t++; 
  }
  can -> Print(Name + ".pdf]");  
  delete can; 
}

void WriteToFile(TFile* F, std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> FitResult, TString Name)
{
  
  F -> mkdir(Name); 
  typedef std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>>::iterator it;  

  for (it i = FitResult.begin(); i != FitResult.end(); i++)
  {
    TString trk_title = i -> first; 
    F -> mkdir(Name + "/" + trk_title); 
    F -> cd(Name + "/" + trk_title); 

    std::vector<TH1F*> Hists = i -> second.second;
    std::vector<float> Error = i -> second.first; 
    BulkWrite(Hists); 
    
    TH1F* H = new TH1F("Error", "Error", Error.size(), 0, Error.size()); 
    for (int p(0); p < H -> GetNbinsX(); p++){H -> SetBinContent(p+1, Error[p]); }
    H -> Write(); 

    F -> cd(); 
    delete H; 
  }
}




void AnalysisLoop()
{

  float m = 0.1;
  std::map<TString, std::vector<float>> Params_SmallWidth; 
  Params_SmallWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_SmallWidth["m_e"] = {m, m, m, m}; 
  Params_SmallWidth["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params_SmallWidth["s_e"] = {0.006, 0.006, 0.006, 0.006};  
  Params_SmallWidth["iterations"] = {1}; 
  Params_SmallWidth["LR_iterations"] = {50}; 
  Params_SmallWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_SmallWidth["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params_SmallWidth["cache"] = {10000}; 
  Params_SmallWidth["x_range"] = {0.1, 9.0}; 

  std::map<TString, std::vector<float>> Params_NormalWidth; 
  Params_NormalWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_NormalWidth["m_e"] = {m, m, m, m}; 
  Params_NormalWidth["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params_NormalWidth["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params_NormalWidth["iterations"] = {1}; 
  Params_NormalWidth["LR_iterations"] = {50}; 
  Params_NormalWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_NormalWidth["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params_NormalWidth["cache"] = {10000}; 
  Params_NormalWidth["x_range"] = {0.1, 9.0}; 

  std::map<TString, std::vector<float>> Params;
  
  std::map<TString, std::vector<TH1F*>> Hists = MC_Reader_All("Merged_MC.root");
  typedef std::map<TString, std::vector<TH1F*>>::iterator it;
 
  TFile* F = new TFile("Fits_Results_New.root", "RECREATE"); 
  for (it i = Hists.begin(); i != Hists.end(); i++)
  {
    TString Name = i -> first; 

    std::vector<TH1F*> Hist = i -> second; 
    if (Name.Contains("radius") || Name.Contains("Truth")){continue;}
    
    bool skip = false; 
    std::vector<TH1F*> Inside = Hists[Name + "_radius_Less"]; 
    std::vector<TH1F*> Inside_Truth = Hists[Name + "_radius_Less_Truth"]; 
    TH1F* Outside = Hists[Name + "_radius_Greater"][0]; 
    for (int z(0); z < 2; z++)
    {
      int E = Inside[z] -> GetEntries(); 
      if (E < 1000){skip = true;} 
    }
    if (Outside -> GetEntries() < 1000){skip = true;}
    if (skip == true){continue;}

    std::vector<TH1F*> Inside_Using; 
    std::vector<std::vector<TH1F*>> Inside_Using_Truth; 
    for (int z(0); z < Inside.size(); z++)
    {
      int E = Inside[z] -> GetEntries(); 
      if (E < 1000){continue;}

      Inside_Using.push_back(Inside[z]); 
      std::vector<TH1F*> Temp; 
      for (int y(z*4); y < (z+1)*4; y++){Temp.push_back(Inside_Truth.at(y));}
      Inside_Using_Truth.push_back(Temp); 
    }

    std::cout << Name << std::endl;

    // Doing the different fit techniques 
    // Normalization fit 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Normal = OnlyNormal(Inside_Using, Outside, Params, Name); 
    AnalysisPlots(Normal, Inside_Using_Truth, Name + "_Normal");  
    WriteToFile(F, Normal, Name + "_Normal"); 

    // Shifting with Normalization
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> ShiftNormal = NormalShift(Inside_Using, Outside, Params, Name);
    AnalysisPlots(ShiftNormal, Inside_Using_Truth, Name + "_ShiftNormal");  
    WriteToFile(F, ShiftNormal, Name + "_ShiftNormal"); 


    // Shifting Normalization with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShiftFFT = ShiftNormalFFT(Inside_Using, Outside, Params_SmallWidth, Name); 
    AnalysisPlots(NormalShiftFFT, Inside_Using_Truth, Name + "_NormalShiftFFT");  
    WriteToFile(F, NormalShiftFFT, Name + "_ShiftNormalFFT"); 

    // Shifting, Normalization, Width with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShiftWidthFFT = NormalWidthDeconvShiftFFT(Inside_Using, Outside, Params_NormalWidth, Name);
    AnalysisPlots(NormalShiftFFT, Inside_Using_Truth, Name + "_NormalShiftWidthFFT");  
    WriteToFile(F, NormalShiftWidthFFT, Name + "_NormalShiftWidthFFT"); 
   
  }
  F -> Close(); 
  delete F;
};

std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> ReadResults()
{
  TFile* F = new TFile("Fits_Results.root", "READ");
  
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

void Reading_Stats()
{
  
  std::map<TString, std::vector<TH1F*>> Truth = MC_Reader_All("Merged_MC.root"); 
  typedef std::map<TString, std::vector<TH1F*>>::iterator xi;  
  
  std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> F = ReadResults();  
  typedef std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>>::iterator it; 
 
  std::vector<TString> Layer = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"};  


  std::map<TString, std::map<TString, std::vector<float>>> Results; 

  for (TString Li : Layer)
  {
    for (TString Ji : JetEnergy)
    {
      TString base = Li + "_" + Ji; 
      TString truth_name = base + "_radius_Less_Truth"; 
      std::vector<TH1F*> Trks_TN = Truth[truth_name];
      
      std::vector<TH1F*> trk_1_truth; 
      std::vector<TH1F*> trk_2_truth; 
      std::vector<TH1F*> trk_3_truth; 
      std::vector<TH1F*> trk_4_truth; 
      
      int sdo = 4; 
      for (int tru(0); tru < Trks_TN.size(); tru++)
      {
        if (tru < sdo){trk_1_truth.push_back(Trks_TN[tru]); }
        if (tru > sdo-1 && tru < 2*sdo){trk_2_truth.push_back(Trks_TN[tru]); }
        if (tru > 2*sdo-1 && tru < 3*sdo){trk_3_truth.push_back(Trks_TN[tru]); }
        if (tru > 3*sdo-1 && tru < 4*sdo){trk_4_truth.push_back(Trks_TN[tru]); }
      }
     
      for (it i = F.begin(); i != F.end(); i++)
      {
        TString Energy_Layer = i -> first; 
        if (!Energy_Layer.Contains(base)){continue;}
        
        std::vector<std::pair<TString, std::vector<TH1F*>>> map = i -> second; 
        std::map<TString, std::vector<TH1F*>> Stats_Normal;  
        std::map<TString, std::vector<TH1F*>> Stats_ShiftNormal;  
        std::map<TString, std::vector<TH1F*>> Stats_ShiftNormalFFT;  
        std::map<TString, std::vector<TH1F*>> Stats_NormalShiftWidthFFT;  
    
        int opt = -1; 
        for (int p(0); p < map.size(); p++)
        {
          auto Fill =[] (std::map<TString, std::vector<TH1F*>> Map, TString key, TString Incoming, std::vector<TH1F*> Pred, std::vector<TH1F*> Truth)
          {
            if (Incoming.Contains(key))
            {
              Map[key] = Pred; 
              Map[key + "_T"] = Truth; 
            }
            return Map; 
          };
          TString sub = map[p].first;
          std::vector<TH1F*> pred = map[p].second;

          if (Energy_Layer.Contains("NormalShiftWidthFFT"))
          {
            Stats_NormalShiftWidthFFT = Fill(Stats_NormalShiftWidthFFT, "ntrk_1", sub, pred, trk_1_truth); 
            Stats_NormalShiftWidthFFT = Fill(Stats_NormalShiftWidthFFT, "ntrk_2", sub, pred, trk_2_truth); 
            Stats_NormalShiftWidthFFT = Fill(Stats_NormalShiftWidthFFT, "ntrk_3", sub, pred, trk_3_truth); 
            Stats_NormalShiftWidthFFT = Fill(Stats_NormalShiftWidthFFT, "ntrk_4", sub, pred, trk_4_truth); 
            opt = 1;  
          }
          else if (Energy_Layer.Contains("ShiftNormalFFT"))
          {
            Stats_ShiftNormalFFT = Fill(Stats_ShiftNormalFFT, "ntrk_1", sub, pred, trk_1_truth); 
            Stats_ShiftNormalFFT = Fill(Stats_ShiftNormalFFT, "ntrk_2", sub, pred, trk_2_truth); 
            Stats_ShiftNormalFFT = Fill(Stats_ShiftNormalFFT, "ntrk_3", sub, pred, trk_3_truth); 
            Stats_ShiftNormalFFT = Fill(Stats_ShiftNormalFFT, "ntrk_4", sub, pred, trk_4_truth); 
            opt = 2; 
          }         
          else if (Energy_Layer.Contains("ShiftNormal"))
          {
            Stats_ShiftNormal = Fill(Stats_ShiftNormal, "ntrk_1", sub, pred, trk_1_truth); 
            Stats_ShiftNormal = Fill(Stats_ShiftNormal, "ntrk_2", sub, pred, trk_2_truth); 
            Stats_ShiftNormal = Fill(Stats_ShiftNormal, "ntrk_3", sub, pred, trk_3_truth); 
            Stats_ShiftNormal = Fill(Stats_ShiftNormal, "ntrk_4", sub, pred, trk_4_truth); 
            opt = 3; 
          }        
          else if (Energy_Layer.Contains("Normal"))
          {
            Stats_Normal = Fill(Stats_Normal, "ntrk_1", sub, pred, trk_1_truth); 
            Stats_Normal = Fill(Stats_Normal, "ntrk_2", sub, pred, trk_2_truth); 
            Stats_Normal = Fill(Stats_Normal, "ntrk_3", sub, pred, trk_3_truth); 
            Stats_Normal = Fill(Stats_Normal, "ntrk_4", sub, pred, trk_4_truth); 
            opt = 4; 
          }        
        }
       
        if (opt == 1)
        {
          std::map<TString, std::vector<float>> NormalShiftWidthFFT_R = AnalysisOutput(Stats_NormalShiftWidthFFT, base + "NormalShiftWidthFFT"); 
          Results[base + "NormalShiftWidthFFT"] = NormalShiftWidthFFT_R; 
        }
        if (opt == 2)
        {
          std::map<TString, std::vector<float>> ShiftNormalFFT_R = AnalysisOutput(Stats_ShiftNormalFFT, base + "ShiftNormalFFT"); 
          Results[base + "ShiftNormalFFT"] = ShiftNormalFFT_R; 
        }
        if (opt == 3)
        {
          std::map<TString, std::vector<float>> ShiftNormal_R = AnalysisOutput(Stats_ShiftNormal, base + "ShiftNormal"); 
          Results[base + "ShiftNormal"] = ShiftNormal_R; 
        }
        if (opt == 4)
        {
          std::map<TString, std::vector<float>> Normal_R = AnalysisOutput(Stats_Normal, base + "Normal"); 
          Results[base + "Normal"] = Normal_R; 
        }
      }
    }
  }
  
  typedef std::map<TString, std::map<TString, std::vector<float>>>::iterator zi;
  for (zi p = Results.begin(); p != Results.end(); p++)
  {
    auto Print =[] (std::map<TString, std::vector<float>> R)
    {
      typedef std::map<TString, std::vector<float>>::iterator x; 
      for (x i = R.begin(); i != R.end(); i++)
      {
        TString key = i -> first; 
        if (key.Contains("FLost")){continue;}
        std::cout << i -> first << "| ";
        int trk = 0;  
        for (float p : i -> second)
        {
          trk++; 
          std::cout << "trk_" << trk << ": "; 
          std::cout<< p << "  |  ";
        }
        std::cout << std::endl;
      }
      std::cout << "FLost2 Predicted: " << R["FLost2_P"][0] << " +- " << R["FLost2_P"][1] << std::endl;
      std::cout << "FLost2 Truth: " << R["FLost2_T"][0] <<  std::endl;


    };
      
    TString name = p -> first; 

    std::cout << "____ " << name << " ____" << std::endl;
    Print(p -> second); 
    std::cout << "" << std::endl;




     






  }











};

void Evaluation()
{
  //AnalysisLoop();  
  Reading_Stats(); 
}


