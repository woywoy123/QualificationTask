#include<PostAnalysis/Fits.h>
#include<PostAnalysis/Evaluate_Fits.h>

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
    std::map<TString, std::vector<float>> Result = ScalingFit(trk, Conv[i]); 

    typedef std::map<TString, std::vector<float>>::iterator it; 
    std::vector<float> Error;  
    for (it p = Result.begin(); p != Result.end(); p++)
    {
      std::vector<float> R = p -> second; 
      Error.push_back(R[1]);  
    }
    Output[trk -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Error, Conv[i]); 
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
    int z(0);  
    for (std::pair<TH1F*, std::vector<float>> H : Result)
    {
      Temp_H.push_back(H.first); 
      Error.push_back(H.second[5]); 
      ntrk_Conv[z] -> Reset(); 
      ntrk_Conv[z] -> Add(H.first, 1); 
      z++;
    }

    BulkDelete(Temp_H);  
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

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Experimental_FFT(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  auto Subtract =[&] (std::vector<TH1F*> ntrk_Conv, TH1F* Data, int trk)
  {
    for (int i(0); i < ntrk_Conv.size(); i++){ if (i != trk){ Data -> Add(ntrk_Conv[i], -1); } }
    Average(Data); 
    Flush({Data}, {ntrk_Conv[trk]});  
  };

  auto Clone =[&] (TH1F* H, TH1F* O)
  {
    float HL = H -> Integral(); 
    float OL = O -> Integral(); 
  
    H -> Scale(1. / HL); 
    O -> Scale(1. / OL); 
    for (int x(0); x < H -> GetNbinsX(); x++)
    {
      float e = H -> GetBinContent(x+1); 
      float f = O -> GetBinContent(x+1); 
      O -> SetBinContent(x+1, e);
    }
    float T = O -> Integral(); 
    O -> Scale(1. / T); 
    H -> Scale(HL); 
    O -> Scale(OL); 
  };


  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  std::vector<std::vector<TH1F*>> Conv; 
  for (int i(0); i < Data.size(); i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, Data.size(), N + Name + "_ExperimentalFFT"); 
    Conv.push_back(ntrk_Conv);
  }
  
  for (int x(0); x < Params["iterations"][0]; x++)
  {
    for (int t(0); t < Data.size(); t++)
    {
      std::vector<TH1F*> conv_trks = ConvolveNTimes(Conv[0][0], Data.size(), "_Con"); 
      
      TString name_D = Data[t] -> GetTitle(); 
      TH1F* trk_Data = (TH1F*)Data[t] -> Clone(name_D + "_C"); 
      
      std::vector<std::pair<TH1F*, std::vector<float>>> F_C = LoopGenAll(Conv[t], trk_Data, Params); 
      BulkDelete(conv_trks);  

      std::vector<TH1F*> ntrk_Conv; 
      for (int p(0); p < F_C.size(); p++){ntrk_Conv.push_back(F_C[p].first);}
      
      Subtract(ntrk_Conv, trk_Data, t); 
      Flush(ntrk_Conv, Conv[t]); 
      
      if (x != Params["iterations"][0] - 1){continue;}
      std::vector<float> Error; 
      for (std::pair<TH1F*, std::vector<float>> H : F_C)
      {
        Error.push_back(H.second[5]); 
      }
      Output[Data[t] -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Error, Conv[t]); 
    }
    
    if (x <= Params["iterations"][0] -2)
    {
      for (int t(0); t < Conv.size(); t++)
      {
        TH1F* pure = Conv[t][t]; 
        for (int z(0); z < Conv.size(); z++){Clone(pure, Conv[z][t]);}
      }
    }
  }

  return Output; 
}

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Experimental_Shift(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  auto Subtract =[&] (std::vector<TH1F*> ntrk_Conv, TH1F* Data, int trk)
  {
    for (int i(0); i < ntrk_Conv.size(); i++){ if (i != trk){ Data -> Add(ntrk_Conv[i], -1); } }
    Average(Data); 
    Flush({Data}, {ntrk_Conv[trk]});  
  };

  auto Clone =[&] (TH1F* H, TH1F* O)
  {
    float HL = H -> Integral(); 
    float OL = O -> Integral(); 
  
    H -> Scale(1. / HL); 
    O -> Scale(1. / OL); 
    for (int x(0); x < H -> GetNbinsX(); x++)
    {
      float e = H -> GetBinContent(x+1); 
      float f = O -> GetBinContent(x+1); 
      O -> SetBinContent(x+1, e);
    }
    float T = O -> Integral(); 
    O -> Scale(1. / T); 
    H -> Scale(HL); 
    O -> Scale(OL); 
  };

  std::vector<std::vector<TH1F*>> Conv; 
  for (int i(0); i < Data.size(); i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, Data.size(), N + Name + "_ExperimentalShift"); 
    Conv.push_back(ntrk_Conv);
  }
 
  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  typedef std::map<TString, std::vector<float>>::iterator xi;   
  for (int i(0); i < Params["iterations"][0]; i++)
  {

    // Clone the data 
    std::vector<std::vector<float>> Errors; 
    for (int x(0); x < Data.size(); x++)
    {
      TString n = Data[x] -> GetTitle(); n += ("_"); 
      TH1F* D = (TH1F*)Data[x] -> Clone(n); 
      
      // Do the fit 
      std::map<TString, std::vector<float>> Fit_res = ScalingShift(D, Conv[x]);

      // Subtract the fit from the data
      Subtract(Conv[x], D, x); 
      
      // Collect error 
      std::vector<float> Temp; 
      for (xi h = Fit_res.begin(); h != Fit_res.end(); h++)
      {
        std::vector<float> e = h -> second; 
        Temp.push_back(e[1]);
      }
      Errors.push_back(Temp); 
    }
      
    
    for (int x(0); x < Conv.size(); x++)
    {
      TH1F* Template = Conv[x][x]; 

      for (int p(0); p < Conv.size(); p++)
      {
        if (x == p) { continue; }
        Clone(Template, Conv[p][x]);  
      }
    }
   
    if (i < Params["iterations"][0] -1){continue;}
    for (int x(0); x < Data.size(); x++)
    {
      Output[Data[x] -> GetTitle()] = std::pair<std::vector<float>, std::vector<TH1F*>>(Errors[x], Conv[x]); 
    }
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
    PlotHists(ntru_t, ntru_p, can); 
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

void AnalysisLoop(int sdo, TString SaveName)
{

  float m = 0.25;
  std::map<TString, std::vector<float>> Params_SmallWidth; 
  Params_SmallWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_SmallWidth["m_e"] = {m, m, m, m}; 
  Params_SmallWidth["s_s"] = {-0.0005, -0.0005, -0.0005, -0.0005};
  Params_SmallWidth["s_e"] = {0.001, 0.006, 0.006, 0.006};  
  Params_SmallWidth["iterations"] = {1}; 
  Params_SmallWidth["LR_iterations"] = {50}; 
  Params_SmallWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_SmallWidth["G_Stdev"] = {0.001, 0.001, 0.001, 0.001}; 
  Params_SmallWidth["cache"] = {10000}; 
  Params_SmallWidth["x_range"] = {0.01, 9.0}; 

  std::map<TString, std::vector<float>> Params_NormalWidth; 
  Params_NormalWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_NormalWidth["m_e"] = {m, m, m, m}; 
  Params_NormalWidth["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params_NormalWidth["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params_NormalWidth["iterations"] = {3}; 
  Params_NormalWidth["LR_iterations"] = {50}; 
  Params_NormalWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_NormalWidth["G_Stdev"] = {0.025, 0.025, 0.025, 0.025}; 
  Params_NormalWidth["cache"] = {10000}; 
  Params_NormalWidth["x_range"] = {0., 7.0}; 

  
  std::map<TString, std::vector<float>> Params;
  
  std::map<TString, std::vector<TH1F*>> Hists = MC_Reader_All("Merged_MC.root");
  typedef std::map<TString, std::vector<TH1F*>>::iterator it;
 
  TFile* F = new TFile(SaveName, "RECREATE"); 
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

     
      std::vector<TH1F*> Temp; 
      for (int y(z*sdo); y < (z+1)*sdo; y++){Temp.push_back(Inside_Truth.at(y));}
      Inside_Using_Truth.push_back(Temp); 

      if (sdo != 1){Inside_Using.push_back(Inside[z]); }
      else 
      {
        Inside_Using.push_back(Inside[z]);
        break;
      }
    }

    typedef std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>>::iterator ix; 

    // ====== Experimental Stuff ======== //

    // Shifting, Normalization, Width with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> ExperimentalShift = Experimental_FFT(Inside_Using, Outside, Params_NormalWidth, Name);
    AnalysisPlots(ExperimentalShift, Inside_Using_Truth, Name + "_AlgorithmExperimental");  
    WriteToFile(F, ExperimentalShift, Name + "_ExperimentalAlgorithm"); 
 
    for (ix x = ExperimentalShift.begin(); x != ExperimentalShift.end(); x++)
    {
      std::pair<std::vector<float>, std::vector<TH1F*>> p = x -> second; 
      BulkDelete(p.second); 
    }   
    // ====== End Experimental Stuff ======== //

    // Doing the different fit techniques 
    // Normalization fit 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Normal = OnlyNormal(Inside_Using, Outside, Params, Name); 
    AnalysisPlots(Normal, Inside_Using_Truth, Name + "_Normal");  
    WriteToFile(F, Normal, Name + "_Normal"); 

    for (ix x = Normal.begin(); x != Normal.end(); x++)
    {
      std::pair<std::vector<float>, std::vector<TH1F*>> p = x -> second; 
      BulkDelete(p.second); 
    }

    // Shifting with Normalization
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> ShiftNormal = NormalShift(Inside_Using, Outside, Params, Name);
    AnalysisPlots(ShiftNormal, Inside_Using_Truth, Name + "_ShiftNormal");  
    WriteToFile(F, ShiftNormal, Name + "_ShiftNormal"); 
    for (ix x = ShiftNormal.begin(); x != ShiftNormal.end(); x++)
    {
      std::pair<std::vector<float>, std::vector<TH1F*>> p = x -> second; 
      BulkDelete(p.second); 
    }
    
    // Shifting Normalization with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShiftFFT = ShiftNormalFFT(Inside_Using, Outside, Params_SmallWidth, Name); 
    AnalysisPlots(NormalShiftFFT, Inside_Using_Truth, Name + "_NormalShiftFFT");  
    WriteToFile(F, NormalShiftFFT, Name + "_ShiftNormalFFT"); 

    for (ix x = NormalShiftFFT.begin(); x != NormalShiftFFT.end(); x++)
    {
      std::pair<std::vector<float>, std::vector<TH1F*>> p = x -> second; 
      BulkDelete(p.second); 
    }
    
    // Shifting, Normalization, Width with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShiftWidthFFT = NormalWidthDeconvShiftFFT(Inside_Using, Outside, Params_NormalWidth, Name);
    AnalysisPlots(NormalShiftWidthFFT, Inside_Using_Truth, Name + "_NormalShiftWidthFFT");  
    WriteToFile(F, NormalShiftWidthFFT, Name + "_NormalShiftWidthFFT"); 
 
    for (ix x = NormalShiftWidthFFT.begin(); x != NormalShiftWidthFFT.end(); x++)
    {
      std::pair<std::vector<float>, std::vector<TH1F*>> p = x -> second; 
      BulkDelete(p.second); 
    }   
  
  }
  F -> Close(); 
  delete F;
}


void SimpleNTrackPlot()
{ 
  std::map<TString, std::vector<TH1F*>> Hists = MC_Reader_All("Merged_MC.root"); 
  typedef std::map<TString, std::vector<TH1F*>>::iterator it;
  
  std::vector<TH1F*> Trk1; 
 
  std::vector<TH1F*> Trk1_Tru1; 
  std::vector<TH1F*> Trk2_Tru2;
  std::vector<TH1F*> Trk3_Tru3; 
  std::vector<TH1F*> Trk4_Tru4; 

  for (it i = Hists.begin(); i != Hists.end(); i++)
  {
    TString Name = i -> first; 

    if (Name.Contains("radius") || Name.Contains("Truth")){continue;}
    
    std::vector<TH1F*> Inside = Hists[Name + "_radius_Less"]; 
    std::vector<TH1F*> Inside_Truth = Hists[Name + "_radius_Less_Truth"]; 
    TH1F* trk1_input = Hists[Name + "_radius_Greater"][0]; 
    Trk1.push_back(trk1_input); 
      
    Trk1_Tru1.push_back(Inside_Truth[0]); 
    Trk2_Tru2.push_back(Inside_Truth[5]); 
    Trk3_Tru3.push_back(Inside_Truth[10]); 
    Trk4_Tru4.push_back(Inside_Truth[15]); 
  }
  
  TH1F* trk1 = SumHists(Trk1, "Trk1_Inside"); 
  std::vector<TH1F*> ntrk_conv = ConvolveNTimes(trk1, 4, "TRK");  

  TH1F* trk1_tru1 = SumHists(Trk1_Tru1, "dEdx_trk1_tru1"); 
  TH1F* trk2_tru2 = SumHists(Trk2_Tru2, "dEdx_trk2_tru2"); 
  TH1F* trk3_tru3 = SumHists(Trk3_Tru3, "dEdx_trk3_tru3"); 
  TH1F* trk4_tru4 = SumHists(Trk4_Tru4, "dEdx_trk4_tru4"); 
  std::vector<TH1F*> trk_tru = {trk1_tru1, trk2_tru2, trk3_tru3, trk4_tru4}; 


  float m = 0.5;
  std::map<TString, std::vector<float>> Params_SmallWidth; 
  Params_SmallWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_SmallWidth["m_e"] = {m, m, m, m}; 
  Params_SmallWidth["s_s"] = {-0.0005, -0.0005, -0.0005, -0.0005}; // the minus sign disables the width 
  Params_SmallWidth["s_e"] = {0.001, 0.006, 0.006, 0.006};  
  Params_SmallWidth["LR_iterations"] = {50}; 
  Params_SmallWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_SmallWidth["G_Stdev"] = {0.001, 0.001, 0.001, 0.001}; 
  Params_SmallWidth["cache"] = {10000}; 
  Params_SmallWidth["x_range"] = {0., 7.0}; 

  std::map<TString, std::vector<float>> Params_NormalWidth; 
  Params_NormalWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_NormalWidth["m_e"] = {m, m, m, m}; 
  Params_NormalWidth["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params_NormalWidth["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params_NormalWidth["LR_iterations"] = {50}; 
  Params_NormalWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_NormalWidth["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params_NormalWidth["cache"] = {10000}; 
  Params_NormalWidth["x_range"] = {0., 6.0}; 

  std::map<TString, std::vector<float>> Params;
 
  // Normalization only tested:
  TString Name = "All_Layer_All_Energy"; 
  std::vector<TH1F*> ntrk_Conv_Normal = ConvolveNTimes(trk1, trk_tru.size(), Name + "_Normal"); 
  for (int g(0); g < trk_tru.size(); g++)
  {
    std::map<TString, std::vector<float>> NormalFit_R = ScalingFit(trk_tru[g], {ntrk_Conv_Normal[g]}); 
  }
  
  // Normalization + Shifting d - x method test: 
  std::vector<TH1F*> ntrk_Conv_ShiftNormal = ConvolveNTimes(trk1, trk_tru.size(), Name + "_ShiftNormal"); 
  for (int g(0); g < trk_tru.size(); g++)
  {
    std::map<TString, std::vector<float>> ShiftNormalFit_R = ScalingShift(trk_tru[g], {ntrk_Conv_ShiftNormal[g]}); 
  }
 
  // Normalization + Shifting using FFT test: 
  std::vector<TH1F*> ntrk_Conv_ShiftNormalFFT = ConvolveNTimes(trk1, trk_tru.size(), Name + "_ShiftNormalFFT"); 
  for (int g(0); g < trk_tru.size(); g++)
  {
    std::vector<std::pair<TH1F*, std::vector<float>>> ShiftNormalFFT_R = FitDeconvolutionPerformance(trk_tru[g], {ntrk_Conv_ShiftNormalFFT[g]}, Params_SmallWidth, Params_SmallWidth["cache"][0], Params_SmallWidth["cache"][0]); 
      
    ntrk_Conv_ShiftNormalFFT[g] -> Reset();
    ntrk_Conv_ShiftNormalFFT[g] -> Add(ShiftNormalFFT_R[0].first, 1); 
    delete ShiftNormalFFT_R[0].first; 
  }
  
  // Normalization + Shifting + Width using FFT test: 
  std::vector<TH1F*> ntrk_Conv_ShiftNormalWidthFFT = ConvolveNTimes(trk1, trk_tru.size(), Name + "_ShiftNormalWidthFFT"); 
  for (int g(0); g < trk_tru.size(); g++)
  {
    std::vector<TH1F*> ShiftNormalWidthFFT_R = LoopGen({ntrk_Conv_ShiftNormalWidthFFT[g]}, trk_tru[g], Params_NormalWidth); 
      
    ntrk_Conv_ShiftNormalWidthFFT[g] -> Reset();
    ntrk_Conv_ShiftNormalWidthFFT[g] -> Add(ShiftNormalWidthFFT_R[0], 1); 
    delete ShiftNormalWidthFFT_R[0]; 
  }

  std::vector<std::vector<float>> Inte_Error; 
  for (int i(0); i < trk_tru.size(); i++)
  {
    TString name = "dEdx_trk"; name += (i+1); name += ("_tru_"); name += (i+1); name += (".pdf"); 
    TCanvas* can = new TCanvas();
    can -> SetLogy(); 
    can -> Print(name + "["); 

    std::vector<float> Temp; 
    TH1F* trk_tru_T = trk_tru[i]; 
    
    TH1F* Normal_trk_tru = ntrk_Conv_Normal[i];
    RatioPlot(Normal_trk_tru, trk_tru_T, can); 
    can -> Print(name);  

    TH1F* NormalShift_trk_tru = ntrk_Conv_ShiftNormal[i]; 
    RatioPlot(NormalShift_trk_tru, trk_tru_T, can); 
    can -> Print(name);  

    TH1F* NormalShiftFFT_trk_tru = ntrk_Conv_ShiftNormalFFT[i]; 
    RatioPlot(NormalShiftFFT_trk_tru, trk_tru_T, can); 
    can -> Print(name);  

    TH1F* ShiftNormalWidthFFT_trk_tru = ntrk_Conv_ShiftNormalWidthFFT[i]; 
    RatioPlot(ShiftNormalWidthFFT_trk_tru, trk_tru_T, can); 
    can -> Print(name);  

    Temp.push_back(ErrorByIntegral(trk_tru_T, Normal_trk_tru)*100); 
    Temp.push_back(ErrorByIntegral(trk_tru_T, NormalShift_trk_tru)*100); 
    Temp.push_back(ErrorByIntegral(trk_tru_T, NormalShiftFFT_trk_tru)*100); 
    Temp.push_back(ErrorByIntegral(trk_tru_T, ShiftNormalWidthFFT_trk_tru)*100); 
    Inte_Error.push_back(Temp); 

    can -> Print(name + "]"); 
    delete can; 
  }
  
  std::vector<TString> Algos = {"Normal", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT"};  
  
  for (int i(0); i < Algos.size(); i++)
  {
    std::vector<float> f = Inte_Error[i]; 
  
    TString out =""; 
    out += (Algos[i]); 
    out +=(" | "); 
    for (int x(0); x < f.size(); x++)
    {
      out += (PrecisionString(f[x], 5)); 
      out += (" | "); 
    }
    std::cout << out << std::endl;
  }


  TH1F* Clear = (TH1F*)trk1 -> Clone("plotting"); 
  Normalize(Clear); 
  Clear -> SetTitle("The 1-Track dE/dx Measurement being Convolved n-Times"); 
  Clear -> SetLineColor(kRed); 

  ntrk_conv[0] -> SetTitle("Track1"); 
  ntrk_conv[1] -> SetTitle("Track2"); 
  ntrk_conv[2] -> SetTitle("Track3"); 
  ntrk_conv[3] -> SetTitle("Track4"); 

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  PlotHists(Clear, ntrk_conv, can); 
  can -> Print("ntrk_conv.pdf"); 

  

}


void Evaluation()
{
  //AnalysisLoop(1, "1_trk_1_tru.root"); 
  //AnalysisLoop(4, "All_Fits.root");
  //Evaluate_Fits("1_trk_1_tru.root");  
  //Evaluate_Fits("All_Fits.root");  
  //EvaluateErrorImpact("1_trk_1_tru.root"); 
  
  //EvaluatePureNTracksIndividually(); 
  //Evaluate_nTrackFits("Individual_Ntracks.root"); 
  SimpleNTrackPlot(); 
}

void EvaluatePureNTracksIndividually()
{
  float m = 0.5;
  std::map<TString, std::vector<float>> Params_SmallWidth; 
  Params_SmallWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_SmallWidth["m_e"] = {m, m, m, m}; 
  Params_SmallWidth["s_s"] = {-0.0005, -0.0005, -0.0005, -0.0005}; // the minus sign disables the width 
  Params_SmallWidth["s_e"] = {0.001, 0.006, 0.006, 0.006};  
  Params_SmallWidth["LR_iterations"] = {50}; 
  Params_SmallWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_SmallWidth["G_Stdev"] = {0.001, 0.001, 0.001, 0.001}; 
  Params_SmallWidth["cache"] = {10000}; 
  Params_SmallWidth["x_range"] = {0.01, 9.0}; 

  std::map<TString, std::vector<float>> Params_NormalWidth; 
  Params_NormalWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_NormalWidth["m_e"] = {m, m, m, m}; 
  Params_NormalWidth["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params_NormalWidth["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params_NormalWidth["LR_iterations"] = {50}; 
  Params_NormalWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_NormalWidth["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params_NormalWidth["cache"] = {10000}; 
  Params_NormalWidth["x_range"] = {0.01, 9.0}; 

  std::map<TString, std::vector<float>> Params;
  
  std::map<TString, std::vector<TH1F*>> Hists = MC_Reader_All("Merged_MC.root");
  typedef std::map<TString, std::vector<TH1F*>>::iterator it;
 
  TFile* F = new TFile("Individual_Ntracks.root", "RECREATE"); 
  int sdo = 4;
  int z = 0; 
  for (it i = Hists.begin(); i != Hists.end(); i++)
  {
    TString Name = i -> first; 

    std::vector<TH1F*> Hist = i -> second; 
    if (Name.Contains("radius") || Name.Contains("Truth")){continue;}
    
    bool skip = false; 
    std::vector<TH1F*> Inside = Hists[Name + "_radius_Less"]; 
    std::vector<TH1F*> Inside_Truth = Hists[Name + "_radius_Less_Truth"]; 
    TH1F* trk1_input = Hists[Name + "_radius_Greater"][0]; 
    for (int z(0); z < 2; z++)
    {
      int E = Inside[z] -> GetEntries(); 
      if (E < 1000){skip = true;} 
    }
    if (trk1_input -> GetEntries() < 1000){skip = true;}
    if (skip == true){continue;}

    std::vector<TH1F*> Inside_Using; 
    std::vector<std::vector<TH1F*>> Inside_Using_Truth; 
    for (int z(0); z < Inside.size(); z++)
    {
      int E = Inside[z] -> GetEntries(); 
      if (E < 1000){continue;}

     
      std::vector<TH1F*> Temp; 
      for (int y(z*sdo); y < (z+1)*sdo; y++){Temp.push_back(Inside_Truth.at(y));}
      Inside_Using_Truth.push_back(Temp); 

      if (sdo != 1){Inside_Using.push_back(Inside[z]); }
      else 
      {
        Inside_Using.push_back(Inside[z]);
        break;
      }
    }
    
    // Truth for n-track, n-true 
    std::vector<TH1F*> Truth_Temp;
    for (int g(0); g < Inside_Using_Truth.size(); g++){Truth_Temp.push_back(Inside_Using_Truth[g][g]);}

    // Normalization only tested:
    std::vector<TH1F*> ntrk_Conv_Normal = ConvolveNTimes(trk1_input, Inside_Using_Truth.size(), Name + "_Normal"); 
    for (int g(0); g < Inside_Using_Truth.size(); g++)
    {
      std::map<TString, std::vector<float>> NormalFit_R = ScalingFit(Inside_Using_Truth[g][g], {ntrk_Conv_Normal[g]}); 
    }
    
    // Normalization + Shifting d - x method test: 
    std::vector<TH1F*> ntrk_Conv_ShiftNormal = ConvolveNTimes(trk1_input, Inside_Using_Truth.size(), Name + "_ShiftNormal"); 
    for (int g(0); g < Inside_Using_Truth.size(); g++)
    {
      std::map<TString, std::vector<float>> ShiftNormalFit_R = ScalingShift(Inside_Using_Truth[g][g], {ntrk_Conv_ShiftNormal[g]}); 
    }
 
    // Normalization + Shifting using FFT test: 
    std::vector<TH1F*> ntrk_Conv_ShiftNormalFFT = ConvolveNTimes(trk1_input, Inside_Using_Truth.size(), Name + "_ShiftNormalFFT"); 
    for (int g(0); g < Inside_Using_Truth.size(); g++)
    {
      std::vector<std::pair<TH1F*, std::vector<float>>> ShiftNormalFFT_R = FitDeconvolutionPerformance(Inside_Using_Truth[g][g], {ntrk_Conv_ShiftNormalFFT[g]}, Params_SmallWidth, Params_SmallWidth["cache"][0], Params_SmallWidth["cache"][0]); 
        
      ntrk_Conv_ShiftNormalFFT[g] -> Reset();
      ntrk_Conv_ShiftNormalFFT[g] -> Add(ShiftNormalFFT_R[0].first, 1); 
      delete ShiftNormalFFT_R[0].first; 
    }
    
    // Normalization + Shifting + Width using FFT test: 
    std::vector<TH1F*> ntrk_Conv_ShiftNormalWidthFFT = ConvolveNTimes(trk1_input, Inside_Using_Truth.size(), Name + "_ShiftNormalWidthFFT"); 
    for (int g(0); g < Inside_Using_Truth.size(); g++)
    {
      std::vector<TH1F*> ShiftNormalWidthFFT_R = LoopGen({ntrk_Conv_ShiftNormalWidthFFT[g]}, Inside_Using_Truth[g][g], Params_NormalWidth); 
        
      ntrk_Conv_ShiftNormalWidthFFT[g] -> Reset();
      ntrk_Conv_ShiftNormalWidthFFT[g] -> Add(ShiftNormalWidthFFT_R[0], 1); 
      delete ShiftNormalWidthFFT_R[0]; 
    }
    
    TString Root_dir = Name + "/"; 
    F -> mkdir(Root_dir); 
    F -> cd(Root_dir); 
    
    TString Normal_dir = Root_dir + "Normal"; 
    F -> mkdir(Normal_dir); 
    F -> cd(Normal_dir); 

    F -> mkdir(Normal_dir + "/Fits"); 
    F -> cd(Normal_dir + "/Fits");  
    BulkWrite(ntrk_Conv_Normal); 
    
    F -> cd(Normal_dir); 
    
    F -> mkdir(Normal_dir + "/Truth"); 
    F -> cd(Normal_dir + "/Truth"); 
    BulkWrite(Truth_Temp);  

    F -> cd(Root_dir); 


    TString NormalShift_dir = Root_dir + "ShiftNormal"; 
    F -> mkdir(NormalShift_dir); 
    F -> cd(NormalShift_dir); 

    F -> mkdir(NormalShift_dir + "/Fits"); 
    F -> cd(NormalShift_dir + "/Fits");  
    BulkWrite(ntrk_Conv_ShiftNormal); 
    
    F -> cd(NormalShift_dir); 
    
    F -> mkdir(NormalShift_dir + "/Truth"); 
    F -> cd(NormalShift_dir + "/Truth"); 
    BulkWrite(Truth_Temp);  

    F -> cd(Root_dir); 


    TString ShiftNormalFFT_dir = Root_dir + "ShiftNormalFFT"; 
    F -> mkdir(ShiftNormalFFT_dir); 
    F -> cd(ShiftNormalFFT_dir); 

    F -> mkdir(ShiftNormalFFT_dir + "/Fits"); 
    F -> cd(ShiftNormalFFT_dir + "/Fits");  
    BulkWrite(ntrk_Conv_ShiftNormalFFT); 
    
    F -> cd(ShiftNormalFFT_dir); 
    
    F -> mkdir(ShiftNormalFFT_dir + "/Truth"); 
    F -> cd(ShiftNormalFFT_dir + "/Truth"); 
    BulkWrite(Truth_Temp);  

    F -> cd(Root_dir);  


    TString ShiftNormalWidthFFT_dir = Root_dir + "ShiftNormalWidthFFT"; 
    F -> mkdir(ShiftNormalWidthFFT_dir); 
    F -> cd(ShiftNormalWidthFFT_dir); 

    F -> mkdir(ShiftNormalWidthFFT_dir + "/Fits"); 
    F -> cd(ShiftNormalWidthFFT_dir + "/Fits");  
    BulkWrite(ntrk_Conv_ShiftNormalWidthFFT); 
    
    F -> cd(ShiftNormalWidthFFT_dir); 
    
    F -> mkdir(ShiftNormalWidthFFT_dir + "/Truth"); 
    F -> cd(ShiftNormalWidthFFT_dir + "/Truth"); 
    BulkWrite(Truth_Temp);  


    F -> cd(); 

    BulkDelete(ntrk_Conv_Normal); 
    BulkDelete(ntrk_Conv_ShiftNormal); 
    BulkDelete(ntrk_Conv_ShiftNormalFFT); 
    BulkDelete(ntrk_Conv_ShiftNormalWidthFFT);
  } 
  F -> Close(); 
}




