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
    int z(0);  
    for (std::pair<TH1F*, std::vector<float>> H : Result)
    {
      Temp_H.push_back(H.first); 
      Error.push_back(H.second[5]); 
      ntrk_Conv[z] -> Reset(); 
      ntrk_Conv[z] -> Add(H.first, 1); 
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

  float m = 0.2;
  std::map<TString, std::vector<float>> Params_SmallWidth; 
  Params_SmallWidth["m_s"] = {-m, -m, -m, -m}; 
  Params_SmallWidth["m_e"] = {m, m, m, m}; 
  Params_SmallWidth["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params_SmallWidth["s_e"] = {0.006, 0.006, 0.006, 0.006};  
  Params_SmallWidth["iterations"] = {1}; 
  Params_SmallWidth["LR_iterations"] = {50}; 
  Params_SmallWidth["G_Mean"] = {0, 0, 0, 0}; 
  Params_SmallWidth["G_Stdev"] = {0.001, 0.001, 0.001, 0.001}; 
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

    // Doing the different fit techniques 
    // Normalization fit 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Normal = OnlyNormal(Inside_Using, Outside, Params, Name); 
    //AnalysisPlots(Normal, Inside_Using_Truth, Name + "_Normal");  
    WriteToFile(F, Normal, Name + "_Normal"); 

    // Shifting with Normalization
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> ShiftNormal = NormalShift(Inside_Using, Outside, Params, Name);
    //AnalysisPlots(ShiftNormal, Inside_Using_Truth, Name + "_ShiftNormal");  
    WriteToFile(F, ShiftNormal, Name + "_ShiftNormal"); 


    // Shifting Normalization with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShiftFFT = ShiftNormalFFT(Inside_Using, Outside, Params_SmallWidth, Name); 
    //AnalysisPlots(NormalShiftFFT, Inside_Using_Truth, Name + "_NormalShiftFFT");  
    WriteToFile(F, NormalShiftFFT, Name + "_ShiftNormalFFT"); 

    // Shifting, Normalization, Width with FFT 
    std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> NormalShiftWidthFFT = NormalWidthDeconvShiftFFT(Inside_Using, Outside, Params_NormalWidth, Name);
    //AnalysisPlots(NormalShiftFFT, Inside_Using_Truth, Name + "_NormalShiftWidthFFT");  
    WriteToFile(F, NormalShiftWidthFFT, Name + "_NormalShiftWidthFFT"); 
  }
  F -> Close(); 
  delete F;
}

void Evaluation()
{
  AnalysisLoop(1, "1_trk_1_tru.root"); 
  AnalysisLoop(4, "All_Fits.root");
  //Evaluate_Fits("Fits_Results.root");  
}




