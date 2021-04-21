#include<PostAnalysis/AlgorithmFunctions.h>

std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params)
{
  std::vector<std::pair<TH1F*, std::vector<float>>> trks = LoopGenAll(ntrk_Conv, Data, Params);  
  
  std::vector<TH1F*> F_C; 
  for (int i(0); i < trks.size(); i++)
  {
    F_C.push_back(trks[i].first); 
  }
  return F_C;
}

std::vector<std::pair<TH1F*, std::vector<float>>> LoopGenAll(std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params)
{
  std::vector<TString> Names_Dec; 
  float r = 0.2; 
  int bins = Data -> GetNbinsX(); 
  float min = Data -> GetXaxis() -> GetXmin(); 
  float max = Data -> GetXaxis() -> GetXmax(); 
  float w = (max - min) / float(bins); 
  float new_min = min - w*bins*r; 
  float new_max = max + w*bins*r; 

  std::vector<TH1F*> PDF_D;
  std::vector<TH1F*> PSF;  
  std::vector<TH1F*> PDF; 
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString nameG = "Gx_"; nameG+=(i+1); 
    TH1F* Gaus = Gaussian(Params["G_Mean"][i], Params["G_Stdev"][i], bins+2*bins*r, new_min, new_max, nameG); 
    PSF.push_back(Gaus);  

    TString name = "Dec_"; name += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* H = new TH1F(name, name, bins+2*bins*r, new_min, new_max); 
    PDF_D.push_back(H); 
  
    TString name_L = "L_"; name_L += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* X = new TH1F(name_L, name_L, bins+2*bins*r, new_min, new_max); 
    PDF.push_back(X); 

    for (int j(0); j < bins; j++)
    {
      X -> SetBinContent(j+1+r*bins, ntrk_Conv[i] -> GetBinContent(j+1)); 
    }
  }
  Normalize(PDF); 
  Normalize(PSF); 
  MultiThreadingDeconvolutionExperimental(PDF, PSF, PDF_D, Params["LR_iterations"][0]); 
  
  TString t_ = Data -> GetTitle(); 
  TH1F* Data_D = new TH1F(t_ +"D", t_ +"D", bins+2*bins*r, new_min, new_max); 
  for (int j(0); j < bins; j++)
  {
    Data_D -> SetBinContent(j+1+r*bins, Data -> GetBinContent(j+1)); 
    Data_D -> SetBinError(j+1+r*bins, Data -> GetBinError(j+1)); 
  }

  std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data_D, PDF_D, Params, Params["cache"][0], Params["cache"][0]);
  delete Data_D; 
  std::vector<std::pair<TH1F*, std::vector<float>>> F_C;

  float bin_min = PDF[0] -> GetXaxis() -> FindBin(min) -1; 
  for (int i(0); i < trk_Fit.size(); i++)
  {
    TString name_N = "N_"; name_N += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* H = new TH1F(name_N, name_N, bins, min, max); 
   
    for (int j(0); j < bins; j++)
    {
      H -> SetBinContent(j+1, trk_Fit[i].first -> GetBinContent(j+1 + bin_min)); 
    }
    
    F_C.push_back(std::pair<TH1F*, std::vector<float>>(H, trk_Fit[i].second));
    delete trk_Fit[i].first;
  } 
  
  BulkDelete(PDF_D); 
  BulkDelete(PSF); 
  BulkDelete(PDF); 
  return F_C; 
}

std::vector<std::vector<TH1F*>> BaseTemplates(TH1F* trk1, int ntrks, TString Name_Algo)
{
  std::vector<std::vector<TH1F*>> Conv; 
  for (int i(0); i < ntrks; i++)
  {
    TString N = "NTRK"; N +=(i+1); N += ("_"); 
    std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(trk1, ntrks, N + Name_Algo); 
    Conv.push_back(ntrk_Conv);
  }
  return Conv; 
}

std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Normalization(std::vector<TH1F*> Data, TH1F* trk1, std::map<TString, std::vector<float>> Params, TString Name)
{
  std::vector<std::vector<TH1F*>> Conv = BaseTemplates(trk1, Data.size(), Name); 

  std::map<TString, std::pair<std::vector<float>, std::vector<TH1F*>>> Output; 
  for (int i(0); i < Data.size(); i++)
  {
    TH1F* trk = Data[i]; 
    std::map<TString, std::vector<float>> Result = ScalingFit(trk, Conv[i]); 

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

































