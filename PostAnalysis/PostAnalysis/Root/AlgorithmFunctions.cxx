#include<PostAnalysis/AlgorithmFunctions.h>

void ShapeSigmoid(TH1F* trk_Fit, TH1F* ntrk_Conv)
{
  for (int z(0); z < trk_Fit -> GetNbinsX(); z++)
  {
    float e = trk_Fit -> GetBinContent(z+1); 
    float f = ntrk_Conv -> GetBinContent(z+1); 
    float log = e; 
    if (std::isnan(log)){log = 0;}
    if (std::isinf(log)){log = 0;}
    ntrk_Conv -> SetBinContent(z+1, log); 
  } 
}

void Scale(TH1F* Data, std::vector<TH1F*> ntrk)
{
  int excess_i = 0;
  float excess_sum;  
  for (int i(0); i < Data ->GetNbinsX(); i++)
  {
    float e = Data -> GetBinContent(i+1); 
    float sum(0);  
    for (int c(0); c < ntrk.size(); c++)
    {
      sum += ntrk[c] -> GetBinContent(i+1);  
    }
    
    if ( sum > e && e > 0)
    {
      excess_sum += sum/e; 
      excess_i++;  
    }
  }
   
  if ( excess_i > 1)
  {
    float average = excess_sum / (float)excess_i; 
    for (int c(0); c < ntrk.size(); c++){ntrk[c] -> Scale((float)1/average);}
  }
}

void ScaleShape(TH1F* Data, std::vector<TH1F*> ntrk)
{
  for (int i(0); i < Data ->GetNbinsX(); i++)
  {
    float e = Data -> GetBinContent(i+1); 
    float sum(0);  
    for (int c(0); c < ntrk.size(); c++)
    {
      sum += ntrk[c] -> GetBinContent(i+1);  
    }
    for (int z(0); z < ntrk.size(); z++)
    {
      float point = ntrk[z] -> GetBinContent(i+1); 
      float ed = point*(e/sum);  
      float log = ed;  
      if (std::isnan(log)){log = 0;}
      if (std::isinf(log)){log = 0;} 
      ntrk[z] -> SetBinContent(i+1, log); 
    }     
  }
}

float FitError(TH1F* Data, std::vector<TH1F*> ntrk)
{
  int bins = Data -> GetNbinsX(); 
  float sum_D = 0; 
  float sum_E = 0;  
  for (int i(0); i < bins; i++)
  {
    float e = Data -> GetBinContent(i+1); 
    float sum_H = 0; 
    for (int p(0); p < ntrk.size(); p++){sum_H += std::abs(ntrk[p] -> GetBinContent(i+1));}
    sum_E += std::abs(sum_H - e);  
  }
  sum_E = sum_E/Data -> Integral();
  return sum_E;
}



void Flush(std::vector<TH1F*> F_C, std::vector<TH1F*> ntrk_Conv, bool sig)
{
  if (F_C.size() != ntrk_Conv.size()) {std::cout << "!!!!!!!!!!!!!!!!!" << std::endl; }
  for (int i(0); i < F_C.size(); i++)
  {
    if (sig == false)
    {
      for (int x(0); x < F_C[i] -> GetNbinsX(); x++)
      {
        ntrk_Conv[i] -> SetBinContent(x+1, F_C[i] -> GetBinContent(x+1));   
      }
    }
    else
    {
      ShapeSigmoid(F_C[i], ntrk_Conv[i]); 
    }
    delete F_C[i];  
  }
}

void Average(TH1F* Data)
{
  std::vector<float> D; 
  for (int i(0); i < Data -> GetNbinsX(); i++)
  {
    float y = Data -> GetBinContent(i+1); 
    float sum = 0; 
    int p = 0; 
    for (int x(0); x < 1; x++)
    {
      float e = Data -> GetBinContent(i + x + 1); 
      if (e < 0){ e = 0; }
      sum += e;
      p++; 
    }
    if ( sum < 1e-10 ){sum = 0;}
    sum = sum / float(p); 
    D.push_back(sum);  
  }
  for (int i(0); i < Data -> GetNbinsX(); i++)
  {
    Data -> SetBinContent(i+1, D[i]); 
  }
}


std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params)
{
  std::vector<TString> Names_Dec; 
  float r = 0.1; 
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
  MultiThreadingDeconvolutionExperimental(PDF, PSF, PDF_D, Params["LR_iterations"][0]); 
  
  std::vector<TH1F*> PDF_D_N; 
  float bin_min = PDF[0] -> GetXaxis() -> FindBin(min) -1; 
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString name_N = "N_"; name_N += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* H = new TH1F(name_N, name_N, bins, min, max); 
    PDF_D_N.push_back(H);  
   
    for (int j(0); j < bins; j++)
    {
      H -> SetBinContent(j+1, PDF_D[i] -> GetBinContent(j+1 + bin_min)); 
    }
  } 
  
  std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data, PDF_D_N, Params, Params["cache"][0], Params["cache"][0]);
  std::vector<TH1F*> F_C;
  for (int i(0); i < trk_Fit.size(); i++){F_C.push_back(trk_Fit[i].first); }

  BulkDelete(PDF_D); 
  BulkDelete(PSF); 
  BulkDelete(PDF); 
  BulkDelete(PDF_D_N); 
  return F_C; 
}


std::map<TString, std::vector<TH1F*>> MainAlgorithm(TH1F* InputTrk1, std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, int trk_Data)
{
  
  auto Smear =[](TH1F* Data, float stdev)
  {
    float lumi = Data -> Integral(); 
    int bins = Data -> GetNbinsX(); 
    float min = Data -> GetXaxis() -> GetXmin(); 
    float max = Data -> GetXaxis() -> GetXmax(); 
    TH1F* Gaus = Gaussian(0 ,stdev , bins, min, max, "Tempo"); 
    Convolution(Data, Gaus, Data); 
    Normalize(Data);
    Data -> Scale(lumi); 
    delete Gaus;  
  }; 
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int iterations = Params["iterations"][0]; 

  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(InputTrk1, 4, "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
  std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, Data_Copy, Params); 
  Flush(F_C, ntrk_Conv, false); 

  TString name = Data[trk_Data] -> GetTitle(); name += (".pdf"); 
  PlotHists(Data_Copy, ntrk_Conv, can); 
  can -> Print(name); 

  std::map<TString, std::vector<TH1F*>> Out; 
  for (int i(0); i < iterations; i++)
  {
    
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> trk; 
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p == trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else{ trk.push_back(ntrk_Conv[p]); }
    }
    Average(Data_Copy); 
    Smear(Data_Copy, 0.05); 
    F_C = LoopGen(trk, Data_Copy, Params); 
    Flush(F_C, trk, true);
    
    PlotHists(Data_Copy, trk, can); 
    can -> Print(name); 
    delete Data_Copy; 


    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> not_trk; 
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p != trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else{ not_trk.push_back(ntrk_Conv[p]); }
    }
    F_C = LoopGen(not_trk, Data_Copy, Params); 
    ScaleShape(Data_Copy, not_trk);  
    Flush(F_C, not_trk, true); 

    PlotHists(Data_Copy, ntrk_Conv, can); 
    can -> Print(name); 
    delete Data_Copy;    

    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, Data_Copy, Params); 
    ScaleShape(Data_Copy, F_C);  
    Flush(F_C, ntrk_Conv, true); 

    PlotHists(Data_Copy, ntrk_Conv, can); 
    can -> Print(name); 
    delete Data_Copy;    
  

  } 

  for (int x(0); x < ntrk_Conv.size(); x++)
  {
    TString name2 = Data[trk_Data]->GetTitle();   
    name2.Remove(name2.Length() - 5, name2.Length()); 
    name2.Remove(0, 12); 
    name = "dEdx_ntrk_"; name += (trk_Data+1); name += ("_ntru_"); name += (x+1); name += ("_"); 
    name2 += ("_R");
    name = name + name2;
    TH1F* H = (TH1F*)ntrk_Conv[x] -> Clone(name); 
    H -> SetTitle(name); 
      
    Out["Result"].push_back(H); 
  }
  Out["Data"].push_back(Data[trk_Data]);
  BulkDelete(ntrk_Conv);
  delete can;

  return Out; 
}


