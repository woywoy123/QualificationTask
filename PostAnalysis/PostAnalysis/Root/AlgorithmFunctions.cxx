#include<PostAnalysis/AlgorithmFunctions.h>

void ShapeSigmoid(TH1F* trk_Fit, TH1F* ntrk_Conv)
{
  for (int z(0); z < trk_Fit -> GetNbinsX(); z++)
  {
    float e = trk_Fit -> GetBinContent(z+1); 
    float f = ntrk_Conv -> GetBinContent(z+1); 
    
    float sig = std::exp(-std::abs((e/f)));
    float log = 1./ (1 + sig); 
   
    if (f <= 0 && e <= 0)
    { 
      f = 0; 
      log = 0; 
      e = std::abs((ntrk_Conv -> GetBinContent(z-1) + e + ntrk_Conv -> GetBinContent(z))/3.); 
    } 
    //std::cout << log << std::endl;
    ntrk_Conv -> SetBinContent(z+1, e*(1-log)+f*log); 
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
      float ed = point*(e/(sum));  
      float sig = std::exp(-std::abs(sum/e));
      float log = 1./ (1 + sig); 
      if (sum <= 0 || ed <= 0)
      {
        ed = std::abs((ntrk[z] -> GetBinContent(i-1) + point + ntrk[z] -> GetBinContent(i))/3.);
        log = 0; 
        point = 0; 
      }
      
      //std::cout << ed << " ::: " << log << " ::: " << point << std::endl;
      ntrk[z] -> SetBinContent(i+1, ed*(1-log)+point*log); 
    }     
  }
}

void Flush(std::vector<TH1F*> F_C, std::vector<TH1F*> ntrk_Conv, bool sig)
{
  for (int i(0); i < F_C.size(); i++)
  {
    if (sig == false)
    {
      ntrk_Conv[i] -> Reset(); 
      ntrk_Conv[i] -> Add(F_C[i], 1); 
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
  for (int i(0); i < Data -> GetNbinsX(); i++)
  {
    float y = Data -> GetBinContent(i+1); 
    float sum = 0; 
    int p = 0; 
    if (y > 0) {continue;}
    for (int x(-1); x < 2; x++)
    {
      float e = Data -> GetBinContent(i + x + 1); 
      if (e < 0){continue;}
      sum += e;
      p++; 
    }
    if (p == 0){ p = 1; }
    sum = sum / float(p); 
    Data -> SetBinContent(i+1, sum); 
  }
}


std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, std::vector<TH1F*> PSF, TH1F* Data, std::map<TString, std::vector<float>> Params)
{
  std::vector<TString> Names_Dec; 
  int index = 0;  
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString name = "Dec_"; name += (ntrk_Conv[i] -> GetTitle()); 
    Names_Dec.push_back(name);
  }
  std::vector<TH1F*> PDF_D = CloneTH1F(Data, Names_Dec);
  
  MultiThreadingDeconvolutionExperimental(ntrk_Conv, PSF, PDF_D, Params["LR_iterations"][0]); 
  std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data, PDF_D, Params, Params["cache"][0], Params["cache"][0]);
  
  std::vector<TH1F*> F_C;
  for (int i(0); i < trk_Fit.size(); i++){F_C.push_back(trk_Fit[i].first); }
  
  for (int i(0); i < PDF_D.size(); i++){delete PDF_D[i];}
  return F_C; 
}


std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, int trk_Data)
{
  
  int iterations = Params["iterations"][0]; 
  int LucyRichardson = Params["LR_iterations"][0]; 
  int bins = Data[0] -> GetNbinsX(); 
  float min = Data[0] -> GetXaxis() -> GetXmin(); 
  float max = Data[0] -> GetXaxis() -> GetXmax(); 
  float width = (max - min) / float(bins); 
  min += width / 2.; 
  max += width / 2.; 
  
  std::vector<TH1F*> PSF; 
  for (int x(0); x < Data.size(); x++)
  {
    TString name = "Gaussian_"; name += (x+1); 
    float m = Params["G_Mean"][x]; 
    float s = Params["G_Stdev"][x];
    TH1F* Gaus = Gaussian(m ,s, bins, min, max, name); 
    PSF.push_back(Gaus);
  }
  std::map<TString, std::vector<TH1F*>> Out; 
  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(Data[0], Data.size(), "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 

  std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
  Flush(F_C, ntrk_Conv, false); 

  TCanvas* can = new TCanvas(); 
  TString name = Data[trk_Data] -> GetTitle(); name += (".pdf"); 
  can -> SetLogy();
  for (int i(0); i < iterations; i++)
  {
    delete Data_Copy; 
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 

    std::vector<TH1F*> Not_trk; 
    std::vector<TH1F*> Not_PSF; 
    std::vector<TH1F*> Truth_Not;  
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x == trk_Data){ Data_Copy -> Add(ntrk_Conv[x], -1); }
      else
      {
        Not_trk.push_back(ntrk_Conv[x]); 
        Not_PSF.push_back(PSF[x]);
      }
    }

    F_C = LoopGen(Not_trk, Not_PSF, Data_Copy, Params);     
    ScaleShape(Data_Copy, F_C); 
    Flush(F_C, Not_trk, true);
    PlotHists(Data_Copy, ntrk_Conv, can); 
    can -> Print(name); 

    delete Data_Copy; 
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 

    std::vector<TH1F*> Delta;
    std::vector<TH1F*> Delta_PSF; 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x != trk_Data){Data_Copy -> Add(ntrk_Conv[x], -1);}
     
      if ( x == trk_Data )
      {
        Delta.push_back(ntrk_Conv[x]); 
        Delta_PSF.push_back(PSF[x]); 
      }
    }
    
    F_C = LoopGen(Delta, Delta_PSF, Data_Copy, Params);     
    ScaleShape(Data_Copy, F_C); 
    Flush(F_C, Delta, true);  

    PlotHists(Data_Copy, ntrk_Conv, can); 
    can -> Print(name); 


    delete Data_Copy; 
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
     
    TString base = "_ntrk_"; base += (trk_Data+1); base += ("_iter_"); base += (i); 
    TString iter = "Iteration_"; iter += (i); 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      TString name = ntrk_Conv[x] -> GetTitle(); name += (base); 
      TH1F* H = (TH1F*)ntrk_Conv[x] -> Clone(name); 
      H -> SetTitle(name); 
      
      for (int p(0); p < H -> GetNbinsX(); p++)
      {
        H -> SetBinError(p+1, 1e-9); 
      }
      Out[iter].push_back(H); 
    }
  }
  BulkDelete(ntrk_Conv);
  BulkDelete(PSF); 
  return Out; 
}


