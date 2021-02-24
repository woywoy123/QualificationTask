#include<PostAnalysis/AlgorithmFunctions.h>

void ShapeSigmoid(TH1F* trk_Fit, TH1F* ntrk_Conv)
{
  for (int z(0); z < trk_Fit -> GetNbinsX(); z++)
  {
    float e = trk_Fit -> GetBinContent(z+1); 
    float f = ntrk_Conv -> GetBinContent(z+1); 
    float log = std::exp((std::log(e)*(0.5)+std::log(f)*(0.5))); 
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
      float log = std::exp((std::log(ed)*0.9 + (0.1)*std::log(point))); 
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


std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, std::vector<TH1F*> PSF, TH1F* Data, std::map<TString, std::vector<float>> Params)
{
  std::vector<TString> Names_Dec; 
  int index = 0;  
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString name = "Dec_"; name += (ntrk_Conv[i] -> GetTitle()); 
    Names_Dec.push_back(name);
    std::cout << ntrk_Conv[i] -> GetTitle() << std::endl;
  }
  std::vector<TH1F*> PDF_D = CloneTH1F(Data, Names_Dec);
  MultiThreadingDeconvolutionExperimental(ntrk_Conv, PSF, PDF_D, Params["LR_iterations"][0]); 
  std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data, PDF_D, Params, Params["cache"][0], Params["cache"][0]);
  for (int i(0); i < Names_Dec.size(); i++){delete PDF_D[i];}
 
  std::vector<TH1F*> F_C;
  for (int i(0); i < trk_Fit.size(); i++){F_C.push_back(trk_Fit[i].first); }
 
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
  int bins = Data[0] -> GetNbinsX(); 
  float min = Data[0] -> GetXaxis() -> GetXmin(); 
  float max = Data[0] -> GetXaxis() -> GetXmax(); 

  std::vector<TH1F*> PSF; 
  for (int x(0); x < Data.size(); x++)
  {
    TString name = "Gaussian_"; name += (x+1); 
    float m = Params["G_Mean"][x]; 
    float s = Params["G_Stdev"][x];
    TH1F* Gaus = Gaussian(m ,s, bins, min, max, name); 
    PSF.push_back(Gaus);
  }
  
  Data[trk_Data] -> SetLineColor(kBlack); 
  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(InputTrk1, Data.size(), "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
  std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
  Flush(F_C, ntrk_Conv, false); 
  
  TString name = Data[trk_Data] -> GetTitle(); name += ("_Experimental.pdf"); 
  std::vector<TString> Names_Dec; 
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString name = "Temp_"; name += (ntrk_Conv[i] -> GetTitle()); 
    Names_Dec.push_back(name);
  }

  PlotHists(Data_Copy, ntrk_Conv, can); 
  can -> Print(name); 
  std::map<TString, std::vector<TH1F*>> Out; 
  for (int i(0); i < iterations; i++)
  {
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> trk; 
    std::vector<TH1F*> psf;
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p != trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else
      { 
        trk.push_back(ntrk_Conv[p]); 
        psf.push_back(PSF[p]);  
      }
    }
    F_C = LoopGen(trk, psf, Data_Copy, Params); 
    ScaleShape(Data_Copy, trk);  
    Flush(F_C, trk, true); 
    PlotHists(Data_Copy, trk, can); 
    can -> Print(name); 
    delete Data_Copy;    
  
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> not_trk; 
    std::vector<TH1F*> not_psf;
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p == trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else
      { 
        not_trk.push_back(ntrk_Conv[p]); 
        not_psf.push_back(PSF[p]);  
      }
    }
    Smear(Data_Copy, 0.05);
    F_C = LoopGen(not_trk, not_psf, Data_Copy, Params); 
    Flush(F_C, not_trk, true); 
    delete Data_Copy;    
    
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    ScalingFit(Data_Copy, ntrk_Conv);  
    ScaleShape(Data_Copy, ntrk_Conv); 
    PlotHists(Data_Copy, ntrk_Conv, can); 
    can -> Print(name); 

    delete Data_Copy;   


    TString base = "_ntrk_"; base += (trk_Data+1); base += ("_iter_"); base += (i); 
    TString iter = "Iteration_"; iter += (i); 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      TString name = ntrk_Conv[x] -> GetTitle(); name += (base); name += (Data[trk_Data]->GetTitle());  
      TH1F* H = (TH1F*)ntrk_Conv[x] -> Clone(name); 
      H -> SetTitle(name); 
      
      for (int p(0); p < H -> GetNbinsX(); p++)
      {
        H -> SetBinError(p+1, 1e-9); 
      }
      Out[iter].push_back(H); 
    }
  }
  Out["Data"].push_back(Data[trk_Data]);
  BulkDelete(ntrk_Conv);
  BulkDelete(PSF); 
  delete can;
  return Out; 
}


