#include<PostAnalysis/BaseFunctions.h>


// Creates Histograms from vector names - Used for bulk generation 
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension)
{
  std::vector<TH1F*> Histograms;
  for (TString name : Names)
  {
    if (Extension != ""){ name+=(Extension); }
    TH1F* hist = new TH1F(name, name, bins, min, max);
    Histograms.push_back(hist);
  }

  return Histograms;
}

// Clone and Reset the histograms 
std::vector<TH1F*> CloneTH1F(TH1F* Hist, std::vector<TString> Names)
{
  float min = Hist -> GetXaxis() -> GetXmin(); 
  float max = Hist -> GetXaxis() -> GetXmax(); 
  int bins = Hist -> GetNbinsX();  
  
  std::vector<TH1F*> Output; 
  for (int i(0); i < Names.size(); i++)
  {
    TH1F* H = new TH1F(Names[i], Names[i], bins, min, max);  
    H -> Reset(); 
    H -> SetTitle(Names[i]); 
    Output.push_back(H); 
  }
  return Output; 
}

// Normalize Histogram
void Normalize(TH1F* Hist)
{
  float Lumi = Hist -> Integral(); 
  Hist -> Scale(1/Lumi); 
}

// Normalize Histograms in bulk 
void Normalize(std::vector<TH1F*> Hists)
{
  for (TH1F* H : Hists){ Normalize(H); }
}

// Normalize vectors 
std::vector<float> Normalize(std::vector<float> V1) 
{
  float sum = 0; 
  for (int i(0); i < V1.size(); i++){sum = sum + V1[i];}
  for (int i(0); i < V1.size(); i++)
  {
    V1[i] = V1[i]/sum;
  }
  return V1; 
}

// TH1F to vector
std::vector<float> ToVector(TH1F* Hist)
{
  int bins = Hist -> GetNbinsX(); 
  std::vector<float> Output; 
   
  for (int i(0); i < bins; i++)
  {
    Output.push_back(Hist -> GetBinContent(i+1));   
  }
  return Output; 
}

// Vector to TH1F
void ToTH1F(std::vector<float> Input, TH1F* Hist)
{
  for (int i(0); i < Input.size(); i++)
  {
    Hist -> SetBinContent(i+1, Input[i]); 
  }
}

// Shift a histogram using bins 
void Shift(TH1F* Hist, int shift)
{
  std::vector<float> Content; 
  for (int i(0); i < Hist -> GetNbinsX(); i++)
  {
    float e = Hist -> GetBinContent(i+1); 
    Content.push_back(e); 
  }

  Hist -> Reset(); 
  for (int i(0); i < Hist -> GetNbinsX(); i++)
  {
    Hist -> SetBinContent(i+1+shift, Content[i]);  
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
    
    if ( sum > e && e > 1)
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

void ScaleShape(TH1F* Data, std::vector<TH1F*> ntrk, float gamma)
{
  for (int i(0); i < Data ->GetNbinsX(); i++)
  {
    float e = Data -> GetBinContent(i+1); 
    float sum(0);  
    for (int c(0); c < ntrk.size(); c++)
    {
      sum += ntrk[c] -> GetBinContent(i+1);  
    }
    if (e != 0)
    {
      float p = std::abs(sum/e);
      
      for (int c(0); c < ntrk.size(); c++)
      {
        float en = ntrk[c] -> GetBinContent(i+1);
        
        ntrk[c] -> SetBinContent(i+1, 0.5*(1./p)*en + 0.5*en);  
      }
    }
    else
    {
      for (int c(0); c < ntrk.size(); c++)
      {
        ntrk[c] -> SetBinContent(i+1, e); 
      }
    }
  }
}

//Variable Name Generator 
std::vector<TString> NameGenerator(int number, TString shorty)
{
  std::vector<TString> Output; 
  
  for (int i(0); i < number; i++)
  {
    TString name = shorty; name += (i+1); 
    Output.push_back(name); 
  }
  return Output; 
}

// Benchmark function 
float Pythagoras(std::vector<float> v1, std::vector<float> v2)
{
  float sum = 0; 
  for (int i(0); i < v1.size(); i++)
  {
    float diff = pow(v1[i] - v2[i], 2); 
    sum += diff;  
  }
  return pow(sum, 0.5); 
}

float SquareError(TH1F* Hist1, TH1F* Hist2, float x_min, float x_max)
{
  int bins = Hist1 -> GetNbinsX(); 
  int bin_min = Hist1 -> GetXaxis() -> FindBin(x_min);
  int bin_max = Hist1 -> GetXaxis() -> FindBin(x_max);
  
  if (x_min == x_max)
  {
    bin_min = 0; 
    bin_max = bins; 
  }
 
  std::vector<float> v1, v2; 
  for (int i(bin_min); i < bin_max; i++)
  {
    v1.push_back(Hist1 -> GetBinContent(i+1)); 
    v2.push_back(Hist2 -> GetBinContent(i+1)); 
  }
  return Pythagoras(v1, v2); 
}

float ErrorByIntegral(TH1F* Hist1, TH1F* Hist2, float x_min, float x_max)
{ 
  int bins = Hist1 -> GetNbinsX(); 
  int bin_min = Hist1 -> GetXaxis() -> FindBin(x_min); 
  int bin_max = Hist1 -> GetXaxis() -> FindBin(x_max); 

  if (x_min == x_max)
  {
    bin_min = 0; 
    bin_max = bins; 
  }
  
  // Calculate the total error from Hist1, compared to Hist2
  float sum = 0; 
  for (int i(bin_min); i < bin_max; i++)
  {
    float h1 = Hist1 -> GetBinContent(i+1); 
    float h2 = Hist2 -> GetBinContent(i+1); 
    sum += std::abs(h1 - h2);
  } 
  float A = Hist2 -> Integral(bin_min, bin_max); 
  return sum/A; 

}

void Statistics(TH1F* H1, TH1F* H2, float x_min, float x_max)
{
  float max = H1 -> GetXaxis() -> GetXmax(); 
  float min = H1 -> GetXaxis() -> GetXmin(); 

  if (x_min == x_max) 
  {
    x_max = max; 
    x_min = min;   
  }

  std::cout << "H1: " << H1 -> GetTitle() << " ::::: H2: " << H2 -> GetTitle() << std::endl;
  std::cout << "Domain selected: " << x_min << " -> " << x_max << std::endl;
  std::cout << "- Absolute Error / Integral of H2: " << ErrorByIntegral(H1, H2, x_min, x_max) << std::endl;
  std::cout << "- KolmogorovTest: " << H2 -> KolmogorovTest(H1) << std::endl;
  std::cout << "- Normalization: H1: " << H1 -> Integral() << " :::: H2: " << H2 -> Integral() << " :::: H1/H2: " << float(H1 -> Integral()) / float(H2 -> Integral()) << std::endl; 
}

float GetMaxValue(TH1F* H)
{
  int bin_m = H -> GetMaximumBin(); 
  float v = H -> GetBinContent(bin_m+1); 
  return v; 
}

void BulkWrite(std::vector<TH1F*> Hist_V)
{
  for (TH1F* H : Hist_V){H -> Write();}
}

void BulkDelete(std::vector<TH1F*> Hist_V)
{
  for (int i(0); i < Hist_V.size(); i++)
  {
    delete Hist_V[i]; 
  }
}

std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, int trk_Data)
{

  auto ShapeSigmoid =[] (TH1F* trk_Fit, TH1F* ntrk_Conv)
  {
    for (int z(0); z < trk_Fit -> GetNbinsX(); z++)
    {
      float e = trk_Fit -> GetBinContent(z+1); 
      float f = ntrk_Conv -> GetBinContent(z+1); 
      
      if ( e == 0 && f == 0) {continue;}
      float dif = std::pow(-1 * std::abs(e-f)/abs(e+f), 1); 
      float sig = 1./(1+std::exp(dif)); 
      ntrk_Conv -> SetBinContent(z+1, e*(1-sig)+f*sig); 
    } 
  };

  auto LoopGen =[&](std::vector<TH1F*> ntrk_Conv, std::vector<TH1F*> PSF, TH1F* Data, std::map<TString, std::vector<float>> Params)
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
    for (int i(0); i < trk_Fit.size(); i++)
    {
      F_C.push_back(trk_Fit[i].first); 
    }
    
    if (ntrk_Conv.size() != 1){ScaleShape(Data, F_C, 1); }
    for (int i(0); i < ntrk_Conv.size(); i++)
    {
      //Scale(Data, F_C); 
      ShapeSigmoid(F_C[i], ntrk_Conv[i]); 
    } 

    for (int i(0); i < PDF_D.size(); i++){delete PDF_D[i];}
    for (int i(0); i < F_C.size(); i++){delete F_C[i];}
  };

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

  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(Data[0], Data.size(), "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
  LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
  Scale(Data_Copy, ntrk_Conv); 

  std::map<TString, std::vector<TH1F*>> Out; 
  for (int i(0); i < iterations; i++)
  {
    LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 

    std::vector<TH1F*> Delta;
    std::vector<TH1F*> Delta_PSF; 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x != trk_Data){Data_Copy -> Add(ntrk_Conv[x], -1);}
      else 
      { 
        Delta.push_back(ntrk_Conv[x]); 
        Delta_PSF.push_back(PSF[x]); 
      }
    }
    
    LoopGen(Delta, Delta_PSF, Data_Copy, Params); 
    
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

