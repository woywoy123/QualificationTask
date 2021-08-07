#include<PostAnalysis/Convolution.h>


// =========================== Convolution ================================= //
std::vector<float> ConvolutionSum(std::vector<float> V1, std::vector<float> V2, int ZeroPointBin)
{
  int n = V1.size(); 
  std::vector<float> conv(2*n, 0); 

  for (int i(0); i < 2*n; i++)
  {
    for (int k(0); i - k >= 0; k++)
    {
      if (k >= n || i - k >= n){continue;}
      conv[i] += V1[i-k] * V2[k]; 
    }
  }
  std::vector<float> result( conv.begin() + ZeroPointBin, conv.begin() + n + ZeroPointBin); 
  return result; 
}

std::vector<TH1F*> ConvolveNTimes(TH1F* Start, int n, std::vector<TString> base, TString extension)
{
  std::vector<TString> Names;  
  for (int i(0); i < n; i++)
  {
    base[i] += (extension);  
  }
  std::vector<TH1F*> Hist_V = CloneTH1F(Start, base);
  for (int i(0); i < Start -> GetNbinsX(); i++)
  {
    Hist_V[0] -> SetBinContent(i+1, Start -> GetBinContent(i+1)); 
  }

  for (int i(0); i < n-1; i++)
  {
    Convolution(Hist_V[i], Start, Hist_V[i+1]); 
  }
  
  Normalize(Hist_V); 
  return Hist_V; 
}

void Convolution(TH1F* Hist1, TH1F* Hist2, TH1F* conv)
{
  int bins = Hist2 -> GetNbinsX(); 
  int ZeroPointBin = Hist2 -> GetXaxis() -> FindBin(0.) - 1;    
  int Padding = 2;   
  
  std::vector<float> H1(bins + 2*Padding, 0), H2(bins + 2*Padding, 0);  
  for ( int i(0); i < bins; i++)
  {
    H1[i+Padding] = Hist1 -> GetBinContent(i+1); 
    H2[i+Padding] = Hist2 -> GetBinContent(i+1);   
  }

  std::vector<float> Conv = ConvolutionSum(H1, H2, ZeroPointBin + Padding); 

  for ( int i(0); i < bins; i++)
  {
    conv -> SetBinContent(i+1, Conv.at(Padding+i)); 
  }
}

// ============================ Deconvolution Stuff ======================== //
std::vector<float> LucyRichardson(std::vector<float> data, std::vector<float> psf, std::vector<float> current, float damp)
{
  float offset = (data.size() - current.size())/2; 
  std::vector<float> next(current.size(), 0); 
  for (int i(0); i < current.size(); i++)
  {
    float sum_j = 0; 
    int upperLimitJ = psf.size(); 
    for (int j(i); j < upperLimitJ; j++)
    {
      float c_j = 0; 
      int upperLimitK = j; 
      for (int k(0); k <= upperLimitK; k++)
      {
        c_j += psf[j-k]*current[k];    
      }
      if (c_j != 0)
      {
        sum_j += data[j] / c_j * psf[j-i];
      }
    }
    next[i] = current[i] * (1+damp*(sum_j-1)); 
    if (next[i] < 0. || std::isnan(next[i]) || std::isinf(next[i]) ){next[i] = 0.;}
  }
  return next;
}

void Deconvolution(TH1F* Signal, TH1F* PSF, TH1F* Out, int iter)
{
  auto ToVector =[](TH1F* H, int bins)
  {
    int bins_H = H -> GetNbinsX(); 

    std::vector<float> Output; 
    for (int i(0); i < bins_H; i++)
    {
      float e = H -> GetBinContent(i+1);
      Output.push_back(e); 
    }
    return Output; 
  };

  int bins = Signal -> GetNbinsX(); 
  std::vector<float> Signal_V = ToVector(Signal, bins); 
  std::vector<float> PSF_V = ToVector(PSF, bins); 
  std::vector<float> Deconv_V = std::vector<float>(bins, 1.); 

  for (int i(0); i < iter; i++)
  {
    Deconv_V = LucyRichardson(Signal_V, PSF_V, Deconv_V, 0.95); 
  }

  int bin_0 = Signal -> GetXaxis() -> FindBin(0.) -1; 
  for (int i(0); i < bins; i++)
  {
    Out -> SetBinContent(i+bin_0+1, Deconv_V[i]);
  }
}

void MultiThreadingDeconvolution(std::vector<TH1F*> Data, std::vector<TH1F*> PSF, std::vector<TH1F*> Result, int Iter)
{

  std::vector<std::thread> th; 
  for (int i(0); i < Result.size(); i++){th.push_back(std::thread(Deconvolution, Data[i], PSF[i], Result[i], Iter));}
  for (std::thread &t : th){t.join();}

}

void Smooth1Trk(TH1F* trk1_start, int iter)
{
  float L = trk1_start -> Integral(); 

  int bins = trk1_start -> GetNbinsX(); 
  float min = trk1_start -> GetXaxis() -> GetXmin(); 
  float max = trk1_start -> GetXaxis() -> GetXmax(); 
  float w = (max - min) / float(bins); 
  float new_min = min - w*bins; 
  float new_max = max + w*bins; 

  TString name = "Dec_"; 
  TH1F* G = new TH1F(name, name, bins+2*bins, new_min, new_max); 
  
  TString name_L = "L_";
  TH1F* X = new TH1F(name_L, name_L, bins+2*bins, new_min, new_max); 

  for (int j(0); j < bins; j++)
  {
    X -> SetBinContent(j+1+bins, trk1_start -> GetBinContent(j+1)); 
  }

  Convolution(X, X, G); 
  Deconvolution(G, X, G, iter);

  for (int j(0); j < bins; j++)
  {
    float e = G -> GetBinContent(j+1+bins); 
    trk1_start -> SetBinContent(j+1, e);    
  }
  Normalize(trk1_start); 
  trk1_start -> Scale(L); 
  
  delete G, X; 
}