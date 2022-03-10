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


