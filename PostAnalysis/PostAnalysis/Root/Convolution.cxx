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

std::vector<TH1F*> ConvolveNTimes(TH1F* Start, int n, TString extension)
{
  std::vector<TString> Names;  
  for (int i(0); i < n; i++)
  {
    TString name = "TRK_"; name +=(i+1); name +=("_"); name += (extension);  
    Names.push_back(name);  
  }
  std::vector<TH1F*> Hist_V = CloneTH1F(Start, Names);
  Hist_V[0] -> Add(Start, 1);  
  
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

std::vector<float> Deconvolution(TH1F* PDF, TH1F* PSF, TH1F* Output, int Max_Iter)
{
  int pdf_bins = PDF -> GetNbinsX(); 
  int psf_bins = PSF -> GetNbinsX(); 
  int out_bins = Output -> GetNbinsX();

  // Domain of PDF and PSF 
  float pdf_min = PDF -> GetXaxis() -> GetXmin(); 
  float pdf_max = PDF -> GetXaxis() -> GetXmax(); 
  float psf_min = PSF -> GetXaxis() -> GetXmin();  
  float psf_max = PSF -> GetXaxis() -> GetXmax(); 
  float out_min = Output -> GetXaxis() -> GetXmin();  
  float out_max = Output -> GetXaxis() -> GetXmax(); 
  
  // Get out of the function - Cant deconvolve 
  std::vector<float> Converge;

  if ((pdf_bins != psf_bins || pdf_min != psf_min || pdf_max != psf_max) || (out_bins != psf_bins || out_min != psf_min || out_max != psf_max))
  {
    Converge.push_back(0);
    std::cout << "###################################" << std::endl;
    std::cout << "Check the centering of the bins...." << std::endl;
    std::cout << pdf_bins << " psf ->: " << psf_bins << " \n " << 
    pdf_min << " psf ->: " << psf_min << " \n " << 
    pdf_max << " :: psf ->: " << psf_max << std::endl;
    std::cout << "###################################" << std::endl;
    return Converge;
  }
  int bins = pdf_bins; 

  // Find the zero bin of the x-axis 
  int bin_0 = PDF -> GetXaxis() -> FindBin(0.) -1; 
  
  // Convert histograms to vectors
  std::vector<float> Temp1 = ToVector(PSF); 
  std::vector<float> Temp2 = ToVector(PDF);
  std::vector<float> PSF_V;  
  std::vector<float> PDF_V; 
  
  for ( int i(0); i < bins*2; i++)
  { 
    if (i < bins)
    { 
      PSF_V.push_back(Temp1.at(i)); 
    }
    else { PSF_V.push_back(0.); }

    if (i >= bins )
    { 
      PDF_V.push_back(Temp2.at(i-bins)); 
    }
    else { PDF_V.push_back(0.); }
  } 

  std::vector<float> Deconv_V(bins*2, 0.5);
  float d_old = 100; 
  float d = 100; 

  for (int i(0); i < Max_Iter; i++)
  {
    d_old = d; 
    std::vector<float> Deconv_Vold = Deconv_V; 
    Deconv_V = LucyRichardson(PDF_V, PSF_V, Deconv_V, 0.95); 

    d = Pythagoras(Deconv_Vold, Deconv_V); 
    Converge.push_back(d);  
    //if (d_old - d == 1e-8){break;} 
  }
  
  // Find the bin where the X axis is 0 
  int out_0 = Output -> GetXaxis() -> FindBin(0.)-1;  
  
  // Populate the output  
  for (int i(0); i < Deconv_V.size(); i++)
  {
    Output -> SetBinContent(-bins + out_0 + i+1, Deconv_V[i]); 
  }
 
  return Converge;
}

void DeconvolutionExperimental(TH1F* Signal, TH1F* PSF, TH1F* Out, int iter)
{
  auto ToVector =[](TH1F* H, int bins)
  {
    int bins_H = H -> GetNbinsX(); 
    int Padding = (bins - bins_H); 

    std::vector<float> Output; 
    for (int i(0); i < bins_H; i++)
    {
      float e = H -> GetBinContent(i+1);
      Output.push_back(e); 
    }
    for (int i(0); i < Padding; i++)
    {
      float e = H -> GetBinContent(bins_H - i - 1); 
      Output.push_back(e); 
    }
    return Output; 
  };

  float r = 1.1; 
  int bins = Signal -> GetNbinsX(); 
  int bin_0 = Signal -> GetXaxis() -> FindBin(0.) -1; 
  std::vector<float> Signal_V = ToVector(Signal, bins*r); 
  std::vector<float> PSF_V = ToVector(PSF, bins*r); 
  std::vector<float> Deconv_V = std::vector<float>(bins*r, 0.5); 

  for (int i(0); i < iter; i++)
  {
    Deconv_V = LucyRichardson(Signal_V, PSF_V, Deconv_V, 0.95); 
  }

  int Padding = ((r - 1)*bins)/2; 
  for (int i(0); i < bins; i++)
  {
    Out -> SetBinContent(i+1 + bin_0, Deconv_V[i]);
  }
}

void MultiThreadingDeconvolution(std::vector<TH1F*> Data, std::vector<TH1F*> PSF, std::vector<TH1F*> Result, int Iter)
{

  std::vector<std::thread> th; 
  for (int i(0); i < Result.size(); i++){th.push_back(std::thread(Deconvolution, Data[i], PSF[i], Result[i], Iter));}
  for (std::thread &t : th){t.join();}

}

void MultiThreadingDeconvolutionExperimental(std::vector<TH1F*> Data, std::vector<TH1F*> PSF, std::vector<TH1F*> Result, int Iter)
{

  std::vector<std::thread> th; 
  for (int i(0); i < Result.size(); i++){th.push_back(std::thread(DeconvolutionExperimental, Data[i], PSF[i], Result[i], Iter));}
  for (std::thread &t : th){t.join();}

}



