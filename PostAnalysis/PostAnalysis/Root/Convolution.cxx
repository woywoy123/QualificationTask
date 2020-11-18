#include<PostAnalysis/Convolution.h>


// =========================== Convolution ================================= //
std::vector<float> ConvolutionFFT(const std::vector<float> V1, const std::vector<float> V2, int ZeroPointBin)
{
  int n = V1.size(); 
  std::vector<float> conv(n, 0); 

  // Create two FFT objects 
  TVirtualFFT* fft1 = TVirtualFFT::FFT(1, &n, "R2C K P"); 
  TVirtualFFT* fft2 = TVirtualFFT::FFT(1, &n, "R2C K P"); 

  // Shift the vector such that the convolution vector starts at the zero point on the domain 
  // First half: Points that are past the Zero point of the domain  
  for (Int_t i(0); i < n - ZeroPointBin; i++)
  {
    fft1 -> SetPoint(i, V1.at(i+ZeroPointBin), 0);
    fft2 -> SetPoint(i, V2.at(i+ZeroPointBin), 0);  
  }

  // Second half: Points that are below the Zero point of the domain 
  for (Int_t i(0); i < ZeroPointBin; i++)
  {
    fft1 -> SetPoint(n - ZeroPointBin + i, V1.at(i), 0);
    fft2 -> SetPoint(n - ZeroPointBin + i, V2.at(i), 0);  
  }
 
  // Change points into Frequency domain  
  fft1 -> Transform(); 
  fft2 -> Transform(); 
  
  // Perform the multiplication of the FFT
  TVirtualFFT* fft2r =TVirtualFFT::FFT(1, &n, "C2R K P "); 
  for ( int i(0); i < n/2 +1; i++)
  {
    Double_t r1, r2, i1, i2; // Define the real and imaginary points 
   
    // Populate the points above  
    fft1 -> GetPointComplex(i, r1, i1);
    fft2 -> GetPointComplex(i, r2, i2); 

    // Perform the multiplication A*A = (r1 + i1)*(r2 - i2) = r1*r2 - i1*i2 + i1*r2 - i2*r1
    Double_t re = r1*r2 - i1*i2; 
    Double_t im = i1*r2 + i2*r1; 
    TComplex t(re, im); 
    fft2r -> SetPointComplex(i, t); 
  }

  // Change fft to real space
  fft2r -> Transform(); 
  for ( int i(0); i < n; i++)
  {
    Double_t r1, i1; 
    fft2r -> GetPointComplex(i, r1, i1); 
    
    // Make sure to revert the domain back
    if (i < n - ZeroPointBin){ conv[i+ZeroPointBin] = r1; }
    else { conv[i-ZeroPointBin] = r1; }
  }

  delete fft1; 
  delete fft2; 
  delete fft2r;

  return conv; 
}

void ResidualRemove(TH1F* Hist)
{
  float bin_m = Hist -> GetMaximumBin(); 
  float T = Hist -> GetBinContent(bin_m); 
  int iter(0); 
  int breaker = 0; 
  for (int i(0); i < bin_m; i++)
  {
    float e = Hist -> GetBinContent(bin_m -1 -i); 
    if (e < T) 
    {
      T = e; 
      iter = bin_m -i; 
      breaker = 0; 
    }
    else { breaker++; }
    if (breaker == 4) {break;}
  }
  
  for (int i(0); i < iter; i++)
  {
    Hist -> SetBinContent(i+1, 0); 
  }
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

  std::vector<float> Conv = ConvolutionFFT(H1, H2, ZeroPointBin + Padding); 

  for ( int i(0); i < bins; i++)
  {
    conv -> SetBinContent(i+1, Conv.at(Padding+i)); 
  }
}



// ============================ Deconvolution Stuff ======================== //
std::vector<float> LucyRichardson(std::vector<float> data, std::vector<float> current, std::vector<float> psf, float damp)
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
        sum_j += data[j] / (c_j * psf[j-i]);
      }
    }
    next[i] = current[i] + (1+damp*(sum_j-1)); 
    if (next[i] < 0. || std::isnan(next[i]) || std::isinf(next[i]) ){next[i] = 0.;}
  }
  return next;
}

std::vector<float> Deconvolution(TH1F* PDF, TH1F* PSF, TH1F* Output, int Max_Iter)
{
  int pdf_bins = PDF -> GetNbinsX(); 
  int psf_bins = PSF -> GetNbinsX(); 

  // Domain of PDF and PSF 
  float pdf_min = PDF -> GetXaxis() -> GetXmin(); 
  float pdf_max = PDF -> GetXaxis() -> GetXmax(); 
  float psf_min = PSF -> GetXaxis() -> GetXmin();  
  float psf_max = PSF -> GetXaxis() -> GetXmax(); 

  // Make sure the bin widths are equal before proceeding 
  float width_pdf = (pdf_max - pdf_min) / float(pdf_bins);   
  float width_psf = (psf_max - psf_min) / float(psf_bins); 

  // Get out of the function - Cant deconvolve 
  std::vector<float> Converge;
  if (width_pdf != width_psf)
  {
    Converge.push_back(0);
    return Converge;
  }

  // Unify their domains by finding the max and mins of the two distributions 
  float domain_min; 
  float domain_max; 
  if (pdf_min < psf_min){domain_min = pdf_min;}
  else { domain_min = psf_min; }

  if (pdf_max > psf_max){domain_max = pdf_max;}
  else { domain_max = psf_max; }

  // Create the new histograms for the new domain definition 
  int bins = (domain_max - domain_min)/width_psf; 

  // Create unique names for these hists for multithreading  
  TString unique = PDF -> GetTitle(); unique += ("_"); unique += (PSF -> GetTitle()); unique += ("_"); 
  TH1F* H1 = new TH1F(unique + "H1", unique + "H1", bins, domain_min, domain_max); 
  TH1F* H2 = new TH1F(unique + "H2", unique + "H2", bins, domain_min, domain_max); 

  // Find the zero bin position of both histograms 
  int psf_0 = PSF -> GetXaxis() -> FindBin(0.) -1;  
  int pdf_0 = PDF -> GetXaxis() -> FindBin(0.) -1;  
  int bin_0 = H1 -> GetXaxis() -> FindBin(0.) - 1; 

  // Create iterators 
  int i_psf = 0; 
  int i_pdf = 0; 
  for (int i(0); i < bins; i++)
  {
    int d_bin_0 = bin_0 - i; 

    // Fill histogram based on PSF
    if (d_bin_0 <= psf_0)
    {
      H1 -> SetBinContent(i+1, PSF -> GetBinContent(i_psf+1)); 
      i_psf++;
    } 

    // Fill histogram based on PDF
    if (d_bin_0 <= pdf_0)
    {
      H2 -> SetBinContent(i+1, PDF -> GetBinContent(i_pdf+1)); 
      i_pdf++;
    } 
  }
  
  // Convert histograms to vectors
  std::vector<float> PSF_V = ToVector(H1);
  std::vector<float> PDF_V = ToVector(H2); 
  std::vector<float> Deconv_V(bins, 0.5);
  float d_old = 100; 
  float d = 100; 

  for (int i(0); i < Max_Iter; i++)
  {
    d_old = d; 
    std::vector<float> Deconv_Vold = Deconv_V; 
    Deconv_V = LucyRichardson(PDF_V, PSF_V, Deconv_V, 0.75); 

    d = Pythagoras(Deconv_Vold, Deconv_V); 
    Converge.push_back(d);  
    if (d_old - d == 1e-8){break;} 
  }
  
  TH1F* Deconv_H = new TH1F(unique + "Deconv", unique + "Deconv", bins, domain_min, domain_max);
  ToTH1F(Deconv_V, Deconv_H); 
  Output -> Add(Deconv_H);   

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  //Deconv_H -> Draw("HIST*"); 
  H1 -> Draw("SAMEHIST-"); 
  H2 -> Draw("SAMEHIST"); 
  can -> Print("debug.pdf"); 


  
  delete H1; 
  delete H2;  
  delete Deconv_H;  
  return Converge;
}






























