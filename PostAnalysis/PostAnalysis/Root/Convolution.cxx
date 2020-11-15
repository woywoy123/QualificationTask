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
std::vector<float> LucyRichardson(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y)
{
  // G - Measured Signal
  // F - Current estimate of PSF
  // H - Estimate of truth signal P(G|F) = H
  // Alg: F(j+1) = F(j)*Sum(i) [ H(ij) * G(i) / Sum(k) [H(jk) * F(k)] ]
  // y - Dampening 
  std::vector<float> PSF(F.size(), 0);

  // Each bin is calculated separately  
  for ( unsigned int i(0); i < H.size(); i++)
  {
    float sum_i(0); 
    for (unsigned int j(i); j < F.size(); j++)
    {
      //Sum(k) [H(jk) * F(k)]
      float sum_k(0); 
      for (int k(0); k <= j; k++)
      {
        float H_jk = H[j-k];
        float F_k = F[k];
        sum_k += H_jk*F_k;
      }
      if (sum_k != 0) 
      { 
        sum_i += (G[j]*H[j-i]) / sum_k; 
      } 

      PSF[i] = F[i] * ( 1 + y * (sum_i -1) );
      if (PSF[i] < 0. || std::isnan(PSF[i]) || std::isinf(PSF[i]))  {PSF[i] = 0.;}
    }  
  }
  return PSF; 
}

