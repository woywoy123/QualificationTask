#include<PostAnalysis/BaseFunctions.h>


std::vector<TH1F*> BaseFunctions::MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension)
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

std::vector<TH1F*> BaseFunctions::MakeTH1F(std::vector<TString> Names, TH1F* Input)
{
  std::vector<TH1F*> Output(Names.size()); 
  for (int i(0); i < Names.size(); i++)
  {
    TH1F* H = (TH1F*)Input -> Clone(Names[i]);
    H -> Reset();
    H -> SetTitle(Names[i]);
    Output[i] = H;  
  }
  return Output;
}

std::vector<TH1F*> BaseFunctions::CopyTH1F(std::vector<TH1F*> Hists, TString Extension)
{
  std::vector<TH1F*> Output(Hists.size()); 
  for (int i(0); i < Hists.size(); i++)
  {
    TString Name = Hists[i] -> GetTitle(); 
    TH1F* H = (TH1F*)Hists[i] -> Clone(Name + Extension);
    H -> SetTitle(Name + Extension);
    Output[i] = H;  
  }
  return Output;
}


std::vector<float> BaseFunctions::Ratio(std::vector<TH1F*> Hists, TH1F* Data)
{
  float Lumi = Data -> Integral(); 
  std::vector<float> Ratios;  
  for (TH1F* H : Hists)
  {
    float e = H -> Integral(); 
    Ratios.push_back(e/Lumi);  
  }
  return Ratios; 
}

std::vector<float> BaseFunctions::Ratio(std::vector<RooRealVar*> Vars, TH1F* Data)
{
  float Lumi = Data -> Integral(); 
  std::vector<float> Output; 
  for (RooRealVar* v : Vars)
  {
    Output.push_back((v -> getVal())/Lumi);   
    delete v; 
  }
  return Output;
}

std::vector<float> BaseFunctions::ClosureAndData(std::vector<TH1F*> Hists, TH1F* Data)
{
  for (TH1F* H : Hists){Data -> Add(H);}
  return Ratio(Hists, Data); 
}

void BaseFunctions::Normalize(TH1F* Hist)
{
  float sum(0); 
  for (int i(0); i < Hist -> GetNbinsX(); i++)
  {
    sum = sum + Hist -> GetBinContent(i+1); 
  }
  Hist -> Scale(1/sum);
}

void BaseFunctions::Normalize(std::vector<TH1F*> Hist)
{
  for (TH1F* H : Hist)
  {
    Normalize(H); 
  }

}

void BaseFunctions::Subtraction(std::vector<TH1F*> ntrk, TH1F* Data, int Exclude, std::vector<float> Ratios)
{
  float lumi = Data -> Integral(); 
  for (int i(0); i < ntrk.size(); i++)
  {
    Normalize(ntrk[i]);
    ntrk[i] -> Scale(Ratios[i]*lumi);
    if ( i != Exclude -1 )
    {  
      Data -> Add(ntrk[i], -1); 
    } 
  }
}

void BaseFunctions::Scale(std::vector<TH1F*> PDFs, std::vector<RooRealVar*> Vars)
{
  for (int i(0); i < PDFs.size(); i++)
  {
    PDFs[i] -> Scale(Vars[i] -> getVal()); 
    delete Vars[i];  
  }
}

std::vector<RooRealVar*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables(Names.size()); 
  for (int i(0); i < Names.size(); i++)
  {
    Variables[i] = new RooRealVar(Names[i]+"_Fit", Names[i], Begin[i], End[i]);  
  }
  return Variables;
}

std::vector<RooRealVar*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<float> Var1, std::vector<float> Var2, std::vector<float> Var3)
{
  std::vector<RooRealVar*> Variables(Names.size());
  for (int i(0); i < Names.size(); i++)
  {
    Variables[i] = new RooRealVar(Names[i], Names[i], Var1[i], Var2[i], Var3[i]); 
  }
  return Variables; 
}

std::vector<RooGaussian*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<RooRealVar*> Mean, std::vector<RooRealVar*> Stdev, RooRealVar* Domain)
{
  std::vector<RooGaussian*> Gaussian(Names.size());
  for (int i(0); i < Names.size(); i++)
  {
    Gaussian[i] = new RooGaussian(Names[i], Names[i], *Domain, *Mean[i], *Stdev[i]); 
  }
  return Gaussian;
}

std::vector<RooGaussian*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<RooRealVar*> V1, std::vector<float> V2, std::vector<float> V3)
{
  std::vector<RooGaussian*> Gaussian(Names.size());
  for (int i(0); i < Names.size(); i++)
  {
    Gaussian[i] = new RooGaussian(Names[i], Names[i], *V1[i], RooFit::RooConst(V2[i]), RooFit::RooConst(V3[i])); 
  }
  return Gaussian;
  
}

std::vector<RooFFTConvPdf*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<RooHistPdf*> PDFs, std::vector<RooGaussian*> Gaus, RooRealVar* Domain)
{
  std::vector<RooFFTConvPdf*> FFT(Names.size());
  for (int i(0); i < Names.size(); i++)
  {
    FFT[i] = new RooFFTConvPdf(Names[i], Names[i], *Domain, *PDFs[i], *Gaus[i]);
  }
  return FFT; 
}

std::vector<RooDataHist*> BaseFunctions::RooData(std::vector<TH1F*> Hist, RooRealVar* Domain)
{
  std::vector<RooDataHist*> DataHist(Hist.size());
  for (int i(0); i < Hist.size(); i++)
  {
    TString name = Hist[i] -> GetTitle(); 
    DataHist[i] = new RooDataHist(name, name, *Domain, Hist[i]);
  }
  return DataHist;
}

RooDataHist* BaseFunctions::RooData(TH1F* Hist, RooRealVar* Domain)
{
  std::vector<TH1F*> H = {Hist};
  std::vector<RooDataHist*> D = RooData(H, Domain);
  return D[0];
}

std::vector<RooHistPdf*> BaseFunctions::RooPDF(std::vector<TH1F*> Hist, RooRealVar* Domain)
{
  std::vector<RooDataHist*> Data = RooData(Hist, Domain);
  std::vector<RooHistPdf*> HistPdf(Hist.size());
  for (int i(0); i < Data.size(); i++)
  {
    TString name = Hist[i] -> GetTitle(); 
    HistPdf[i] = new RooHistPdf(name, name, *Domain, *Data[i]);
  }
  return HistPdf; 
}

RooArgList BaseFunctions::RooList(std::vector<RooRealVar*> Vector)
{
  RooArgList L;
  for (RooRealVar* a : Vector){L.add(*a);}
  return L; 
}

RooArgList BaseFunctions::RooList(std::vector<RooHistPdf*> Vector)
{
  RooArgList L;
  for (RooHistPdf* a : Vector){L.add(*a);}
  return L; 
}

float BaseFunctions::ChiSquare(std::vector<float> V1, std::vector<float> V2)
{
  float di(0);
  for (int i(0); i < V1.size(); i++)
  {
    di += pow(V1[i] - V2[i], 2);
  }
  return di;
}

void BaseFunctions::PredictionTruthPrint(std::vector<float> Truth, std::vector<float> Prediction)
{
  for (int i(0); i < Truth.size(); i++)
  {
    std::cout << "trk-" << i+1 << " Truth: " << Truth[i] << " Prediction: " << Prediction[i] << std::endl;
  }
}

std::vector<float> BaseFunctions::LucyRichardson(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y)
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

void BaseFunctions::ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv)
{
  int nBins_2 = Hist2 -> GetNbinsX();
  int nBins_1 = Hist1 -> GetNbinsX();
   
  // Make sure the bins are equally sized  
  if ( nBins_2 != nBins_1 ) { return; }

  // Convert TH1 to vector 
  std::vector<float> H1, H2;
  for ( unsigned int i(0); i < nBins_2; i++ )
  {
    H1.push_back(Hist1 -> GetBinContent(i+1));
    H2.push_back(Hist2 -> GetBinContent(i+1));
  }

  // Set bin content of conv histogram 
  std::vector<float> Conv = ConvolveHists(H1, H2);
  for ( unsigned int i(0); i < nBins_2; i++ )
  {
    conv -> SetBinContent(i+1, Conv.at(i));
  }
}

std::vector<float> BaseFunctions::ConvolveHists(std::vector<float> Hist1, std::vector<float> Hist2)
{
  int n = Hist1.size();
  std::vector<float> conv(n, 0);

  // Initialize the FFT method by giving it the data points  
  TVirtualFFT* fft1 = TVirtualFFT::FFT(1, &n, "R2C K P");
  TVirtualFFT* fft2 = TVirtualFFT::FFT(1, &n, "R2C K P");
  
  for ( Int_t i(0); i < n; i++ )
  {
    fft1 -> SetPoint(i, Hist1.at(i), 0); 
    fft2 -> SetPoint(i, Hist2.at(i), 0); 
  }
  fft1 -> Transform();
  fft2 -> Transform();

  // Main part of the FFT 
  TVirtualFFT* fft2r = TVirtualFFT::FFT(1, &n, "C2R K P");  
  for ( Int_t i(0); i < n/2 +1; i++ )
  {
    Double_t r1, r2, i1, i2;
    fft1 -> GetPointComplex(i, r1, i1);   
    fft2 -> GetPointComplex(i, r2, i2); 
    
    Double_t re = r1*r2 - i1*i2;
    Double_t im = r1*i2 + r2*i1;
    
    TComplex t(re, im);
    fft2r -> SetPointComplex(i, t);  
  }

  // Reverse FFT to real space
  fft2r -> Transform();
  for ( Int_t i(0); i < n; i++ )
  {
    Double_t r1, i1;
    fft2r -> GetPointComplex(i, r1, i1);
    conv[i] = r1;
  }
  delete fft1; 
  delete fft2;
  delete fft2r;
  
  return conv;
}

std::vector<float> BaseFunctions::TH1FDataVector(TH1F* Data, float Offset)
{
  int bins = Data -> GetNbinsX();
  std::vector<float> Data_Vector(bins + bins*Offset, 0);
  for (int i(0); i < bins -1; i++)
  {
    Data_Vector[i] = Data -> GetBinContent(i+1);
    if (i < bins*Offset+1) { Data_Vector[i+bins-1] = Data -> GetBinContent(bins - i - 1); }
  }

  return Data_Vector;
}

void BaseFunctions::ToTH1F(std::vector<float> Vector, TH1F* Hist)
{
  for (int i(0); i < Vector.size(); i++)
  {
    Hist -> SetBinContent(i+1, Vector[i]);
  }
}

void BaseFunctions::ShiftExpandTH1F(TH1F* In, TH1F* Out, int start)
{
  int binI = In -> GetNbinsX();
  int binO = Out -> GetNbinsX();
  int Padding = (binO - binI)/2;

  for (int i(0); i < binI; i++)
  {
    float e = In -> GetBinContent(i+1);
    Out -> SetBinContent(Padding + start + i + 1, e);
  }
}

void BaseFunctions::ShiftExpandTH1F(std::vector<TH1F*> In, std::vector<TH1F*> Out, int start)
{
  for (int i(0); i < In.size(); i++)
  {
    ShiftExpandTH1F(In[i], Out[i], start);
  }
}

//void BaseFunctions::ResidualRemove(TH1F* Hist)
//{
//  int max = Hist -> GetMaximumBin();
// 
//  float T = 1e10; 
//  int iter = 0; 
//  int breaker = 0; 
//  for (int i(1); i < max; i++)
//  {
//    float e = Hist -> GetBinContent(max - i - 1);
//    if (e < T)
//    {
//      T = e;
//      iter = max - i;
//      breaker=0;
//    }
//    if ( e > T ){breaker++;}
//    if (breaker == 2){break;} 
//  }
//  for (int i(0); i < iter; i++)
//  {
//    Hist -> SetBinContent(i+1, 1e-9);
//  }
//}

void BaseFunctions::ResidualRemove(TH1F* Hist)
{
  float bin_m = Hist -> GetMaximumBin();
  float sum = Hist -> Integral(1, bin_m);
  float Average = sum/bin_m;
 
  float T = sum;  
  int iter(0); 
  int breaker = 0; 
  for (int i(0); i < bin_m; i++)
  {
    float e = Hist -> GetBinContent(bin_m-i-1);
    if (e < T && e <= Average)
    {
      T = e;
      iter = bin_m - i;
      breaker = 0;  
    }
    else{breaker++;}
    if (breaker == 2) {break;}
  }
 
  for (int i(0); i < iter; i++)
  {
    Hist -> SetBinContent(i+1, 1e-20);
  }

}

