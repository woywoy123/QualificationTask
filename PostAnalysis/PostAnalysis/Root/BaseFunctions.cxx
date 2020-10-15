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
  float Lumi = 1;//Data -> Integral(); 
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

void BaseFunctions::Scale(std::vector<TH1F*> PDFs, std::vector<RooRealVar*> Vars)
{
  for (int i(0); i < PDFs.size(); i++)
  {
    PDFs[i] -> Scale(Vars[i] -> getVal()); 
    delete Vars[i];  
  }
}

std::vector<RooRealVar*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<float> Guess)
{
  std::vector<RooRealVar*> Variables(Guess.size());
  for (int i(0); i < Guess.size(); i++)
  {
    Variables[i] = new RooRealVar(Names[i], Names[i], Guess[i]); 
  }
  return Variables;
}

std::vector<RooRealVar*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables(Begin.size()); 
  for (int i(0); i < Begin.size(); i++)
  {
    Variables[i] = new RooRealVar(Names[i], Names[i], Begin[i], End[i]);  
  }
  return Variables;
}

std::vector<RooRealVar*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<float> Var1, std::vector<float> Var2, std::vector<float> Var3)
{
  std::vector<RooRealVar*> Variables(Var1.size());
  for (int i(0); i < Var1.size(); i++)
  {
    Variables[i] = new RooRealVar(Names[i], Names[i], Var1[i], Var2[i], Var3[i]); 
  }
  return Variables; 
}

std::vector<RooGaussian*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<RooRealVar*> Mean, std::vector<RooRealVar*> Stdev, RooRealVar* Domain)
{
  std::vector<RooGaussian*> Gaussian(Mean.size());
  for (int i(0); i < Mean.size(); i++)
  {
    Gaussian[i] = new RooGaussian(Names[i], Names[i], *Domain, *Mean[i], *Stdev[i]); 
  }
  return Gaussian;
}

std::vector<RooGaussian*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<RooRealVar*> V1, std::vector<float> V2, std::vector<float> V3)
{
  std::vector<RooGaussian*> Gaussian(V1.size());
  for (int i(0); i < V1.size(); i++)
  {
    Gaussian[i] = new RooGaussian(Names[i], Names[i], *V1[i], RooFit::RooConst(V2[i]), RooFit::RooConst(V3[i])); 
  }
  return Gaussian;
  
}

std::vector<RooFFTConvPdf*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<RooGaussian*> Gaus, std::vector<RooHistPdf*> PDFs, RooRealVar* Domain)
{
  std::vector<RooFFTConvPdf*> FFT(PDFs.size());
  for (int i(0); i < PDFs.size(); i++)
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
    DataHist[i] = new RooDataHist(name, name, *Domain, RooFit::Import(*Hist[i]));
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
    TString name = Hist[i] -> GetTitle(); name += ("_Fit"); 
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
	float sum_V1(0); 
	float sum_V2(0); 
	for (int i(0); i < V1.size(); i++)
	{
		sum_V1 += V1[i]; 
		sum_V2 += V2[i]; 
	}

  for (int i(0); i < V1.size(); i++)
  {
    di += pow(V1[i]*1/sum_V1 - V2[i]*1/sum_V2, 2);
  }
  return di;
}

float BaseFunctions::ChiSquare(TH1F* H1, TH1F* H2)
{
	auto vector =[](TH1F* H1)
	{
		std::vector<float> V1;
		for (int i(0); i < H1 -> GetNbinsX(); i++)
		{
			V1.push_back(H1 -> GetBinContent(i+1)); 
		}
		return V1; 
	};
	
	std::vector<float> v1 = vector(H1);
	std::vector<float> v2 = vector(H2); 
	return ChiSquare(v1, v2); 		

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
  for ( int i(0); i < nBins_2; i++ )
  {
    H1.push_back(Hist1 -> GetBinContent(i+1));
    H2.push_back(Hist2 -> GetBinContent(i+1));
  }

  // Set bin content of conv histogram 
  std::vector<float> Conv = ConvolveHists(H1, H2);
  for ( int i(0); i < nBins_2; i++ )
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
    float e = In -> GetBinContent(i +1);
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

void BaseFunctions::ResidualRemove(TH1F* Hist)
{
  float bin_m = Hist -> GetMaximumBin();
  float T = Hist -> GetBinContent(bin_m);  
  int iter(0); 
  int breaker = 0; 
  for (int i(0); i < bin_m; i++)
  {
    float e = Hist -> GetBinContent(bin_m-i-1);
    if (e < T)
    {
      T = e;
      iter = bin_m - i;
      breaker = 0;  
    }
    else{breaker++;}
    if (breaker == 4) {break;}
  }
 
  for (int i(0); i < iter; i++)
  {
    Hist -> SetBinContent(i+1, 1e-8);
  }

}

float BaseFunctions::FLost(std::vector<TH1F*> ntrk_Data, std::vector<std::vector<TH1F*>> Truth_Sets)
{

  TH1F* trk1_D = ntrk_Data[0]; 
  TH1F* trk2_D = ntrk_Data[1];
  TH1F* trk3_D = ntrk_Data[2];
  TH1F* trk4_D = ntrk_Data[3];

  std::vector<TH1F*> trk1_PDF = Truth_Sets[0];
  std::vector<TH1F*> trk2_PDF = Truth_Sets[1]; 
  std::vector<TH1F*> trk3_PDF = Truth_Sets[2]; 
  std::vector<TH1F*> trk4_PDF = Truth_Sets[3]; 

  std::vector<TString> Names = {"ntrk1", "ntrk2", "ntrk3", "ntrk4"};

  // Lost Tracks in 1-Track 
  // === Sum up the total number of tracks 
  float Lumi_D1 = trk1_D -> Integral(); 
  float e1_D1 = trk1_PDF[0] -> Integral();
  float e2_D1 = trk1_PDF[1] -> Integral();
  float e3_D1 = trk1_PDF[2] -> Integral();
  float e4_D1 = trk1_PDF[3] -> Integral();
  
  // === Sum up the total number of tracks 
  float Lumi_D2 = trk2_D -> Integral(); 
  float e1_D2 = trk2_PDF[0] -> Integral();
  float e2_D2 = trk2_PDF[1] -> Integral();
  float e3_D2 = trk2_PDF[2] -> Integral();
  float e4_D2 = trk2_PDF[3] -> Integral();
 
  // === Sum up the total number of tracks 
  float Lumi_D3 = trk3_D -> Integral(); 
  float e1_D3 = trk3_PDF[0] -> Integral();
  float e2_D3 = trk3_PDF[1] -> Integral();
  float e3_D3 = trk3_PDF[2] -> Integral();
  float e4_D3 = trk3_PDF[3] -> Integral();
   
  // === Sum up the total number of tracks 
  float Lumi_D4 = trk4_D -> Integral(); 
  float e1_D4 = trk4_PDF[0] -> Integral();
  float e2_D4 = trk4_PDF[1] -> Integral();
  float e3_D4 = trk4_PDF[2] -> Integral();
  float e4_D4 = trk4_PDF[3] -> Integral();
   
  // FLost 1: The number of tracks lost in the 1-track measurement 
  float Nom_1 = 1*e2_D1 + 2*e3_D1 + 3*e4_D1;
  float Den_1 = e1_D1 + 2*e2_D1 + 3*e3_D1 + 4*e4_D1;  
 
  // FLost 2: The number of tracks lost in the 1-track measurement 
  float Nom_2 = 1*e3_D2 + 2*e4_D2;
  float Den_2 = 2*e2_D2 + 3*e3_D2 + 4*e4_D2; 
 
  // FLost 3: The number of tracks lost in the 1-track measurement 
  float Nom_3 = 1*e4_D3;
  float Den_3 = 3*e3_D3 + 4*e4_D3; 
 
  float FLost = (Nom_1 + 2*Nom_2 + 3*Nom_3) / (Den_1 + Den_2 + Den_3);
  
  return FLost;
}

void BaseFunctions::SetBinError(TH1F* Hist, double Error)
{
  for (int i(0); i < Hist -> GetNbinsX(); i++)
  {
    Hist -> SetBinError(i+1, Error); 
  }
}

void BaseFunctions::SetBinError(std::vector<TH1F*> Hist, double Error)
{
  for (TH1F* H : Hist)
  {
    SetBinError(H, Error);  
  }
}

void BaseFunctions::CopyBinErrors(TH1F* Source, TH1F* Target)
{
  int bins = Source -> GetNbinsX(); 
  for (int i(0); i < bins; i++)
  {
    float es = Source -> GetBinError(i+1); 
    if (es == 0) { es = 1e-5; }
    Target -> SetBinError(i+1, es); 
  }
}

void BaseFunctions::CopyBinErrors(std::vector<TH1F*> Source, std::vector<TH1F*> Target)
{
  for (int i(0); i < Source.size(); i++)
  {
    CopyBinErrors(Source[i], Target[i]); 
  }
}

void BaseFunctions::CopyBinErrors(TH1F* Source, std::vector<TH1F*> Target)
{
  for (int i(0); i < Target.size(); i++)
  {
    CopyBinErrors(Source, Target[i]);     
  }
}

void BaseFunctions::SetPercentError(TH1F* Hist, float percent)
{
  int bins = Hist -> GetNbinsX(); 
  for (int i(0); i < bins; i++)
  {
    float e = Hist -> GetBinContent(i+1); 
    Hist -> SetBinError(i+1, e*percent);
  }
}

void BaseFunctions::SetPercentError(std::vector<TH1F*> Hists, float percent)
{
  for (int i(0); i < Hists.size(); i++)
  {
    SetPercentError(Hists[i], percent); 
  }
}

