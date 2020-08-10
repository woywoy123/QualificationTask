#include<PostAnalysis/Functions.h>

std::vector<TH1F*> Functions::MakeTH1F(std::vector<TString> Names, int bins, int min, int max, TString Extension)
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

void Functions::FillTH1F_From_File(std::vector<TH1F*> Histograms, TFile* File, TString DetectorLayer, TString Extension)
{
  File -> cd(DetectorLayer);
  for (TH1F* hist : Histograms)
  {
    TString hist_name = hist -> GetName();
    if ( Extension != "") { hist_name = RemoveExtension(hist_name, Extension); }
    hist -> Add((TH1F*)gDirectory -> Get(hist_name));
  } 
  File -> cd();
}

void Functions::FillTH1F_From_File(TH1F* Histogram, TFile* File, TString DetectorLayer, std::vector<TString> List, TString Extension)
{
  File -> cd(DetectorLayer);
  for (TString HistName : List)
  {
    if ( Extension != "") { HistName = RemoveExtension(HistName, Extension); }
    Histogram -> Add((TH1F*)gDirectory -> Get(HistName));
  }
   File -> cd();
}

std::vector<TString> Functions::SplitString(TString token, TString Split)
{
  std::vector<TString> Output;
  TObjArray *array = Split.Tokenize(token);
  for (unsigned int x = 0; x < array -> GetEntries(); x++)
  {
    Output.push_back(((TObjString*)(array -> At(x))) -> String());
  }
  return Output;
}

TString Functions::RemoveExtension(TString Name, TString Extension)
{
  int l = Name.Sizeof();
  int l_e = Extension.Sizeof();
  return Name.Remove(l-l_e, l);
}

TH1F* Functions::VectorToTH1F(std::vector<float> Vec, TString name, int bins, int min, int max)
{
  TH1F* hist = new TH1F(name, name, bins, min, max);
  for (int i(0); i < bins; i++)
  {
    hist -> SetBinContent(i+1, Vec.at(i));  
  } 
  return hist;
}

void Functions::VectorToTH1F(std::vector<float> Vec, TH1F* hist, int StartBin)
{ 
  for (int i(0); i < Vec.size(); i++)
  {
    hist -> SetBinContent(i+1 + StartBin, Vec.at(i));  
  } 
}

std::vector<float> Functions::TH1FToVector(TH1F* hist, int CustomLength, int start)
{

  std::vector<float> Out;
  for (int i(0); i < start; i++)
  {
    Out.push_back(0);
  }

  for (int i(0); i < hist -> GetNbinsX(); i++)
  {
    Out.push_back(hist -> GetBinContent(i+1)); 
  }
  
  int delta = CustomLength - hist -> GetNbinsX() - start;
  for (int i(0); i < delta; i++)
  {
    Out.push_back(0);
  }

  return Out;
}

void Functions::CutTH1F(TH1F* Input, TH1F* Output)
{
  int bins = Output -> GetNbinsX();
  for (int i(0); i < bins; i++)
  {
    Output -> SetBinContent(i+1, Input -> GetBinContent(i+1));
  }
}

void Functions::ExpandTH1F(TH1F* Input, TH1F* Output, int start)
{
  int bins = Input -> GetNbinsX();
  for (int i(0); i < bins; i++)
  {
    Output -> SetBinContent(start + i +1, Input -> GetBinContent(i+1)); 
  }
}
// ================== Fitting Base Classes ==============================//
RooHistPdf* Fit_Functions::ConvertTH1FtoPDF(RooDataHist* Histogram, TString Name, RooRealVar* domain)
{
  Name+=("_PDF");
  RooHistPdf* PDF = new RooHistPdf(Name, Name, *domain, *Histogram);
  return PDF;
}

RooDataHist* Fit_Functions::ConvertTH1toDataHist(TH1F* Hist, RooRealVar* domain)
{
  TString name = Hist -> GetName();
  RooDataHist* Histo = new RooDataHist(name, name, *domain, Hist);
  return Histo;
}

RooDataHist* Fit_Functions::ConvertTH1toDataHist(TH1* Hist, RooRealVar* domain)
{
  TString name = Hist -> GetName();
  RooDataHist* Histo = new RooDataHist(name, name, *domain, Hist);
  return Histo;
}

RooRealVar* Fit_Functions::GenerateVariable(TString name, float begin, float end)
{
  RooRealVar* var = new RooRealVar(name, name, begin, end);
  return var;
}

RooGaussian* Fit_Functions::GenerateGaussian(TString name, RooRealVar* mean, RooRealVar* stdev, RooRealVar* domain)
{
  RooGaussian* g = new RooGaussian(name, name, *domain, *mean, *stdev);
  return g;
}

RooFFTConvPdf* Fit_Functions::GenerateConvolve(TString name, RooHistPdf* PDF, RooGaussian* Gaus, RooRealVar* domain)
{
  RooFFTConvPdf* Conv = new RooFFTConvPdf(name, name, *domain, *PDF, *Gaus);
  return Conv;
}


RooArgList Fit_Functions::VectorToArgList(std::vector<RooRealVar*> Vector)
{
  RooArgList List; 
  for (RooRealVar* arg : Vector)
  {
    List.add(*arg);
  }
  return List;
}

RooArgList Fit_Functions::VectorToArgList(std::vector<RooHistPdf*> Vector)
{
  RooArgList List; 
  for (RooHistPdf* arg : Vector)
  {
    List.add(*arg);
  }
  return List;
}


// ================= Fitting Derived Classes ===========================//
std::vector<RooDataHist*> Fit_Functions::ConvertTH1toDataHist(std::vector<TH1F*> Histograms, RooRealVar* domain)
{
  std::vector<RooDataHist*> DataHists;
  for (TH1F* hist : Histograms)
  {
    RooDataHist* Histo = ConvertTH1toDataHist(hist, domain);
    DataHists.push_back(Histo);
  }
  return DataHists;
}

std::vector<RooHistPdf*> Fit_Functions::ConvertTH1FtoPDF(std::vector<TH1F*> Histograms, RooRealVar* domain)
{
  std::vector<RooHistPdf*> PDFs;
  std::vector<RooDataHist*> DataHists = ConvertTH1toDataHist(Histograms, domain);
  for (unsigned int x = 0; x < DataHists.size(); x++)
  {
    RooDataHist* Histo = DataHists.at(x);
    PDFs.push_back(ConvertTH1FtoPDF(Histo, Histograms.at(x) -> GetName(), domain)); 
  }
  return PDFs;
}

std::vector<RooRealVar*> Fit_Functions::GenerateVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables(Names.size());
  for ( unsigned int x = 0; x < Names.size(); x++)
  {
         Variables[x] = GenerateVariable(Names.at(x), Begin.at(x), End.at(x));
  }
  return Variables;
}

std::vector<RooGaussian*> Fit_Functions::GaussianVariables(std::vector<TString> Names, std::vector<RooRealVar*> Mean, std::vector<RooRealVar*> Stdev, RooRealVar* Domain)
{
  std::vector<RooGaussian*> Gaussian(Names.size());
  for ( unsigned int i(0); i < Names.size(); i++)
  {
    Gaussian[i] = GenerateGaussian(Names[i], Mean[i], Stdev[i], Domain); 
  }
  return Gaussian;
}

std::vector<RooFFTConvPdf*> Fit_Functions::ConvolveVariables(std::vector<TString> Names, std::vector<RooHistPdf*> PDFs, std::vector<RooGaussian*> Gaus, RooRealVar* domain)
{
  std::vector<RooFFTConvPdf*> Conv(Names.size());
  for (unsigned int i(0); i < Names.size(); i++)
  {
    Conv[i] = GenerateConvolve(Names[i], PDFs[i], Gaus[i], domain); 
  }
  return Conv;
}

std::vector<float> Fit_Functions::LRDeconvolution(std::vector<float> G, std::vector<float> H, std::vector<float> F, float y)
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

void Fit_Functions::ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv, int StartBin)
{
  int nBins_2 = Hist2 -> GetNbinsX();
  int nBins_1 = Hist1 -> GetNbinsX();
   
  // Make sure the bins are equally sized  
  if ( nBins_2 != nBins_1 ) { return; }

  // Convert TH1 to vector 
  std::vector<float> H1, H2;
  for ( unsigned int i(StartBin); i < nBins_2; i++ )
  {
    H1.push_back(Hist1 -> GetBinContent(i+1));
    H2.push_back(Hist2 -> GetBinContent(i+1));
  }

  // Set bin content of conv histogram 
  std::vector<float> Conv = ConvolveHists(H1, H2);
  for ( unsigned int i(0); i < nBins_2 - StartBin; i++ )
  {
    conv -> SetBinContent(i+1 + StartBin, Conv.at(i));
  }
}

std::vector<float> Fit_Functions::ConvolveHists(std::vector<float> Hist1, std::vector<float> Hist2)
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

// version using RooFit
std::vector<float> Fit_Functions::TailReplace(TH1F* hist, std::vector<float> deconv)
{
  // Bin length 
  int d_bins = deconv.size();
  int h_bins = hist -> GetNbinsX();
  int delta_bins = d_bins - h_bins;
  int start = delta_bins;

  // Create the histograms 
  TH1F* Deconv = new TH1F("Deconv", "Deconv", d_bins, 0, d_bins);
  TH1F* Hist = (TH1F*)Deconv -> Clone("Hist");

  // Fill the histograms with the data entries 
  for (int i(0); i < d_bins; i++)
  {
    Deconv -> SetBinContent(i+1, deconv[i]);
    if ( i < delta_bins ){ Hist -> SetBinContent(h_bins + i + 1, hist -> GetBinContent(h_bins - i - 1) + 1); }   
    if ( i < h_bins ){ Hist -> SetBinContent(i+1, hist -> GetBinContent(i+1)); }
  }
 
  Normalizer(Deconv);
  Normalizer(Hist);
  
  int peak = Deconv -> GetMaximumBin() + 5;
  // ======= Perform a fit of the tail 

  RooRealVar x("x", "x", 1, d_bins -1);
  x.setRange("Tail", peak, h_bins - 400);

  RooRealVar s("s", "s", -50, 50);

  RooFormulaVar del("del", "del", "x-s", RooArgSet(x, s));

  RooDataHist De("De", "De", x, Deconv); // <---- PDF
  RooDataHist Hi("hi", "hi", x, Hist); // <---- we fit this to Deconv
    
  RooHistPdf model("model", "model", del, x, De, 1);
  RooFitResult* Result = model.fitTo(De, RooFit::Range("Tail"), RooFit::PrintLevel(0));
  
  // ===== Replace the Deconv with the hist ==== //
  float sum_d = Deconv -> Integral();
  float sum_h = Hist -> Integral();  
  int shift = s.getVal();
  
  for (int i(peak); i < d_bins; i++)
  {
    float h = Hist -> GetBinContent(i+1);
    float d = Deconv -> GetBinContent(i+1);
    Deconv -> SetBinContent(i+1 + shift, h*(sum_h/sum_d));
  }

  std::vector<float> Dec(d_bins, 0);
  for (int i(0); i < d_bins; i++)
  {
    Dec[i] = Deconv -> GetBinContent(i+1); 
  }

  delete Deconv;
  delete Hist;
  delete Result;
  return Dec;
}

// This is the closure test of the tail replace 
std::vector<float> Fit_Functions::TailReplaceClosure(TH1F* hist, std::vector<float> deconv)
{
   
  Functions F;
  Fit_Functions f;

  // ========= Constants and Inputs =================== //
  const float bins = hist -> GetNbinsX();
  const float v_size = deconv.size();
  const float offset = (v_size - bins)/bins;
 
  // ========= Add Padding to the histograms ========== //
  TH1F* Source = new TH1F("Source", "Source", deconv.size(), 0, 20);
  TH1F* Target = new TH1F("Target", "Target", deconv.size(), 0, 20);

  std::vector<float> Replacement = f.TH1FDataVector(hist, offset);
   
  F.VectorToTH1F(Replacement, Source);
  F.VectorToTH1F(deconv, Target);

  float TLumi = Target -> Integral();
  float SLumi = Source -> Integral();
 
  const int UpTo = std::round(2*(Target -> GetMaximumBin()));
  if (Target -> GetBinContent(UpTo) == deconv[0]) { return deconv; }

  for (int i(UpTo); i < v_size; i++) 
  {
    float t_e = Source -> GetBinContent(i + 1);
    float d_e = Target -> GetBinContent(i + 1);
    Target -> SetBinContent(i+1, t_e * (TLumi/SLumi) );    
  }

  std::vector<float> De = F.TH1FToVector(Target);

  delete Target;
  delete Source;
  return De;

}

std::vector<RooRealVar*> Fit_Functions::FitPDFtoData(std::vector<TH1F*> PDFs, TH1F* Data, float min, float max)
{
  Fit_Functions f;

  std::vector<TString> Names;
  for ( int i(0); i < PDFs.size(); i++)
  {
    Names.push_back(PDFs[i] -> GetName());
  }

  // Define required variables
  std::vector<float> Parms_S(Names.size(), 0);
  std::vector<float> Parms_E(Names.size(), Data -> Integral());

  // Variable Generation
  RooRealVar* x = new RooRealVar("x", "x", min, max);
  std::vector<RooRealVar*> vars = f.GenerateVariables(Names, Parms_S, Parms_E);
  std::vector<RooHistPdf*> pdfs = f.ConvertTH1FtoPDF(PDFs, x); 
  RooArgList var_List = f.VectorToArgList(vars); 
  RooArgList pdfs_List = f.VectorToArgList(pdfs);  
  RooDataHist* D = f.ConvertTH1toDataHist(Data, x);

  // Model 
  RooAddPdf model("model", "model", pdfs_List, var_List);
  model.fitTo(*D);

  delete x; 
  delete D;
  for (int i(0); i < pdfs.size(); i++)
  {
    delete pdfs[i];
  }
  return vars;
}

void Fit_Functions::Normalizer(TH1F* Hist)
{
  int nbins = Hist -> GetNbinsX();
  float sum(0);
  for ( int i(0); i < nbins; i++)
  {
    float e = Hist -> GetBinContent(i+1);
    sum += e; 
  }
  Hist -> Scale(1/sum);
}

std::vector<float> Fit_Functions::Fractionalizer(std::vector<RooRealVar*> vars, TH1F* Data)
{
  float lumi = Data -> Integral();
  float p1 = vars.at(0) -> getVal();
  float p2 = vars.at(1) -> getVal();
  float p3 = vars.at(2) -> getVal();
  float p4 = vars.at(3) -> getVal();  
  std::vector<float> frac = {p1/lumi, p2/lumi, p3/lumi, p4/lumi}; 
  return frac;
}

void Fit_Functions::Subtraction(std::vector<TH1F*> nTrk, TH1F* Target, int ntrk, std::vector<RooRealVar*> var)
{
  std::vector<float> frac = Fractionalizer(var, Target);
  float lumi = Target -> Integral();
   
  for (int i(0); i < nTrk.size(); i++)
  {
    Normalizer(nTrk.at(i));
    float f = frac[i]*lumi;
    if ( i != 1 )
    { 
      nTrk.at(i) -> Scale(f);
      Target -> Add(nTrk.at(i), -1); 
    }
    delete var[i];
    Normalizer(nTrk.at(i));
  }
  
  ArtifactRemove(Target); 
}

void Fit_Functions::ArtifactRemove(TH1F* Hist, TString mode)
{
  int index(0);
  if (mode == "b")
  { 
    index = ArtifactRemoveBackward(Hist); 
    for (int i(0); i < index; i++){ Hist -> SetBinContent(i+1, 0); }
  }
  else if (mode == "f")
  { 
    index = ArtifactRemoveForward(Hist);  
    for (int i(index); i < Hist -> GetNbinsX(); i++){ Hist -> SetBinContent(i + 1, 0); }      
  }
}

int Fit_Functions::ArtifactRemoveForward(TH1F* Hist)
{
  int max = Hist -> GetMaximumBin(); 
  float f = Hist -> GetBinContent(max);
  int index = 0; 
  for (int i(max); i < Hist -> GetNbinsX(); i++)
  {
    float e = Hist -> GetBinContent(i + 1);
    
    if (e < f) 
    { 
      f = e; 
      index++;
    }
    else {break;}
  }
  return max + index;
}

int Fit_Functions::ArtifactRemoveBackward(TH1F* Hist)
{
  int max = Hist -> GetMaximumBin(); 
  if (max < 2) 
  { 
    Hist -> SetBinContent(max, 0);
    max = Hist -> GetMaximumBin();
  }

  float f = Hist -> GetBinContent(max);
  int index = 0; 
  for (int i(0); i < max; i++)
  {
    float e = Hist -> GetBinContent(max - i - 1);
    
    if ( e < f ) 
    { 
      f = e; 
    }
    else { index++; }
  }
  return index;
}

std::vector<float> Fit_Functions::TH1FDataVector(TH1F* Data, float Offset)
{
  int bins = Data -> GetNbinsX();
  std::vector<float> Data_Vector(bins + bins*Offset, 0);
  for (int i(0); i < bins; i++)
  {
    Data_Vector[i] = Data -> GetBinContent(i+1);
    if (i < bins*Offset) { Data_Vector[i+bins] = Data -> GetBinContent(bins - i - 1); }
  }

  return Data_Vector;
}

void Fit_Functions::GaussianGenerator(float mean, float stdev, int N, TH1F* Hist)
{
  gRandom = new TRandom();

  for (int i(0); i < N; i++)
  {
    Hist -> Fill(gRandom -> Gaus(mean, stdev));
  }
  
  delete gRandom; 
}

std::vector<RooRealVar*> Fit_Functions::GaussianConvolutionFit(std::vector<TH1F*> PDFs, TH1F* trk2, float min, float max, float offset, float stdev_s, float stdev_e, float mean_s, float mean_e)
{
  // Namespace
  Functions F;
 
  // Constants in the code 
  const float Padding = max / 2;
  const int bins = trk2 -> GetNbinsX(); 
  const float ss = (max-min) / bins;
  const int Pad = std::round(std::abs(Padding) / ss);
  const float GaussianMean = 0; /// need to add later
  const float STDev = 0.1; // Need to add later
  const float Damp = 0.75;
  const int iteration = 50;

  // Define the static Gaussian PSF
  TH1F* PSF_HL = new TH1F("PSF_HL", "PSF_HL", bins + 2*Pad, min - Padding, max + Padding);   
  GaussianGenerator(GaussianMean, STDev, 500000, PSF_HL);
  Normalizer(PSF_HL); 
  std::vector<float> PSF_V = F.TH1FToVector(PSF_HL);
  delete PSF_HL;

  std::vector<TH1F*> PDF_HLs;
  for (TH1F* trk : PDFs)
  {
    // Configure the histograms 
    TString name = trk -> GetName(); name+=("_PDF");
    TH1F* PDF = new TH1F(name, name, bins + 2*Pad, min - Padding, max + Padding);
    
    // Fill the histograms with an offset 
    std::vector<float> temp = TH1FDataVector(trk, offset);
    F.VectorToTH1F(temp, PDF, Pad);

    // Normalize the histograms 
    Normalizer(PDF);
    
    // Convert to vector 
    std::vector<float> PDF_V = F.TH1FToVector(PDF);

    // Start the deconvolution 
    std::vector<float> deconv(PSF_V.size(), 0.5);
    for (int i(0); i < iteration; i++)
    {
      deconv = LRDeconvolution(PDF_V, PSF_V, deconv, Damp);  
    }
    
    // Write deconv to a histogram 
    PDF -> Reset(); 
    F.VectorToTH1F(deconv, PDF, Pad);
    Normalizer(PDF);
    PDF_HLs.push_back(PDF);
  }

  // Define the range of the dEdx
  RooRealVar* x = new RooRealVar("x", "x", min, max); 

  // Define the Gaussian Parameter: Mean
  std::vector<TString> Means_String = { "m1", "m2", "m3", "m4" };
  std::vector<float> Means_Begin(Means_String.size(), mean_s);
  std::vector<float> Means_End(Means_String.size(), mean_e);
  std::vector<RooRealVar*> Means = GenerateVariables(Means_String, Means_Begin, Means_End);

  // Define the Gaussian Parameter: Standard Deviation
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4" };
  std::vector<float> Stdev_Begin(Stdev_String.size(), stdev_s);
  std::vector<float> Stdev_End(Stdev_String.size(), stdev_e);
  std::vector<RooRealVar*> Stdev = GenerateVariables(Stdev_String, Stdev_Begin, Stdev_End);

  // Define the Gaussian Variables
  std::vector<TString> Gaus_String = { "g1", "g2", "g3", "g4"};
  std::vector<RooGaussian*> G_Vars = GaussianVariables(Gaus_String, Means, Stdev, x);
 
  // Import the PDFs as a RooDataHist
  std::vector<RooHistPdf*> PDF_Vars = ConvertTH1FtoPDF(PDF_HLs, x);
   
  // Define the ntrack coefficients:
  float Lumi = trk2 -> Integral();
  std::vector<TString> C_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4" };
  std::vector<float> C_Begin = { 0, 0, 0, 0 };
  std::vector<float> C_End = { Lumi, Lumi, Lumi, Lumi };
  std::vector<RooRealVar*> C_Vars = GenerateVariables(C_String, C_Begin, C_End);

  // Convolve the PDFs with the Gaussians
  std::vector<TString> Conv_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4" };
  std::vector<RooFFTConvPdf*> Conv_Vars = ConvolveVariables(Conv_String, PDF_Vars, G_Vars, x);
  
  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgList(*Conv_Vars[0], *Conv_Vars[1], *Conv_Vars[2], *Conv_Vars[3]), RooArgList(*C_Vars[0], *C_Vars[1], *C_Vars[2], *C_Vars[3]));   
  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk2_D = new RooDataHist("trk2_D", "trk2_D", *x, trk2); 
  model.fitTo(*trk2_D, Constrain(*Means[0]), Constrain(*Means[1]), Constrain(*Means[2]), Constrain(*Means[3]), Constrain(*Stdev[0]), Constrain(*Stdev[1]), Constrain(*Stdev[2]), Constrain(*Stdev[3]));

  for (unsigned int i(0); i < Means.size(); i++)
  {
    delete Means[i];
    delete Stdev[i];
    delete G_Vars[i];
    delete PDF_Vars[i];
    delete Conv_Vars[i];
    delete PDF_HLs[i];
  }

  delete trk2_D; 
  return C_Vars; 
}

// ============= Benchmarking Class ========== //
float Benchmark::WeightedEuclidean(std::vector<float> v1, std::vector<float> v2)
{
  // First check that the vectors have the same length
  if ( v1.size() != v2.size() ) 
  { 
    std::cout << v1.size() << v2.size() << std::endl;
    std::cout << "Invalid data vectors passed!!!!" << std::endl; 
    return 0;
  }
  float dist(0);
  float normV1(0);
  float normV2(0);
  
  // Calculating the sum over the vector entries
  for ( int i(0); i < v1.size(); i++)
  {
    normV1 += v1.at(i);
    normV2 += v2.at(i);
  }

  normV1 = 1/normV1;
  normV2 = 1/normV2;

  float sum(0); 
  for (int i(0); i < v1.size(); i++)
  {
    sum += pow( v1.at(i)*normV1 - v2.at(i)*normV2, 2); 
  }
  return sum;
}

float Benchmark::PythagoreanDistance(std::vector<float> v1, std::vector<float> v2)
{
  float di(0);
  for ( int i(0); i < v1.size(); i++)
  {
    di += pow(v1.at(i) - v2.at(i), 2);
  }
  return pow(di, 0.5);

}
