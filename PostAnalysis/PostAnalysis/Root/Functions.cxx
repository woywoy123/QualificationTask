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

void Functions::VectorToTH1F(std::vector<float> Vec, TH1F* hist)
{ 
  for (int i(0); i < hist -> GetNbinsX(); i++)
  {
    hist -> SetBinContent(i+1, Vec.at(i));  
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

RooRealVar* Fit_Functions::GenerateVariable(TString name, double begin, double end)
{
  RooRealVar* var = new RooRealVar(name, name, begin, 0., end);
  return var;
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

std::vector<RooRealVar*> Fit_Functions::GenerateVariables(std::vector<TString> Names, std::vector<double> Begin, std::vector<double> End)
{
  std::vector<RooRealVar*> Variables;
  for ( unsigned int x = 0; x < Names.size(); x++)
  {
    RooRealVar* Var = GenerateVariable(Names.at(x), Begin.at(x), End.at(x));
    Variables.push_back(Var);
  }
  return Variables;
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

void Fit_Functions::ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv, int offset)
{
  int nBins_2 = Hist2 -> GetNbinsX();
  int nBins_1 = Hist1 -> GetNbinsX();
   
  // Make sure the bins are equally sized  
  if ( nBins_2 != nBins_1 ) { return; }

  // Convert TH1 to vector 
  std::vector<float> H1, H2, Con;
  for ( unsigned int i(0); i < nBins_2 - offset; i++ )
  {
    H1.push_back(Hist1 -> GetBinContent(i+1));
    H2.push_back(Hist1 -> GetBinContent(i+1));
    Con.push_back(Hist1 -> GetBinContent(i+1));
  }
 
  // Set bin content of conv histogram 
  std::vector<float> Conv = ConvolveHists(H1, H2);
  std::cout << Conv.size() << std::endl;
  for ( unsigned int i(0); i < nBins_2 - offset; i++ )
  {
    conv -> SetBinContent(i+1, Conv.at(i));
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

// ================== Plotting Class ================================= //
TCanvas* Plot_Functions::GeneratePlot(TString Title, RooRealVar* range, RooDataHist* Data, RooAddPdf Model, std::vector<RooHistPdf*> PDFs, std::vector<TString> pdf_titles)
{
  // Initialize the frame for Plotting 
  RooPlot* xframe = range -> frame(RooFit::Title(Title)); 
  Data -> plotOn(xframe, RooFit::Name("Measured"));
  Model.plotOn(xframe, RooFit::Name("Model Fit"));
  
  // Loop over the histograms 
  for (unsigned i = 0; i < PDFs.size(); i++) 
  {
    Model.plotOn(xframe, RooFit::Name(pdf_titles[i]), RooFit::Components(*PDFs[i]), RooFit::LineStyle(kDotted), RooFit::LineColor(Constants::Colors[i]));
  }
 
  // Create the canvas  
  auto f = new TCanvas();
  gPad -> SetLogy();
  gStyle -> SetOptStat(0);
  
  xframe -> SetXTitle("dEdx [MeV g^{-1} cm^{2}]");
  xframe -> SetMinimum(1);
  xframe -> Draw();

  // Create the legend 
  TLegend* Legend = new TLegend(0.9, 0.9, 0.75, 0.75);
  Legend -> AddEntry("Measured", "Measured Distribution"); 
  Legend -> AddEntry("Model Fit", "Model Fit");
  for ( TString name : pdf_titles )
  {
    Legend -> AddEntry(name, name);
  }
  Legend -> Draw();

  return f;
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

