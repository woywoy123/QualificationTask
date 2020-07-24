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



// Work in progress... 
void Fit_Functions::ConvolveHists(TH1F* Hist1, TH1F* Hist2, TH1F* conv, float min, float max)
{
  RooRealVar* x = new RooRealVar("x", "x", min, max);
  RooDataHist* HistD1 = ConvertTH1toDataHist(Hist1, x);
  RooDataHist* HistD2 = ConvertTH1toDataHist(Hist2, x);

  RooHistPdf* HistPDF1 = ConvertTH1FtoPDF(HistD1, "Hist1", x);
  RooHistPdf* HistPDF2 = ConvertTH1FtoPDF(HistD2, "Hist2", x);

  RooNumConvPdf* H12 = new RooNumConvPdf("ConvH1H2", "ConvH1H2", *x, *HistPDF1, *HistPDF2);
  TF1* tf1 = new TF1(*H12 -> asTF(RooArgList(*x)));
  conv -> GetListOfFunctions() -> Add(tf1); 

  // Create the canvas  
  auto f = new TCanvas();
  gPad -> SetLogy();
  gStyle -> SetOptStat(0);
 
  // Initialize the frame for Plotting 
  RooPlot* xframe = x -> frame(RooFit::Title("Hello")); 
  H12 -> plotOn(xframe, RooFit::Name("Measured"));
  
  xframe -> SetMinimum(1);
  xframe -> Draw();




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

