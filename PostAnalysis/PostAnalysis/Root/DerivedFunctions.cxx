#include<PostAnalysis/DerivedFunctions.h>

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, 
                                                     TH1F* Data, 
                                                     RooRealVar* Domain, 
                                                     std::vector<float> Begin, 
                                                     std::vector<float> End, 
                                                     std::vector<TString> Names)
{
  BaseFunctions B;
  B.Normalize(Hists);
  std::vector<RooRealVar*> Var = B.RooVariables(Names, Begin, End);
  std::vector<RooHistPdf*> Hist = B.RooPDF(Hists, Domain);
  RooDataHist* data = B.RooData(Data, Domain);

  RooAddPdf model("Model", "Model", B.RooList(Hist), B.RooList(Var));
  model.fitTo(*data);
  
  for (int i(0); i < Var.size(); i++)
  {
    delete Hist[i];
  }
  delete data;
  
  return Var;
} 

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, TH1F* Data, float min, float max)
{
  RooRealVar* x = new RooRealVar("x", "x", min, max);
  std::vector<float> Begin(Hists.size(), 0);
  std::vector<float> End(Hists.size(), Data -> Integral());
  std::vector<RooRealVar*> Vars = FitToData(Hists, Data, x, Begin, End, Constants::FitNames);
  return Vars; 
}

int DerivedFunctions::NumericalShift(TH1F* H1, TH1F* H2)
{
  BaseFunctions B; 

  int bin1 = H1 -> GetNbinsX();
  int bin2 = H2 -> GetNbinsX();
  int delta = bin1 - bin2;
  TH1F* Hist;
  TH1F* Clone;
  int max;
  
  // Means bin1 is bigger  
  if (delta > 0)
  {
    max = bin1;
    Hist = (TH1F*)H1 -> Clone("Hist");
    Hist -> Reset(); 
    B.ShiftExpandTH1F(H2, Hist, 0);
    Clone = (TH1F*)H1 -> Clone("Clone");
  }
  else
  {
    max = bin2;
    Hist = (TH1F*)H2 -> Clone("Hist"); 
    Hist -> Reset();
    B.ShiftExpandTH1F(H1, Hist, 0);
    Clone = (TH1F*)H2 -> Clone("Clone");
  }

  B.Normalize({Hist, Clone});
  TH1F* Temp = (TH1F*)Clone -> Clone("Temp");
  Temp -> Reset(); 
  
  float T = 100; 
  int shifted;
  for (int i(-max); i < max; i++)
  {
    B.ShiftExpandTH1F(Hist, Temp, i);
   
    float e(0); 
    for (int x(0); x < Temp -> GetNbinsX(); x++)
    {
      e += std::abs(Temp -> GetBinContent(x+1) - Clone -> GetBinContent(x+1));
    }
    
    if (T > e)
    {
      T = e; 
      shifted = i+1;
    }
  } 
  if (delta < 0){return -shifted;}

  delete Temp; 
  delete Clone; 
  delete Hist;  
  return shifted;
}

void DerivedFunctions::ReplaceShiftTail(TH1F* Source, TH1F* Target)
{
  BaseFunctions B;

  int bin_s = Target -> GetNbinsX(); 
  float Lumi = Target -> Integral();
  float Lumi_S = Source -> Integral();
  B.Normalize({Source, Target});

  // Extend the source and target hist to avoid nan 
  std::vector<float> Source_V = B.TH1FDataVector(Source, 0.1);
  std::vector<float> Target_V = B.TH1FDataVector(Target, 0.1);  
  TH1F* Source_E = new TH1F("Source_E", "Source_E", Source_V.size(), 0, Source_V.size());
  TH1F* Target_E = new TH1F("Target_E", "Target_E", Target_V.size(), 0, Target_V.size()); 
  B.ToTH1F(Source_V, Source_E); 
  B.ToTH1F(Target_V, Target_E);

  // Define the tail to be replaced and move away from peak 
  int max_bin = Target -> GetMaximumBin() + bin_s*0.05;
  int max_bin_s = Source -> GetMaximumBin();
  TH1F* Temp = (TH1F*)Source -> Clone("Temp"); 
  
  for (int i(0); i < max_bin_s; i++)
  {
    Temp -> SetBinContent(i+1, 0);
  }
  int shift = NumericalShift(Temp, Target); 
   
  float t_S = Source_E -> GetBinContent(max_bin + 1 - shift -1);
  float t_T = Target_E -> GetBinContent(max_bin + 1);
  for (int i(max_bin); i < bin_s; i++) 
  {
    float S_e = Source_E -> GetBinContent(i + 1 - shift -1);
    float T_e = Target_E -> GetBinContent(i + 1); 
    Target -> SetBinContent(i+1, S_e*(t_T/t_S));
  }
  

  Target -> Scale(Lumi);
  Source -> Scale(Lumi_S);
  delete Temp; 
  delete Source_E;
  delete Target_E;
}

std::vector<TH1F*> DerivedFunctions::nTRKGenerator(TH1F* trk1, TH1F* trk2, float offset, int iter)
{
  BaseFunctions B;
  std::vector<float> trk2_V = B.TH1FDataVector(trk2, offset);
  std::vector<float> deconv(trk2_V.size(), 0.5);

  // Generate the histograms as output
  std::vector<TString> Names = {"TRK1", "TRK2", "TRK3", "TRK4"};  
  std::vector<TH1F*> PDFs = B.MakeTH1F(Names, trk2);
 
  // Deconvolve the 2trk data into 1 trk  
  for (int i(0); i < iter; i++)
  {
    deconv = B.LucyRichardson(trk2_V, deconv, deconv); 
  }
  
  // Replace the tail of the deconv vector 
  B.ToTH1F(deconv, PDFs[0]);

  // === TRK1
  ReplaceShiftTail(trk1, PDFs[0]);
  
  // === TRK2
  B.ConvolveHists(PDFs[0], PDFs[0], PDFs[1]);  
  
  // === TRK3
  B.ConvolveHists(PDFs[1], PDFs[0], PDFs[2]);  
 
  // === TRK4 
  B.ConvolveHists(PDFs[1], PDFs[1], PDFs[3]);  

  B.Normalize(PDFs);
  return PDFs;
}

TH1F* DerivedFunctions::GaussianConvolve(TH1F* Hist, float mean, float stdev, int Toys)
{
  DistributionGenerators DG; 
  BaseFunctions B;

  // Define a padding
  int bins = Hist -> GetNbinsX();
  int Padding = bins/2;
  TString name = Hist -> GetTitle();
  TH1F* Hist_L = new TH1F("Hist_L", "Hist_L", 2*bins, -Padding, bins + Padding);
  TH1F* Gaus_L = new TH1F("Gaus_L", "Gaus_L", 2*bins, -Padding, bins + Padding);
  TH1F* GxH_L = new TH1F("GxH_L" + name, "GxH_L" + name, 2*bins, -Padding, bins + Padding);

  // Convert the Hist to a longer version for deconvolution 
  std::vector<float> Hist_V = B.TH1FDataVector(Hist, 0);
  B.ShiftExpandTH1F(Hist, Hist_L, -Padding);

  // Generate gaussian distribution 
  DG.Gaussian(mean, stdev, Toys, Gaus_L);
  B.Normalize(Gaus_L);

  // Convolve the Hist with Gaussian 
  B.ConvolveHists(Gaus_L, Hist_L, GxH_L);
  B.Normalize(GxH_L);

  delete Gaus_L; 
  delete Hist_L; 
  return GxH_L; 
}

void DerivedFunctions::FitGaussian(std::vector<TH1F*> GxTrk, std::vector<TH1F*> PDFs, float mean, float stdev, float m_s, float m_e, float s_s, float s_e, float offset, int iter)
{
  DistributionGenerators DG; 
  BaseFunctions B; 
  Plotting P; 
  
  // Define Padding 
  int bins_pdf = PDFs[0] -> GetNbinsX(); 
  int bins_gxt = GxTrk[0] -> GetNbinsX(); 
  int Padding = (bins_gxt - bins_pdf)/2; 

  std::vector<TString> names; 
  for (int i(0); i < PDFs.size(); i++)
  {
    TString name = PDFs[i] -> GetTitle(); name += (i+1);  
    names.push_back(name); 
  } 
  std::vector<TH1F*> Deconv = B.MakeTH1F(names, PDFs[0]);  
  
  for (int i(0); i < PDFs.size(); i++)
  {
    TH1F* PDF = PDFs[i];  
    std::vector<float> PDF_V = B.TH1FDataVector(PDF, offset);  
    std::vector<float> deconv(PDF_V.size(), 0.5);
    
    TH1F* PSF = new TH1F("Gaus", "Gaus", deconv.size(), -deconv.size()/2, deconv.size()); 
    DG.Gaussian(mean, stdev, Constants::GaussianToys, PSF);
    std::vector<float> PSF_V = B.TH1FDataVector(PSF, 0);
    delete PSF; 
  
  
    for (int x(0); x < iter; x++) 
    {
      deconv = B.LucyRichardson(PDF_V, PSF_V, deconv);
    }
    
    B.ToTH1F(deconv, Deconv[i]);   
  }

  P.PlotHists(Deconv, Deconv); 
  
  // Need to fix up the deconvolution. Seems to be shifted......  
  
}




















// Doesnt work....
void DerivedFunctions::RooShift(TH1F* H1, TH1F* H2)
{
  BaseFunctions B; 

  int bin1 = H1 -> GetNbinsX();
  int bin2 = H2 -> GetNbinsX();
  int delta = bin1 - bin2;
  TH1F* Hist;
  TH1F* Clone;
  int max;
  
  // Means bin1 is bigger  
  if (delta > 0)
  {
    max = bin1;
    Hist = (TH1F*)H1 -> Clone("Hist");
    Hist -> Reset(); 
    B.ShiftExpandTH1F(H2, Hist, 0);
    Clone = (TH1F*)H1 -> Clone("Clone");
  }
  else
  {
    max = bin2;
    Hist = (TH1F*)H2 -> Clone("Hist"); 
    Hist -> Reset();
    B.ShiftExpandTH1F(H1, Hist, 0);
    Clone = (TH1F*)H2 -> Clone("Clone");
  }
  
 // Plotting P;
 // TCanvas* can = P.PlotHists(Clone, Hist);

  // Normalize the two hists 
  //B.Normalize({Hist, Clone});

  int peak = Clone -> GetMaximumBin();

  // ======= RooFit Section ====== //
  RooRealVar x("x", "x", -250, 500);
  RooRealVar s("s", "s", 0, -5, 5);

  RooFormulaVar del("del", "del", "x-s", RooArgSet(x, s));

  RooDataHist H("De", "De", x, Hist); 
  RooDataHist C("hi", "hi", x, Clone);
    
  RooHistPdf model("model", "model", del, x, H, 1);
  model.fitTo(C);

  Plotting P;
  TCanvas* can = P.PlotHists(model, x, C); 
  can -> Draw();
}
