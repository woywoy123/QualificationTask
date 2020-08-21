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


void DerivedFunctions::RemoveArtifact(TH1F* Conv)
{
  int max = Conv -> GetMaximumBin();
  float T = 0; 
  int iter = 0;  
  for (int i(0); i < max; i++)
  {
    float e = Conv -> GetBinContent(i+1); 
    float f = Conv -> GetBinContent(i+2);
    float delta = f/e; 
    if (delta > T)
    {
      T = delta;
      iter = i; 
    }
  }

  for (int i(0); i < iter; i++)
  {
    float e = Conv -> GetBinContent(iter+1);
    float f = Conv -> GetBinContent(i+1); 
    Conv -> SetBinContent(i+1, std::pow(f, T*(iter - i))); 
  }
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
  
  for (TH1F* H : PDFs)
  {
    RemoveArtifact(H); 
  }

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

std::map<TString, float> DerivedFunctions::FitGaussian(TH1F* GxTrk, std::vector<TH1F*> PDFs, float mean, float stdev, float m_s, float m_e, float s_s, float s_e, float offset, int iter)
{
  DistributionGenerators DG; 
  BaseFunctions B; 
  
  // Define Padding 
  int bins_pdf = PDFs[0] -> GetNbinsX(); 

  // Create temorary histograms 
  std::vector<TString> names; 
  for (int i(0); i < PDFs.size(); i++)
  {
    TString name = PDFs[i] -> GetTitle(); name += (i+1);  
    names.push_back(name); 
  } 

  std::vector<TH1F*> Hists = B.MakeTH1F(names, GxTrk);   
  for (int i(0); i < PDFs.size(); i++)
  {
    // Import and extend the TH1F* data
    std::vector<float> PDF_V = B.TH1FDataVector(PDFs[i], offset);  
   
    // Rescale the data to have the same length has the GxTrk 
    TString v = "Temp"; v+=(i+1);
    TH1F* Temp = new TH1F(v, v, PDF_V.size(), 0, PDF_V.size());
    B.ToTH1F(PDF_V, Temp);
    B.ShiftExpandTH1F(Temp, Hists[i], bins_pdf*offset/2); 
    delete Temp;

    // Convert into a vector for the deconv 
    PDF_V = B.TH1FDataVector(Hists[i], 0);
    Hists[i] -> Reset();

    // Define the prior of the same length as the input data
    std::vector<float> deconv(PDF_V.size(), 0.5);    

    // Create A gaussian distribution    
    DG.Gaussian(mean, stdev, Constants::GaussianToys, Hists[i]);
    B.Normalize(Hists[i]); 

    // Convert the PSF histogram into a vector 
    std::vector<float> PSF_V = B.TH1FDataVector(Hists[i], 0);
    Hists[i] -> Reset(); 
   
    for (int x(0); x < iter; x++) 
    {
      deconv = B.LucyRichardson(PDF_V, PSF_V, deconv);
    }
    
    B.ToTH1F(deconv, Hists[i]);   
  }
  
  // Convert the histograms back to the original format 
  std::vector<TH1F*> DeconvPDFs;
  for (int i(0); i < PDFs.size(); i++)
  {
    std::vector<float> Data = B.TH1FDataVector(Hists[i]);  
    delete Hists[i];
  
    TString name = PDFs[i] -> GetTitle(); name += (i+1);
    TH1F* H = new TH1F(name, name, bins_pdf, 0, bins_pdf);  
    H -> Reset();
    B.ToTH1F(Data, H);
    DeconvPDFs.push_back(H); 
  }
  TH1F* Data = new TH1F("Data", "Data", bins_pdf, 0, bins_pdf);
  B.ShiftExpandTH1F(GxTrk, Data, 0);
 
  // Define the range of the dEdx
  RooRealVar* x = new RooRealVar("x", "x", 0, bins_pdf); 

  // Define the Gaussian Parameter: Mean
  std::vector<TString> Means_String = { "m1", "m2", "m3", "m4" };
  std::vector<float> Means_Begin(Means_String.size(), m_s);
  std::vector<float> Means_End(Means_String.size(), m_e);
  std::vector<RooRealVar*> Means = B.RooVariables(Means_String, Means_Begin, Means_End);

  // Define the Gaussian Parameter: Standard Deviation
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4" };
  std::vector<float> Stdev_Begin(Stdev_String.size(), s_s);
  std::vector<float> Stdev_End(Stdev_String.size(), s_e);
  std::vector<RooRealVar*> Stdev = B.RooVariables(Stdev_String, Stdev_Begin, Stdev_End);

  // Define the Gaussian Variables
  std::vector<TString> Gaus_String = { "g1", "g2", "g3", "g4"};
  std::vector<RooGaussian*> G_Vars = B.RooVariables(Gaus_String, Means, Stdev, x);
 
  // Import the PDFs as a RooDataHist
  std::vector<RooHistPdf*> PDF_Vars = B.RooPDF(DeconvPDFs, x);
   
  // Define the ntrack coefficients:
  float Lumi = GxTrk -> Integral();
  std::vector<TString> C_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4" };
  std::vector<float> C_Begin = { -1, -1, -1, -1 };
  std::vector<float> C_End = { Lumi, Lumi, Lumi, Lumi };
  std::vector<RooRealVar*> C_Vars = B.RooVariables(C_String, C_Begin, C_End);

  // Convolve the PDFs with the Gaussians
  std::vector<TString> Conv_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4" };
  std::vector<RooFFTConvPdf*> Conv_Vars = B.RooVariables(Conv_String, PDF_Vars, G_Vars, x);
  
  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgList(*Conv_Vars[0], *Conv_Vars[1], *Conv_Vars[2], *Conv_Vars[3]), RooArgList(*C_Vars[0], *C_Vars[1], *C_Vars[2], *C_Vars[3]));   

  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk2_D = new RooDataHist("trk2_D", "trk2_D", *x, Data); 
  model.fitTo(*trk2_D, RooFit::Constrain(*Means[0]), 
                       RooFit::Constrain(*Means[1]), 
                       RooFit::Constrain(*Means[2]), 
                       RooFit::Constrain(*Means[3]), 
                       RooFit::Constrain(*Stdev[0]), 
                       RooFit::Constrain(*Stdev[1]), 
                       RooFit::Constrain(*Stdev[2]), 
                       RooFit::Constrain(*Stdev[3]));
 
  std::map<TString, float> Out;
  for (int i(0); i < Stdev.size(); i++)
  {
    Out[Means_String[i]] = Means[i] -> getVal();
    Out[Stdev_String[i]] = Stdev[i] -> getVal();
    Out[C_String[i]] = C_Vars[i] -> getVal();

    delete Means[i]; 
    delete Stdev[i]; 
    delete C_Vars[i]; 
    delete G_Vars[i]; 
    delete Conv_Vars[i]; 
    delete Hists[i];
    delete PDF_Vars[i];
  }
  delete trk2_D; 
  delete x;
  delete Data;

  return Out;   
}

std::vector<TH1F*> DerivedFunctions::MainAlgorithm(TH1F* trk1, TH1F* trk2, std::vector<float> Params, float offset, float Gamma, int iter, int cor_loop)
{
  BaseFunctions B;
  Plotting P; 

  float mean = Params[0]; 
  float stdev = Params[1]; 
  float m_s = Params[2]; 
  float m_e = Params[3]; 
  float s_s = Params[4];
  float s_e = Params[5];
  std::vector<TH1F*> PDFs;
  const std::vector<TString> Names = {"n_trk1", "n_trk2", "n_trk3", "n_trk4"}; 
  const std::vector<TString> Stdev = {"s1", "s2", "s3", "s4"}; 
  const std::vector<TString> Mean = {"m1", "m2", "m3", "m4"};

  for (int x(0); x < cor_loop; x++)
  { 
   
    PDFs = nTRKGenerator(trk1, trk2, offset, iter);
    std::map<TString, float> Parameters = FitGaussian(trk2, PDFs, mean, stdev, m_s, m_e, s_s, s_e, offset, iter);  
    
    for (int i(0); i < Names.size(); i++)
    {
      TH1F* H = GaussianConvolve(PDFs[i], Parameters[Mean[i]], Parameters[Stdev[i]]);
      PDFs[i] -> Reset();
      B.ShiftExpandTH1F(H, PDFs[i]);
      delete H; 
      float e = Parameters[Names[i]];  
      PDFs[i] -> Scale(Gamma*e);   
    }

    trk2 -> Add(PDFs[0], -1);
    trk2 -> Add(PDFs[2], -1);
    trk2 -> Add(PDFs[3], -1);

    if (x == cor_loop){return PDFs;}     
    for (int i(0); i < PDFs.size(); i++){delete PDFs[i];}
  
  }
  return PDFs;
}

// Note to self: 
// - Try doing a subsequent data fit to the post GaussianConvolve

















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
