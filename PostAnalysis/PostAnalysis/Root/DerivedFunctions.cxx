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
  RooFitResult* stat = model.fitTo(*data, RooFit::Save());
  std::cout << "#########" << stat -> status() << std::endl; 
  
   
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
  std::vector<float> Begin(Hists.size(), 0.);
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

// Likely useless function. Need to see if I will remove this 
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
  B.Normalize(PDFs);
  RemoveArtifact(PDFs[0]);

  // === TRK2
  B.ConvolveHists(PDFs[0], PDFs[0], PDFs[1]); 
  B.Normalize(PDFs[1]);
  B.ResidualRemove(PDFs[1]); 
  
  // === TRK3
  B.ConvolveHists(PDFs[1], PDFs[0], PDFs[2]);  
  B.Normalize(PDFs[2]);
  B.ResidualRemove(PDFs[2]); 
 
  // === TRK4 
  B.ConvolveHists(PDFs[1], PDFs[1], PDFs[3]);  
  B.Normalize(PDFs[3]);
  B.ResidualRemove(PDFs[3]); 
 
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

    delete Temp;
  }
  
  // Convert the histograms back to the original format 
  std::vector<TH1F*> DeconvPDFs;
  for (int i(0); i < PDFs.size(); i++)
  {
    std::vector<float> Data = B.TH1FDataVector(Hists[i]);  
  
    TString name = PDFs[i] -> GetTitle(); name += ("_"); name += (i+1);
    TH1F* H = new TH1F(name, name, bins_pdf, 0, bins_pdf);  
    H -> Reset();
    B.ToTH1F(Data, H);
    DeconvPDFs.push_back(H); 
  }
  TH1F* Data_G = new TH1F("Data_Gaus_Fit", "Data_Gaus_Fit", bins_pdf, 0, bins_pdf);
  B.ShiftExpandTH1F(GxTrk, Data_G, 0);
  B.Normalize(Data_G);
 
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
  float Lumi = Data_G -> Integral();
  std::vector<TString> C_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4" };
  std::vector<float> C_Begin = { 0., 0., 0., 0. };
  std::vector<float> C_End = { Lumi, Lumi, Lumi, Lumi };
  std::vector<RooRealVar*> C_Vars = B.RooVariables(C_String, C_Begin, C_End);

  // Convolve the PDFs with the Gaussians
  std::vector<TString> Conv_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4" };
  std::vector<RooFFTConvPdf*> Conv_Vars = B.RooVariables(Conv_String, PDF_Vars, G_Vars, x);
  
  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgList(*Conv_Vars[0], *Conv_Vars[1], *Conv_Vars[2], *Conv_Vars[3]), RooArgList(*C_Vars[0], *C_Vars[1], *C_Vars[2], *C_Vars[3]));   

  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk2_D = new RooDataHist("trk2_D", "trk2_D", *x, Data_G); 

  model.fitTo(*trk2_D, RooFit::Constrain(*Means[0]), 
                                            RooFit::Constrain(*Means[1]), 
                                            RooFit::Constrain(*Means[2]), 
                                            RooFit::Constrain(*Means[3]), 
                                            RooFit::Constrain(*Stdev[0]), 
                                            RooFit::Constrain(*Stdev[1]), 
                                            RooFit::Constrain(*Stdev[2]), 
                                            RooFit::Constrain(*Stdev[3]));
  
  RooFitResult* stat = model.fitTo(*trk2_D, RooFit::Save(), RooFit::SumW2Error(true), RooFit::NumCPU(4)); 
  
  
   
  // Plotting 
  //RooPlot* xframe = x -> frame(RooFit::Title("Hello"));
  //trk2_D -> plotOn(xframe, RooFit::Name("Data"), RooFit::DataError(RooAbsData::None), RooFit::XErrorSize(0));
  //model.plotOn(xframe, RooFit::Name("1trk"), RooFit::Components(*Conv_Vars[0]), RooFit::LineStyle(kDotted), RooFit::LineColor(kBlue));
  //model.plotOn(xframe, RooFit::Name("2trk"), RooFit::Components(*Conv_Vars[1]), RooFit::LineStyle(kDotted), RooFit::LineColor(kCyan)); 
  //model.plotOn(xframe, RooFit::Name("3trk"), RooFit::Components(*Conv_Vars[2]), RooFit::LineStyle(kDotted), RooFit::LineColor(kOrange));
  //model.plotOn(xframe, RooFit::Name("4trk"), RooFit::Components(*Conv_Vars[3]), RooFit::LineStyle(kDotted), RooFit::LineColor(kGreen));

  //TCanvas* z = new TCanvas();
  //gPad -> SetLogy();
  //xframe -> SetMinimum(1e-5);
  //xframe -> Draw();
  //z -> Update();
 
  std::map<TString, float> Out;

  if (!stat){Out["Status"] = -1;}
  else{Out["Status"] = stat -> status();}

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
    delete DeconvPDFs[i];
  }
  delete trk2_D; 
  delete x;
  delete Data_G;
  delete stat;

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

  RooRealVar* x = new RooRealVar("x", "x", 0, trk1 -> GetNbinsX());
  
  TCanvas* can = new TCanvas("Outside", "Outside");

  TH1F* trk2_L = new TH1F("trk2_L", "trk2_L", trk2 -> GetNbinsX(), 0, trk2 -> GetNbinsX());
  TH1F* trk1_L = new TH1F("trk1_L", "trk1_L", trk1 -> GetNbinsX(), 0, trk1 -> GetNbinsX());
  B.ShiftExpandTH1F(trk2, trk2_L); 
  B.ShiftExpandTH1F(trk1, trk1_L); 

  for (int x(0); x < cor_loop; x++)
  { 
    TH1F* trk2_Copy = (TH1F*)trk2_L -> Clone("trk2_Copy");
    for (int v(0); v < 5; v++)
    { 
      PDFs = nTRKGenerator(trk1, trk2_Copy, offset, iter);
      std::map<TString, float> Parameters = FitGaussian(trk2_Copy, PDFs, mean, stdev, m_s, m_e, s_s, s_e, offset, iter);  
  
      for (int i(0); i < Names.size(); i++)
      {
        TH1F* H = GaussianConvolve(PDFs[i], Parameters[Mean[i]], Parameters[Stdev[i]]);
        PDFs[i] -> Reset();
        B.ShiftExpandTH1F(H, PDFs[i]);
        delete H; 
     
        float Norm(0); 
        for (TString N : Names)
        { 
          float e = Parameters[N]; 
          Norm = Norm + e; 
        }
        
          
        float z = Parameters[Names[i]];  
        if (Parameters["Status"] == 0 || Parameters["Status"] == 4)
        {
          PDFs[i] -> Scale(0.05*z);
        }
      }
      if (Parameters["Status"] == 0)
      {
        trk2_Copy -> Add(PDFs[0], -1);
        trk2_Copy -> Add(PDFs[2], -1);
        trk2_Copy -> Add(PDFs[3], -1);
      }
      else 
      {
        trk2_Copy -> Reset();
        trk2_Copy -> Add(trk2_L);
      }
      
      std::cout << "###################" << v << std::endl; 
      std::cout << "###################" << Parameters["Status"] << std::endl; 
     
      if( v < 4)
      {
        for (int i(0); i < PDFs.size(); i++){delete PDFs[i];}
      }
    }
    
    B.Normalize(PDFs);  
    std::vector<RooRealVar*> vars_trk1 = FitToData(PDFs, trk1_L, 0, trk1_L -> GetNbinsX()); 
    float t1_t1 = vars_trk1[0] -> getVal(); 
    float t2_t1 = vars_trk1[1] -> getVal();
    float t3_t1 = vars_trk1[2] -> getVal(); 
    float t4_t1 = vars_trk1[3] -> getVal();
  
    PDFs[0] -> Scale(0.25*(t1_t1)); 
    PDFs[1] -> Scale(0.25*(t2_t1));
    PDFs[2] -> Scale(0.25*(t3_t1));
    PDFs[3] -> Scale(0.25*(t4_t1)); 
    
    trk1_L -> Add(PDFs[1], -1);
    trk1_L -> Add(PDFs[2], -1);
    trk1_L -> Add(PDFs[3], -1); 
    
    B.Normalize(PDFs);
    std::vector<RooRealVar*> vars = FitToData(PDFs, trk2_L, 0, trk2_L -> GetNbinsX());
    float t1 = vars[0] -> getVal(); 
    float t2 = vars[1] -> getVal();
    float t3 = vars[2] -> getVal(); 
    float t4 = vars[3] -> getVal();
  
    PDFs[0] -> Scale(Gamma*(t1)); 
    PDFs[1] -> Scale(Gamma*(t2));
    PDFs[2] -> Scale(Gamma*(t3));
    PDFs[3] -> Scale(Gamma*(t4)); 
    
    trk2_L -> Add(PDFs[0], -1);
    trk2_L -> Add(PDFs[2], -1);
    trk2_L -> Add(PDFs[3], -1);
   
    if (x > 2){B.ResidualRemove(trk2_L);} 
     
    P.PlotHists(PDFs, trk2_L, can);  
    can -> Update(); 



    //iter = iter + 25; 
  
    std::cout << "################### " << x << std::endl;

    delete trk2_Copy;
    if (x == cor_loop){return PDFs;}     
    for (int i(0); i < PDFs.size(); i++){delete PDFs[i];} 
    std::vector<float> scales = B.Ratio(vars, trk2_L);  
    std::vector<float> scales1 = B.Ratio(vars_trk1, trk1_L); 
  }
  return PDFs;
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
