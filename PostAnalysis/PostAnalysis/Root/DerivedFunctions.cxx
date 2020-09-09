#include<PostAnalysis/DerivedFunctions.h>

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, 
                                                     TH1F* Data, 
                                                     RooRealVar* Domain, 
                                                     std::vector<float> Begin, 
                                                     std::vector<float> End, 
                                                     std::vector<TString> Names)
{
  BaseFunctions B;
  std::vector<RooRealVar*> Var = B.RooVariables(Names, Begin, End);
  std::vector<RooHistPdf*> Hist = B.RooPDF(Hists, Domain);
  RooDataHist* data = B.RooData(Data, Domain);

  RooAddPdf model("Model", "Model", B.RooList(Hist), B.RooList(Var));
  RooFitResult* stat = model.fitTo(*data, RooFit::Save());
   
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

void DerivedFunctions::SafeScale(std::vector<TH1F*> PDFs, TH1F* Data)
{
  
  std::vector<RooRealVar*> var = FitToData(PDFs, Data, 0, Data -> GetNbinsX()); 
  for (int i(0); i < var.size(); i++)
  {
    PDFs[i] -> Scale(var[i] -> getVal()); 
    delete var[i];
  }

  int bins = Data -> GetNbinsX(); 
  for (int i(0); i < bins; i++)
  {
    float sum_pdfs = 0; 
    for (int x(0); x < PDFs.size(); x++)
    {
      sum_pdfs += PDFs[x] -> GetBinContent(i+1); 
    }
    
    float data_e = Data -> GetBinContent(i+1); 
    float ratio = data_e/sum_pdfs;
  
    if (sum_pdfs == 0) {ratio = 0;}
    for (int x(0); x < PDFs.size(); x++)
    {
      float e = PDFs[x] -> GetBinContent(i+1); 
      PDFs[x] -> SetBinContent(i+1, ratio*e); 
    }
  }
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

void DerivedFunctions::ReplaceShiftTail(TH1F* Source, TH1F* Target, float offset)
{
  BaseFunctions B;

  int bin_s = Target -> GetNbinsX(); 
  float Lumi = Target -> Integral();
  float Lumi_S = Source -> Integral();
  B.Normalize({Source, Target});

  // Extend the source and target hist to avoid nan 
  std::vector<float> Source_V = B.TH1FDataVector(Source, offset);
  std::vector<float> Target_V = B.TH1FDataVector(Target, offset);  
  TH1F* Source_E = new TH1F("Source_E", "Source_E", Source_V.size(), 0, Source_V.size());
  TH1F* Target_E = new TH1F("Target_E", "Target_E", Target_V.size(), 0, Target_V.size()); 
  B.ToTH1F(Source_V, Source_E); 
  B.ToTH1F(Target_V, Target_E);

  // Define the tail to be replaced and move away from peak 
  int max_bin = Target -> GetMaximumBin() + bin_s*0.01;
  int max_bin_s = Source -> GetMaximumBin();
  TH1F* Temp = (TH1F*)Source -> Clone("Temp"); 
  
  for (int i(0); i < max_bin_s; i++)
  {
    Temp -> SetBinContent(i+1, 0);
  }
  int shift = NumericalShift(Temp, Target); 
  
  float t_S = Source_E -> GetBinContent(max_bin + 1 - shift -1);
  float t_T = Target_E -> GetBinContent(max_bin + 1);
  float sum_S = 0;
  float sum_T = 0;  
  for (int i(max_bin); i < bin_s; i++) 
  {
    float S_e = Source_E -> GetBinContent(i + 1 - shift -1);
    float T_e = Target_E -> GetBinContent(i + 1); 
    sum_S += S_e;
    sum_T += T_e;  
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
  std::vector<TString> Names = {"TRK1", "TRK2", "TRK3", "TRK4", "TRK5"};  
  std::vector<TH1F*> PDFs = B.MakeTH1F(Names, trk2);
 
  // Deconvolve the 2trk data into 1 trk  
  for (int i(0); i < iter; i++)
  {
    deconv = B.LucyRichardson(trk2_V, deconv, deconv); 
  }
  
  // Replace the tail of the deconv vector 
  B.ToTH1F(deconv, PDFs[0]);

  // === TRK1
  ReplaceShiftTail(trk1, PDFs[0], offset); 
  B.Normalize(PDFs);
  B.ResidualRemove(PDFs[0]);

  // === TRK2
  B.ConvolveHists(PDFs[0], PDFs[0], PDFs[1]); 
  B.Normalize(PDFs[1]);
  B.ResidualRemove(PDFs[1]); 

  // === TRK3
  B.ConvolveHists(PDFs[1], PDFs[0], PDFs[2]);  
  B.Normalize(PDFs[2]);
  B.ResidualRemove(PDFs[2]); 

  // === TRK4 
  B.ConvolveHists(PDFs[2], PDFs[0], PDFs[3]);  
  B.Normalize(PDFs[3]);
  B.ResidualRemove(PDFs[3]); 

  // === TRK5
  B.ConvolveHists(PDFs[3], PDFs[0], PDFs[4]);  
  B.Normalize(PDFs[4]);
  B.ResidualRemove(PDFs[4]); 

 
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
  TH1F* GxH_L = new TH1F("GxH_L", "GxH_L", 2*bins, -Padding, bins + Padding);

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

std::map<TString, float> DerivedFunctions::FitGaussian(TH1F* GxTrk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> Params, float offset, int iter)
{
  auto GausConvolute = [](TH1F* PDF, float offset, int iter, TH1F* Hist, TH1F* Gaus)
  {
    BaseFunctions B;
    int bins = PDF -> GetNbinsX();
    std::vector<float> PDF_V = B.TH1FDataVector(PDF, offset);
    TString v = PDF -> GetTitle(); v += ("_Temp");
    TH1F* Temp = new TH1F(v, v, PDF_V.size(), 0, PDF_V.size()); 
    B.ToTH1F(PDF_V, Temp);
    
    B.ShiftExpandTH1F(Temp, Hist, bins*offset/2);
    Temp -> Reset();
    
    PDF_V = B.TH1FDataVector(Hist, 0);
    Hist -> Reset();

    std::vector<float> deconv(PDF_V.size(), 0.5);
    B.ShiftExpandTH1F(Gaus, Hist);

    std::vector<float> PSF_V = B.TH1FDataVector(Hist, 0); 
    Hist -> Reset(); 

    for (int x(0); x < iter; x++) 
    {
      deconv = B.LucyRichardson(PDF_V, PSF_V, deconv);
    }
    
    B.ToTH1F(deconv, Hist); 
       
    delete Temp;    
  }; 
  
  DistributionGenerators DG; 
  BaseFunctions B;  
  Plotting P; 


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
  TH1F* Gaussian = (TH1F*)Hists[0] -> Clone("Gaussian");
  std::vector<float> Gaus_Setting = Params["Gaussian"];
  DG.Gaussian(Gaus_Setting[0], Gaus_Setting[1], Constants::GaussianToys, Gaussian);
  B.Normalize(Gaussian); 
  
  std::vector<std::thread> th; 
  for (int i(0); i < PDFs.size(); i++)
  {
    th.push_back(std::thread(GausConvolute, PDFs[i], offset, iter, Hists[i], Gaussian));
  }

  for (std::thread &t : th){ t.join(); }
   
  // Convert the histograms back to the original format 
  std::vector<TH1F*> DeconvPDFs;
  for (int i(0); i < PDFs.size(); i++)
  {
    std::vector<float> Data = B.TH1FDataVector(Hists[i]);  
  
    TString name = PDFs[i] -> GetTitle(); name += ("_"); name += (i+1);
    TH1F* H = new TH1F(name, name, bins_pdf, 0, bins_pdf);  
    H -> Reset();
    B.ToTH1F(Data, H);
    B.ResidualRemove(H);
    DeconvPDFs.push_back(H); 
  }
  TH1F* Data_G = new TH1F("Data_Gaus_Fit", "Data_Gaus_Fit", bins_pdf, 0, bins_pdf);
  B.ShiftExpandTH1F(GxTrk, Data_G, 0);
  Data_G -> SetBinContent(bins_pdf, Data_G -> GetBinContent(bins_pdf-1));
  B.Normalize(Data_G);

  // Define the range of the dEdx
  RooRealVar* x = new RooRealVar("x", "x", 0, bins_pdf); 

  // Define the Gaussian Parameter: Mean
  std::vector<TString> Means_String = { "m1", "m2", "m3", "m4", "m5"};
  std::vector<float> Means_Begin = Params["m_s"];
  std::vector<float> Means_End = Params["m_e"];
  std::vector<float> Means_Initial = Params["m_i"];
  std::vector<RooRealVar*> Means = B.RooVariables(Means_String, Means_Begin, Means_End);

  // Define the Gaussian Parameter: Standard Deviation
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4", "s5"};
  std::vector<float> Stdev_Begin = Params["s_s"];
  std::vector<float> Stdev_End = Params["s_e"]; 
  std::vector<float> Stdev_Initial = Params["s_i"];
  std::vector<RooRealVar*> Stdev = B.RooVariables(Stdev_String, Stdev_Begin, Stdev_End);

  // Define the Gaussian Variables
  std::vector<TString> Gaus_String = { "g1", "g2", "g3", "g4", "g5"};
  std::vector<RooGaussian*> G_Vars = B.RooVariables(Gaus_String, Means, Stdev, x);

  // Import the PDFs as a RooDataHist
  B.Normalize(DeconvPDFs);
  std::vector<RooHistPdf*> PDF_Vars = B.RooPDF(DeconvPDFs, x);

  // Define the ntrack coefficients:
  float Lumi = Data_G -> Integral(); 
  std::vector<TString> C_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4", "n_trk5" };
  std::vector<float> C_Begin = { 0., 0., 0., 0., 0.};
  std::vector<float> C_End = { Lumi, Lumi, Lumi, Lumi, Lumi };
  std::vector<RooRealVar*> C_Vars = B.RooVariables(C_String, C_Begin, C_End);

  // Convolve the PDFs with the Gaussians
  std::vector<TString> Conv_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4", "P5xG5"};
  std::vector<RooFFTConvPdf*> Conv_Vars = B.RooVariables(Conv_String, PDF_Vars, G_Vars, x);
  
//  //=== Constraint variables: Standard Deviation 
//  std::vector<TString> sGaus_Constrain_String = {"sg1", "sg2", "sg3", "sg4", "sg5"};
//  std::vector<float> Stdev_Constant_V1 = Params["s_c_v1"];
//  std::vector<float> Stdev_Constant_V2 = Params["s_c_v2"];
//  std::vector<RooGaussian*> sGaus_Constrain = B.RooVariables(sGaus_Constrain_String, Stdev, Stdev_Constant_V1, Stdev_Constant_V2);
//
//  //=== Constraint variables: Mean 
//  std::vector<TString> mGaus_Constrain_String = {"mg1", "mg2", "mg3", "mg4", "mg5"};
//  std::vector<float> Mean_Constant_V1 = Params["m_c_v1"];
//  std::vector<float> Mean_Constant_V2 = Params["m_c_v2"];
//  std::vector<RooGaussian*> mGaus_Constrain = B.RooVariables(mGaus_Constrain_String, Means, Mean_Constant_V1, Mean_Constant_V2);

  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk2_D = new RooDataHist("trk2_D", "trk2_D", *x, Data_G); 
  
  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgSet(*Conv_Vars[0], *Conv_Vars[1], *Conv_Vars[2], *Conv_Vars[3], *Conv_Vars[4]), 
                                    RooArgSet(*C_Vars[0], *C_Vars[1], *C_Vars[2], *C_Vars[3], *C_Vars[4])); 

  //model.fitTo(*trk2_D, RooFit::ExternalConstraints(*mGaus_Constrain[0]),
  //                     RooFit::ExternalConstraints(*mGaus_Constrain[1]), 
  //                     RooFit::ExternalConstraints(*mGaus_Constrain[2]),
  //                     RooFit::ExternalConstraints(*mGaus_Constrain[3]),
  //                     RooFit::ExternalConstraints(*mGaus_Constrain[4]));

  RooFitResult* stat = model.fitTo(*trk2_D, RooFit::SumW2Error(true), RooFit::Save());//, RooFit::ExternalConstraints(*sGaus_Constrain[0]),
                                                            //RooFit::ExternalConstraints(*sGaus_Constrain[1]),
                                                            //RooFit::ExternalConstraints(*sGaus_Constrain[2]),  
                                                            //RooFit::ExternalConstraints(*sGaus_Constrain[3]),  
                                                            //RooFit::ExternalConstraints(*sGaus_Constrain[4])); 
   
    
  //TCanvas* can = P.PlotHists(model, x, Conv_Vars, trk2_D);
  //can -> Print("out.pdf"); 
 
  
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
   // delete sGaus_Constrain[i];
   // delete mGaus_Constrain[i];
  }
  delete trk2_D; 
  delete x;
  delete Data_G;
  delete stat;

  return Out;   
}

std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> DerivedFunctions::MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, float Gamma, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure)
{
  auto MakeGaussianConvoluted = [](TString n, std::vector<TH1F*> PDFs, std::map<TString, float> trkn_Params, TH1F* trkn_L, float Gamma)
  {
    BaseFunctions B;
    DerivedFunctions DF;
    const std::vector<TString> Names = {"n_trk1", "n_trk2", "n_trk3", "n_trk4", "n_trk5"}; 
    const std::vector<TString> Stdev = {"s1", "s2", "s3", "s4", "s5"}; 
    const std::vector<TString> Mean = {"m1", "m2", "m3", "m4", "m5"};
    
    std::vector<TH1F*> GxTrkN(Names.size());
    for (int i(0); i < Names.size(); i++)
    {
      TString Name = n; Name += Names[i];
      TH1F* trknXg = (TH1F*)trkn_L -> Clone(Name);
      trknXg -> SetTitle(Name);
      trknXg -> Reset();
      TH1F* trkn_H = DF.GaussianConvolve(PDFs[i], trkn_Params[Mean[i]], trkn_Params[Stdev[i]]);
      B.ShiftExpandTH1F(trkn_H, trknXg);
      delete trkn_H;

      float lumi = trkn_L -> Integral();
      
      // Check for nan
      if (trkn_Params[Names[i]] < 0){Gamma = 0;}
      trknXg -> Scale(Gamma*trkn_Params[Names[i]]*lumi);
      GxTrkN[i] = trknXg;
    }
    return GxTrkN;
  };

  auto Parallel = [&MakeGaussianConvoluted](TH1F* trk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> Params, float offset, int iter, std::vector<TH1F*> GxTrk, TString name)
  {
    DerivedFunctions DF; 
    std::map<TString, float> trk_Param = DF.FitGaussian(trk, PDFs, Params, offset, iter);
    std::vector<TH1F*> Out = MakeGaussianConvoluted(name, PDFs, trk_Param, trk, 1);
    for (int i(0); i < Out.size(); i++)
    {
      GxTrk.push_back(Out[i]);     
    }
    return GxTrk;
  };

  BaseFunctions B;
  Plotting P; 
  TH1::AddDirectory(false);

  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> Output; 
  RooRealVar* x = new RooRealVar("x", "x", 0, ntrk[0] -> GetNbinsX());  
  TH1F* trk1_L = new TH1F("trk1_L", "trk1_L", ntrk[0] -> GetNbinsX(), 0, ntrk[0] -> GetNbinsX());
  TH1F* trk2_L = new TH1F("trk2_L", "trk2_L", ntrk[1] -> GetNbinsX(), 0, ntrk[1] -> GetNbinsX());
  TH1F* trk3_L = new TH1F("trk3_L", "trk3_L", ntrk[2] -> GetNbinsX(), 0, ntrk[2] -> GetNbinsX());
  TH1F* trk4_L = new TH1F("trk4_L", "trk4_L", ntrk[3] -> GetNbinsX(), 0, ntrk[3] -> GetNbinsX());
  TH1F* trk5_L = new TH1F("trk5_L", "trk5_L", ntrk[4] -> GetNbinsX(), 0, ntrk[4] -> GetNbinsX());

  B.ShiftExpandTH1F(ntrk[0], trk1_L); 
  B.ShiftExpandTH1F(ntrk[1], trk2_L);  
  B.ShiftExpandTH1F(ntrk[2], trk3_L); 
  B.ShiftExpandTH1F(ntrk[3], trk4_L); 
  B.ShiftExpandTH1F(ntrk[4], trk5_L); 

  TH1F* trk1_L_C = (TH1F*)trk1_L -> Clone("trk1_L_C");
  TH1F* trk2_L_C = (TH1F*)trk2_L -> Clone("trk2_L_C");
  TH1F* trk3_L_C = (TH1F*)trk3_L -> Clone("trk3_L_C");
  TH1F* trk4_L_C = (TH1F*)trk4_L -> Clone("trk4_L_C");
  TH1F* trk5_L_C = (TH1F*)trk5_L -> Clone("trk5_L_C");
    
  // =========================== Plotting variables 
  TString Title_Params = "GlobalSettings-Gamma-"; Title_Params += (Gamma); Title_Params += ("-iter-"); Title_Params += (iter); Title_Params += ("-cor_loop-"); Title_Params += (cor_loop); 
   
  std::vector<TString> Name = {"m_s", "m_e", "s_s", "s_e"}; 
  for ( TString n : Name )
  {
    std::vector<float> t = Params[n]; 
    Title_Params += (n + ":");
    for (float x : t)
    {
      double h = (int)(10000*x);
      h = h/10000;
      Title_Params += (h); Title_Params += (",");
    }
    Title_Params += (":");
  }
  Title_Params += (".pdf");

  TCanvas* can = new TCanvas();
  can -> Print(Title_Params + "[");
  gStyle -> SetOptStat(0);
  // ======================================================= 

  for (int x(0); x < cor_loop; x++)
  { 
    std::vector<TH1F*> GxTrk1;   
    std::vector<TH1F*> GxTrk2;
    std::vector<TH1F*> GxTrk3;
    std::vector<TH1F*> GxTrk4; 
    std::vector<TH1F*> GxTrk5; 
    std::vector<TH1F*> PDFs;  

    // Generate the minimal version of the PDFs
    PDFs = nTRKGenerator(trk1_L, trk2_L, offset, iter);

    // Generate the Gaussian Parameters for the convolution from a fit 
    GxTrk1 = Parallel(trk1_L, PDFs, Params, offset, iter, GxTrk1, "GxTrk1");  
    //GxTrk2 = Parallel(trk2_L, PDFs, Params, offset, iter, GxTrk2, "GxTrk2"); 
    GxTrk3 = Parallel(trk3_L, PDFs, Params, offset, iter, GxTrk3, "GxTrk3");
    GxTrk4 = Parallel(trk4_L, PDFs, Params, offset, iter, GxTrk4, "GxTrk4"); 
    GxTrk5 = Parallel(trk5_L, PDFs, Params, offset, iter, GxTrk5, "GxTrk5"); 

    float sc = cor_loop -1; 
    float it = x; 
    sc = (1/sc)*it;    
     
    const float del = 0.05; 
    int sub = 5; 
    trk2_L -> Reset(); 
    trk2_L -> Add(trk2_L_C); 
    for (int v(0); v < sub ; v++)
    {
      std::vector<TH1F*> PDFx = nTRKGenerator(trk1_L, trk2_L, offset, iter);
      std::vector<TH1F*> Temp; 
      Temp = Parallel(trk2_L, PDFs, Params, offset, iter, Temp, "GxTrk2");

      iter = iter +1;
      trk2_L -> Add(Temp[0], -0.25);  
      trk2_L -> Add(Temp[2], -del);     
      trk2_L -> Add(Temp[3], -del);   
      trk2_L -> Add(Temp[4], -del);   
      iter++;
    
      if ( v < sub -1)
      { 
        for (int z(0); z < Temp.size(); z++)
        {
          delete Temp[z];
          delete PDFx[z];   
        } 
      }
      else 
      {
        for (int z(0); z < Temp.size(); z++)
        {
          GxTrk2.push_back(Temp[z]);
          delete PDFx[z];
        }
      }
    }
     
    trk2_L -> Reset();
    trk2_L -> Add(trk2_L_C); 
    trk3_L -> Reset();
    trk3_L -> Add(trk3_L_C);
    trk4_L -> Reset();
    trk4_L -> Add(trk4_L_C);
    trk5_L -> Reset();
    trk5_L -> Add(trk4_L_C);

    SafeScale(GxTrk2, trk2_L); 
    SafeScale(GxTrk3, trk3_L); 
    SafeScale(GxTrk4, trk4_L); 
    iter = iter +10;
    
    trk2_L -> Add(GxTrk2[0], -1);
    trk2_L -> Add(GxTrk2[2], -sc);
    trk2_L -> Add(GxTrk2[3], -sc);
    trk2_L -> Add(GxTrk2[4], -sc);

    trk3_L -> Add(GxTrk3[0], -sc);
    trk3_L -> Add(GxTrk3[1], -sc);
    trk3_L -> Add(GxTrk3[3], -sc);
    trk3_L -> Add(GxTrk3[4], -sc);

    trk4_L -> Add(GxTrk4[0], -sc);
    trk4_L -> Add(GxTrk4[1], -sc);
    trk4_L -> Add(GxTrk4[2], -sc);
    trk4_L -> Add(GxTrk4[4], -sc);

    trk5_L -> Add(GxTrk5[0], -sc);
    trk5_L -> Add(GxTrk5[1], -sc);
    trk5_L -> Add(GxTrk5[2], -sc);                            
    trk5_L -> Add(GxTrk5[3], -sc);

    std::cout << "################### " << x << std::endl;
   
    // ==== Trk1 
    std::vector<TString> Names1 = {"trk1_F1", "trk2_F1", "trk3_F1", "trk4_F1", "trk5_F1"}; 
    std::vector<TH1F*> Trk1_PDFs = B.MakeTH1F(Names1, ntrk[0]); 
    B.ShiftExpandTH1F(GxTrk1, Trk1_PDFs);
    
    // ==== Trk2
    std::vector<TString> Names2 = {"trk1_F2", "trk2_F2", "trk3_F2", "trk4_F2", "trk5_F2"}; 
    std::vector<TH1F*> Trk2_PDFs = B.MakeTH1F(Names2, ntrk[1]); 
    B.ShiftExpandTH1F(GxTrk2, Trk2_PDFs);

    // ==== Trk3 
    std::vector<TString> Names3 = {"trk1_F3", "trk2_F3", "trk3_F3", "trk4_F3", "trk5_F3"}; 
    std::vector<TH1F*> Trk3_PDFs = B.MakeTH1F(Names3, ntrk[2]); 
    B.ShiftExpandTH1F(GxTrk3, Trk3_PDFs);

    // ==== Trk4 
    std::vector<TString> Names4 = {"trk1_F4", "trk2_F4", "trk3_F4", "trk4_F4", "trk5_F4"}; 
    std::vector<TH1F*> Trk4_PDFs = B.MakeTH1F(Names4, ntrk[3]); 
    B.ShiftExpandTH1F(GxTrk4, Trk4_PDFs);

    //// ==== Plotting output 
    //// Trk1 
    //can -> Clear();
    //can -> Divide(1); 
    //P.PlotHists({Trk1_PDFs}, Closure[0], ntrk[0], can); 
    //can -> Print(Title_Params);

    //// Trk2 
    //can -> Clear();
    //can -> Divide(1);
    //P.PlotHists({Trk2_PDFs}, Closure[1], ntrk[1], can); 
    //can -> Print(Title_Params);

    //// Trk3 
    //can -> Clear();
    //can -> Divide(1);
    //P.PlotHists({Trk3_PDFs}, Closure[2], ntrk[2], can); 
    //can -> Print(Title_Params);

    //// Trk4 
    //can -> Clear();
    //can -> Divide(1);
    //P.PlotHists({Trk4_PDFs}, Closure[3], ntrk[3], can); 
    //can -> Print(Title_Params);
    
    // Final Plot where we compare subtracted with the predicted PDF and the real distribution 
    std::vector<TString> Names_Sub = {"trk1_Sub", "trk2_Sub", "trk3_Sub", "trk4_Sub"}; 
    std::vector<TH1F*> Tracks = B.MakeTH1F(Names_Sub, ntrk[2]); 
    B.ShiftExpandTH1F({trk1_L, trk2_L, trk3_L, trk4_L}, Tracks);
    std::vector<TH1F*> Truth = {Closure[0][0], Closure[1][1], Closure[2][2], Closure[3][3]};
 
    // All combined
    can -> Clear();
    can -> Divide(2, 2);
    P.PlotHists({Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs}, Closure, ntrk, can);   
    can -> Print(Title_Params);    

    can -> Clear();
    can -> Divide(2,2);
    P.PlotHists({Trk1_PDFs[0], Trk2_PDFs[1], Trk3_PDFs[2], Trk4_PDFs[3]}, Tracks, Truth, can);
    can -> Print(Title_Params);




    // Output
    Output[0] = std::make_pair(trk1_L, GxTrk1);
    Output[1] = std::make_pair(trk2_L, GxTrk2);
    Output[2] = std::make_pair(trk3_L, GxTrk3); 
    Output[3] = std::make_pair(trk4_L, GxTrk4);

    if (x == cor_loop -1){break;}     
    
    // Clean up memory  
    for (int i(0); i < GxTrk2.size(); i++)
    {
      delete GxTrk2[i];      
      delete Trk2_PDFs[i];
      delete PDFs[i];
     
      delete GxTrk1[i];
      delete GxTrk3[i];
      delete GxTrk4[i];
      delete GxTrk5[i];

      delete Trk1_PDFs[i];
      delete Trk3_PDFs[i];
      delete Trk4_PDFs[i];   
    }
  }
  can -> Print(Title_Params + ")");


  // Clear up some memory
  delete can; 
  //for (int i(0); i < Output.size(); i++)
  //{
  //  for (int x(0); x < Output[i].second.size(); x++)
  //  { 
  //    delete Output[i].second.at(x);
  //  }
  //}
  return Output;
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
