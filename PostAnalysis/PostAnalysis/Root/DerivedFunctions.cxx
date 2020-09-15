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
  
  int bins = Data -> GetNbinsX(); 
  for (int i(0); i < bins; i++)
  {
    // Get the sum of the PDFs for the bin we are looking at 
    float sum_pdfs = 0; 
    for (int x(0); x < PDFs.size(); x++)
    {
      sum_pdfs += PDFs[x] -> GetBinContent(i+1); 
    }
   
    // Get the number of entries within the data for that bin 
    float data_e = Data -> GetBinContent(i+1); 
    
    // Get the ratio between the data and the PDFs predictions 
    float ratio = data_e/sum_pdfs;  

    // Some fail safe conditions 
    if (sum_pdfs == 0) {ratio = 1e-21;}
    if (std::isnan(ratio)) {ratio = 1e-21;}

    // Update the bins of the PDFs
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
    Target -> SetBinContent(i+1, S_e); //*(t_T/t_S));
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

  // === TRK2
  B.ConvolveHists(PDFs[0], PDFs[0], PDFs[1]); 
  B.Normalize(PDFs[1]);

  // === TRK3
  B.ConvolveHists(PDFs[1], PDFs[0], PDFs[2]);  
  B.Normalize(PDFs[2]);

  // === TRK4 
  B.ConvolveHists(PDFs[2], PDFs[0], PDFs[3]);  
  B.Normalize(PDFs[3]);

  // === TRK5
  B.ConvolveHists(PDFs[3], PDFs[0], PDFs[4]);  
  B.Normalize(PDFs[4]);
 
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

  // Get the number of bins of the PDFs
  int bins_pdf = PDFs[0] -> GetNbinsX(); 

  // Generate the names for the temporary histograms and generate them 
  std::vector<TString> names; 
  for (int i(0); i < PDFs.size(); i++)
  {
    TString name = PDFs[i] -> GetTitle(); name += (i+1);  
    names.push_back(name); 
  } 
  std::vector<TH1F*> Hists = B.MakeTH1F(names, GxTrk); 


  // Generate a Gaussian distribution and normalize it 
  TH1F* Gaussian = (TH1F*)Hists[0] -> Clone("Gaussian");
  std::vector<float> Gaus_Setting = Params["Gaussian"];
  DG.Gaussian(Gaus_Setting[0], Gaus_Setting[1], Constants::GaussianToys, Gaussian);
  B.Normalize(Gaussian); 
 
  // Deconvolve the above Gaussian with the track PDF using a multithreaded approach 
  std::vector<std::thread> th; 
  for (int i(0); i < PDFs.size(); i++){th.push_back(std::thread(GausConvolute, PDFs[i], offset, iter, Hists[i], Gaussian));}
  for (std::thread &t : th){ t.join(); }
 
  // Convert the histograms back to the original format (the original length of the histogram
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

  // Transfer the Data to a new histogram to remove the need to know min max values 
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
  std::vector<RooRealVar*> Means = B.RooVariables(Means_String, Means_Begin, Means_End);

  // Define the Gaussian Parameter: Standard Deviation
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4", "s5"};
  std::vector<float> Stdev_Begin = Params["s_s"];
  std::vector<float> Stdev_End = Params["s_e"]; 
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
  
  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk2_D = new RooDataHist("trk2_D", "trk2_D", *x, Data_G); 
  
  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgSet(*Conv_Vars[0], *Conv_Vars[1], *Conv_Vars[2], *Conv_Vars[3], *Conv_Vars[4]), 
                                    RooArgSet(*C_Vars[0], *C_Vars[1], *C_Vars[2], *C_Vars[3], *C_Vars[4])); 
  RooFitResult* stat = model.fitTo(*trk2_D, RooFit::SumW2Error(true), RooFit::Save());
   
   
  // Store the output in a map and delete pointers 
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

std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> DerivedFunctions::MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure)
{
  // Lambda Function: This does the convolution and takes care of the labels 
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
      if (std::isnan(trkn_Params[Names[i]])){Gamma = 0;}
      trknXg -> Scale(Gamma*trkn_Params[Names[i]]*lumi);
      GxTrkN[i] = trknXg;
    }
    return GxTrkN;
  };

  auto DynamicVariables = [] (std::map<TString, std::vector<float>> *Params, std::map<TString, float> FitPredictions)
  {  
    std::vector<TString> Means_String = {"m1", "m2", "m3", "m4", "m5"};
    std::vector<TString> Stdev_String = {"s1", "s2", "s3", "s4", "s5"};
   
    std::map<TString, std::vector<float>> &P = *Params;
    for ( int i(0); i < Means_String.size(); i++ )
    {
      TString M = Means_String[i]; 
      TString S = Stdev_String[i]; 
      
      float P_M = FitPredictions[M];  
      float P_S = FitPredictions[S]; 
 
      if (P["m_s"][i] == P_M){P["m_s"][i] = P_M++;}
      if (P["m_e"][i] == P_M){P["m_e"][i] = P_M++;}
      if (P["s_s"][i] == P_S){P["s_s"][i] = P_S++;}
      if (P["s_e"][i] == P_S){P["s_e"][i] = P_S*0.5;}
    } 
  };

  // Lambda Function: Condensed version of the above. Wanted to use this as a multithread function but RooFit doesnt seem to be multithreading safe 
  auto Parallel = [&MakeGaussianConvoluted, &DynamicVariables](TH1F* trk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> *Params, float offset, int iter, std::vector<TH1F*> GxTrk, TString name)
  {
    DerivedFunctions DF; 
    std::map<TString, float> trk_Param = DF.FitGaussian(trk, PDFs, *Params, offset, iter);
    std::vector<TH1F*> Out = MakeGaussianConvoluted(name, PDFs, trk_Param, trk, 1);
    DynamicVariables(Params, trk_Param);
    for (int i(0); i < Out.size(); i++){GxTrk.push_back(Out[i]);}
    return GxTrk;
  };

  BaseFunctions B;
  Plotting P; 
  TH1::AddDirectory(false);
  
  // Create histograms for the ntrack data
  TH1F* trk1_L = new TH1F("trk1_L", "trk1_L", ntrk[0] -> GetNbinsX(), 0, ntrk[0] -> GetNbinsX());
  TH1F* trk2_L = new TH1F("trk2_L", "trk2_L", ntrk[1] -> GetNbinsX(), 0, ntrk[1] -> GetNbinsX());
  TH1F* trk3_L = new TH1F("trk3_L", "trk3_L", ntrk[2] -> GetNbinsX(), 0, ntrk[2] -> GetNbinsX());
  TH1F* trk4_L = new TH1F("trk4_L", "trk4_L", ntrk[3] -> GetNbinsX(), 0, ntrk[3] -> GetNbinsX());
  TH1F* trk5_L = new TH1F("trk5_L", "trk5_L", ntrk[4] -> GetNbinsX(), 0, ntrk[4] -> GetNbinsX());

  // Fill those histograms with the ntrack data
  B.ShiftExpandTH1F(ntrk[0], trk1_L); 
  B.ShiftExpandTH1F(ntrk[1], trk2_L);  
  B.ShiftExpandTH1F(ntrk[2], trk3_L); 
  B.ShiftExpandTH1F(ntrk[3], trk4_L); 
  B.ShiftExpandTH1F(ntrk[4], trk5_L); 

  // Make clones of these for restoring the state after each iteration 
  TH1F* trk1_L_C = (TH1F*)trk1_L -> Clone("trk1_L_C");
  TH1F* trk2_L_C = (TH1F*)trk2_L -> Clone("trk2_L_C");
  TH1F* trk3_L_C = (TH1F*)trk3_L -> Clone("trk3_L_C");
  TH1F* trk4_L_C = (TH1F*)trk4_L -> Clone("trk4_L_C");
  TH1F* trk5_L_C = (TH1F*)trk5_L -> Clone("trk5_L_C");
    
  // Forward declaration 
  std::map<TString, std::vector<float>> Params1;
  std::map<TString, std::vector<float>> Params2;
  std::map<TString, std::vector<float>> Params3;
  std::map<TString, std::vector<float>> Params4;
  std::map<TString, std::vector<float>> Params5;

  // ============= Output Variables ============= //
  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> Output; 
  TH1F* FLOST_Prediction = new TH1F("FLost_Pred", "FLost_Pred", cor_loop, 0, cor_loop); 
  TH1F* FLOST_Truth = new TH1F("FLost_Truth", "FLost_Truth", cor_loop, 0, cor_loop); 

  //TCanvas* can_HD = new TCanvas();
  //TString Title_HD = "out.pdf";
  for (int x(0); x < cor_loop; x++)
  { 
    // Forward declaration 
    std::vector<TH1F*> GxTrk1;   
    std::vector<TH1F*> GxTrk2;
    std::vector<TH1F*> GxTrk3;
    std::vector<TH1F*> GxTrk4; 
    std::vector<TH1F*> GxTrk5; 
    std::vector<TH1F*> PDFs;  

    // Generate the minimal version of the PDFs without any Gaussian convolution 
    PDFs = nTRKGenerator(trk1_L, trk2_L, offset, iter);

    // Create separate Params containers so that we can add some dynamics 
    if ( x == 0)
    {
      Params1 = Params; 
      Params2 = Params; 
      Params3 = Params; 
      Params4 = Params; 
      Params5 = Params; 
    }

    // Generate the Gaussian Parameters for the convolution from a fit
    GxTrk1 = Parallel(trk1_L, PDFs, &Params1, offset, iter, GxTrk1, "GxTrk1");  
    GxTrk2 = Parallel(trk2_L, PDFs, &Params2, offset, iter, GxTrk2, "GxTrk2"); 
    GxTrk3 = Parallel(trk3_L, PDFs, &Params3, offset, iter, GxTrk3, "GxTrk3");
    GxTrk4 = Parallel(trk4_L, PDFs, &Params4, offset, iter, GxTrk4, "GxTrk4"); 
    GxTrk5 = Parallel(trk5_L, PDFs, &Params5, offset, iter, GxTrk5, "GxTrk5"); 

    // Increase the iteration after each correction loop 
    iter = iter + 2;

    // Reset the state of the data histograms 
    trk1_L -> Reset();
    trk1_L -> Add(trk1_L_C);     
    trk2_L -> Reset();
    trk2_L -> Add(trk2_L_C); 
    trk3_L -> Reset();
    trk3_L -> Add(trk3_L_C);
    trk4_L -> Reset();
    trk4_L -> Add(trk4_L_C);
    trk5_L -> Reset();
    trk5_L -> Add(trk5_L_C);

    // Scale the PDFs according to the number of events in a bin 
    SafeScale(GxTrk1, trk1_L); 
    SafeScale(GxTrk2, trk2_L); 
    SafeScale(GxTrk3, trk3_L); 
    SafeScale(GxTrk4, trk4_L); 
    SafeScale(GxTrk5, trk5_L); 

    // Do the subtraction  
    trk1_L -> Add(GxTrk1[1], -1); 
    trk1_L -> Add(GxTrk1[2], -1); 
    trk1_L -> Add(GxTrk1[3], -1); 
    trk1_L -> Add(GxTrk1[4], -1); 
 
    trk2_L -> Add(GxTrk2[0], -1);
    trk2_L -> Add(GxTrk2[2], -1);
    trk2_L -> Add(GxTrk2[3], -1);
    trk2_L -> Add(GxTrk2[4], -1);

    trk3_L -> Add(GxTrk3[0], -1);
    trk3_L -> Add(GxTrk3[1], -1);
    trk3_L -> Add(GxTrk3[3], -1);
    trk3_L -> Add(GxTrk3[4], -1);

    trk4_L -> Add(GxTrk4[0], -1);
    trk4_L -> Add(GxTrk4[1], -1);
    trk4_L -> Add(GxTrk4[2], -1);
    trk4_L -> Add(GxTrk4[4], -1);

    trk5_L -> Add(GxTrk5[0], -1);
    trk5_L -> Add(GxTrk5[1], -1);
    trk5_L -> Add(GxTrk5[2], -1);                            
    trk5_L -> Add(GxTrk5[3], -1);

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

    // ==== Trk5 
    std::vector<TString> Names5 = {"trk1_F5", "trk2_F5", "trk3_F5", "trk4_F5", "trk5_F5"}; 
    std::vector<TH1F*> Trk5_PDFs = B.MakeTH1F(Names5, ntrk[4]); 
    B.ShiftExpandTH1F(GxTrk5, Trk5_PDFs);

    // Create the variabls needed for the Plotting below  
    std::vector<TString> Names_Sub = {"trk1_Sub", "trk2_Sub", "trk3_Sub", "trk4_Sub"}; 
    std::vector<TH1F*> Tracks = B.MakeTH1F(Names_Sub, ntrk[2]); 
    B.ShiftExpandTH1F({trk1_L, trk2_L, trk3_L, trk4_L}, Tracks);
    std::vector<TH1F*> Truth = {Closure[0][0], Closure[1][1], Closure[2][2], Closure[3][3]};
    // ======== Section for the output ========== //
    // === Save the PDFs and the subtracted Hists

    // === Save the subtracted Hists
    std::vector<TString> Names_Pure = {"trk1_Pure", "trk2_Pure", "trk3_Pure", "trk4_Pure", "trk5_Pure"};    
    std::vector<TH1F*> Track_Out = B.MakeTH1F(Names_Pure, ntrk[2]); 
    B.ShiftExpandTH1F({trk1_L, trk2_L, trk3_L, trk4_L, trk5_L}, Track_Out);

    std::vector<TH1F*> Pred_PDF;
    std::vector<std::vector<TH1F*>> Set = {Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs, Trk5_PDFs, Track_Out, ntrk, Closure[0], Closure[1], Closure[2], Closure[3], Closure[4]}; 

    for (std::vector<TH1F*> PDF_S : Set)
    {
      for (int g(0); g < PDF_S.size(); g++)
      {
        TString name_pdf = PDF_S[g] -> GetName(); name_pdf += (x); 
        Pred_PDF.push_back((TH1F*)PDF_S[g] -> Clone(name_pdf));
      }
    }

    // === Save the FLost progress
    float FLost_Pred = B.FLost(ntrk, {Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs, Trk5_PDFs}); 
    float FLost_MC = B.FLost(ntrk, Closure);   
    FLOST_Prediction -> SetBinContent(x+1, FLost_Pred); 
    FLOST_Truth -> SetBinContent(x+1, FLost_MC);   
    
    TString name_out =  "FLost_P.at."; name_out += (x); 
    TH1F* FL_P_Copy = (TH1F*)FLOST_Prediction -> Clone(name_out);

    name_out =  "FLost_T.at."; name_out += (x); 
    TH1F* FL_T_Copy = (TH1F*)FLOST_Truth -> Clone(name_out);
    Pred_PDF.push_back(FL_T_Copy); 

    Output[x] = std::make_pair(FL_P_Copy, Pred_PDF);
 
    // Clean up memory  
    for (int i(0); i < GxTrk2.size(); i++)
    {
      delete PDFs[i];
     
      delete GxTrk1[i]; 
      delete GxTrk2[i];      
      delete GxTrk3[i];
      delete GxTrk4[i];
      delete GxTrk5[i];

      delete Trk1_PDFs[i];
      delete Trk2_PDFs[i];
      delete Trk3_PDFs[i];
      delete Trk4_PDFs[i];   
      delete Trk5_PDFs[i];
    }
  }
  //can -> Print(Title_Params + ")");
  //can_HD -> Print(Title_HD + ")"); 
 
  //delete can; 
  //delete can_HD; 
  return Output;
}


