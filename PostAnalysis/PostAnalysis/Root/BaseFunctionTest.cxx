#include<PostAnalysis/BaseFunctionTest.h>

void PlotLandau()
{
  int bins = 1000; 
  int min = 0; 
  int max = 20; 

  std::vector<TString> Names = {"Landau"}; 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<TH1F*> Lan = Landau(Names, {1.0}, LandauParams, 500000, bins, min, max); 

  TCanvas* can = new TCanvas(); 
  can -> Print("Landau.pdf["); 
  PlotHists(Lan[0], can);
  can -> Print("Landau.pdf"); 
  can -> Print("Landau.pdf]"); 
}

void PlotGaussian()
{
  TH1F* Gaus = Gaussian(0, 1, 1000, -5, 5);
  
  TCanvas* can = new TCanvas();
  can -> Print("Gaussian.pdf["); 
  Gaus -> Draw("SAMEHIST"); 
  can -> Update();  
  can -> Print("Gaussian.pdf"); 
  can -> Print("Gaussian.pdf]"); 
}

void PlotGaussianXGaussian()
{
  int bins = 100; 
  float stdev = 0.5; 
  float mean = -1; 
  float min = -10; 
  float max = 10; 

  // Generate Gaussian histograms 
  TH1F* Gaus1 = Gaussian(mean, stdev, bins, min, max, "1"); 
  TH1F* Gaus2 = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "2");
  
  TCanvas* can = new TCanvas(); 
  can -> Print("Gaussian_Convolution_Test.pdf["); 
  Gaus1 -> Draw("HIST");
  Gaus2 -> Draw("SAMEHIST"); 
  can -> Print("Gaussian_Convolution_Test.pdf"); 
  can -> Clear();

  // Create the Gaussian for the convolution test and the output 
  TH1F* Gaus_Test = Gaussian(mean, stdev, bins, min, max, "Convolution"); 
  TH1F* Conv_Gaus = Gaussian(0, 0.1, bins, min, max, "ConvolutionResult"); 

  // Perform the convolution 
  Convolution(Gaus_Test, Gaus_Test, Conv_Gaus); 

  // Normalize the histograms 
  Normalize(Conv_Gaus);
  Normalize(Gaus2); 

  Gaus2 -> Draw("HIST");
  Conv_Gaus -> Draw("SAMEHIST"); 
  can -> Print("Gaussian_Convolution_Test.pdf"); 
  
  std::cout << "########## Analytical Truth Parameters #########" << std::endl;
  std::cout << "# Parameters of Gaussian 1-> Mean: " << mean << " Standard Deviation: " << stdev << std::endl; 
  std::cout << "# Parameters of Gaussian 2 (Our Target)-> Mean: " << mean + mean << " Standard Deviation: " << stdev*sqrt(2) << std::endl; 
  std::cout << "################################################" << std::endl;
  std::cout << std::endl;
  std::cout << "######### Output Convolution Parameters ########" << std::endl;
  std::cout << "# Parameters of Convolution-> Mean: " << Conv_Gaus -> GetMean() << " Standard Deviation: " << Conv_Gaus -> GetRMS() << std::endl; 
  std::cout << "# Integral: Analytical -> "<< Gaus2 -> Integral() << " Convolved -> " << Conv_Gaus -> Integral() << std::endl;
  can -> Print("Gaussian_Convolution_Test.pdf]"); 
}

void PlotLandauXLandau()
{
  int bins = 1000; 
  float min = 0; 
  float max = 20; 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 
   
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 500000, bins, min, max); 

  TCanvas* can = new TCanvas(); 
  can -> Print("Landau_Convolution.pdf["); 
  can -> SetLogy();
  PlotHists(Gen_Landau, can);  
  can -> Print("Landau_Convolution.pdf"); 
  can -> Clear(); 
  Normalize(Gen_Landau); 

  std::vector<TString> Names_Result = {"Result2", "Result3", "Result4"};  
  std::vector<TH1F*> Results = CloneTH1F(Gen_Landau[0], Names_Result); 
  Convolution(Gen_Landau[0], Gen_Landau[0], Results[0]); 
  Convolution(Results[0], Gen_Landau[0], Results[1]); 
  Convolution(Results[1], Gen_Landau[0], Results[2]); 
 
  // Landau2  
  can -> SetLogy(); 
  RatioPlot(Gen_Landau[1], Results[0], can); 
  can -> Print("Landau_Convolution.pdf"); 
  can -> Clear(); 

  // Landau3
  RatioPlot(Gen_Landau[2], Results[1], can); 
  can -> Print("Landau_Convolution.pdf"); 
  can -> Clear(); 

  // Landau4
  RatioPlot(Gen_Landau[3], Results[2], can); 
  can -> Print("Landau_Convolution.pdf"); 
  can -> Clear(); 
  can -> Print("Landau_Convolution.pdf]"); 
 
  // Compare output with Landau2, Landau3, Landau4 - NOT Landau1!! 
  Stats({Gen_Landau[1], Gen_Landau[2], Gen_Landau[3]}, Results); 
}

void PlotLandauXGaussian()
{
  float min = 0; 
  float max = 20; 
  int bins = 100; 
  
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1}; 
   
  std::vector<TString> Names = {"Landau1"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 500000, bins, min, max); 
  Normalize(Gen_Landau);

  TH1F* Gaus = Gaussian(0, 0.1, bins, min-5, max-5, "1"); 
   
  std::vector<TString> Con_Names = {"GaussianXLandau1"};
  std::vector<TH1F*> Gaus_Landau = Landau(Con_Names, COMP, LandauParams, 500000, bins, min, max); 

  ConvolutionExperimental(Gen_Landau[0], Gaus, Gaus_Landau[0]); 
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  Gaus_Landau[0] -> Draw("HIST");
  Gaus -> Draw("HISTSAME"); 
  Gen_Landau[0] -> Draw("SAMEHIST");  
  can -> Print("LandauXGaussian.pdf"); 
}

void PlotDeconvGausXGaus()
{
  float mean = 0; 
  float stdev = 0.5;
  int bins = 250; 
  float min = -8; 
  float max = 8; 

  TH1F* Gaus_Target = Gaussian(mean, stdev, bins, min, max, "Deconvolution_Target"); 
  TH1F* Gaus_Start = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "Source_Distribution");
  
  TH1F* Gaus_Solution = (TH1F*)Gaus_Target -> Clone("Gaus_Solution"); 
  Gaus_Solution -> Reset();   
  Gaus_Solution -> SetTitle("Gaus_Solution"); 
  
  std::cout << "###################### Analytical Gaussians #######################" << std::endl;
  std::cout << "Source Gaussian -> Mean: " << mean+mean << " Standard Deviation: " << stdev*sqrt(2) << std::endl;
  std::cout << "Target Gaussian -> Mean: " << mean << " Standard Deviation: " << stdev << std::endl;
  std::cout << std::endl;   

  Deconvolution(Gaus_Start, Gaus_Target, Gaus_Solution, 100); 

  // Normalize Solution
  int m = Gaus_Solution -> GetMaximumBin(); 
  Gaus_Solution -> Scale(1 / Gaus_Solution -> GetBinContent(m)); 
  float mean_g = Gaus_Solution -> GetMean(); 
  float stdev_g = Gaus_Solution -> GetRMS(); 
  
  std::cout << "###################### Deconvolution Gaussians #####################" << std::endl;
  std::cout << "Deconvolved Gaussian -> Mean: " << mean_g << " Standard Deviation: " << stdev_g << std::endl;
  std::cout << std::endl;

  Stats({Gaus_Solution}, {Gaus_Target}); 
  
  TCanvas* can = new TCanvas();
  Gaus_Solution -> Draw("HIST*"); 
  //Gaus_Start -> Draw("SAMEHIST"); 
  Gaus_Target -> Draw("SAMEHIST"); 
  can -> Print("DeconvGausXGaus.pdf"); 
}

void PlotDeconvLandauXLandau()
{

  std::cout << "###################### Deconvolution Landaus #######################" << std::endl;

  int bins = 200; 
  float min = 0; 
  float max = 20;
  float centering = (max-min)/float(bins);
  int Iters = 200;
  float width = (max - min)/float(bins); 
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 
   
  std::vector<TString> Names = {"Landau1_Dec", "Landau2_Dec", "Landau3_Dec", "Landau4_Dec"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 

  std::vector<TString> Name_Results = {"Result1", "Result2", "Result3"};
  std::vector<TH1F*> Results = MakeTH1F(Name_Results, bins, min-width/2, max-width/2, "_Deconvolution");
  
  std::vector<TString> Name_Converge = {"Convergence1", "Convergence2", "Convergence3"};
  std::vector<TH1F*> Converge_H = MakeTH1F(Name_Converge, Iters, 0, Iters);

  // Start with Landau 4 and deconvolve with Landau1
  std::vector<float> Converge3 = Deconvolution(Gen_Landau[3], Gen_Landau[0], Results[2], Iters);
  std::vector<float> Converge2 = Deconvolution(Gen_Landau[2], Gen_Landau[0], Results[1], Iters);
  std::vector<float> Converge1 = Deconvolution(Gen_Landau[1], Gen_Landau[0], Results[0], Iters);  
 
  // Convert the convergence vectors to TH1Fs
  ToTH1F(Converge1, Converge_H[0]);  
  ToTH1F(Converge2, Converge_H[1]);
  ToTH1F(Converge3, Converge_H[2]);

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print("DeconvLandauXLandau.pdf["); 

  PlotHists(Gen_Landau, can);
  can -> Print("DeconvLandauXLandau.pdf"); 
  can -> Clear(); 

  for (int i(0); i < 3; i++)
  {
    Gen_Landau[i] -> SetLineStyle(kDashed);
    RatioPlot(Results[i], Gen_Landau[i], can);  
    can -> Print("DeconvLandauXLandau.pdf"); 
    can -> Clear(); 

    Converge_H[i] -> Draw("HIST*"); 
    can -> Print("DeconvLandauXLandau.pdf"); 
    can -> Clear(); 
  }

  can -> Print("DeconvLandauXLandau.pdf]"); 

}

void PlotDeconvLandauXGaussian()
{
  
  int bins = 300; 
  float min = -2; 
  float max = 24;
  int Iters = 400; 
  float centering = (max-min)/float(bins);
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 
  Normalize(Gen_Landau); 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(0, 0.25, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(0, 0.25, bins, min, max, "Original2"); 
  TH1F* Gaussian3 = Gaussian(0, 0.25, bins, min, max, "Original3"); 
  TH1F* Gaussian4 = Gaussian(0, 0.25, bins, min, max, "Original4"); 
  std::vector<TH1F*> PSF_Original = {Gaussian1, Gaussian2, Gaussian3, Gaussian4}; 
  Normalize(PSF_Original); 

  // Convolved histograms 
  std::vector<TString> Name_Conv = {"Landau1XGaussian1", "Landau2XGaussian2", "Landau3XGaussian3", "Landau4XGaussian4"};
  std::vector<TH1F*> Convs = MakeTH1F(Name_Conv, bins, min-centering/2, max-centering/2);

  // Convolve the Landaus with the Gaussians 
  Convolution(Gen_Landau[0], Gaussian1, Convs[0]); 
  Convolution(Gen_Landau[1], Gaussian2, Convs[1]); 
  Convolution(Gen_Landau[2], Gaussian3, Convs[2]); 
  Convolution(Gen_Landau[3], Gaussian4, Convs[3]); 

  // Reconstructed PSF
  std::vector<TString> Name_PSF = {"Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4"};
  std::vector<TH1F*> PSF = MakeTH1F(Name_PSF, bins, min-centering/2, max-centering/2, "Reconstructed");

  // Reconstructed Landau
  std::vector<TString> Name_Landau = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Landau_Recon = MakeTH1F(Name_Landau, bins, min-centering/2, max-centering/2, "Reconstructed");

  // Start with reverting the gaussian smearing 
  MultiThreadingDeconvolution(Convs, PSF_Original, Landau_Recon, Iters); 

  // Start return the Gaussian smearing PSF
  MultiThreadingDeconvolution(Convs, Gen_Landau, PSF, Iters);  

  TCanvas* can = new TCanvas(); 
  can -> Print("DeconvLandauXGaussian.pdf[");
  can -> SetLogy(); 

  for (int i(0); i < Gen_Landau.size(); i++)
  {
    // Plot the reverting of the Landau
    //Normalize(Landau_Recon[0]); 
    RatioPlot(Landau_Recon[i], Gen_Landau[i],  can); 
    can -> Print("DeconvLandauXGaussian.pdf");
    can -> Clear(); 
    
    // Plot the Gaussian PSF
    RatioPlot(PSF[i], PSF_Original[i], can); 
    can -> Print("DeconvLandauXGaussian.pdf");
    can -> Clear(); 
  }
  can -> Print("DeconvLandauXGaussian.pdf]");
}

void PlotGaussianDeconvolutionFit()
{
  int bins = 100; 
  float min = -8; 
  float max = 8; 
  float mean = 0; 
  float stdev = 1; 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(mean, stdev, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "Target"); 
  for (int i(0); i < bins; i++)
  {
    Gaussian2 -> SetBinError(i+1, 1e-9); 
  }

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-1}; 
  Params["m_e"] = {1}; 
  Params["s_s"] = {0.1}; 
  Params["s_e"] = {2}; 
  Params["x_range"] = {-8, 8}; 
  std::vector<TH1F*> H = FitDeconvolution(Gaussian2, {Gaussian1}, Params, 1000, 10000); 
  TCanvas* can = new TCanvas(); 
  RatioPlot(Gaussian2, H[0], can); 
  can -> Print("GaussianDeconvolutionFit.pdf"); 
  Stats({Gaussian2}, {H[0]}); 
}

void PlotLandauXGausFit()
{
  int bins = 1000; 
  float min = -2; 
  float max = 20; 
  float mean = 0; 
  float stdev = 0.5; 

  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 
  TH1F* Gaussian1 = Gaussian(mean, stdev, bins, min, max, "Original1"); 
  
  // Data 
  TH1F* FakeData = (TH1F*)Gen_Landau[0] -> Clone("Data"); 
  FakeData -> SetTitle("Data"); 
  FakeData -> Reset(); 
  FakeData -> Add(Gen_Landau[0]);
  Convolution(FakeData, Gaussian1, FakeData);  
   
  for (int i(0); i < bins; i++)
  {
    FakeData -> SetBinError(i+1, 1e-9); 
  }

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-1, -1, -1, -1}; 
  Params["m_e"] = {1, 1, 1, 1}; 
  Params["s_s"] = {0.001, 0.001, 0.001, 0.001}; 
  Params["s_e"] = {2, 2, 2, 2}; 
  Params["x_range"] = {0, 18}; 
  std::vector<TH1F*> H = FitDeconvolution(FakeData, Gen_Landau, Params, 1000, 10000); 
  TCanvas* can = new TCanvas(); 
  RatioPlot(FakeData, H[0], can); 
  can -> Print("LandauXGaussianFit.pdf"); 
  Stats({FakeData}, {H[0]}); 
}

void PlotNLandauXNGausFit()
{
  
  int bins = 1000; 
  float min = -4; 
  float max = 22;
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {0.6, 0.2, 0.1, 0.1}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(0, 0.1, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(0, 0.1, bins, min, max, "Original2"); 
  TH1F* Gaussian3 = Gaussian(0, 0.1, bins, min, max, "Original3"); 
  TH1F* Gaussian4 = Gaussian(0, 0.1, bins, min, max, "Original4"); 
  std::vector<TH1F*> PSF_Original = {Gaussian1, Gaussian2, Gaussian3, Gaussian4}; 

  // Create the fake data set histograms 
  std::vector<TString> T_N = {"Trk1", "Trk2", "Trk3", "Trk4"}; 
  std::vector<TH1F*> D = CloneTH1F(Gaussian1, T_N); 

  // Convolve the Generated Landaus
  Convolution(Gen_Landau[0], PSF_Original[0], D[0]); 
  Convolution(Gen_Landau[1], PSF_Original[1], D[1]); 
  Convolution(Gen_Landau[2], PSF_Original[2], D[2]); 
  Convolution(Gen_Landau[3], PSF_Original[3], D[3]); 

  // Stack theses conolved histgrams 
  std::vector<TString> N = {"FakeData"}; 
  TH1F* FakeData = CloneTH1F(Gaussian1, N)[0]; 
  FakeData -> Add(D[0]); 
  FakeData -> Add(D[1]); 
  FakeData -> Add(D[2]); 
  FakeData -> Add(D[3]); 

  for (int i(0); i < bins; i++)
  {
    FakeData -> SetBinError(i+1, 1e-9); 
    Gen_Landau[0] -> SetBinError(i+1, 1e-9); 
    Gen_Landau[1] -> SetBinError(i+1, 1e-9); 
    Gen_Landau[2] -> SetBinError(i+1, 1e-9); 
    Gen_Landau[3] -> SetBinError(i+1, 1e-9); 
  }

  // Plot the truth 
  TCanvas* can = new TCanvas();  
  can -> Print("NLandauXNGausFit.pdf["); 
  can -> SetLogy(); 
  PlotHists(FakeData, D, can); 
  can -> Print("NLandauXNGausFit.pdf");
  
  // Perform the fit  
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.1, -0.1, -0.1, -0.1}; 
  Params["m_e"] = {0.1, 0.1, 0.1, 0.1}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {1, 1, 1, 1};  
  Params["x_range"] = {-0.4, 20}; 
  std::vector<TH1F*> Results = FitDeconvolution(FakeData, Gen_Landau, Params, 10000, 100);
  
  // Plot the output and the results
  PlotHists(FakeData, D, Results, can); 
  can -> Print("NLandauXNGausFit.pdf");

  // Create ratio plots of the truth vs the prediction
  PlotAndBenchmark(D, Results, "NLandauXNGausFit.pdf", can);
  can -> Print("NLandauXNGausFit.pdf]"); 
}

void PlotDeconvolutionFit()
{
  
  int bins  = 1000; 
  float min = -2;  //-4.02; 
  float max = 18;  //22;
  int Iters = 100; 
  float centering = (max-min)/float(bins);
  std::cout << centering << std::endl; 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {0.6, 0.3, 0.05, 0.05}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(0, 0.075, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(0, 0.075, bins, min, max, "Original2"); 
  TH1F* Gaussian3 = Gaussian(0, 0.075, bins, min, max, "Original3"); 
  TH1F* Gaussian4 = Gaussian(0, 0.075, bins, min, max, "Original4"); 
  std::vector<TH1F*> PSF_Original = {Gaussian1, Gaussian2, Gaussian3, Gaussian4}; 

  // Create a fake dataset by overlaying multiple Landaus 
  TH1F* Data = new TH1F("Data", "Data", bins, min-centering/2, max - centering/2); 
  Data -> Add(Gen_Landau[0]); 
  Data -> Add(Gen_Landau[1]); 
  Data -> Add(Gen_Landau[2]); 
  Data -> Add(Gen_Landau[3]); 

  for (int i(0); i < bins; i++)
  {
    Data -> SetBinError(i+1, 1e-9); 
    //Gen_Landau[0] -> SetBinError(i+1, 1e-9); 
    //Gen_Landau[1] -> SetBinError(i+1, 1e-9); 
    //Gen_Landau[2] -> SetBinError(i+1, 1e-9); 
    //Gen_Landau[3] -> SetBinError(i+1, 1e-9); 
  }

  // Deconvolved PDFs
  std::vector<TString> Name_Conv = {"Landau1_D", "Landau2_D", "Landau3_D", "Landau4_D"};
  std::vector<TH1F*> D_PDFs = CloneTH1F(Gen_Landau[0], Name_Conv);
  
  // Deconvolve the "PDFs" with the Gaussian 
  MultiThreadingDeconvolution(Gen_Landau, PSF_Original, D_PDFs, Iters); 

  // Plot the initial deconvolved PDFs with the Gaussians that were used for the deconvolution   
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print("DeconvolutionFit.pdf["); 
  for (int i(0); i < Gen_Landau.size(); i++)
  {
    Gen_Landau[i] -> SetLineColor(kRed); 
    Gen_Landau[i] -> GetYaxis() -> SetRangeUser(0.1, Gen_Landau[i] -> Integral()); 
    Gen_Landau[i] -> Draw("HIST"); 
    D_PDFs[i] -> Draw("SAMEHIST");   
    PSF_Original[i] -> SetLineColor(kRed); 
    PSF_Original[i] -> SetLineStyle(kRed); 
    PSF_Original[i] -> Draw("SAMEHIST"); 
    can -> Print("DeconvolutionFit.pdf"); 
    can -> Clear();  
  }

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.1, -0.1, -0.1, -0.1}; 
  Params["m_e"] = {0.1, 0.1, 0.1, 0.1}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params["x_range"] = {-0.4, 20}; 
  std::vector<TH1F*> Result = FitDeconvolution(Data, D_PDFs, Params, 100000, 1000000);
   
  // Plot the outs and do some statistical analysis
  PlotHists(Data, Gen_Landau, Result, can); 
  can -> Print("DeconvolutionFit.pdf"); 
 
  can -> Clear(); 
  PlotAndBenchmark(Gen_Landau, Result, "DeconvolutionFit.pdf", can, 0, 20); 
  
  can -> Print("DeconvolutionFit.pdf]"); 
}
