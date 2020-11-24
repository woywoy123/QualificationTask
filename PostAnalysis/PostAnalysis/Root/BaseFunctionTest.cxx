#include<PostAnalysis/BaseFunctionTest.h>

void TestLandau(TFile* F)
{
  F -> mkdir("TestLandau"); 
  F -> cd("TestLandau"); 
  int bins = 1000; 
  int min = 0; 
  int max = 20; 

  std::vector<TString> Names = {"Landau"}; 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<TH1F*> Lan = Landau(Names, {1.0}, LandauParams, 500000, bins, min, max); 
  
  // Write to file
  Lan[0] -> Write(); 

}

void TestGaussian(TFile* F)
{
  F -> mkdir("TestGaussian");
  F -> cd("TestGaussian"); 
  TH1F* Gaus = Gaussian(0, 1, 1000, -5, 5);

  // Write to file
  Gaus -> Write();  
}

void TestGaussianXGaussian(TFile* F)
{
  F -> mkdir("TestGaussianXGaussian"); 
  F -> cd("TestGaussianXGaussian"); 

  int bins = 100; 
  float stdev = 0.5; 
  float mean = -1; 
  float min = -10; 
  float max = 10; 

  // Generate Gaussian histograms 
  TH1F* Gaus1 = Gaussian(mean, stdev, bins, min, max, "1"); 
  TH1F* Gaus2 = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "2");

  // Define the titles of the histograms  
  TString n1 = "Analytical - M: "; n1 += (mean); n1 += (" STDEV: "); n1 += (stdev);  
  TString n2 = "Analytical - M: "; n2 += (mean+mean); n2 += (" STDEV: "); n2 += (stdev*sqrt(2));  
  Gaus1 -> SetTitle(n1); 
  Gaus2 -> SetTitle(n2); 
  Normalize(Gaus1); 
  Normalize(Gaus2); 

  // Write to file
  Gaus1 -> Write(); 
  Gaus2 -> Write(); 

  // Create the Gaussian for the convolution test and the output 
  TH1F* Conv_Gaus = (TH1F*)Gaus1 -> Clone("Solution"); 
  Conv_Gaus -> SetTitle("Solution"); 
  
  // Perform the convolution 
  Convolution(Gaus1, Gaus1, Conv_Gaus); 

  // Normalize the histograms 
  Normalize(Conv_Gaus);

  // write to file
  Conv_Gaus -> Write(); 
  
  std::cout << "########## Analytical Truth Parameters #########" << std::endl;
  std::cout << "# Parameters of Gaussian 1-> Mean: " << mean << " Standard Deviation: " << stdev << std::endl; 
  std::cout << "# Parameters of Gaussian 2 (Our Target)-> Mean: " << mean + mean << " Standard Deviation: " << stdev*sqrt(2) << std::endl; 
  std::cout << "################################################" << std::endl;
  std::cout << std::endl;
  std::cout << "######### Output Convolution Parameters ########" << std::endl;
  std::cout << "# Parameters of Convolution-> Mean: " << Conv_Gaus -> GetMean() << " Standard Deviation: " << Conv_Gaus -> GetRMS() << std::endl; 
  std::cout << "# Integral: Analytical -> "<< Gaus2 -> Integral() << " Convolved -> " << Conv_Gaus -> Integral() << std::endl;
}

void TestLandauXLandau(TFile* F)
{
  F -> mkdir("TestLandauXLandau"); 
  F -> cd("TestLandauXLandau"); 
  
  int bins = 1000; 
  float min = 0; 
  float max = 20; 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 
   
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 
  
  std::vector<TString> Names_Result = {"Result2", "Result3", "Result4"};  
  std::vector<TH1F*> Results = CloneTH1F(Gen_Landau[0], Names_Result); 
  Convolution(Gen_Landau[0], Gen_Landau[0], Results[0]); 
  Convolution(Results[0], Gen_Landau[0], Results[1]); 
  Convolution(Results[1], Gen_Landau[0], Results[2]); 

  // Normalize the distributions
  Normalize(Results); 
  Normalize(Gen_Landau); 
 
  // Compare output with Landau2, Landau3, Landau4 - NOT Landau1!! 
  Stats({Gen_Landau[1], Gen_Landau[2], Gen_Landau[3]}, Results); 

  // Write to file
  BulkWrite(Gen_Landau); 
  BulkWrite(Results); 
}

void TestLandauXGaussian(TFile* F)
{
  F -> mkdir("TestLandauXGaussian"); 
  F -> cd("TestLandauXGaussian"); 

  float min = -2; 
  float max = 22; 
  int bins = 1000; 
  float mean = 0; 
  float stdev = 0.1; 
  
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1}; 
  
  // Generate the Landau 
  std::vector<TString> Names = {"Landau1"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 500000, bins, min, max); 
  Normalize(Gen_Landau);

  // Generate the Gaussian
  TH1F* Gaus = Gaussian(mean, stdev, bins, min, max, "1"); 
  TString name = "M: "; name +=(mean); name +=(" StDev: "); name +=(stdev);
  Gaus -> SetTitle(name);  
  
  // Convolve the Gaussian with the Landau 
  std::vector<TString> Con_Names = {"GaussianXLandau1"};
  std::vector<TH1F*> Gaus_Landau = Landau(Con_Names, COMP, LandauParams, 500000, bins, min, max); 
  Convolution(Gen_Landau[0], Gaus, Gaus_Landau[0]); 

  // Write to file 
  Gen_Landau[0] -> Write(); 
  Gaus -> Write(); 
  Gaus_Landau[0] -> Write();  
}

void TestDeconvGausXGaus(TFile* F)
{
  F -> mkdir("TestDeconvGausXGaus");  
  F -> cd("TestDeconvGausXGaus");  
  
  float mean = 0; 
  float stdev = 0.5;
  int bins = 250; 
  float min = -8; 
  float max = 8; 

  TH1F* Gaus_Target = Gaussian(mean, stdev, bins, min, max, "Target"); 
  TH1F* Gaus_Start = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "Source");
  
  TH1F* Gaus_Solution = (TH1F*)Gaus_Target -> Clone("Solution"); 
  Gaus_Solution -> Reset();   
  Gaus_Solution -> SetTitle("Solution"); 
  
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
  
  Gaus_Target -> Write(); 
  Gaus_Start -> Write(); 
  Gaus_Solution -> Write(); 
}

void TestDeconvLandauXLandau(TFile* F)
{
  F -> mkdir("TestDeconvLandauXLandau"); 
  F -> cd("TestDeconvLandauXLandau"); 

  std::cout << "###################### Deconvolution Landaus #######################" << std::endl;
  int bins = 200; 
  float min = 0; 
  float max = 20;
  float centering = (max-min)/float(bins);
  int Iters = 200;
  float width = (max - min)/float(bins); 
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 
   
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 

  std::vector<TString> Name_Results = {"Result1", "Result2", "Result3"};
  std::vector<TH1F*> Results = CloneTH1F(Gen_Landau[0], Name_Results);
  
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

  // Write to file 
  BulkWrite(Gen_Landau); 
  BulkWrite(Results);
  BulkWrite(Converge_H); 
  
  // Evalutate the deconvolution 
  Stats({Gen_Landau[0], Gen_Landau[1], Gen_Landau[2]}, Results); 
}

void TestDeconvLandauXGaussian(TFile* F)
{
  F -> mkdir("TestDeconvLandauXGaussian"); 
  F -> cd("TestDeconvLandauXGaussian"); 
  
  int bins = 200; 
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
  TH1F* Gaussian1 = Gaussian(0, 0.25, bins, min, max, "Gaussian1"); 
  TH1F* Gaussian2 = Gaussian(0, 0.25, bins, min, max, "Gaussian2"); 
  TH1F* Gaussian3 = Gaussian(0, 0.25, bins, min, max, "Gaussian3"); 
  TH1F* Gaussian4 = Gaussian(0, 0.25, bins, min, max, "Gaussian4"); 
  std::vector<TH1F*> PSF_Original = {Gaussian1, Gaussian2, Gaussian3, Gaussian4}; 
  Normalize(PSF_Original); 

  // Convolved histograms 
  std::vector<TString> Name_Conv = {"Landau1XGaussian1", "Landau2XGaussian2", "Landau3XGaussian3", "Landau4XGaussian4"};
  std::vector<TH1F*> Convs = CloneTH1F(Gen_Landau[0], Name_Conv);

  // Convolve the Landaus with the Gaussians 
  Convolution(Gen_Landau[0], Gaussian1, Convs[0]); 
  Convolution(Gen_Landau[1], Gaussian2, Convs[1]); 
  Convolution(Gen_Landau[2], Gaussian3, Convs[2]); 
  Convolution(Gen_Landau[3], Gaussian4, Convs[3]); 

  // Wrtie to file 
  BulkWrite(Gen_Landau); 
  BulkWrite(PSF_Original); 
  BulkWrite(Convs); 

  // Reconstructed PSF
  std::vector<TString> Name_PSF = {"R_Gaus_1", "R_Gaus_2", "R_Gaus_3", "R_Gaus_4"};
  std::vector<TH1F*> PSF = CloneTH1F(Gen_Landau[0], Name_PSF);

  // Reconstructed Landau
  std::vector<TString> Name_Landau = {"R_Landau_1", "R_Landau_2", "R_Landau_3", "R_Landau_4"};
  std::vector<TH1F*> Landau_Recon = CloneTH1F(Gen_Landau[0], Name_Landau);

  // Start with reverting the gaussian smearing 
  MultiThreadingDeconvolution(Convs, PSF_Original, Landau_Recon, Iters); 

  // Start return the Gaussian smearing PSF
  MultiThreadingDeconvolution(Convs, Gen_Landau, PSF, Iters);  

  // Write to file
  BulkWrite(Landau_Recon);  
  BulkWrite(PSF); 
}

void TestGaussianDeconvolutionFit(TFile* F)
{
  F -> mkdir("TestGaussianDeconvolutionFit"); 
  F -> cd("TestGaussianDeconvolutionFit"); 

  int bins = 100; 
  float min = -8; 
  float max = 8; 
  float mean = 1; 
  float stdev = 1; 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(mean, stdev, bins, min, max, "Source"); 
  TH1F* Gaussian2 = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "Target"); 
  for (int i(0); i < bins; i++)
  {
    Gaussian2 -> SetBinError(i+1, 1e-9); 
  }

  // Define the fitting parameters
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-5}; 
  Params["m_e"] = {5}; 
  Params["s_s"] = {0.1}; 
  Params["s_e"] = {4}; 
  Params["x_range"] = {-8, 8}; 
  std::vector<TH1F*> H = FitDeconvolution(Gaussian2, {Gaussian1}, Params, 1000, 10000); 
 
  // Check the statistcs of the fit  
  Stats({Gaussian2}, {H[0]}); 

  Gaussian1 -> Write(); 
  Gaussian2 -> Write(); 
  H[0] -> Write(); 
}

void TestLandauXGausFit(TFile* F)
{
  F -> mkdir("TestLandauXGausFit"); 
  F -> cd("TestLandauXGausFit"); 

  int bins = 500; 
  float min = -2; 
  float max = 18; 
  float mean = 0; 
  float stdev = 0.1; 

  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 
  TH1F* Gaussian1 = Gaussian(mean, stdev, bins, min, max, "1"); 
  
  // Fake Data 
  TH1F* FakeData = (TH1F*)Gen_Landau[0] -> Clone("Data"); 
  FakeData -> SetTitle("Data"); 
  FakeData -> Reset(); 
  FakeData -> Add(Gen_Landau[0]);
  Convolution(FakeData, Gaussian1, FakeData);  
  for (int i(0); i < bins; i++){FakeData -> SetBinError(i+1, 1e-9);}

  // Perform the fit to the Fake data and see if we can revert the Gaussian convolution 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-1, -1, -1, -1}; 
  Params["m_e"] = {1, 1, 1, 1}; 
  Params["s_s"] = {0.001, 0.001, 0.001, 0.001}; 
  Params["s_e"] = {2, 2, 2, 2}; 
  Params["x_range"] = {0, 18}; 
  std::vector<TH1F*> H = FitDeconvolution(FakeData, Gen_Landau, Params, 10000, 100000); 
  Stats({FakeData}, {H[0]}); 

  // Write to file 
  BulkWrite(Gen_Landau); 
  Gaussian1 -> Write();
  FakeData -> Write(); 
  BulkWrite(H); 

}

void TestNLandauXNGausFit(TFile* F)
{
  F -> mkdir("TestNLandauXNGausFit"); 
  F -> cd("TestNLandauXNGausFit"); 
  
  int bins = 500; 
  float min = -2; 
  float max = 18;
 
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

  for (int i(0); i < bins; i++){FakeData -> SetBinError(i+1, 1e-9);}

  // Perform the fit  
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.2, -0.2, -0.2, -0.2}; 
  Params["m_e"] = {0.2, 0.2, 0.2, 0.2}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.5, 0.5, 0.5, 0.5};  
  Params["x_range"] = {-0.4, 16}; 
  std::vector<TH1F*> Results = FitDeconvolution(FakeData, Gen_Landau, Params, 10000, 100000);

  // Write to file 
  BulkWrite(Gen_Landau); 
  BulkWrite(PSF_Original); 
  BulkWrite(D); 
  FakeData -> Write(); 
  BulkWrite(Results); 
}

void TestDeconvolutionFit(TFile* F)
{
  F -> mkdir("TestDeconvolutionFit"); 
  F -> cd("TestDeconvolutionFit"); 

  int bins  = 500; 
  float min = -2;   
  float max = 18;  
  int Iters = 250; 
  float centering = (max-min)/float(bins);
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {0.6, 0.3, 0.05, 0.05}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(0, 0.075, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(0, 0.08, bins, min, max, "Original2"); 
  TH1F* Gaussian3 = Gaussian(0, 0.09, bins, min, max, "Original3"); 
  TH1F* Gaussian4 = Gaussian(0, 0.1, bins, min, max, "Original4"); 
  std::vector<TH1F*> PSF_Original = {Gaussian1, Gaussian2, Gaussian3, Gaussian4}; 

  // Create a fake dataset by overlaying multiple Landaus 
  TH1F* Data = new TH1F("Data", "Data", bins, min-centering/2, max - centering/2); 
  Data -> Add(Gen_Landau[0]); 
  Data -> Add(Gen_Landau[1]); 
  Data -> Add(Gen_Landau[2]); 
  Data -> Add(Gen_Landau[3]); 
  for (int i(0); i < bins; i++){Data -> SetBinError(i+1, 1e-9);}

  // Deconvolved PDFs
  std::vector<TString> Name_Conv = {"Landau1_D", "Landau2_D", "Landau3_D", "Landau4_D"};
  std::vector<TH1F*> D_PDFs = CloneTH1F(Gen_Landau[0], Name_Conv);
  
  // Deconvolve the "PDFs" with the Gaussian 
  MultiThreadingDeconvolution(Gen_Landau, PSF_Original, D_PDFs, Iters); 

  // Perform the fit 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.1, -0.1, -0.1, -0.1}; 
  Params["m_e"] = {0.1, 0.1, 0.1, 0.1}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params["x_range"] = {-0.4, 20}; 
  std::vector<TH1F*> Result = FitDeconvolution(Data, D_PDFs, Params, 10000, 100000);
  Stats(Result, Gen_Landau, 0.4, 16); 

  // Write to file
  BulkWrite(Gen_Landau); 
  BulkWrite(PSF_Original); 
  BulkWrite(D_PDFs); 
  Data -> Write(); 
  BulkWrite(Result); 


}
