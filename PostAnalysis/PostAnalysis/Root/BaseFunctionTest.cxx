#include<PostAnalysis/BaseFunctionTest.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/AlgorithmFunctions.h>

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
  TH1F* Gaus1 = Gaussian(mean, stdev, bins, min, max, "Conv_1"); 
  TH1F* Gaus2 = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "Conv_2");

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

  Deconvolution(Gaus_Start, Gaus_Target, Gaus_Solution, 200); 

  // Normalize Solution
  int m = Gaus_Solution -> GetMaximumBin(); 
  Gaus_Solution -> Scale(1 / Gaus_Solution -> GetBinContent(m)); 
  float mean_g = Gaus_Solution -> GetMean(); 
  float stdev_g = Gaus_Solution -> GetRMS(); 
  
  std::cout << "###################### Deconvolution Gaussians #####################" << std::endl;
  std::cout << "Deconvolved Gaussian -> Mean: " << mean_g << " Standard Deviation: " << stdev_g << std::endl;
  std::cout << std::endl;
  
  Gaus_Target -> Write(); 
  Gaus_Start -> Write(); 
  Gaus_Solution -> Write(); 
}

void TestDeconvLandauXLandau(TFile* F)
{
  F -> mkdir("TestDeconvLandauXLandau"); 
  F -> cd("TestDeconvLandauXLandau"); 

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
  float stdev = 0.1; 
  float mean = 0; 
  TH1F* Gaussian1 = Gaussian(mean, stdev, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(mean, stdev, bins, min, max, "Original2"); 
  TH1F* Gaussian3 = Gaussian(mean, stdev, bins, min, max, "Original3"); 
  TH1F* Gaussian4 = Gaussian(mean, stdev, bins, min, max, "Original4"); 
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
  Params["s_e"] = {0.2, 0.2, 0.2, 0.2};  
  Params["x_range"] = {0, 16}; 
  std::vector<std::pair<TH1F*, std::vector<float>>> Results_Perf = FitDeconvolutionPerformance(FakeData, Gen_Landau, Params, 100000, 100000);

  std::vector<TH1F*> Results; 
  for (int i(0); i < Results_Perf.size(); i++)
  {
    std::pair<TH1F*, std::vector<float>> O = Results_Perf[i]; 
    TH1F* H = O.first; 
    std::vector<float> Pre = O.second; 
    TString name = H -> GetTitle(); 
   
    std::cout << "Truth vs Prediction of histogram: " << name << std::endl;  
    float Inte = H -> Integral();   
    float Lumi = FakeData -> Integral(); 
    std::cout << "Fraction Truth: " << COMP[i] << " Prediction: " << Inte / Lumi << " Error: " << Pre[5] /Lumi << std::endl; 
    std::cout << "Mean Truth: " << mean << " Prediction: " << Pre[0] << " Error: " << Pre[3] << std::endl; 
    std::cout << "Stdev Truth: " << stdev << " Prediction: " << Pre[1] << " Error: " << Pre[4] << std::endl; 
    Results.push_back(H); 
  }


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
  std::vector<TH1F*> Result = FitDeconvolution(Data, D_PDFs, Params, 100000, 1000000);

  // Write to file
  BulkWrite(Gen_Landau); 
  BulkWrite(PSF_Original); 
  BulkWrite(D_PDFs); 
  Data -> Write(); 
  BulkWrite(Result); 
}

void TestComparisonBinCenteringLandauXLandau(TFile* F)
{
  F -> mkdir("TestComparisonBinCenteringLandauXLandau"); 
  F -> cd("TestComparisonBinCenteringLandauXLandau"); 
  
  int bins = 100; 
  float min = -2; 
  float max = 18;
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 
  
  // Create the names for the two Landaus. One being the non bin centered and the other being centered. These are the numerically generated Landaus 
  // Wrong  
  std::vector<TString> Name_W = {"Landau1_W", "Landau2_W", "Landau3_W", "Landau4_W"};
  std::vector<TH1F*> Landau_W = WrongLandau(Name_W, COMP, LandauParams, 500000, bins, min, max); 

  // Normal 
  std::vector<TString> Name_N = {"Landau1", "Landau2", "Landau3", "Landau4"}; 
  std::vector<TH1F*> Landau_N = Landau(Name_N, COMP, LandauParams, 500000, bins, min, max); 
  
  // Create Result Histograms 
  // Wrong  
  std::vector<TString> Name_RW = {"Landau1_RW", "Landau2_RW", "Landau3_RW", "Landau4_RW"};
  std::vector<TH1F*> Landau_RW = CloneTH1F(Landau_W[0], Name_RW); 

  // Normal
  std::vector<TString> Name_RN = {"Landau1_RN", "Landau2_RN", "Landau3_RN", "Landau4_RN"};
  std::vector<TH1F*> Landau_RN = CloneTH1F(Landau_N[0], Name_RN); 

  // Start the Convolution with Landau1
  // Wrong 
  Convolution(Landau_W[0], Landau_W[0], Landau_RW[1]);
  Convolution(Landau_RW[1], Landau_W[0], Landau_RW[2]);
  Convolution(Landau_RW[2], Landau_W[0], Landau_RW[3]);

  // Normal
  Convolution(Landau_N[0], Landau_N[0], Landau_RN[1]);
  Convolution(Landau_RN[1], Landau_N[0], Landau_RN[2]);
  Convolution(Landau_RN[2], Landau_N[0], Landau_RN[3]);

  // Normalize all the histograms
  Normalize(Landau_W); 
  Normalize(Landau_RW); 
  Normalize(Landau_N); 
  Normalize(Landau_RN); 

  // Write to File 
  BulkWrite(Landau_W); 
  BulkWrite(Landau_RW); 
  BulkWrite(Landau_N); 
  BulkWrite(Landau_RN); 
}

void TestOscillationLucyRichardson(TFile* F)
{
  F -> mkdir("TestOscillationLucyRichardson"); 
  F -> cd("TestOscillationLucyRichardson"); 

  float min = -2; 
  float max = 18;
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1}; 
  
  // Normal 
  TH1F* Landau_500 = Landau({"B_500"}, COMP, LandauParams, 500000, 500, min, max)[0]; 
  TH1F* Landau_200 = Landau({"B_200"}, COMP, LandauParams, 500000, 200, min, max)[0]; 
  TH1F* Landau_100 = Landau({"B_100"}, COMP, LandauParams, 500000, 100, min, max)[0]; 
  Landau_500 -> SetTitle("Landau with 500 bins"); 
  Landau_200 -> SetTitle("Landau with 200 bins");  
  Landau_100 -> SetTitle("Landau with 100 bins");
  Landau_500 -> Write(); 
  Landau_200 -> Write(); 
  Landau_100 -> Write(); 
 
  // Gaussians being defined 
  TH1F* G_500_05 = Gaussian(0, 0.3, 500, min, max, "G_1");  
  TH1F* G_500_075 = Gaussian(0, 0.2, 500, min, max, "G_2");
  TH1F* G_500_1 = Gaussian(0, 0.05, 500, min, max, "G_3");
  G_500_05 -> SetTitle("Gaussian with 500 bins - M: 0, S: 0.3");   
  G_500_075 -> SetTitle("Gaussian with 500 bins - M: 0, S: 0.2"); 
  G_500_1 -> SetTitle("Gaussian with 500 bins - M: 0, S: 0.05"); 
  G_500_05 -> Write(); 
  G_500_075 -> Write(); 
  G_500_1 -> Write(); 

  TH1F* G_200_05 = Gaussian(0, 0.2, 200, min, max, "G_4"); 
  TH1F* G_100_05 = Gaussian(0, 0.2, 100, min, max, "G_5"); 
  G_200_05 -> SetTitle("Gaussian with 200 bins - M: 0, S: 0.2"); 
  G_100_05 -> SetTitle("Gaussian with 100 bins - M: 0, S: 0.2"); 
  G_200_05 -> Write(); 
  G_100_05 -> Write(); 

  // Perform the Deconvolution with different Gaussian Widths 
  std::vector<TString> Result_500_G_N = {"B_500_G_05_I_100_x", "B_500_G_075_I_100", "B_500_G_1_I_100"};
  std::vector<TH1F*> Results_500_G = CloneTH1F(Landau_500, Result_500_G_N);  
  MultiThreadingDeconvolution({Landau_500, Landau_500, Landau_500}, {G_500_05, G_500_075, G_500_1}, Results_500_G, 100); 
  Results_500_G[0] -> SetTitle("G 500 bins - M: 0, S: 0.3, 100 I");
  Results_500_G[1] -> SetTitle("G 500 bins - M: 0, S: 0.2, 100 I");
  Results_500_G[2] -> SetTitle("G 500 bins - M: 0, S: 0.05, 100 I");
  Results_500_G[0] -> Write(); 
  Results_500_G[1] -> Write(); 
  Results_500_G[2] -> Write(); 

  // Perform the Deconvolution with different number of iterations 
  std::vector<TString> Result_500_I_N = {"B_500_G_05_I_100", "B_500_G_05_I_200", "B_500_G_05_I_300"};
  std::vector<TH1F*> Results_500_I = CloneTH1F(Landau_500, Result_500_I_N);  
  Deconvolution(Landau_500, G_500_05, Results_500_I[0], 100); 
  Deconvolution(Landau_500, G_500_05, Results_500_I[1], 200); 
  Deconvolution(Landau_500, G_500_05, Results_500_I[2], 300); 
  Results_500_I[0] -> SetTitle("G 500 bins - M: 0, S: 0.3, 100 I");  
  Results_500_I[1] -> SetTitle("G 500 bins - M: 0, S: 0.3, 200 I"); 
  Results_500_I[2] -> SetTitle("G 500 bins - M: 0, S: 0.3, 300 I"); 
  Results_500_I[0] -> Write(); 
  Results_500_I[1] -> Write(); 
  Results_500_I[2] -> Write(); 

  // Perform the Deconvolution with different number of bins 
  TH1F* Results_B_100 = CloneTH1F(Landau_100, {"B_100_G_05_I_100_b"})[0];  
  TH1F* Results_B_200 = CloneTH1F(Landau_200, {"B_200_G_05_I_100_b"})[0];  
  TH1F* Results_B_500 = CloneTH1F(Landau_500, {"B_500_G_05_I_100_b"})[0];  
  MultiThreadingDeconvolution({Landau_100, Landau_200, Landau_500}, {G_100_05, G_200_05, G_500_075}, {Results_B_100, Results_B_200, Results_B_500}, 100);  

  Results_B_100 -> SetTitle("G 100 bins - M: 0, S: 0.2, 100 I");
  Results_B_200 -> SetTitle("G 200 bins - M: 0, S: 0.2, 100 I");
  Results_B_500 -> SetTitle("G 500 bins - M: 0, S: 0.2, 100 I");
  Results_B_100 -> Write(); 
  Results_B_200 -> Write(); 
  Results_B_500 -> Write(); 
}

void TestAlgorithm(TFile* F)
{

  F -> Write(); 

  auto Sum_Hist =[] (std::vector<TH1F*> Hists, TString Name)
  {
    TH1F* Hist = (TH1F*)Hists[0] -> Clone(Name); 
    Hist -> Reset();  
    Hist -> SetTitle(Name); 
    for (int i(0); i < Hists.size(); i++)
    {
      Hist -> Add(Hists[i]); 
    }
    return Hist; 
  };

  auto Write_To_File =[&] (std::map<TString, std::vector<TH1F*>> Map)
  {
    typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
    for (it i = Map.begin(); i != Map.end(); i++)
    {
      std::vector<TH1F*> H = i ->second; 
      BulkWrite(H); 
      BulkDelete(H); 
    }
  }; 

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.001, -0.001, -0.001, -0.001}; 
  Params["m_e"] = {0.001, 0.001, 0.001, 0.001}; 
  Params["s_s"] = {0.04, 0.04, 0.04, 0.04};
  Params["s_e"] = {0.15, 0.15, 0.15, 0.15};  
  Params["x_range"] = {0.01, 13.}; 
  Params["iterations"] = {100}; 
  Params["LR_iterations"] = {100}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.075, 0.075, 0.075, 0.075}; 
  Params["cache"] = {10000}; 

  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarloLayerEnergy(Dir); 
  F -> ReOpen("UPDATE"); 
  F -> mkdir("TestAlgorithm"); 
  F -> cd("TestAlgorithm"); 

  std::vector<TH1F*> Track1 = MC["Track_1_All"];
  std::vector<TH1F*> Track2 = MC["Track_2_All"];
  std::vector<TH1F*> Track3 = MC["Track_3_All"];
  std::vector<TH1F*> Track4 = MC["Track_4_All"];
  TH1F* Trk1 = Sum_Hist(Track1, "trk1_data"); 
  TH1F* Trk2 = Sum_Hist(Track2, "trk2_data"); 
  TH1F* Trk3 = Sum_Hist(Track3, "trk3_data"); 
  TH1F* Trk4 = Sum_Hist(Track4, "trk4_data"); 
  std::vector<TH1F*> Data = {Trk1, Trk2, Trk3, Trk4}; 
  
  std::map<TString, std::vector<TH1F*>> Trk1_Measurements = MainAlgorithm(Data, Params, 0); 
  Write_To_File(Trk1_Measurements);  
  
  std::map<TString, std::vector<TH1F*>> Trk2_Measurements = MainAlgorithm(Data, Params, 1); 
  Write_To_File(Trk2_Measurements);  
  
  std::map<TString, std::vector<TH1F*>> Trk3_Measurements = MainAlgorithm(Data, Params, 2); 
  Write_To_File(Trk3_Measurements);  

  std::map<TString, std::vector<TH1F*>> Trk4_Measurements = MainAlgorithm(Data, Params, 3); 
  Write_To_File(Trk4_Measurements);  
}

void TestAlgorithmJetEnergy(TFile* F)
{
  auto Sum_Hist =[] (std::vector<TH1F*> Hists, TString Name)
  {
    TH1F* Hist = (TH1F*)Hists[0] -> Clone(Name); 
    Hist -> Reset();  
    Hist -> SetTitle(Name); 
    for (int i(0); i < Hists.size(); i++)
    {
      Hist -> Add(Hists[i]); 
    }
    return Hist; 
  };

  F -> Write(); 
  std::vector<TString> Layers = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.001, -0.001, -0.001, -0.001}; 
  Params["m_e"] = {0.001, 0.001, 0.001, 0.001}; 
  Params["s_s"] = {0.025, 0.025, 0.025, 0.025};
  Params["s_e"] = {0.075, 0.075, 0.075, 0.075};  
  Params["x_range"] = {0.01, 9.8}; 
  Params["iterations"] = {100}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 

  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarloLayerEnergy(Dir); 
  F -> ReOpen("UPDATE"); 
  F -> mkdir("TestAlgorithmJetEnergy"); 
  F -> cd("TestAlgorithmJetEnergy"); 

  typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
  for (int i(0); i < JetEnergy.size(); i++)
  {
    TString Jet = JetEnergy[i]; 
    TString b1 = "Track_1_" + Jet; 
    TString b2 = "Track_2_" + Jet; 
    TString b3 = "Track_3_" + Jet; 
    TString b4 = "Track_4_" + Jet; 
    TH1F* Track1 = Sum_Hist(MC[b1], b1); 
    TH1F* Track2 = Sum_Hist(MC[b2], b2); 
    TH1F* Track3 = Sum_Hist(MC[b3], b3); 
    TH1F* Track4 = Sum_Hist(MC[b4], b4); 
    std::vector<TH1F*> Data = {Track1, Track2, Track3, Track4};  
    
    // Check if the Track data hists actually have data 
    float L1 = Track1 -> Integral(); 
    float L2 = Track2 -> Integral(); 
    float L3 = Track3 -> Integral(); 
    float L4 = Track4 -> Integral(); 
    
    if (L1 == 0 || L2 == 0 || L3 == 0 || L4 == 0){ continue; }
    TString dir = "TestAlgorithmJetEnergy/"+ Jet; 
    F -> mkdir(dir); 
    F -> cd(dir); 
    std::map<TString, std::vector<TH1F*>> Trk1_Measurements = MainAlgorithm(Data, Params, 0); 
    it x1 = Trk1_Measurements.end(); 
    std::vector<TH1F*> H1 = Trk1_Measurements[x1 -> first]; 
    BulkWrite(H1); 
    BulkWrite(MC[b1]); 

    std::map<TString, std::vector<TH1F*>> Trk2_Measurements = MainAlgorithm(Data, Params, 1); 
    it x2 = Trk2_Measurements.end(); 
    std::vector<TH1F*> H2 = Trk2_Measurements[x2 -> first]; 
    BulkWrite(H2); 
    BulkWrite(MC[b2]); 

    std::map<TString, std::vector<TH1F*>> Trk3_Measurements = MainAlgorithm(Data, Params, 2); 
    it x3 = Trk3_Measurements.end(); 
    std::vector<TH1F*> H3 = Trk3_Measurements[x3 -> first]; 
    BulkWrite(H3); 
    BulkWrite(MC[b3]); 

    std::map<TString, std::vector<TH1F*>> Trk4_Measurements = MainAlgorithm(Data, Params, 3); 
    it x4 = Trk4_Measurements.end(); 
    std::vector<TH1F*> H4 = Trk4_Measurements[x4 -> first]; 
    BulkWrite(H4); 
    BulkWrite(MC[b4]); 

    F -> cd("TestAlgorithmJetEnergy"); 
  }

}

void TestReadFile(TFile* F)
{
  F -> mkdir("TestReadFile"); 
  F -> cd("TestReadFile"); 

  TH1F* H = new TH1F("", "", 100, 0, 1); 
  H -> Write(); 
}

void TestReadFileTrackEnergy(TFile* F)
{
  F -> mkdir("TestReadFileTrackEnergy"); 
  F -> cd ("TestReadFileTrackEnergy"); 
  TH1F* H = new TH1F("", "", 100, 0, 1); 
  H -> Write(); 
}


void TestMonteCarloMatchConvolution(TFile* F)
{ 
  F -> Write();
  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarlo(Dir); 

  F -> ReOpen("UPDATE"); 
  F -> mkdir("TestMonteCarloMatchConvolution"); 
  F -> cd("TestMonteCarloMatchConvolution"); 

  // Get the 1 Track distribution
  std::vector<TH1F*> trk1 = MC["trk1_All"]; 

  // Get the 2 Track distribution
  std::vector<TH1F*> trk2 = MC["trk2_All"]; 

  // Get the 3 Track distribution
  std::vector<TH1F*> trk3 = MC["trk3_All"]; 

  // Get the 4 Track distribution
  std::vector<TH1F*> trk4 = MC["trk4_All"]; 

  // Test in Two cases:
  // - How well does the convolution of 1 track match the other tracks 
  // - How well does the convolution + Fitting of 1 track match other tracks 
  // Is there an improvement?
  
  TH1F* trk1_tru1 = trk1[0]; 
  TH1F* trk2_tru2 = trk2[1]; 
  TH1F* trk3_tru3 = trk3[2]; 
  TH1F* trk4_tru4 = trk4[3]; 

  // Parameters of the histograms 
  int bins = trk1_tru1 -> GetNbinsX(); 
  float min = trk1_tru1 -> GetXaxis() -> GetXmin(); 
  float max = trk1_tru1 -> GetXaxis() -> GetXmax(); 
  float width = (max - min)/float(bins); 
  min+=width/2.;
  max+=width/2.; 

  // Just convolution 
  std::vector<TH1F*> ntrks_Conv = ConvolveNTimes(trk1_tru1, 4, "Conv"); 
  Normalize(ntrks_Conv); 

  // Convolution + Deconvolution + Fit 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.01, -0.01, -0.01, -0.01}; 
  Params["m_e"] = {0.01, 0.01, 0.01, 0.01}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params["x_range"] = {0.1, 9.8}; 

  TH1F* Gaus = Gaussian(0, 0.4, bins, min, max, "Original1");  
  std::vector<TH1F*> PSF = {Gaus, Gaus, Gaus, Gaus};  
  std::vector<TH1F*> ntrks_F = ConvolveNTimes(trk1_tru1, 4, "Fit");
  
  std::vector<TString> Names_De ={"trk1_Deconv", "trk2_Deconv", "trk3_Deconv", "trk4_Deconv"}; 
  std::vector<TH1F*> PDF_D = CloneTH1F(trk1_tru1, Names_De); 
  MultiThreadingDeconvolution(ntrks_F, PSF, PDF_D, 200); 



  std::vector<std::pair<TH1F*, std::vector<float>>> trk1_Fit = FitDeconvolutionPerformance(trk1_tru1, {PDF_D[0]}, Params, 100000, 100000);
  std::pair<TH1F*, std::vector<float>> trk1_pair = trk1_Fit[0]; 
  TH1F* trk1_F = trk1_pair.first; 
  std::vector<float> Perf_trk1 = trk1_pair.second; 
  TH1F* trk1_Fit_ = (TH1F*)trk1_F -> Clone("trk1_Fit"); 
  trk1_Fit_ -> SetTitle("TRK_1_FIT"); 
  delete trk1_F; 
 
  std::vector<std::pair<TH1F*, std::vector<float>>> trk2_Fit = FitDeconvolutionPerformance(trk2_tru2, {PDF_D[1]}, Params, 100000, 100000);
  std::pair<TH1F*, std::vector<float>> trk2_pair = trk2_Fit[0]; 
  TH1F* trk2_F = trk2_pair.first; 
  std::vector<float> Perf_trk2 = trk2_pair.second; 
  TH1F* trk2_Fit_ = (TH1F*)trk2_F -> Clone("trk2_Fit"); 
  trk2_Fit_ -> SetTitle("TRK_2_FIT"); 
  delete trk2_F; 
 

  std::vector<std::pair<TH1F*, std::vector<float>>> trk3_Fit = FitDeconvolutionPerformance(trk3_tru3, {PDF_D[2]}, Params, 100000, 100000);
  std::pair<TH1F*, std::vector<float>> trk3_pair = trk3_Fit[0]; 
  TH1F* trk3_F = trk3_pair.first; 
  std::vector<float> Perf_trk3 = trk3_pair.second; 
  TH1F* trk3_Fit_ = (TH1F*)trk3_F -> Clone("trk3_Fit"); 
  trk3_Fit_ -> SetTitle("TRK_3_FIT"); 
  delete trk3_F; 
 

  std::vector<std::pair<TH1F*, std::vector<float>>> trk4_Fit = FitDeconvolutionPerformance(trk4_tru4, {PDF_D[3]}, Params, 100000, 100000);
  std::pair<TH1F*, std::vector<float>> trk4_pair = trk4_Fit[0]; 
  TH1F* trk4_F = trk4_pair.first; 
  std::vector<float> Perf_trk4 = trk4_pair.second; 
  TH1F* trk4_Fit_ = (TH1F*)trk4_F -> Clone("trk4_Fit"); 
  trk4_Fit_ -> SetTitle("TRK_4_FIT"); 
  delete trk4_F; 


  trk1_tru1 -> Write();
  trk2_tru2 -> Write(); 
  trk3_tru3 -> Write();
  trk4_tru4 -> Write(); 

  BulkWrite(ntrks_Conv);  
 
  Gaus -> Write();
  BulkWrite(ntrks_F); 
  trk1_Fit_ -> Write();  
  trk2_Fit_ -> Write(); 
  trk3_Fit_ -> Write(); 
  trk4_Fit_ -> Write(); 


  std::cout << "############################################" << std::endl;
  std::cout << "/// Track 1 " << std::endl; 
  std::cout << "Normalization: " << Perf_trk1[2] << " Error: " << Perf_trk1[5] / trk1_Fit_ -> Integral() << std::endl;

  std::cout << "/// Track 2 " << std::endl; 
  std::cout << "Normalization: " << Perf_trk2[2] << " Error: " << Perf_trk2[5] / trk2_Fit_ -> Integral() << std::endl;
 
  std::cout << "/// Track 3 " << std::endl; 
  std::cout << "Normalization: " << Perf_trk3[2] << " Error: " << Perf_trk3[5] / trk3_Fit_ -> Integral() << std::endl;

  std::cout << "/// Track 4 " << std::endl; 
  std::cout << "Normalization: " << Perf_trk4[2] << " Error: " << Perf_trk4[5] / trk4_Fit_ -> Integral() << std::endl;



}

void TestMonteCarloFit(TFile* F)
{
  auto Fake_Data = [](std::vector<TH1F*> Data, TString Name)
  {
    TH1F* D = (TH1F*)Data[0] -> Clone(Name); 
    D -> Reset(); 
    D -> SetTitle(Name); 
    for (int i(0); i < Data.size(); i++){D -> Add(Data[i]);}
    return D;  
  }; 

  auto CopyAll = [](std::vector<TH1F*> Fits, std::vector<TString> Names)
  {
    std::vector<TH1F*> Out; 
    for ( int i(0); i < Fits.size(); i++)
    {
      TH1F* H = (TH1F*)Fits[i] -> Clone(Names[i]); 
      H -> SetTitle(Names[i]); 
      delete Fits[i]; 
      Out.push_back(H); 
    }
    return Out;
  };
 

  F -> Write();
  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarlo(Dir); 

  F -> ReOpen("UPDATE"); 
  F -> mkdir("TestMonteCarloFit"); 
  F -> cd("TestMonteCarloFit"); 

  // =================== Preparation of data ======================== //
  // Get the 1 Track distribution
  std::vector<TH1F*> trk1 = MC["trk1_All"]; 

  // Get the 2 Track distribution
  std::vector<TH1F*> trk2 = MC["trk2_All"]; 

  // Get the 3 Track distribution
  std::vector<TH1F*> trk3 = MC["trk3_All"]; 

  // Get the 4 Track distribution
  std::vector<TH1F*> trk4 = MC["trk4_All"]; 

  // Create the dataset we will fit:
  TH1F* trk1_D = Fake_Data(trk1, "Track-1"); 
  TH1F* trk2_D = Fake_Data(trk2, "Track-2"); 
  TH1F* trk3_D = Fake_Data(trk3, "Track-3"); 
  TH1F* trk4_D = Fake_Data(trk4, "Track-4"); 

  // Start Hist 
  TH1F* trk1_tru1 = trk1[0]; 

  // Parameters of the histograms 
  int bins = trk1_tru1 -> GetNbinsX(); 
  float min = trk1_tru1 -> GetXaxis() -> GetXmin(); 
  float max = trk1_tru1 -> GetXaxis() -> GetXmax(); 
  float width = (max - min)/float(bins); 
  min+=width/2.;
  max+=width/2.; 

  // Set Error small so the fit performs better
  //for (int i(0); i < trk1_D -> GetNbinsX(); i++){trk1_D -> SetBinError(i+1, 1e-9);}
  //for (int i(0); i < trk2_D -> GetNbinsX(); i++){trk2_D -> SetBinError(i+1, 1e-9);}
  //for (int i(0); i < trk3_D -> GetNbinsX(); i++){trk3_D -> SetBinError(i+1, 1e-9);}
  //for (int i(0); i < trk4_D -> GetNbinsX(); i++){trk4_D -> SetBinError(i+1, 1e-9);}
  // ======================= End of Data Preparation ====================== //

  // Just convolution 
  std::vector<TH1F*> ntrks_Conv = ConvolveNTimes(trk1_tru1, 4, "Conv"); 
  Normalize(ntrks_Conv); 

  // Convolution + Deconvolution + Fit 
  std::map<TString, std::vector<float>> Params; 
  float m = 0.1; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params["s_e"] = {0.5, 0.5, 0.5, 0.5};  
  Params["x_range"] = {-1, 11}; 
  float stdev = 0.01;
  float mean = 0; 

  // Create the Raw Convolved PDFs used for the entire algo 
  TH1F* Gaus = Gaussian(mean, stdev, bins, min, max, "Original1");  
  std::vector<TH1F*> PSF = {Gaus, Gaus, Gaus, Gaus};  
  std::vector<TH1F*> ntrks_F = ConvolveNTimes(trk1_tru1, 4, "Fit");
 
  // Perform the deconvolution 
  std::vector<TString> Names_De ={"trk1_Deconv", "trk2_Deconv", "trk3_Deconv", "trk4_Deconv"}; 
  std::vector<TH1F*> PDF_D = CloneTH1F(trk1_tru1, Names_De); 
  MultiThreadingDeconvolution(ntrks_F, PSF, PDF_D, 250); 
  
  // ====== Apply the fit to the data ========= //
  std::vector<std::pair<TH1F*, std::vector<float>>> trk1_fit = FitDeconvolutionPerformance(trk1_D, PDF_D, Params, 1000000, 100000);
  TH1F* trk1_tru1_f = trk1_fit[0].first; 
  TH1F* trk1_tru2_f = trk1_fit[1].first; 
  TH1F* trk1_tru3_f = trk1_fit[2].first; 
  TH1F* trk1_tru4_f = trk1_fit[3].first; 
  std::vector<TH1F*> trk1_Fit = {trk1_tru1_f, trk1_tru2_f, trk1_tru3_f, trk1_tru4_f};  
  
  // Parameter output of fit 
  std::vector<float> trk1_tru1_fout = trk1_fit[0].second; 
  std::vector<float> trk1_tru2_fout = trk1_fit[1].second;
  std::vector<float> trk1_tru3_fout = trk1_fit[2].second;
  std::vector<float> trk1_tru4_fout = trk1_fit[3].second;
  
  std::vector<std::pair<TH1F*, std::vector<float>>> trk2_fit = FitDeconvolutionPerformance(trk2_D, PDF_D, Params, 1000000, 100000);
  TH1F* trk2_tru1_f = trk2_fit[0].first; 
  TH1F* trk2_tru2_f = trk2_fit[1].first; 
  TH1F* trk2_tru3_f = trk2_fit[2].first; 
  TH1F* trk2_tru4_f = trk2_fit[3].first; 
  std::vector<TH1F*> trk2_Fit = {trk2_tru1_f, trk2_tru2_f, trk2_tru3_f, trk2_tru4_f};  
   
  // Parameter output of fit 
  std::vector<float> trk2_tru1_fout = trk2_fit[0].second; 
  std::vector<float> trk2_tru2_fout = trk2_fit[1].second;
  std::vector<float> trk2_tru3_fout = trk2_fit[2].second;
  std::vector<float> trk2_tru4_fout = trk2_fit[3].second;
 
  std::vector<std::pair<TH1F*, std::vector<float>>> trk3_fit = FitDeconvolutionPerformance(trk3_D, PDF_D, Params, 1000000, 100000);
  TH1F* trk3_tru1_f = trk3_fit[0].first; 
  TH1F* trk3_tru2_f = trk3_fit[1].first; 
  TH1F* trk3_tru3_f = trk3_fit[2].first; 
  TH1F* trk3_tru4_f = trk3_fit[3].first; 
  std::vector<TH1F*> trk3_Fit = {trk3_tru1_f, trk3_tru2_f, trk3_tru3_f, trk3_tru4_f}; 

  // Parameter output of fit 
  std::vector<float> trk3_tru1_fout = trk3_fit[0].second; 
  std::vector<float> trk3_tru2_fout = trk3_fit[1].second;
  std::vector<float> trk3_tru3_fout = trk3_fit[2].second;
  std::vector<float> trk3_tru4_fout = trk3_fit[3].second;
 
  std::vector<std::pair<TH1F*, std::vector<float>>> trk4_fit = FitDeconvolutionPerformance(trk4_D, PDF_D, Params, 1000000, 100000);
  TH1F* trk4_tru1_f = trk4_fit[0].first; 
  TH1F* trk4_tru2_f = trk4_fit[1].first; 
  TH1F* trk4_tru3_f = trk4_fit[2].first; 
  TH1F* trk4_tru4_f = trk4_fit[3].first; 
  std::vector<TH1F*> trk4_Fit = {trk4_tru1_f, trk4_tru2_f, trk4_tru3_f, trk4_tru4_f}; 

  // Parameter output of fit 
  std::vector<float> trk4_tru1_fout = trk4_fit[0].second; 
  std::vector<float> trk4_tru2_fout = trk4_fit[1].second;
  std::vector<float> trk4_tru3_fout = trk4_fit[2].second;
  std::vector<float> trk4_tru4_fout = trk4_fit[3].second;

  std::vector<TString> trk1_Names = {"Track_1_Tru_1_Fit", "Track_1_Tru_2_Fit", "Track_1_Tru_3_Fit", "Track_1_Tru_4_Fit"}; 
  std::vector<TH1F*> trk1_F = CopyAll(trk1_Fit, trk1_Names); 
 
  std::vector<TString> trk2_Names = {"Track_2_Tru_1_Fit", "Track_2_Tru_2_Fit", "Track_2_Tru_3_Fit", "Track_2_Tru_4_Fit"}; 
  std::vector<TH1F*> trk2_F = CopyAll(trk2_Fit, trk2_Names); 
 
  std::vector<TString> trk3_Names = {"Track_3_Tru_1_Fit", "Track_3_Tru_2_Fit", "Track_3_Tru_3_Fit", "Track_3_Tru_4_Fit"}; 
  std::vector<TH1F*> trk3_F = CopyAll(trk3_Fit, trk3_Names); 
 
  std::vector<TString> trk4_Names = {"Track_4_Tru_1_Fit", "Track_4_Tru_2_Fit", "Track_4_Tru_3_Fit", "Track_4_Tru_4_Fit"}; 
  std::vector<TH1F*> trk4_F = CopyAll(trk4_Fit, trk4_Names); 

  // Monte Carlo truth
  BulkWrite(trk1); 
  BulkWrite(trk2); 
  BulkWrite(trk3); 
  BulkWrite(trk4); 

  // Stacked data
  trk1_D -> Write(); 
  trk2_D -> Write(); 
  trk3_D -> Write(); 
  trk4_D -> Write(); 
  
  // Convs
  BulkWrite(ntrks_F); 
  BulkWrite(PSF); 
  BulkWrite(PDF_D); 

  // Fits
  BulkWrite(trk1_F); 
  BulkWrite(trk2_F); 
  BulkWrite(trk3_F); 
  BulkWrite(trk4_F); 

  std::cout << "################### Gaussian Fit Parameters ################################## " << std::endl;
  auto PrintStats = [] (std::vector<float> Stats, std::map<TString, std::vector<float>> Params, int index, float mean, float stdev)
  {
    std::cout << "Original Gaussian Parameters -> mean: " << mean << " standard dev: " << stdev << std::endl;
    std::cout << "Range of fit -> mean: " << Params["m_s"][index] << " -> " << Params["m_e"][index] << std::endl;
    std::cout << "Range of fit -> standard dev: " << Params["s_s"][index] << " -> " << Params["s_e"][index] << std::endl;
    std::cout << "Final Fit Output -> mean: " << Stats[0] << " +- " << Stats[3] << " standard dev: " << Stats[1] << " +- " << Stats[4] << std::endl;
    std::cout << "___________________________________________" << std::endl;
    std::cout << std::endl;
  };

  std::cout << "Track 1, Truth 1" << std::endl;
  PrintStats(trk1_tru1_fout, Params, 0, mean, stdev);
  std::cout << "Track 1, Truth 2" << std::endl;
  PrintStats(trk1_tru2_fout, Params, 1, mean, stdev);
  std::cout << "Track 1, Truth 3" << std::endl;
  PrintStats(trk1_tru3_fout, Params, 2, mean, stdev);
  std::cout << "Track 1, Truth 4" << std::endl;
  PrintStats(trk1_tru4_fout, Params, 3, mean, stdev);

  std::cout << "Track 2, Truth 1" << std::endl;
  PrintStats(trk2_tru1_fout, Params, 0, mean, stdev);
  std::cout << "Track 2, Truth 2" << std::endl;
  PrintStats(trk2_tru2_fout, Params, 1, mean, stdev);
  std::cout << "Track 2, Truth 3" << std::endl;
  PrintStats(trk2_tru3_fout, Params, 2, mean, stdev);
  std::cout << "Track 2, Truth 4" << std::endl;
  PrintStats(trk2_tru4_fout, Params, 3, mean, stdev);

  std::cout << "Track 3, Truth 1" << std::endl;
  PrintStats(trk3_tru1_fout, Params, 0, mean, stdev);
  std::cout << "Track 3, Truth 2" << std::endl;
  PrintStats(trk3_tru2_fout, Params, 1, mean, stdev);
  std::cout << "Track 3, Truth 3" << std::endl;
  PrintStats(trk3_tru3_fout, Params, 2, mean, stdev);
  std::cout << "Track 3, Truth 4" << std::endl;
  PrintStats(trk3_tru4_fout, Params, 3, mean, stdev);

  std::cout << "Track 4, Truth 1" << std::endl;
  PrintStats(trk4_tru1_fout, Params, 0, mean, stdev);
  std::cout << "Track 4, Truth 2" << std::endl;
  PrintStats(trk4_tru2_fout, Params, 1, mean, stdev);
  std::cout << "Track 4, Truth 3" << std::endl;
  PrintStats(trk4_tru3_fout, Params, 2, mean, stdev);
  std::cout << "Track 4, Truth 4" << std::endl;
  PrintStats(trk4_tru4_fout, Params, 3, mean, stdev);
}

