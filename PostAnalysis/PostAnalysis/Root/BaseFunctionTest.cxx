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
  std::cout << "# Parameters of Gaussian 2-> Mean: " << mean + mean << " Standard Deviation: " << stdev*sqrt(2) << std::endl; 
  std::cout << "# Peak position of Gaussian 1 (Mean): " << Gaus1 -> GetMean() << " Gaussian 2 (Mean): " << Gaus2 -> GetMean() << std::endl; 
  std::cout << "# Standard Deviation of Gaussian 1 (STDEV): " << Gaus1 -> GetRMS() << " Gaussian 2 (STDEV): " << Gaus2 -> GetRMS() << std::endl; 
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
  Normalize(Results); 
 
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

