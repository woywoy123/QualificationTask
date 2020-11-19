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
  ResidualRemove(Results[0]);  
  ResidualRemove(Results[1]); 
  ResidualRemove(Results[2]); 
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
  float mean = 1; 
  float stdev = 0.5;
  int bins = 100; 
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
  Gaus_Start -> Draw("SAMEHIST"); 
  Gaus_Target -> Draw("SAMEHIST"); 
  can -> Print("DeconvGausXGaus.pdf"); 
}

void PlotDeconvLandauXLandau()
{

  std::cout << "###################### Deconvolution Landaus #######################" << std::endl;

  int bins = 250; 
  float min = 0; 
  float max = 20;
  int Iters = 500; 
  float centering = (max-min)/float(bins);
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 
   
  std::vector<TString> Names = {"Landau1_Dec", "Landau2_Dec", "Landau3_Dec", "Landau4_Dec"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 

  std::vector<TString> Name_Results = {"Result1", "Result2", "Result3"};
  std::vector<TH1F*> Results = MakeTH1F(Name_Results, bins, min-centering/2, max-centering/2, "_Deconvolution");
  
  std::vector<TString> Name_Converge = {"Convergence1", "Convergence2", "Convergence3"};
  std::vector<TH1F*> Converge_H = MakeTH1F(Name_Converge, Iters, 0, Iters);

  // Start with Landau 4 and deconvolve with Landau1
  std::vector<float> Converge3 = Deconvolution(Gen_Landau[3], Gen_Landau[2], Results[1], Iters);
  //std::vector<float> Converge2 = Deconvolution(Gen_Landau[2], Gen_Landau[1], Results[1], Iters);
  //std::vector<float> Converge1 = Deconvolution(Gen_Landau[1], Gen_Landau[0], Results[0], Iters);  
 
  ////// Convert the convergence vectors to TH1Fs
  //ToTH1F(Converge1, Converge_H[0]);  
  //ToTH1F(Converge2, Converge_H[1]);
  ToTH1F(Converge3, Converge_H[1]);

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print("DeconvLandauXLandau.pdf["); 

  PlotHists(Gen_Landau, can);
  can -> Print("DeconvLandauXLandau.pdf"); 
  can -> Clear(); 

  //for (int i(0); i < 3; i++)
  //{
    Gen_Landau[2] -> SetLineStyle(kDashed);
    RatioPlot(Gen_Landau[2], Results[1], can);  
    can -> Print("DeconvLandauXLandau.pdf"); 
    can -> Clear(); 

    Converge_H[1] -> Draw("HIST*"); 
    can -> Print("DeconvLandauXLandau.pdf"); 
    can -> Clear(); 
  //}

  can -> Print("DeconvLandauXLandau.pdf]"); 

  // Note to self: Understand why the deconvolution is behaving unexpectedly. 
}

void PlotDeconvLandauXGaussian()
{
  
  int bins = 300; 
  float min = -2; 
  float max = 20;
  int Iters = 1000; 
  float centering = (max-min)/float(bins);
 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<float> COMP = {1, 1, 1, 1}; 

  // Define the Landaus    
  std::vector<TString> Names = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Gen_Landau = Landau(Names, COMP, LandauParams, 5000000, bins, min, max); 
  Normalize(Gen_Landau); 

  // Define the Gaussians 
  TH1F* Gaussian1 = Gaussian(0, 0.5, bins, min, max, "Original1"); 
  TH1F* Gaussian2 = Gaussian(0, 0.5, bins, min, max, "Original2"); 
  TH1F* Gaussian3 = Gaussian(0, 0.5, bins, min, max, "Original3"); 
  TH1F* Gaussian4 = Gaussian(0, 0.5, bins, min, max, "Original4"); 
  std::vector<TH1F*> PSF_Original = {Gaussian1, Gaussian2, Gaussian3, Gaussian4}; 
  Normalize(PSF_Original); 

  // Convolved histograms 
  std::vector<TString> Name_Conv = {"Landau1XGaussian1", "Landau2XGaussian2", "Landau3XGaussian3", "Landau4XGaussian4"};
  std::vector<TH1F*> Convs = MakeTH1F(Name_Conv, bins, min, max);

  // Convolve the Landaus with the Gaussians 
  Convolution(Gen_Landau[0], Gaussian1, Convs[0]); 
  Convolution(Gen_Landau[1], Gaussian2, Convs[1]); 
  Convolution(Gen_Landau[2], Gaussian3, Convs[2]); 
  Convolution(Gen_Landau[3], Gaussian4, Convs[3]); 

  // Reconstructed PSF
  std::vector<TString> Name_PSF = {"Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4"};
  std::vector<TH1F*> PSF = MakeTH1F(Name_PSF, bins, min, max, "Reconstructed");

  // Reconstructed Landau
  std::vector<TString> Name_Landau = {"Landau1", "Landau2", "Landau3", "Landau4"};
  std::vector<TH1F*> Landau_Recon = MakeTH1F(Name_Landau, bins, min, max, "Reconstructed");

  // Covergence graph 
  std::vector<TString> Name_Converge = {"Convergence1", "Convergence2", "Convergence3", "Convergence4"};
  std::vector<TH1F*> Converge_H = MakeTH1F(Name_Converge, Iters, 0, Iters);


  // Start with reverting the gaussian smearing 
  std::vector<float> Converge1 = Deconvolution(Convs[0], Gaussian1, Landau_Recon[0] , Iters);
  //std::vector<float> Converge2 = Deconvolution(Convs[1], Gaussian2, Landau_Recon[1] , Iters);
  //std::vector<float> Converge3 = Deconvolution(Convs[2], Gaussian3, Landau_Recon[2] , Iters);
  //std::vector<float> Converge4 = Deconvolution(Convs[3], Gaussian4, Landau_Recon[3] , Iters);  
 
  //// Start return the Gaussian smearing PSF
  //Deconvolution(Convs[0], Gen_Landau[0], PSF[0], Iters);
  //Deconvolution(Convs[1], Gen_Landau[1], PSF[1], Iters);
  //Deconvolution(Convs[2], Gen_Landau[2], PSF[2], Iters);
  //Deconvolution(Convs[3], Gen_Landau[3], PSF[3], Iters);  

  //// Convert the convergence vectors to TH1Fs
  ToTH1F(Converge1, Converge_H[0]);  
  //ToTH1F(Converge2, Converge_H[1]);
  //ToTH1F(Converge3, Converge_H[2]);
  //ToTH1F(Converge4, Converge_H[3]);

  TCanvas* can = new TCanvas(); 
  can -> Print("DeconvLandauXGaussian.pdf[");
  can -> SetLogy(); 

  //for (int i(0); i < Gen_Landau.size(); i++)
  //{
    // Plot the reverting of the Landau
    //Normalize(Landau_Recon[0]); 
    RatioPlot(Gen_Landau[0], Landau_Recon[0], can); 
    can -> Print("DeconvLandauXGaussian.pdf");
    can -> Clear(); 
    
    Converge_H[0] -> Draw("HIST"); 
    can -> Print("DeconvLandauXGaussian.pdf");
    can -> Clear(); 

  //  // Plot the Gaussian PSF
  //  RatioPlot(PSF_Original[i], PSF[i], can); 
  //  can -> Print("DeconvLandauXGaussian.pdf");
  //  can -> Clear(); 
  //}
  can -> Print("DeconvLandauXGaussian.pdf]");










}
