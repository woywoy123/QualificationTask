#include<PostAnalysis/PresentationFigures.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/BaseFunctions.h>

void GaussianXGaussian()
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
  PlotHists({Gaus1, Gaus2}, {"Analytical_Gaussian_Mean_-1_Stdev_0.5", "Analytical_Gaussian_Mean_-2_Stdev_0.7071"}, can);
  can -> Print("Gaussian_Convolution_Test.pdf");
  can -> Clear();

  // Create the Gaussian for the convolution test and the output 
  TH1F* Gaus_Test = Gaussian(mean, stdev, bins, min, max, "Convolution"); 
  TH1F* Conv_Gaus = Gaussian(0, 0.1, bins, min, max, "ConvolutionResult"); 

  // Perform the convolution 
  Convolution(Gaus_Test, Gaus_Test, Conv_Gaus); 

  // Normalize the histograms 
  Conv_Gaus -> Draw("HIST*");
  Normalize(Conv_Gaus);
  Normalize(Gaus2); 

  PlotHists({Conv_Gaus, Gaus2}, {"Gaussian_FFTConvolved", "Analytical_Gaussian_Target"}, can);
  can -> Print("Gaussian_Convolution_Test.pdf"); 
  can -> Print("Gaussian_Convolution_Test.pdf]"); 
} 
