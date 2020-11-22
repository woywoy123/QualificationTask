#include<PostAnalysis/PresentationFigures.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/BaseFunctions.h>

void FigureCompiler(TFile* F)
{

  // Get all the Histogram names in the file, including the directories 
  std::map<TString, std::vector<TString>> Map; 
  for (TObject* key : *F -> GetListOfKeys())
  {
    auto k = dynamic_cast<TKey*>(key);
    std::cout << k -> GetName() << "  " << k -> GetClassName() << std::endl;
    Map[(TString)k -> GetName()] = {}; 

  
  
  
   
  }






  //// GaussianXGaussian Figures
  //if ( Directory(F, "GaussianXGaussian") == true )
  //{
  //  F -> cd("GaussianXGaussian"); 
  //  TH1F* G1 = (TH1F*)gDirectory -> Get("Gaussian_1"); 
  //  TH1F*

  //}







}

void GaussianXGaussian(TFile* F)
{
  F -> mkdir("GaussianXGaussian"); 
  F -> cd("GaussianXGaussian");  
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
  
  // Set the titles on the histograms 
  Gaus1 -> SetTitle(n1); 
  Gaus2 -> SetTitle(n2); 

  // Write to file 
  Gaus1 -> Write();
  Gaus2 -> Write(); 
   
  // Create a solution TH1F
  TH1F* Solution = (TH1F*)Gaus1 -> Clone("Convolved_Result"); 
  Solution -> Reset();  
  Solution -> SetTitle("Convolution Output"); 
   
  // Perform the Convolution with Gaus1
  Convolution(Gaus1, Gaus1, Solution); 
  
  // Write solution to file 
  Solution -> Write();  

  F -> cd(); 
} 

void LandauXLandau()
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
}
