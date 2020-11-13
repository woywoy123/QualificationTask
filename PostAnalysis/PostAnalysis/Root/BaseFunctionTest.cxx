#include<PostAnalysis/BaseFunctionTest.h>

void PlotLandau()
{
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<TH1F*> Lan = MakeTH1F({"Landau"}, 1000, 0, 20); 
  Landau(Lan, {1.0}, LandauParams, 500000, 0, 20); 

  TCanvas* can = new TCanvas(); 
  can -> Print("Landau.pdf["); 
  PlotHists(Lan[0], can);
  can -> Print("Landau.pdf"); 
  can -> Print("Landau.pdf]"); 
}

void PlotGaussian()
{
  std::vector<TH1F*> Gaus = MakeTH1F({"Gaussian"}, 1000, -20, 20); 
  Gaussian(0, 1, Gaus[0]);
  
  TCanvas* can = new TCanvas();
  can -> Print("Gaussian.pdf["); 
  Gaus[0] -> Draw("SAMEHIST"); 
  can -> Update();  
  can -> Print("Gaussian.pdf"); 
  can -> Print("Gaussian.pdf]"); 
}

void PlotLandauXLandauConvolved()
{
  int bins = 100; 
  std::vector<float> LandauParams = {1, 0.9, 0.1}; 
  std::vector<TH1F*> Lan = MakeTH1F({"Landau1", "Landau2", "Landau3", "Landau4"}, bins, 0, 20); 
  Landau(Lan, {1, 1, 1, 1}, LandauParams, 50000000, 0, 20); 
  Normalize(Lan); 
  
  std::vector<TH1F*> Results = MakeTH1F({"Landau Convolved One", "Landau Convolved Two", "Landau Convolved Three"}, bins, 0, 20); 
  Convolution(Lan[0], Lan[0], Results[0]); 
  Convolution(Results[0], Lan[0], Results[1]); 
  Convolution(Results[1], Lan[0], Results[2]); 
  //ArtifactRemove(Results);  
  
  TCanvas* can = new TCanvas();
  can -> Print("LandauXLandau.pdf["); 
  RatioPlot(Results[0], Lan[1], can); 
  can -> Print("LandauXLandau.pdf"); 
  can -> Clear(); 

  RatioPlot(Results[1], Lan[2], can); 
  can -> Print("LandauXLandau.pdf"); 
  can -> Clear(); 

  RatioPlot(Results[2], Lan[3], can); 
  can -> Print("LandauXLandau.pdf"); 
  can -> Clear(); 
  
  can -> Print("LandauXLandau.pdf]"); 
}

void PlotGaussianXGaussianConvolved()
{
  int bins = 1001; 
  std::vector<TH1F*> Gaus = MakeTH1F({"Gaus Convolved once with Gaussian mean = 1 and stdev = 0.1", "Gaus Convolved with itself multiple times Mean = 0, Stdev = 0.1"}, bins, -10, 10); 
  std::vector<TH1F*> Results = MakeTH1F({"Gaus1", "Gaus Convolved Two", "Gaus Convolved Three"}, bins, -10, 10); 
  Gaussian(0, 0.1, Gaus[0]); 
  Gaussian(1, 0.1, Gaus[1]); 
 
  Convolution(Gaus[0], Gaus[1], Results[0]); 
  Convolution(Gaus[0], Gaus[0], Results[1]); 
  Convolution(Results[1], Gaus[0], Results[2]); 
    
  TCanvas* can = new TCanvas();
  can -> Print("GaussianXGaussian.pdf["); 
  RatioPlot(Gaus[0], Results[0], can); 
  can -> Print("GaussianXGaussian.pdf"); 
  can -> Clear(); 

  RatioPlot(Gaus[0], Results[1], can); 
  can -> Print("GaussianXGaussian.pdf"); 
  can -> Clear(); 

  RatioPlot(Gaus[0], Results[2], can); 
  can -> Print("GaussianXGaussian.pdf"); 
  can -> Clear(); 
  
  can -> Print("GaussianXGaussian.pdf]"); 
}
