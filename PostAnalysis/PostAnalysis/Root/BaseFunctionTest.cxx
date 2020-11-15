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
  TH1F* Gaus = Gaussian(0, 1, 1000, -5, 5);
  
  TCanvas* can = new TCanvas();
  can -> Print("Gaussian.pdf["); 
  Gaus -> Draw("SAMEHIST"); 
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
  //Convolution(Lan[0], Lan[0], Results[0]); 
  //Convolution(Results[0], Lan[0], Results[1]); 
  //Convolution(Results[1], Lan[0], Results[2]); 
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











}
