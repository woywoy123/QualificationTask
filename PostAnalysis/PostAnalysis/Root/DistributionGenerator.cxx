#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/BaseFunctions.h>

TH1F* Landau(float mean, float stdev, int bins, float min, float max, TString Extension)
{
  float Bin_Width = (max - min)/float(bins); 
   
  // Create a title for the histogram   
  TString name = "Landau_M:"; name +=(mean); name +=("_STDEV: "); name += (stdev); name += ("_"); name += (Extension); 
  TH1F* Lan = new TH1F(name, name, bins, min, max); 

  int mp = Lan -> GetXaxis() -> FindBin(mean); 
  float c = Lan -> GetXaxis() -> GetBinCenter(mp); 
  TF1* g = new TF1("Landau", "TMath::Landau(x, [0], [1], 0)", min, max); 
  g -> SetParameters(1, c, stdev); 
  Lan -> Add(g); 
  
  delete g;  
  return Lan;  
}

TH1F* Gaussian(float mean, float stdev, int bins, float min, float max, TString Extension)
{
  float Bin_Width = (max - min)/float(bins); 
   
  // Make sure that the distribution goes through at the centre of the bins. 
  //min -= Bin_Width*0.5; 
  //max -= Bin_Width*0.5; 
  
  // Create a title for the histogram   
  TString name = "Gaussian_M:"; name +=(mean); name +=("_STDEV: "); name += (stdev); name += ("_"); name += (Extension); 
  TH1F* Gaus = new TH1F(name, name, bins, min, max); 

  int mp = Gaus -> GetXaxis() -> FindBin(mean); 
  float c = Gaus -> GetXaxis() -> GetBinCenter(mp); 
  TF1* g = new TF1("gaus", "gaus", min, max); 
  g -> SetParameters(1, c, stdev); 
  Gaus -> Add(g); 
  
  delete g;  
  return Gaus;  
}
