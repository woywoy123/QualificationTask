#include<PostAnalysis/DistributionGenerator.h>

void Landau(std::vector<TH1F*> Hists, std::vector<float> COMP, std::vector<float> Parameters, int Number, float min, float max)
{
  // Define the generator and initialize the parameters 
  TF1 Lan("Lan", "landau", min, max);
  for (int i(0); i < Parameters.size(); i++){Lan.SetParameter(i, Parameters[i]);}
 
  // Fill the hist vector
  for (int i(0); i < Number; i++)
  {
    for (int y(0); y < Hists.size(); y++)
    {
      float r(0);
      for (int x(0); x < y+1; x++){r = r + Lan.GetRandom();}
      if (COMP[y] > 0){ Hists[y] -> Fill(r, COMP[y]); } 
    }   
  }
}

TH1F* Gaussian(float mean, float stdev, int bins, float min, float max, TString Extension)
{
  float Bin_Width = (max - min)/float(bins); 
   
  // Make sure that the distribution goes through at the centre of the bins. 
  min -= Bin_Width*0.5; 
  max -= Bin_Width*0.5; 
  
  // Create a title for the histogram   
  TString name = "Gaussian_"; name += (mean); name +=("_"); name += (stdev); name += ("_"); name += (Extension); 
  TH1F* Gaus = new TH1F(name, name, bins, min, max); 
  
  TF1* g = new TF1("gaus", "gaus", min, max); 
  g -> SetParameters(1, mean, stdev); 
  Gaus -> Add(g); 
  
  delete g;  
  return Gaus;  
}
