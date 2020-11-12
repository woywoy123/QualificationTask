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

void Gaussian(float mean, float stdev, TH1F* Hist)
{
  float bins = Hist -> GetNbinsX(); 
  TF1* g = new TF1("gaus", "gaus", -float(bins)/2, float(bins)/2); 
  g -> SetParameters(1, mean, stdev);
  Hist -> Add(g); 
}
