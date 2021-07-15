#include<PostAnalysis/Statistics.h>


float ErrorByIntegral(TH1F* Hist1, TH1F* Hist2, float x_min, float x_max)
{ 
  int bins = Hist1 -> GetNbinsX(); 
  int bin_min = Hist1 -> GetXaxis() -> FindBin(x_min); 
  int bin_max = Hist1 -> GetXaxis() -> FindBin(x_max); 

  if (x_min == x_max)
  {
    bin_min = 0; 
    bin_max = bins; 
  }
  
  // Calculate the total error from Hist1, compared to Hist2
  float sum = 0; 
  for (int i(bin_min); i < bin_max; i++)
  {
    float h1 = Hist1 -> GetBinContent(i+1); 
    float h2 = Hist2 -> GetBinContent(i+1); 
    sum += std::abs(h1 - h2);
  } 
  float A = Hist2 -> Integral(); 
  float r = std::abs((sum / A)); 
  if (std::isinf(r)){r = 100;}
  if (std::isnan(r)){r = 100;}
  return r; 

}

std::vector<float> Flost2(std::vector<std::vector<TH1F*>> ntrk, std::vector<std::vector<float>> Err)
{
  std::vector<float> fi = {0, 0}; 
  for (int i(0); i < ntrk.size(); i++)
  {
    if (ntrk[i].size() > 1 && i < fi.size())
    {  
      fi[i] = ntrk[i][1] -> Integral(); 
    }
  }

  // Convention: n_trk_tru
  float n_1_2 = fi[0];
  float n_2_2 = fi[1];
  float FLost2 = float(n_1_2) / float(2.*(n_1_2 + n_2_2)); 
  
  std::vector<std::vector<float>> Error = {{0, 0}, {0, 0}, {0, 0}};
  for (int i(0); i < Err.size(); i++)
  {
    for (int x(0); x < Err[i].size(); x++)
    {
      if (x < Error[i].size() && i < Error.size()){Error[i][x] = Err[i][x];}
    }
  }

  float E_1_2 = Error[0][1];  
  float E_2_2 = Error[1][1];  
  
  float E_flost2 = (1. / std::pow(n_1_2 + n_2_2, 2)) * std::pow( std::pow(E_2_2, 2) + std::pow(n_2_2 * E_1_2, 2), 0.5); 

  return {FLost2, E_flost2};
}

std::vector<float> Flost3(std::vector<std::vector<TH1F*>> ntrk, std::vector<std::vector<float>> Err)
{
  std::vector<float> fi = {0, 0, 0, 0}; 
  for (int i(0); i < ntrk.size(); i++)
  {
    if (ntrk[i].size() > 2 && i < fi.size())
    {  
      fi[i] = ntrk[i][2] -> Integral(); 
    }
  }

  // Convention: n_trk_tru
  float n_1_3 = fi[0];
  float n_2_3 = fi[1];
  float n_3_3 = fi[2];
  float n_4_3 = 0; //fi[3]; 
  float FLost3 = float(2*n_1_3 + n_2_3) / float(3.*(n_1_3 + n_2_3 + n_3_3 + n_4_3)); 
  
  std::vector<std::vector<float>> Error = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  for (int i(0); i < Err.size(); i++)
  {
    for (int x(0); x < Err[i].size(); x++)
    {
      if (x < Error[i].size() && i < Error.size()){Error[i][x] = Err[i][x];}
    }
  }

  float E_1_3 = Error[0][2];  
  float E_2_3 = Error[1][2];  
  float E_3_3 = Error[2][2];  
  float E_4_3 = Error[3][2]; 

  float denom = std::pow((n_1_3 + n_2_3 + n_3_3 + n_4_3), 2); 
  float E_flost3 = (1 / denom) * std::pow(float(1/3)*std::pow((n_2_3 + 2*n_3_3 + n_4_3)*E_1_3, 2) + std::pow( (n_3_3 + n_4_3 - n_1_3) * E_2_3, 2) + std::pow( 2*n_1_3 + n_2_3, 2) * (std::pow(E_3_3, 2) + std::pow(E_4_3, 2)), 0.5); 

  return {FLost3, E_flost3};
}

