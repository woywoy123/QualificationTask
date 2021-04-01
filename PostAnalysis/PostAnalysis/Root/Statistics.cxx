#include<PostAnalysis/Statistics.h>


// Benchmark function 
float Pythagoras(std::vector<float> v1, std::vector<float> v2)
{
  float sum = 0; 
  for (int i(0); i < v1.size(); i++)
  {
    float diff = pow(v1[i] - v2[i], 2); 
    sum += diff;  
  }
  return pow(sum, 0.5); 
}

float SquareError(TH1F* Hist1, TH1F* Hist2, float x_min, float x_max)
{
  int bins = Hist1 -> GetNbinsX(); 
  int bin_min = Hist1 -> GetXaxis() -> FindBin(x_min);
  int bin_max = Hist1 -> GetXaxis() -> FindBin(x_max);
  
  if (x_min == x_max)
  {
    bin_min = 0; 
    bin_max = bins; 
  }
 
  std::vector<float> v1, v2; 
  for (int i(bin_min); i < bin_max; i++)
  {
    v1.push_back(Hist1 -> GetBinContent(i+1)); 
    v2.push_back(Hist2 -> GetBinContent(i+1)); 
  }
  return Pythagoras(v1, v2); 
}

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
  float A = Hist2 -> Integral(bin_min, bin_max); 
  float r = sum / A; 
  if (std::isinf(r)){r = 0;}
  return r; 

}

void Statistics(TH1F* H1, TH1F* H2, float x_min, float x_max)
{
  float max = H1 -> GetXaxis() -> GetXmax(); 
  float min = H1 -> GetXaxis() -> GetXmin(); 

  if (x_min == x_max) 
  {
    x_max = max; 
    x_min = min;   
  }

  std::cout << "H1: " << H1 -> GetTitle() << " ::::: H2: " << H2 -> GetTitle() << std::endl;
  std::cout << "Domain selected: " << x_min << " -> " << x_max << std::endl;
  std::cout << "- Absolute Error / Integral of H2: " << ErrorByIntegral(H1, H2, x_min, x_max)*100 << "%" << std::endl;
  std::cout << "- KolmogorovTest: " << H2 -> KolmogorovTest(H1) << std::endl;
  std::cout << "- Normalization: H1: " << H1 -> Integral() << " :::: H2: " << H2 -> Integral() << " :::: H1/H2: " << float(H1 -> Integral())*100 / float(H2 -> Integral()) << "%" << std::endl; 
  std::cout << std::endl;
}


std::vector<float> Flost2(std::vector<std::vector<TH1F*>> ntrk, std::vector<TH1F*>Err)
{
  std::vector<float> fi = {0, 0, 0}; 
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
  float n_3_2 = fi[2];
  float FLost2 = float(n_1_2 + n_3_2) / float(2.*(n_1_2 + n_2_2 + n_3_2)); 
  
  std::vector<std::vector<float>> Error = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int i(0); i < Err.size(); i++)
  {
    for (int x(0); x < Err[i] -> GetNbinsX(); x++)
    {
      if (x < Error[i].size() && i < Error.size()){Error[i][x] = Err[i] -> GetBinContent(x+1);}
    }
  }

  float E_1_2 = Error[0][1];  
  float E_2_2 = Error[1][1];  
  float E_3_2 = Error[2][1];  
  
  float E_flost2 = (1./(2*std::pow(n_1_2 + n_2_2 + n_3_2, 2))) * std::pow(std::pow(n_2_2 - n_1_2 - n_3_2, 2)*(std::pow(E_1_2, 2) + std::pow(E_3_2, 2)) + std::pow(2*(n_1_2 + n_3_2), 2)*std::pow(E_2_2, 2), 0.5); 

  return {FLost2, E_flost2};
}

std::vector<float> Flost3(std::vector<std::vector<TH1F*>> ntrk, std::vector<TH1F*>Err)
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
  float n_4_3 = fi[3]; 
  float FLost3 = float(2*n_1_3 + n_2_3 + n_4_3) / float(3.*(n_1_3 + n_2_3 + n_3_3 + n_4_3)); 
  
  std::vector<std::vector<float>> Error = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  for (int i(0); i < Err.size(); i++)
  {
    for (int x(0); x < Err[i] -> GetNbinsX(); x++)
    {
      if (x < Error[i].size() && i < Error.size()){Error[i][x] = Err[i] -> GetBinContent(x+1);}
    }
  }

  float E_1_3 = Error[0][2];  
  float E_2_3 = Error[1][2];  
  float E_3_3 = Error[2][2];  
  float E_4_3 = Error[3][2]; 
  
  float E_flost3 = float( float(1.) / float( 3. * std::pow(n_1_3 + n_2_3 + n_3_3 + n_4_3, 2)) ) * std::sqrt( std::pow(n_2_3 + 2*n_3_3 + n_4_3, 2) * std::pow(E_1_3, 2 ) + std::pow(n_3_3 - n_1_3, 2) * (std::pow(E_4_3, 2) + std::pow(E_2_3, 2)) + std::pow( 2*n_1_3 + n_2_3 + n_4_3, 2 ) * std::pow(E_3_3, 2));

  return {FLost3, E_flost3};
}

