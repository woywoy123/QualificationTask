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
  return sum/A; 

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


