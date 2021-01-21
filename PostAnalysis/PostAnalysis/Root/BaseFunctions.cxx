#include<PostAnalysis/BaseFunctions.h>


// Creates Histograms from vector names - Used for bulk generation 
std::vector<TH1F*> MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension)
{
  std::vector<TH1F*> Histograms;
  for (TString name : Names)
  {
    if (Extension != ""){ name+=(Extension); }
    TH1F* hist = new TH1F(name, name, bins, min, max);
    Histograms.push_back(hist);
  }

  return Histograms;
}

// Clone and Reset the histograms 
std::vector<TH1F*> CloneTH1F(TH1F* Hist, std::vector<TString> Names)
{
  float min = Hist -> GetXaxis() -> GetXmin(); 
  float max = Hist -> GetXaxis() -> GetXmax(); 
  int bins = Hist -> GetNbinsX();  
  
  std::vector<TH1F*> Output; 
  for (int i(0); i < Names.size(); i++)
  {
    TH1F* H = new TH1F(Names[i], Names[i], bins, min, max);  
    H -> Reset(); 
    H -> SetTitle(Names[i]); 
    Output.push_back(H); 
  }
  return Output; 
}

// Normalize Histogram
void Normalize(TH1F* Hist)
{
  float Lumi = Hist -> Integral(); 
  Hist -> Scale(1/Lumi); 
}

// Normalize Histograms in bulk 
void Normalize(std::vector<TH1F*> Hists)
{
  for (TH1F* H : Hists){ Normalize(H); }
}

// Normalize vectors 
std::vector<float> Normalize(std::vector<float> V1) 
{
  float sum = 0; 
  for (int i(0); i < V1.size(); i++){sum = sum + V1[i];}
  for (int i(0); i < V1.size(); i++)
  {
    V1[i] = V1[i]/sum;
  }
  return V1; 
}

// TH1F to vector
std::vector<float> ToVector(TH1F* Hist)
{
  int bins = Hist -> GetNbinsX(); 
  std::vector<float> Output; 
   
  for (int i(0); i < bins; i++)
  {
    Output.push_back(Hist -> GetBinContent(i+1));   
  }
  return Output; 
}

// Vector to TH1F
void ToTH1F(std::vector<float> Input, TH1F* Hist)
{
  for (int i(0); i < Input.size(); i++)
  {
    Hist -> SetBinContent(i+1, Input[i]); 
  }
}

// Shift a histogram using bins 
void Shift(TH1F* Hist, int shift)
{
  std::vector<float> Content; 
  for (int i(0); i < Hist -> GetNbinsX(); i++)
  {
    float e = Hist -> GetBinContent(i+1); 
    Content.push_back(e); 
  }

  Hist -> Reset(); 
  for (int i(0); i < Hist -> GetNbinsX(); i++)
  {
    Hist -> SetBinContent(i+1+shift, Content[i]);  
  }
}

//Variable Name Generator 
std::vector<TString> NameGenerator(int number, TString shorty)
{
  std::vector<TString> Output; 
  
  for (int i(0); i < number; i++)
  {
    TString name = shorty; name += (i+1); 
    Output.push_back(name); 
  }
  return Output; 
}

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
  std::cout << "- Absolute Error / Integral of H2: " << ErrorByIntegral(H1, H2, x_min, x_max) << std::endl;
  std::cout << "- KolmogorovTest: " << H2 -> KolmogorovTest(H1) << std::endl;
  std::cout << "- Normalization: H1: " << H1 -> Integral() << " :::: H2: " << H2 -> Integral() << " :::: H1/H2: " << float(H1 -> Integral()) / float(H2 -> Integral()) << std::endl; 
}

float GetMaxValue(TH1F* H)
{
  int bin_m = H -> GetMaximumBin(); 
  float v = H -> GetBinContent(bin_m+1); 
  return v; 
}

void BulkWrite(std::vector<TH1F*> Hist_V)
{
  for (TH1F* H : Hist_V){H -> Write();}
}

void BulkDelete(std::vector<TH1F*> Hist_V)
{
  for (int i(0); i < Hist_V.size(); i++)
  {
    delete Hist_V[i]; 
  }
}



