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

// Translate bin number to x in domain range
float BinToDomain(int Bin_Number, int bins, float min, float max)
{
  float bin_width = (max - min)/float(bins); 
  min -= bin_width; 
  return min + Bin_Number*bin_width; 
}

























