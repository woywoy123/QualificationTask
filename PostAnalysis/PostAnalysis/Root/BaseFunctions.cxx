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

TH1F* SumHists(std::vector<TH1F*> Hists_V, TString name)
{
  TH1F* Data = (TH1F*)Hists_V[0] -> Clone(name); 
  Data -> Reset(); 
  Data -> SetTitle(name); 

  for (int i(0); i < Hists_V.size(); i++)
  {
    Data -> Add(Hists_V[i], 1);
  }
  return Data;
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

std::vector<TString> NameGenerator(std::vector<TH1F*> Hists, TString append)
{
  std::vector<TString> Output; 
  for ( TH1F* H : Hists )
  {
    TString name = H -> GetTitle(); name += ("_" + append); 
    Output.push_back(name); 
  }
  return Output;  
}

std::vector<TH1F*> BulkClone(std::vector<TH1F*> Hists, std::vector<TString> Name)
{
  std::vector<TH1F*> Output; 
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = (TH1F*)Hists[i] -> Clone(Name[i]);
    H -> SetTitle(Name[i]); 
    Output.push_back(H);  
  }
  return Output; 
}


void Flush(std::vector<TH1F*> F_C, std::vector<TH1F*> ntrk_Conv, bool DeletePointer)
{
  if (F_C.size() != ntrk_Conv.size()) {std::cout << "!!!!!!!!!!!!!!!!!" << std::endl; return;}
  for (int i(0); i < F_C.size(); i++)
  {
    float HL = F_C[i] -> Integral(); 
    float OL = ntrk_Conv[i] -> Integral(); 
    
    F_C[i] -> Scale(1. / HL); 
    ntrk_Conv[i] -> Scale(1. / OL); 
    for (int x(0); x < F_C[i] -> GetNbinsX(); x++)
    {
      float e = F_C[i] -> GetBinContent(x+1); 
      float f = ntrk_Conv[i] -> GetBinContent(x+1); 
      ntrk_Conv[i] -> SetBinContent(x+1, e);
    }
    F_C[i] -> Scale(HL); 
    ntrk_Conv[i] -> Scale(OL); 
    
    
    if (DeletePointer){delete F_C[i];}
  }
}

void MatchBins(std::vector<TH1F*> In, TH1F* Data)
{
  int bins = In[0] -> GetNbinsX(); 
  for (int i(0); i < bins; i++)
  {
    float e = Data -> GetBinContent(i+1); 
    
    float s = 0; 
    for (int c(0); c < In.size(); c++){s += In[c] -> GetBinContent(i+1);}
    
    float r = e / s; 
    if (std::isnan(r)){continue;}
    for (int c(0); c < In.size(); c++)
    {
      float g = In[c] -> GetBinContent(c+1); 
      In[c] -> SetBinContent(c+1, g*r); 
    }
  }
}

void SubtractData(std::vector<TH1F*> In, TH1F* Data, int trk, bool trutrk)
{
  for (int i(0); i < In.size(); i++)
  {
    Average(In[i]); 
    if (i != trk && trutrk == false){ Data -> Add(In[i], -1);}
    if (i == trk && trutrk){ Data -> Add(In[i], -1);}
  }
  Average(Data); 
}

void Average(TH1F* Data)
{
  for (int i(0); i < Data -> GetNbinsX(); i++)
  {
    float y = Data -> GetBinContent(i+1); 
    if ( y <= 0){ y = 1e-10; } 
    if ( std::isnan(y) ) { y = 1e-10; }
    Data -> SetBinContent(i+1, y); 
  }
}

void SmoothHist(TH1F* Hist, int iter)
{
  TH1F* ori = (TH1F*)Hist -> Clone("Original"); 
  
  int order = 4; 
  float sigma = 0.01; 
  
  TF1 gaus("mygaus", "gaus", -200, 200); 
  gaus.SetNpx(10000); 
  gaus.SetParameters(1., 0., sigma); 
  
  const int ndata = Hist -> GetNbinsX(); 
  
  std::vector<float> x(ndata, 0); 
  std::vector<float> y(ndata, 0); 
  std::vector<float> sig(ndata, 0); 
  std::vector<float> wt2(ndata, 0); 

  for (int j(0); j < ndata; j++)
  {
    int j1 = j+1; 
    x[j] = Hist -> GetBinCenter(j1); 
    y[j] = Hist -> GetBinContent(j1); 
  
    if (y[j] != 0)
    {
      sig[j] = Hist -> GetBinError(j1); 
      wt2[j] = pow(sig[j], -2); 
    }
    else
    {
      sig[j] = 1.; 
      wt2[j] = 1.; 
    }
  }
  
  for (int j(0); j < ndata; j++)
  {
    TMatrixD A(ndata, order+1); 
    TMatrixDSym W(ndata); 
    TVectorD m(ndata); 

    for (int i(0); i < ndata; i++)
    {
      float dx = x[i] - x[j]; 
      float gausval = gaus.Eval(dx); 
      A(i, 0) = 1; 
      for (int k(1); k < order +1; k++){A(i, k) = A(i, k-1)*dx;}
      
      W(i, i) = gausval*wt2[i]; 
      m(i) = y[i]; 
    }

    TMatrixD ATW = A.T() * W; 
    TMatrixD C = ATW * A.T(); 
    C.Invert(); 

    auto results = C*ATW*m; 

    Hist -> SetBinContent(j+1, results(0)); 
    Hist -> SetBinError(j+1, sqrt(C(0, 0)));
  }

  TH1F* h = (TH1F*)Hist -> Clone("temp"); 
  Convolution(h, h, h); 
  Deconvolution(h, Hist, Hist, 20);
  Normalize(Hist); 
  Hist -> Scale(ori -> Integral());
  delete h; 

  int start = ori -> GetXaxis() -> FindBin(2.0); 
  float t = ori -> GetBinContent(start); 
  float g = Hist -> GetBinContent(start); 
  float s = t/g;  
  Hist -> Scale(s); 

  for (int i(0); i < start; i++){Hist -> SetBinContent(i+1, ori -> GetBinContent(i+1));}
  for (int i(0); i < ndata; i++){Hist -> SetBinError(i+1, 0.);}

  Hist -> SetLineColor(kRed); 
  delete ori; 
}

void Average(std::vector<TH1F*> Data){for (TH1F* H : Data){ Average(H); }}

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

void BulkDelete(std::map<TString, std::vector<TH1F*>> Hists)
{
  typedef std::map<TString, std::vector<TH1F*>>::iterator it; 
  for (it i = Hists.begin(); i != Hists.end(); i++)
  {
    TString n = i -> first; 
    if (n == "Data"){continue;}
    std::vector<TH1F*> Hists_V = Hists[n]; 
    for (TH1F* H : Hists_V){ delete H; }
  }
}

void CoutText(TString *Input, int v, TString Text)
{
  
  for (int c(0); c < v; c++)
  {
    *Input += (Text); 
  }
}

TString PrecisionString(float number, int precision, bool sci)
{
  TString out; 
  std::ostringstream o; 
  o.precision(precision); 
  
  if (sci){o << std::scientific << number;}
  else {o << std::fixed << number;}
  out += (o.str()); 
  return out; 
}







