#include<PostAnalysis/Experimental.h>

std::vector<TGraph*> NumericalLandau(std::vector<float> COMP, std::vector<float> Parameters, float min, float max, int points, int toys)
{
  TF1 Lan("Lan", "landau", min, max); 
  for (int i(0); i < Parameters.size(); i++)
  {
    Lan.SetParameter(i, Parameters[i]);
  } 

  std::vector<std::vector<float>> L(COMP.size()); 
  
  for (int x(0); x < toys; x++)
  {
    float r = 0; 
    for (int i(0); i < COMP.size(); i++)
    {
      r = r + Lan.GetRandom(); 
      L[i].push_back(r); 
    }
  }
  
  std::vector<TGraph*> Output; 
  float delta = (max - min)/float(points); 
  for (int i(0); i < L.size(); i++ )
  {
    std::map<float, std::pair<float, int>> Landau_map; 
    for (int x(1); x <= points; x++)
    {
      float x_range = delta * x; 
      Landau_map[x_range] = std::pair<float, int>(x_range, 1); 
    }
    
    std::map<float, std::pair<float, int>>::iterator m; 
    for (m = Landau_map.begin(); m != Landau_map.end(); m++)
    {
      float c_range = m -> first; 
      std::pair<float, int> e = m -> second; 
      int hits = e.second; 
      float dEdx = e.first; 
  
      for (int x(0); x < L[i].size(); x++)
      {
        float h = L[i][x]; 
  
        if (c_range < h && h <= c_range + delta)
        {
          hits++; 
          dEdx = dEdx + h;  
        }
      }
      dEdx = dEdx/float(hits); 
      Landau_map[c_range] = std::pair<float, int>(dEdx, hits); 
    }

    Double_t x[points], y[points]; 
    Int_t v = 0; 
    for (m = Landau_map.begin(); m != Landau_map.end(); m++)
    {
      std::pair<float, int> p = m -> second; 
      x[v] = p.first; 
      y[v] = p.second; 
      v++; 
    }
    
    TGraph* gr = new TGraph(points, x, y); 
    Output.push_back(gr);  
  }
  return Output; 
}

void GraphicalLandau()
{
  std::vector<TGraph*> Land = NumericalLandau({1, 1}, {1, 0.9, 0.1}, 0, 20, 100, 10000000); 
  
  TCanvas* can = new TCanvas(); 
  can -> Print("Graph.pdf["); 
  for (TGraph* G : Land)
  {
    can -> SetLogy(); 
    G -> Draw("AC"); 
    can -> Print("Graph.pdf"); 
    can -> Clear();
  } 
  can -> Print("Graph.pdf]"); 
}

void FindSymmetryGaussian()
{
  float max = 10.5; 
  float min = 0.5;
  int bins = 101;  
  float width = (max-min)/(float)bins; width = width/2;  
  min -= width; 
  max -= width; 
 
  TH1F* Gaus = Gaussian(1, 0.1, bins, min, max);
  TH1F* G2 = (TH1F*)Gaus -> Clone("G"); 
  G2 -> Reset();
  std::vector<float> Forward; 
  std::vector<float> Backward;  
  for ( int i(0); i < Gaus -> GetNbinsX(); i++)
  {
    Forward.push_back(Gaus -> GetBinContent(i+1)); 
    Backward.push_back(Gaus -> GetBinContent(bins - i -1));  
  }
  
  G2 -> SetLineStyle(kDashed);
  G2 -> SetLineColor(kRed); 
  TCanvas* can = new TCanvas();

  float low = 0;  
  int shift = 0;
  bool symmetric = false; 
  for (int i(0); i < bins; i++)
  {
    float e = Forward[i]; 
    for (int j(0); j < bins; j++)
    {
      float f = Backward[j]; 
      float dif = e + f; 
      if (dif > low)
      {
        low = dif; 
        shift = j;
      } 
    }
    for (int k(0); k < bins; k++)
    {
      G2 -> SetBinContent(k-shift+1, Backward[k]); 
    }

    std::cout << shift << std::endl;
    Gaus -> Draw("HIST"); 
    G2 -> Draw("SAMEHIST"); 
    can -> Print("debug.pdf");
    
    if (shift == 0){symmetric = true;} 
  }
   
  
  
  
  
  
   
 
 







}

