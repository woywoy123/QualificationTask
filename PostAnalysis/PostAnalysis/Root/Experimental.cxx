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

  TGraph* g = Land[0]; 
  int points = g -> GetN(); 
  Double_t x[points], y[points]; 

  std::vector<float> x_v, y_v; 
  for (int v(0); v < points; v++)
  {
    g -> GetPoint(v, x[v], y[v]); 
    x_v.push_back(x[v]); y_v.push_back(y[v]); 
  }
  
  y_v = ConvolutionFFT(y_v, y_v); 
  y_v = Normalize(y_v); 
  for (int v(0); v < points; v++)
  {
    x[v] = x_v[v]; 
    y[v] = y_v[v]; 
  }
  TGraph* gr = new TGraph(points, x, y); 

  for (int v(0); v < points; v++)
  {
    Land[1] -> GetPoint(v, x[v], y[v]); 
    x_v[v] = x[v]; 
    y_v[v] = y[v]; 
  }
  y_v = Normalize(y_v); 
  for (int v(0); v < points; v++)
  {
    x[v] = x_v[v]; 
    y[v] = y_v[v]; 
  }
  TGraph* g_n = new TGraph(points, x, y); 
 
  TMultiGraph *mg = new TMultiGraph();
  mg -> Add(gr, "AC"); 
  mg -> Add(g_n, "AC");
  mg -> Draw("AC"); 
  can -> Print("Graph.pdf"); 

  can -> Print("Graph.pdf]"); 

}

void ConvolutionTest()
{
  int bins = 100; 
  float stdev = 0.5; 
  float mean = -1; 
  float min = -10; 
  float max = 10; 

  // Generate Gaussian histograms 
  TH1F* Gaus1 = Gaussian(mean, stdev, bins, min, max, "1"); 
  TH1F* Gaus2 = Gaussian(mean+mean, stdev*sqrt(2), bins, min, max, "2");
  
  TCanvas* can = new TCanvas(); 
  can -> Print("Gaussian_Convolution_Test.pdf["); 
  Gaus1 -> Draw("HIST");
  Gaus2 -> Draw("SAMEHIST"); 
  can -> Print("Gaussian_Convolution_Test.pdf"); 
  can -> Clear();

  // Create the Gaussian for the convolution test and the output 
  TH1F* Gaus_Test = Gaussian(mean, stdev, bins, min, max, "Convolution"); 
  TH1F* Conv_Gaus = Gaussian(0, 0.1, bins, min, max, "ConvolutionResult"); 

  // Perform the convolution 
  Convolution_Experimental(Gaus_Test, Gaus_Test, Conv_Gaus); 

  // Normalize the histograms 
  Normalize(Conv_Gaus);
  Normalize(Gaus2); 

  Gaus2 -> Draw("HIST");
  Conv_Gaus -> Draw("SAMEHIST"); 
  can -> Print("Gaussian_Convolution_Test.pdf"); 
  
  std::cout << "########## Analytical Truth Parameters #########" << std::endl;
  std::cout << "# Parameters of Gaussian 1-> Mean: " << mean << " Standard Deviation: " << stdev << std::endl; 
  std::cout << "# Parameters of Gaussian 2-> Mean: " << mean + mean << " Standard Deviation: " << stdev*sqrt(2) << std::endl; 
  std::cout << "# Peak position of Gaussian 1 (Mean): " << Gaus1 -> GetMean() << " Gaussian 2 (Mean): " << Gaus2 -> GetMean() << std::endl; 
  std::cout << "# Standard Deviation of Gaussian 1 (STDEV): " << Gaus1 -> GetRMS() << " Gaussian 2 (STDEV): " << Gaus2 -> GetRMS() << std::endl; 
  std::cout << "################################################" << std::endl;
  std::cout << std::endl;
  std::cout << "######### Output Convolution Parameters ########" << std::endl;
  std::cout << "# Parameters of Convolution-> Mean: " << Conv_Gaus -> GetMean() << " Standard Deviation: " << Conv_Gaus -> GetRMS() << std::endl; 
  std::cout << "# Integral: Analytical -> "<< Gaus2 -> Integral() << " Convolved -> " << Conv_Gaus -> Integral() << std::endl;
  can -> Print("Gaussian_Convolution_Test.pdf]"); 
}
