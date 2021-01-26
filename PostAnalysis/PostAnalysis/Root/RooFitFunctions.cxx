#include<PostAnalysis/RooFitFunctions.h>

std::vector<RooRealVar*> RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables; 
  for (int i(0); i < Names.size(); i++)
  {
    RooRealVar* r = new RooRealVar(Names[i], Names[i], (Begin[i] + End[i])/2., Begin[i], End[i]); 
    Variables.push_back(r); 
  }
  return Variables; 
}

std::vector<RooGaussian*> RooGaussianVariables(std::vector<TString> Names, RooRealVar* domain, std::vector<RooRealVar*> mean, std::vector<RooRealVar*> stdev)
{
  std::vector<RooGaussian*> Variables; 
  for (int i(0); i < Names.size(); i++)
  {
    RooGaussian* g = new RooGaussian(Names[i], Names[i], *domain, *mean[i], *stdev[i]);  
    Variables.push_back(g); 
  }
  return Variables; 
}

RooDataHist* RooDataVariable(TString Name, RooRealVar* domain, TH1F* Hist)
{
  RooDataHist* H = new RooDataHist(Name, Name, *domain, RooFit::Import(*Hist)); 
  return H;
}

std::vector<RooHistPdf*> RooPDFVariables(std::vector<TString> names, RooRealVar* domain, std::vector<TH1F*> Hist)
{
  std::vector<RooHistPdf*> Out; 
  for (int i(0); i < names.size(); i++)
  {
    RooDataHist* H = RooDataVariable(names[i], domain, Hist[i]); 
    RooHistPdf* P = new RooHistPdf(names[i]+"p", names[i] + "p", *domain, *H); 
    Out.push_back(P);  
  }
  return Out; 
}

std::vector<RooFFTConvPdf*> RooFFTVariables(std::vector<TString> Names, RooRealVar* domain, std::vector<RooHistPdf*> PDFs, std::vector<RooGaussian*> Gaus) 
{
  std::vector<RooFFTConvPdf*> Out; 
  for (int i(0); i < Names.size(); i++)
  {
    RooFFTConvPdf* p = new RooFFTConvPdf(Names[i] + "pxg", Names[i] + "pxg", *domain, *PDFs[i], *Gaus[i]);
    p -> setBufferFraction(0); 
    Out.push_back(p); 
  }
  return Out; 
}

std::vector<TH1F*> FitDeconvolution(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache, int cache)
{
  // First we get the domain of the Data histogram we are fitting 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  int n_vars = PDF_H.size(); 
 
  // Input variables;  
  // Standard Deviation 
  std::vector<float> s_s = Params["s_s"]; 
  std::vector<float> s_e = Params["s_e"]; 
  std::vector<TString> s_N = NameGenerator(n_vars, "_s"); 
  
  // Mean 
  std::vector<float> m_s = Params["m_s"]; 
  std::vector<float> m_e = Params["m_e"]; 
  std::vector<TString> m_N = NameGenerator(n_vars, "_m"); 

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", Params["x_range"][0], Params["x_range"][1]); 
  if (fft_cache != 0){x -> setBins(fft_cache, "fft");}
  if (cache != 0){ x -> setBins(cache, "cache");}
   
  // Declare the gaussian variables used to conduct the fit 
  std::vector<RooRealVar*> s_vars = RooVariables(s_N, s_s, s_e); 
  std::vector<RooRealVar*> m_vars = RooVariables(m_N, m_s, m_e);
    
  // Create the names for the Gaussian Objects 
  std::vector<TString> g_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<RooGaussian*> g_vars = RooGaussianVariables(g_N, x,  m_vars, s_vars); 
  
  // Create the Luminosity variables for the fit
  std::vector<TString> l_N = NameGenerator(n_vars, "_L");  
  std::vector<float> l_s(n_vars, 0); 
  std::vector<float> l_e(n_vars, 1.1 * Data -> Integral());  
  std::vector<RooRealVar*> l_vars = RooVariables(l_N, l_s, l_e); 
  
  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N = NameGenerator(n_vars, "");   
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, PDF_H); 
    
  // Create the Convolution PDFs 
  std::vector<RooFFTConvPdf*> fft_vars = RooFFTVariables(pdf_N, x, pdf_vars, g_vars);

  // Combine the variables into a single ArgSet and create the model
  RooArgSet Conv;
  RooArgSet N;   
  for (int i(0); i < n_vars; i++)
  {
    Conv.add(*fft_vars[i]); 
    N.add(*l_vars[i]);  
  }
  RooAddPdf model("model", "model", Conv, N); 
  
  // Call the data 
  RooDataHist* D = RooDataVariable("data", x, Data); 
  model.fitTo(*D, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Range("fit")); 
  //PlotRooFit(model, x, D);  

  // Create a histogram vector to store the solution 
  std::vector<TString> out_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<TH1F*> Out_H = CloneTH1F(PDF_H[0], out_N); 
  
  for (int i(0); i < n_vars; i++)
  {
    float s = s_vars[i] -> getVal(); 
    float m = m_vars[i] -> getVal(); 
    float n = l_vars[i] -> getVal(); 

    TH1F* G = Gaussian(m, s, bins, x_min, x_max);  
    Convolution(G, PDF_H[i], Out_H[i]); 
    Normalize(Out_H[i]);
    Out_H[i] -> Scale(n); 
    delete G; 
  }
 
  // Flush all the variables 
  for (int i(0); i < n_vars; i++)
  {  
    delete s_vars[i]; 
    delete m_vars[i]; 
    delete g_vars[i]; 
    delete l_vars[i]; 
    delete pdf_vars[i]; 
    delete fft_vars[i]; 
  }
  delete x; 
  delete D; 
  return Out_H; 
}


std::vector<std::pair<TH1F*, std::vector<float>>> FitDeconvolutionPerformance(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache, int cache)
{
  // First we get the domain of the Data histogram we are fitting 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  int n_vars = PDF_H.size(); 

  // Input variables;  
  // Standard Deviation 
  std::vector<float> s_s = Params["s_s"]; 
  std::vector<float> s_e = Params["s_e"]; 
  std::vector<TString> s_N = NameGenerator(n_vars, "_s"); 
  
  // Mean 
  std::vector<float> m_s = Params["m_s"]; 
  std::vector<float> m_e = Params["m_e"]; 
  std::vector<TString> m_N = NameGenerator(n_vars, "_m"); 

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", Params["x_range"][0], Params["x_range"][1]); 
  if (fft_cache != 0){x -> setBins(fft_cache, "fft");}
  if (cache != 0){ x -> setBins(cache, "cache");}
   
  // Declare the gaussian variables used to conduct the fit 
  std::vector<RooRealVar*> s_vars = RooVariables(s_N, s_s, s_e); 
  std::vector<RooRealVar*> m_vars = RooVariables(m_N, m_s, m_e);
    
  // Create the names for the Gaussian Objects 
  std::vector<TString> g_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<RooGaussian*> g_vars = RooGaussianVariables(g_N, x,  m_vars, s_vars); 
  
  // Create the Luminosity variables for the fit
  std::vector<TString> l_N = NameGenerator(n_vars, "_L");  
  std::vector<float> l_s(n_vars, 0); 
  std::vector<float> l_e(n_vars, Data -> Integral());  
  std::vector<RooRealVar*> l_vars = RooVariables(l_N, l_s, l_e); 
  
  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N = NameGenerator(n_vars, "");   
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, PDF_H); 
    
  // Create the Convolution PDFs 
  std::vector<RooFFTConvPdf*> fft_vars = RooFFTVariables(pdf_N, x, pdf_vars, g_vars);

  // Combine the variables into a single ArgSet and create the model
  RooArgSet Conv;
  RooArgSet N;   
  for (int i(0); i < n_vars; i++)
  {
    Conv.add(*fft_vars[i]); 
    N.add(*l_vars[i]);  
  }
  RooAddPdf model("model", "model", Conv, N); 
  
  // Call the data 
  RooDataHist* D = RooDataVariable("data", x, Data); 
  model.fitTo(*D, RooFit::SumW2Error(true), RooFit::NumCPU(16), RooFit::Range("fit")); 
  //PlotRooFit(model, x, fft_vars, D);  
 
  // Create a histogram vector to store the solution 
  std::vector<TString> out_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<TH1F*> Out_H = CloneTH1F(PDF_H[0], out_N); 
 
  std::vector<std::pair<TH1F*, std::vector<float>>> Out;  
  for (int i(0); i < n_vars; i++)
  {
    float s = s_vars[i] -> getVal(); 
    float m = m_vars[i] -> getVal(); 
    float n = l_vars[i] -> getVal(); 
    float s_e = s_vars[i] -> getError(); 
    float m_e = m_vars[i] -> getError(); 
    float n_e = l_vars[i] -> getError(); 
    
    std::vector<float> P = {m, s, n, m_e, s_e, n_e}; 

    TF1 T = TF1(*fft_vars[i] -> asTF(RooArgList(*x))); 
    T.SetNpx(bins); 
    TH1* H = T.CreateHistogram(); 

    Out_H[i] -> Add(H, 1); 
    Normalize(Out_H[i]);
    Out_H[i] -> Scale(n); 
    Out.push_back(std::pair<TH1F*, std::vector<float>>(Out_H[i], P)); 
  }
 
  // Flush all the variables 
  for (int i(0); i < n_vars; i++)
  {  
    delete s_vars[i]; 
    delete m_vars[i]; 
    delete g_vars[i]; 
    delete l_vars[i]; 
    delete pdf_vars[i]; 
    delete fft_vars[i]; 
  }
  delete x; 
  delete D; 
  return Out; 
}











