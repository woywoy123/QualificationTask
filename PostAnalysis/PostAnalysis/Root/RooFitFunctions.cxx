#include<PostAnalysis/RooFitFunctions.h>
#include<RooMinimizer.h>
using namespace RooFit;

std::vector<RooRealVar*> RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables; 
  for (int i(0); i < Names.size(); i++)
  {
    RooRealVar* r = new RooRealVar(Names[i], Names[i], Begin[i], End[i]); 
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

std::vector<RooDataHist*> RooDataVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<TH1F*> Hist)
{
  std::vector<RooDataHist*> Out; 
  for (int i(0); i < Hist.size(); i++)
  {
    Out.push_back(RooDataVariable(Name[i], domain, Hist[i])); 
  }
  return Out; 
}

std::vector<RooHistPdf*> RooPDFVariables(std::vector<TString> names, RooRealVar* domain, std::vector<RooDataHist*> Hist)
{
  std::vector<RooHistPdf*> Out; 
  for (int i(0); i < names.size(); i++)
  {
    RooHistPdf* P = new RooHistPdf(names[i]+"_p", names[i] + "_p", *domain, *Hist[i]); 
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
  std::vector<RooDataHist*> D_vars = RooDataVariable(pdf_N, x, PDF_H); 
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, D_vars); 
    
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
    delete D_vars[i]; 
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

  Normalize(PDF_H); 
  for (int j(0); j < bins; j++)
  {
    for (int i(0); i < PDF_H.size(); i++)
    {
      PDF_H[i] -> SetBinError(j+1, 1e-12); 
    }
    //Data -> SetBinError(j+1, 1e-12); 

  }

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
  if (cache != 0){ x -> setBins(cache); }
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
  std::vector<TString> pdf_N = NameGenerator(n_vars, "_D");   
  std::vector<RooDataHist*> D_vars = RooDataVariable(pdf_N, x, PDF_H); 
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, D_vars); 
    
  // Create the Convolution PDFs 
  std::vector<RooFFTConvPdf*> fft_vars = RooFFTVariables(pdf_N, x, pdf_vars, g_vars);

  for (int i(0); i < s_vars.size(); i++)
  {
    if (cache != 0){ s_vars[i] -> setBins(cache); }
    if (fft_cache != 0){s_vars[i] -> setBins(fft_cache, "fft");}
    if (cache != 0){ s_vars[i] -> setBins(cache, "cache");}
  
    if (cache != 0){ m_vars[i] -> setBins(cache); }
    if (fft_cache != 0){m_vars[i] -> setBins(fft_cache, "fft");}
    if (cache != 0){ m_vars[i] -> setBins(cache, "cache");}

    if (cache != 0){ l_vars[i] -> setBins(cache); }
    if (fft_cache != 0){l_vars[i] -> setBins(fft_cache, "fft");}
    if (cache != 0){ l_vars[i] -> setBins(cache, "cache");}

    fft_vars[i] -> setInterpolationOrder(9);
    fft_vars[i] -> setBufferFraction(1); 
  }

  // Combine the variables into a single ArgSet and create the model
  RooArgSet Conv;
  RooArgSet N;   
  RooArgSet S; 
  for (int i(0); i < n_vars; i++)
  {
    Conv.add(*fft_vars[i]); 
    N.add(*l_vars[i]); 
    S.add(*s_vars[i]); 
    m_vars[i] -> setConstant(); 
  }
  RooAddPdf model("model", "model", Conv, N); 
  
  // Call the data 
  RooDataHist* D = RooDataVariable("data", x, Data); 
  //model.fitTo(*D, RooFit::SumW2Error(true), RooFit::NumCPU(16), RooFit::Range("fit"), RooFit::Strategy(2)); 
  
  RooAbsReal* nll = model.createNLL(*D, RooFit::Range("fit"), RooFit::NumCPU(16), RooFit::SumW2Error(true)); 
  RooMinimizer* pg = new RooMinimizer(*nll); 
  pg -> setMaxIterations(1000000); 
  pg -> setMaxFunctionCalls(1000000); 
  pg -> setMinimizerType("Minuit"); 
  pg -> setEps(1e-8); 
  pg -> setOffsetting(true);  
  pg -> setStrategy(2); 
  pg -> optimizeConst(true);
  pg -> migrad();
  pg -> minos(); 
  pg -> fit("m"); 


  PlotRooFit(model, x, fft_vars, D);  
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
    
    TH1F* G = Gaussian(m, s, bins, x_min, x_max);  
    Convolution(G, PDF_H[i], Out_H[i]); 
    
    
    //TF1 T = TF1(*fft_vars[i] -> asTF(RooArgList(*x))); 
    //T.SetNpx(bins); 
    //TH1* H = T.CreateHistogram(); 
    //Out_H[i] -> Add(H, 1); 

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
    delete D_vars[i]; 
    delete fft_vars[i]; 
  }
  delete x; 
  delete D; 
  delete nll; 
  delete pg;
  return Out; 
}

std::vector<float> ScalingFit(TH1F* Data, std::vector<TH1F*> PDF_H)
{
  // First we get the domain of the Data histogram we are fitting 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  int n_vars = PDF_H.size(); 

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", 0, x_max);  
  // Create the Luminosity variables for the fit
  std::vector<TString> l_N = NameGenerator(n_vars, "_L");  
  std::vector<float> l_s(n_vars, 0); 
  std::vector<float> l_e(n_vars, Data -> Integral());  
  std::vector<RooRealVar*> l_vars = RooVariables(l_N, l_s, l_e); 

  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N = NameGenerator(n_vars, "_pdf");   
  std::vector<RooDataHist*> D_vars = RooDataVariable(pdf_N, x, PDF_H); 
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, D_vars); 
 
  // Combine the variables into a single ArgSet and create the model
  RooArgList N;   
  RooArgList PDFs;
  for (int i(0); i < l_vars.size(); i++)
  {
    PDFs.add(*pdf_vars[i]); 
    N.add(*l_vars[i]);  
  }
  RooAddPdf model("model", "model", PDFs, N); 
  
  // Call the data 
  RooDataHist* D = RooDataVariable("data", x, Data); 

  RooAbsReal* nll = model.createNLL(*D, RooFit::Range("fit"));  
  RooMinimizer* pg = new RooMinimizer(*nll);   
  pg -> setMaxIterations(100000); 
  pg -> setMaxFunctionCalls(100000); 
  pg -> setEps(1e-4); 
  pg -> setStrategy(2); 
  pg -> setMinimizerType("Minuit2"); 
  pg -> migrad();
  pg -> minos(N); 
  pg -> fit("m"); 

  std::vector<float> Out;  
  Normalize(PDF_H); 
  for (int i(0); i < l_vars.size(); i++)
  {
    float n = l_vars[i] -> getVal(); 
    Out.push_back(n); 
    PDF_H[i] -> Scale(n);   
    delete l_vars[i]; 
    delete pdf_vars[i]; 
    delete D_vars[i]; 
  }
  delete nll; 
  delete pg; 
  delete x; 
  return Out;
}

std::vector<std::pair<TH1F*, std::vector<float>>> FitWithConstraint(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache, int cache)
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
  if (cache != 0){ x -> setBins(cache); }
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
  std::vector<TString> pdf_N = NameGenerator(n_vars, "_D");   
  std::vector<RooDataHist*> D_vars = RooDataVariable(pdf_N, x, PDF_H); 
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, D_vars); 
    
  // Create the Convolution PDFs 
  std::vector<RooFFTConvPdf*> fft_vars = RooFFTVariables(pdf_N, x, pdf_vars, g_vars);
  for (int i(0); i < s_vars.size(); i++){fft_vars[i] -> setBufferFraction(1);}


  // Combine the variables into a single ArgSet and create the model
  RooArgSet Conv;
  RooArgSet N;   
  for (int i(0); i < n_vars; i++)
  {
    Conv.add(*fft_vars[i]); 
    N.add(*l_vars[i]); 
    N.add(*m_vars[i]); 
    N.add(*s_vars[i]); 
  }
  RooAddPdf model("model", "model", Conv, N);
  model.Print("VT");

  // ============================= Adding the Constraints of the Mean and Stdev ========================= // 
  std::vector<RooGaussian*> s_con; 
  std::vector<RooGaussian*> m_con; 
  for (int i(0); i < s_vars.size(); i++)
  {
    RooGaussian* f = new RooGaussian(s_N[i] + "_con", s_N[i] + "_con", *s_vars[i], RooFit::RooConst(0.02), RooFit::RooConst(0.01)); 
    RooGaussian* c = new RooGaussian(m_N[i] + "_con", m_N[i] + "_con", *m_vars[i], RooFit::RooConst(0.02), RooFit::RooConst(0.01));
    s_con.push_back(f); 
    m_con.push_back(c); 
  }
  // =================================================================================================== //
  
  RooArgSet Parameters; 
  RooArgSet Con_Function; 
  for (int i(0); i < s_con.size(); i++)
  {
    Parameters.add(*s_vars[i]); 
    Parameters.add(*m_vars[i]); 
    Con_Function.add(*s_con[i]); 
    Con_Function.add(*m_con[i]); 
  }
  
  RooProdPdf modelc("modelc", "modelc", RooArgSet(model, Con_Function)); 
  modelc.Print("VT");

  // Call the data 
  RooDataHist* D = RooDataVariable("data", x, Data); 
  modelc.fitTo(*D, Constrain(Parameters), RooFit::SumW2Error(true), RooFit::NumCPU(16), RooFit::Range("fit"), RooFit::Strategy(2)); 
  
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
    
    //Params["G_Stdev"][i] = s; 

    std::vector<float> P = {m, s, n, m_e, s_e, n_e}; 
    
    TH1F* G = Gaussian(m, s, bins, x_min, x_max);  
    Convolution(G, PDF_H[i], Out_H[i]); 
    
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
    delete D_vars[i]; 
    delete fft_vars[i]; 
    delete s_con[i]; 
    delete m_con[i]; 
  }
  delete x; 
  delete D; 
  return Out; 
}

TH1F* ExplicitConstraining(TH1F* Data, TH1F* PDF, std::map<TString, std::vector<float>> Params)
{
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  
  // Create the domain: x -> dE/dx 
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", Params["x_range"][0], Params["x_range"][1]); 
  x -> setBins(Params["cache"][0], "cache"); 
  x -> setBins(Params["cache"][0], "fft"); 
  
  // Gaussian Variables for the FFT 
  RooRealVar* m = new RooRealVar("mean", "mean", Params["G_Mean"][0], Params["m_s"][0], Params["m_e"][0]); 
  RooRealVar* s = new RooRealVar("stdev", "stdev", Params["G_Stdev"][0], Params["s_s"][0], Params["s_e"][0]); 
  RooGaussian* g_ = new RooGaussian("gaus", "gaus", *x, *m, *s); 

  // PDF related variables 
  RooRealVar* l = new RooRealVar("lumi", "lumi", 0, Data -> Integral()); 
  RooDataHist* pdf_D = new RooDataHist("PDF_D", "PDF_D", *x, RooFit::Import(*PDF)); 
  RooHistPdf* pdf = new RooHistPdf("PDF", "PDF", *x, *pdf_D); 

  // Perform the Convolution of the PDF with the Gaussian 
  RooFFTConvPdf* PDFxG = new RooFFTConvPdf("PDFxG", "PDFxG", *x, *pdf, *g_); 

  // Define the model that will be used 
  RooAddPdf model("model", "model", RooArgSet(*PDFxG), RooArgSet(*l)); 

  // ===== Define the contraints of variables ===== //
  RooGaussian* g_s = new RooGaussian("stdev_con", "stdev_con", *s, RooFit::RooConst(0.02), RooFit::RooConst(0.01)); 

  // Build the constrained model 
  RooProdPdf modelc("modelc", "modelc", RooArgSet(model, *g_s)); 
 
  RooDataHist* D = new RooDataHist("D", "D", *x, RooFit::Import(*Data));
  modelc.fitTo(*D, Constrain(*s)); 
 
  TF1 T = TF1(*PDFxG -> asTF(RooArgList(*x))); 
  T.SetNpx(bins); 
  TH1* H = T.CreateHistogram(); 
  PDF -> Reset(); 
  PDF -> Add(H, 1); 
  Normalize(PDF); 

  PDF -> Scale(l -> getVal()); 

  delete x; 
  delete m; 
  delete s; 
  delete g_; 
  delete l; 
  delete pdf_D; 
  delete pdf; 
  delete PDFxG; 
  delete g_s; 
  delete D; 
  
  return PDF; 
}

TH1F* ExplicitConstrainingExternal(TH1F* Data, TH1F* PDF, std::map<TString, std::vector<float>> Params)
{
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  
  // Create the domain: x -> dE/dx 
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", Params["x_range"][0], Params["x_range"][1]); 
  x -> setBins(Params["cache"][0], "cache"); 
  x -> setBins(Params["cache"][0], "fft"); 
  
  // Gaussian Variables for the FFT 
  RooRealVar* m = new RooRealVar("mean", "mean", Params["G_Mean"][0], Params["m_s"][0], Params["m_e"][0]); 
  RooRealVar* s = new RooRealVar("stdev", "stdev", Params["G_Stdev"][0], Params["s_s"][0], Params["s_e"][0]); 
  RooGaussian* g_ = new RooGaussian("gaus", "gaus", *x, *m, *s); 

  // PDF related variables 
  RooRealVar* l = new RooRealVar("lumi", "lumi", 0, Data -> Integral()); 
  RooDataHist* pdf_D = new RooDataHist("PDF_D", "PDF_D", *x, RooFit::Import(*PDF)); 
  RooHistPdf* pdf = new RooHistPdf("PDF", "PDF", *x, *pdf_D); 

  // Perform the Convolution of the PDF with the Gaussian 
  RooFFTConvPdf* PDFxG = new RooFFTConvPdf("PDFxG", "PDFxG", *x, *pdf, *g_); 

  // Define the model that will be used 
  RooAddPdf model("model", "model", RooArgSet(*PDFxG), RooArgSet(*l)); 

  // ===== External constraint
  RooGaussian* g_s_external = new RooGaussian("stdev_con_external", "stdev_con_external", *s, RooFit::RooConst(0.02), RooFit::RooConst(0.01)); 
  
  RooDataHist* D = new RooDataHist("D", "D", *x, RooFit::Import(*Data));
  model.fitTo(*D, ExternalConstraints(*g_s_external)); 


  TF1 T = TF1(*PDFxG -> asTF(RooArgList(*x))); 
  T.SetNpx(bins); 
  TH1* H = T.CreateHistogram(); 
  PDF -> Reset(); 
  PDF -> Add(H, 1); 
  Normalize(PDF); 

  PDF -> Scale(l -> getVal()); 

  delete x; 
  delete m; 
  delete s; 
  delete g_; 
  delete l; 
  delete pdf_D; 
  delete pdf; 
  delete PDFxG; 
  delete g_s_external; 
  delete D; 
  
  return PDF; 
}





























