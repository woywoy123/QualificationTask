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

std::vector<RooRealVar*> RooRealVarGuess(std::vector<TString> Names, float Guess, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Output; 
  for (int i(0); i < Names.size(); i++)
  {
    RooRealVar* var = new RooRealVar(Names[i], Names[i], Guess, Begin[i], End[i]);  
    Output.push_back(var); 
  }
  return Output; 
}; 


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
  std::vector<float> l_e(n_vars, Data -> Integral());  
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
  model.fitTo(*D, RooFit::SumW2Error(true), RooFit::Range("fit")); 
  
  TString plot_t = Data -> GetTitle(); plot_t += ("_trks_"); plot_t += (PDF_H.size());  plot_t += ("_FitDeconvolution_PULL.pdf"); 
  RooFitPullPlot(model, x, fft_vars, D, plot_t);  

  // Create a histogram vector to store the solution 
  std::vector<TString> out_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<TH1F*> Out_H = CloneTH1F(PDF_H[0], out_N); 
  
  for (int i(0); i < n_vars; i++)
  {
    float s = s_vars[i] -> getVal(); 
    float m = m_vars[i] -> getVal(); 
    float n = l_vars[i] -> getVal(); 

    Params["G_Mean"][i] = m; 
    Params["G_Stdev"][i] = s; 

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
 
  for (int i(0); i < PDF_H.size(); i++){Average(PDF_H[i]);}
  Normalize(PDF_H); 
  float L = Data -> Integral();
  Normalize(Data); 

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
   
  // Declare the gaussian variables used to conduct the fit 
  std::vector<RooRealVar*> s_vars = RooVariables(s_N, s_s, s_e); 
  std::vector<RooRealVar*> m_vars = RooVariables(m_N, m_s, m_e);

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", Params["x_range"][0], Params["x_range"][1]); 
  if (fft_cache != 0){x -> setBins(fft_cache, "fft");}
  if (cache != 0){ x -> setBins(cache, "cache");}
    
  // Create the names for the Gaussian Objects 
  std::vector<TString> g_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<RooGaussian*> g_vars = RooGaussianVariables(g_N, x,  m_vars, s_vars); 
  
  // Create the Luminosity variables for the fit
  std::vector<TString> l_N = NameGenerator(n_vars, "_L");  
  std::vector<float> l_s(n_vars, 0); 
  std::vector<float> l_e(n_vars, 1.1*Data -> Integral());  
  std::vector<RooRealVar*> l_vars = RooVariables(l_N, l_s, l_e); 
  
  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N = NameGenerator(n_vars, "_pdf");   
  std::vector<RooDataHist*> D_vars = RooDataVariable(pdf_N, x, PDF_H); 
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, D_vars); 
  
  // Create the Convolution PDFs 
  std::vector<RooFFTConvPdf*> fft_vars = RooFFTVariables(pdf_N, x, pdf_vars, g_vars);
  
  for (int i(0); i < s_vars.size(); i++)
  {
    //fft_vars[i] -> setInterpolationOrder(9);
    fft_vars[i] -> setBufferFraction(1); 
    l_vars[i] -> setBins(cache, "cache"); 
    s_vars[i] -> setBins(cache, "cache"); 
    m_vars[i] -> setBins(cache, "cache"); 
  
    if (Params["m_e"][i] < 0){m_vars[i] -> setConstant(true); }
    if (Params["s_s"][i] < 0){s_vars[i] -> setConstant(true); }
    
  }

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
  //model.fitTo(*D, RooFit::SumW2Error(true), Minimizer("Minuit"), NumCPU(8), RooFit::Range("fit")); 
  RooAbsReal* nll = model.createNLL(*D, RooFit::Range("fit"), RooFit::SumW2Error(true));  
  RooMinimizer* pg = new RooMinimizer(*nll);   
  pg -> setMaxIterations(100000); 
  pg -> setMaxFunctionCalls(100000); 
  //pg -> setStrategy(2); 
  pg -> setMinimizerType("Minuit");
  //pg -> optimizeConst(1); 
  //pg -> setEps(1e-1); 
  pg -> migrad();
  pg -> fit("h");  
  pg -> cleanup(); 
 
  TString plot_t = Data -> GetTitle(); plot_t += ("_trks_"); plot_t += (PDF_H.size());  plot_t += ("_FitDeconvolution_PULL.pdf"); 
  RooFitPullPlot(model, x, fft_vars, D, plot_t);  

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
    Out_H[i] -> Scale(L*n); 
    Out.push_back(std::pair<TH1F*, std::vector<float>>(Out_H[i], P));
  }
  Data -> Scale(L);

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

std::map<TString, std::vector<float>> ScalingFit(TH1F* Data, std::vector<TH1F*> PDF_H)
{
  for (int i(0); i < PDF_H.size(); i++){Average(PDF_H[i]);}
  Normalize(PDF_H);

  float L = Data -> Integral(); 
  Normalize(Data); 
  
  // First we get the domain of the Data histogram we are fitting 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  int n_vars = PDF_H.size(); 

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", x_min, x_max);  
  
  // Create the Luminosity variables for the fit
  std::vector<TString> l_N = NameGenerator(n_vars, "_L");  
  std::vector<float> l_s(n_vars, 0); 
  std::vector<float> l_e(n_vars, 1.1*Data -> Integral());  
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
  model.fitTo(*D, RooFit::Range("fit"), RooFit::Strategy(1), RooFit::SumW2Error(true)); 

  TString plot_t = Data -> GetTitle(); plot_t += ("_trks_"); plot_t += (PDF_H.size());  plot_t += ("_Normal_PULL.pdf"); 
  RooFitPullPlot(model, x, pdf_vars, D, plot_t);  

  std::map<TString, std::vector<float>> Out;  
  Normalize(PDF_H); 
  for (int i(0); i < l_vars.size(); i++)
  {
    float n = l_vars[i] -> getVal(); 
    float n_e = l_vars[i] -> getError(); 
    std::vector<float> Temp = {n, n_e};
    TString Name = "L"; Name += (i+1); 
    Out[Name] = Temp; 

    PDF_H[i] -> Scale(n*L);   
    delete l_vars[i]; 
    delete pdf_vars[i]; 
    delete D_vars[i]; 
  }
  Data -> Scale(L); 
  delete x; 
  return Out;
}

std::map<TString, std::vector<float>> ScalingShift(TH1F* Data, std::vector<TH1F*> PDF_H)
{
  for (int i(0); i < PDF_H.size(); i++){Average(PDF_H[i]);}
 
  //float L = Data -> Integral(); 
  Normalize(PDF_H); 
  //TH1F* Data_C = (TH1F*)Data -> Clone("Data_Copy");  
  //Normalize(Data_C); 

  // First we get the domain of the Data histogram we are fitting 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  int n_vars = PDF_H.size();  
    
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  x -> setRange("fit", 0.2, 8);  

  std::vector<TString> L_N = NameGenerator(PDF_H.size(), "_L"); 
  std::vector<float> L_S(PDF_H.size(), 0); 
  std::vector<float> L_E(PDF_H.size(), 2.0*Data -> Integral()); 
  std::vector<RooRealVar*> Lp = RooVariables(L_N, L_S, L_E); 
  
  std::vector<TString> D_N = NameGenerator(PDF_H.size(), "_D"); 
  std::vector<float> D_S(PDF_H.size(), -0.2); 
  std::vector<float> D_E(PDF_H.size(), 0.2); 
  std::vector<RooRealVar*> D = RooVariables(D_N, D_S, D_E);  

  std::vector<TString> df = NameGenerator(PDF_H.size(), "_d"); 
  std::vector<TString> n_d = NameGenerator(PDF_H.size(), "_data"); 
  
  RooArgSet Lumi; 
  RooArgSet PDF; 
  std::vector<RooHistPdf*> PDF_Var; 
  std::vector<RooDataHist*> Data_H; 
  std::vector<RooFormulaVar*> dist; 
  for (int i(0); i < PDF_H.size(); i++)
  {
    RooFormulaVar* d = new RooFormulaVar(df[i], "x - @1", RooArgList(*x, *D[i])); 
    
    RooDataHist* O = new RooDataHist(n_d[i], n_d[i], *x, PDF_H[i]); 
    RooHistPdf* O_PDF = new RooHistPdf(n_d[i]+"_pdf", n_d[i]+"_pdf", *d, *x, *O, 1); 
    
    PDF.add(*O_PDF); 
    Lumi.add(*Lp[i]); 

    PDF_Var.push_back(O_PDF); 
    Data_H.push_back(O);
    dist.push_back(d); 
  }
 
  RooDataHist Dat("Data", "Data", *x, Data);
  RooAddPdf model("model", "model", PDF, Lumi); 
  model.fitTo(Dat, RooFit::SumW2Error(true), RooFit::Range("fit")); 
  //RooAbsReal* nll = model.createNLL(Dat, RooFit::SumW2Error(true));  
  //RooMinimizer* pg = new RooMinimizer(*nll);   
  //pg -> setMaxIterations(100000); 
  //pg -> setMaxFunctionCalls(100000); 
  ////pg -> setStrategy(2); 
  //pg -> setMinimizerType("Minuit");
  ////pg -> optimizeConst(1); 
  //pg -> setEps(1e-6); 
  ////pg -> migrad();
  //pg -> fit("m");  
  //pg -> cleanup(); 
  //delete nll; 
  //delete pg; 


  TString plot_t = Data -> GetTitle(); plot_t += ("_trks_"); plot_t += (PDF_H.size());  plot_t += ("_ScalingShift_PULL.pdf"); 
  RooFitPullPlot(model, x, PDF_Var, &Dat, plot_t);  
  
  std::map<TString, std::vector<float>> Out;  
  for (int i(0); i < PDF_H.size(); i++)
  {
    float l = Lp[i] -> getVal(); 
    float l_e = Lp[i] -> getError(); 
    float d = D[i] -> getVal(); 
    
    std::vector<float> Temp = {l, l_e};
    TString Name = "L"; Name += (i+1); 
    Out[Name] = Temp; 

    TF1 T = TF1(*PDF_Var[i] -> asTF(RooArgList(*x))); 
    T.SetNpx(Data -> GetNbinsX()); 
    TH1* H = T.CreateHistogram(); 

    PDF_H[i] -> Reset(); 
    PDF_H[i] -> Add(H, 1); 

    Normalize(PDF_H[i]); 
    PDF_H[i] -> Scale(l); 
  }
  
  //delete Data_C;  
  delete x; 
  for (int i(0); i < PDF_H.size(); i++)
  {
    delete PDF_Var[i]; 
    delete Data_H[i]; 
    delete Lp[i]; 
    delete D[i]; 
    delete dist[i]; 
  }
  return Out; 
}

std::vector<TH1F*> IterativeFitting(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, int fft_cache, int cache)
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
   
  // Declare the gaussian variables used to conduct the fit 
  std::vector<RooRealVar*> s_vars = RooVariables(s_N, s_s, s_e); 
  std::vector<RooRealVar*> m_vars = RooVariables(m_N, m_s, m_e);

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max);
  x -> setRange("fit", Params["x_range"][0], Params["x_range"][1]); 
  if (fft_cache != 0){x -> setBins(fft_cache, "fft");}
  if (cache != 0){ x -> setBins(cache, "cache");}
    
  // Create the names for the Gaussian Objects 
  std::vector<TString> g_N = NameGenerator(n_vars, "_GxT"); 
  std::vector<RooGaussian*> g_vars = RooGaussianVariables(g_N, x,  m_vars, s_vars); 
  
  // Create the Luminosity variables for the fit
  std::vector<TString> l_N = NameGenerator(n_vars, "_L");  
  std::vector<float> l_s(n_vars, 0); 
  std::vector<float> l_e(n_vars, 2 * Data -> Integral());  
  std::vector<RooRealVar*> l_vars = RooVariables(l_N, l_s, l_e); 
  
  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N = NameGenerator(n_vars, "_pdf");   
  std::vector<RooDataHist*> D_vars = RooDataVariable(pdf_N, x, PDF_H); 
  std::vector<RooHistPdf*> pdf_vars = RooPDFVariables(pdf_N, x, D_vars); 
  
  // Create the Convolution PDFs 
  std::vector<RooFFTConvPdf*> fft_vars = RooFFTVariables(pdf_N, x, pdf_vars, g_vars);

  for (int i(0); i < s_vars.size(); i++)
  {
    fft_vars[i] -> setInterpolationOrder(9);
    fft_vars[i] -> setBufferFraction(1); 
  }

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
 
  RooAbsReal* nll = model.createNLL(*D, RooFit::Range("fit"));  
  RooMinimizer* pg = new RooMinimizer(*nll);   
  pg -> setMaxIterations(100000); 
  pg -> setMaxFunctionCalls(100000); 
  pg -> setEps(1e-14); 
  pg -> setStrategy(2); 
  pg -> setMinimizerType("Minuit2"); 
  pg -> migrad();
  pg -> minos(N); 
  pg -> fit("m");  
  
  // Create a histogram vector to store the solution 
  for (int i(0); i < n_vars; i++)
  {
    float s = s_vars[i] -> getVal(); 
    float m = m_vars[i] -> getVal(); 
    float n = l_vars[i] -> getVal(); 
   
    TF1 T = TF1(*fft_vars[i] -> asTF(RooArgList(*x))); 
    T.SetNpx(bins); 
    TH1* H = T.CreateHistogram(); 
    PDF_H[i] -> Reset(); 
    PDF_H[i] -> Add(H, 1); 
    Normalize(PDF_H[i]); 

    PDF_H[i] -> Scale(n); 
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
  return PDF_H; 
}






























