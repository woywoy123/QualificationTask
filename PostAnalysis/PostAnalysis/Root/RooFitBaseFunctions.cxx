#include<PostAnalysis/RooFitBaseFunctions.h>

std::map<TString, std::vector<float>> Normalization(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{
  Average(PDF_H);
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  std::cout << bins << std::endl;
  
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 

  // Normalization variables
  std::vector<RooRealVar*> l_vars = ProtectionRealVariable("l", PDF_H, Params, 1e-9, (Data -> Integral() +1)); 

  // PDF Data variables 
  std::vector<TString> pdf_N_D = NameGenerator(PDF_H, "_D"); 
  std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, PDF_H); 
  
  // PDF PDF variables
  std::vector<TString> pdf_N_P = NameGenerator(PDF_H, "_P"); 
  std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 

  RooArgList N; 
  RooArgList PDFs; 
  for (int i(0); i < l_vars.size(); i++)
  {
    PDFs.add(*pdf_P[i]); 
    N.add(*l_vars[i]); 
  }
  RooAddPdf model("model", "model", PDFs, N); 

  // Fitting the PDFs to the Data 
  RooDataHist D("data", "data", *x, Import(*Data)); 
  RooFitResult* re = MinimizationCustom(model, &D, Params, x); 

  std::map<TString, std::vector<float>> Output;  
  CaptureResults(re, &Output);  
 
  if (Name != "")
  {
    TString plot_t = Name + "PULL.pdf"; 
    RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  }

  Normalize(PDF_H); 
  for (int i(0); i < l_vars.size(); i++)
  {
    float n = l_vars[i] -> getVal(); 
    float n_e = l_vars[i] -> getError(); 
    Output["Normalization"].push_back(n); 
    Output["Normalization_Error"].push_back(n_e); 
  
    PDF_H[i] -> Scale(n);
  }

  BulkDelete(l_vars); 
  BulkDelete(pdf_D); 
  BulkDelete(pdf_P); 
  delete x; 
  
  return Output; 
}


std::map<TString, std::vector<float>> NormalizationShift(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{
  TH1F* Copy_D = (TH1F*)Data -> Clone("Temp"); 

  Average(PDF_H);

  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 
  
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 

  // Normalization variables
  std::vector<RooRealVar*> l_vars = ProtectionRealVariable("l", PDF_H, Params, 1e-9, (Data -> Integral() +1)); 

  // Shift variables
  std::vector<RooRealVar*> d_vars = ProtectionRealVariable("dx", PDF_H, Params, -0.5, 0.5);   

  RooArgList N; 
  RooArgList PDFs; 
  std::vector<RooRealVar*> All;
  std::vector<RooDataHist*> pdf_D; 
  std::vector<RooHistPdf*> pdf_P; 
  std::vector<RooFormulaVar*> pdf_F; 

  std::vector<TString> data_n = NameGenerator(PDF_H, "_D"); 
  std::vector<TString> pdf_n = NameGenerator(PDF_H, "_P"); 
  for (int i(0); i < PDF_H.size(); i++)
  {
    TString n = "dx_"; n+= (i+1); 
    RooFormulaVar* dx = new RooFormulaVar(n, "@0 - @1", RooArgList(*x, *d_vars[i])); 
    RooDataHist* DH = new RooDataHist(data_n[i], data_n[i], *x, PDF_H[i]); 
    RooHistPdf* PH = new RooHistPdf(pdf_n[i], pdf_n[i], *dx, *x, *DH, 2); 

    pdf_D.push_back(DH); 
    pdf_P.push_back(PH); 
    pdf_F.push_back(dx); 

    PDFs.add(*PH); 
    N.add(*l_vars[i]);
    All.push_back(l_vars[i]); 
    All.push_back(d_vars[i]); 
  }
  RooAddPdf model("model", "model", PDFs, N); 
  
  // Fitting the PDFs to the Data 
  RooDataHist D("data", Data -> GetTitle(), *x, Copy_D); 
  //model.fixAddCoefNormalization(N);
  
  RooFitResult* re = MinimizationCustom(model, &D, Params, x);
  std::map<TString, std::vector<float>> Output;  
  CaptureResults(re, &Output);
  
  if (Name != "")
  {
    TString plot_t = Name + "PULL.pdf"; 
    RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  }

  for (int i(0); i < l_vars.size(); i++)
  {
    float n = l_vars[i] -> getVal(); 
    float n_e = l_vars[i] -> getError(); 
    float d = d_vars[i] -> getVal(); 
    float d_e = d_vars[i] -> getError(); 

    Output["Normalization"].push_back(n); 
    Output["Normalization_Error"].push_back(n_e); 
    Output["Shift"].push_back(d); 
    Output["Shift_Error"].push_back(d_e); 
   
    CopyPDFToTH1F(pdf_P[i], x, PDF_H[i], Data); 
    Normalize(PDF_H[i]); 

    PDF_H[i] -> Scale(n);
  }

  delete Copy_D; 
  BulkDelete(l_vars); 
  BulkDelete(d_vars); 
  BulkDelete(pdf_D); 
  BulkDelete(pdf_P); 
  BulkDelete(pdf_F); 
  delete x; 
  
  return Output; 
}

std::map<TString, std::vector<float>> ConvolutionFFT(TH1F* Data_Org, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{
  auto IntoOutput =[&] (std::map<TString, std::vector<float>>* output, 
                        std::vector<RooRealVar*> N, 
                        std::vector<RooRealVar*> M, 
                        std::vector<RooRealVar*> S, 
                        std::vector<TH1F*> PDF_H, 
                        std::vector<RooFFTConvPdf*> PxG, 
                        RooRealVar* x, 
                        TH1F* Data
  )
  {
    for (int i(0); i < N.size(); i++)
    {
      float n = N[i] -> getVal(); 
      float n_er = N[i] -> getError(); 
      float m = M[i] -> getVal(); 
      float m_er = M[i] -> getError(); 
      float s = S[i] -> getVal(); 
      float s_er = S[i] -> getError(); 
  
      (*output)["Normalization"].push_back(n); 
      (*output)["Normalization_Error"].push_back(n_er); 
      (*output)["Mean"].push_back(m); 
      (*output)["Mean_Error"].push_back(m_er); 
      (*output)["Stdev"].push_back(s); 
      (*output)["Stdev_Error"].push_back(s_er); 
  
      CopyPDFToTH1F(PxG[i], x, PDF_H[i], Data); 
      Normalize(PDF_H[i]); 
      PDF_H[i] -> Scale(n); 
    }
  }; 
  
  Average(PDF_H);
  TH1F* Data = (TH1F*)Data_Org -> Clone("temp"); 
  if (Params["Auto"].size() != 0){AutoRanges(&Params, PDF_H);}

  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 

  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  if (Params["fft_cache"].size() != 0)
  {
    x -> setBins(Params["fft_cache"][0], "fft"); 
    x -> setBins(Params["fft_cache"][0], "cache");  
  }
  
  // Base Variables 
  std::vector<RooRealVar*> l_vars = ProtectionRealVariable("l", PDF_H, Params, 1e-9, (Data -> Integral() +1)); 
  std::vector<RooRealVar*> m_vars = ProtectionRealVariable("m", PDF_H, Params, -0.01, 0.01); 
  std::vector<RooRealVar*> s_vars = ProtectionRealVariable("s", PDF_H, Params, 0.005, 0.0051); 

  // Create the resolution model: Gaussian 
  std::vector<TString> g_N = NameGenerator(PDF_H, "_Gx");
  std::vector<RooGaussian*> g_vars = RooGaussianVariable(g_N, x, m_vars, s_vars); 

  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N_D = NameGenerator(PDF_H, "_D"); 
  std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, PDF_H); 
  std::vector<TString> pdf_N_P = NameGenerator(PDF_H, "_P"); 
  std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 
  
  // Convolve the Gaussian and the PDFs
  std::vector<TString> pxg_N = NameGenerator(PDF_H, "_PxG"); 
  std::vector<RooFFTConvPdf*> PxG_vars = RooFFTVariable(pxg_N, x, g_vars, pdf_P); 

  RooArgList L; 
  RooArgList PxG; 
  for (int i(0); i < PxG_vars.size(); i++)
  {
    PxG_vars[i] -> setCacheObservables( RooArgSet(*x) );
    PxG_vars[i] -> setBufferFraction(1.); 
    PxG.add(*PxG_vars[i]); 
    L.add(*l_vars[i]); 
  }


  // Fitting the PDFs to the Data 
  RooDataHist D("data", "data", *x, Data); 
  RooAddPdf model("model", "model", PxG, L); 
  //model.fixAddCoefNormalization(L);
  RooFitResult* re = MinimizationCustom(model, &D, Params, x);  

  std::map<TString, std::vector<float>> Output;
  CaptureResults(re, &Output);

  TString plot_t = Name; plot_t += "PULL.pdf"; 
  RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  IntoOutput(&Output, l_vars, m_vars, s_vars, PDF_H, PxG_vars, x, Data_Org); 
  
  BulkDelete(l_vars); 
  BulkDelete(m_vars); 
  BulkDelete(s_vars); 
  BulkDelete(g_vars); 
  BulkDelete(pdf_D); 
  BulkDelete(pdf_P); 
  BulkDelete(PxG_vars); 
  delete Data; 
  delete x; 
  return Output; 
}

std::map<TString, std::vector<float>> IncrementalFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{
  auto IntoOutput =[&] (std::map<TString, std::vector<float>>* output, 
                        std::vector<RooRealVar*> N, 
                        std::vector<RooRealVar*> M, 
                        std::vector<RooRealVar*> S, 
                        std::vector<TH1F*> PDF_H, 
                        std::vector<RooFFTConvPdf*> PxG, 
                        RooRealVar* x, 
                        TH1F* Data
  )
  {
    for (int i(0); i < N.size(); i++)
    {
      float n = N[i] -> getVal(); 
      float n_er = N[i] -> getError(); 
      float m = M[i] -> getVal(); 
      float m_er = M[i] -> getError(); 
      float s = S[i] -> getVal(); 
      float s_er = S[i] -> getError(); 
  
      (*output)["Normalization"].push_back(n); 
      (*output)["Normalization_Error"].push_back(n_er); 
      (*output)["Mean"].push_back(m); 
      (*output)["Mean_Error"].push_back(m_er); 
      (*output)["Stdev"].push_back(s); 
      (*output)["Stdev_Error"].push_back(s_er); 
  
      CopyPDFToTH1F(PxG[i], x, PDF_H[i], Data); 
      Normalize(PDF_H[i]); 
      PDF_H[i] -> Scale(n);
    }
  }; 
  
  Average(PDF_H);   
  std::vector<RooRealVar*> l_vars = ProtectionRealVariable("l", PDF_H, Params, 1e-9, (Data -> Integral() +1)); 
  std::vector<RooRealVar*> m_vars = ProtectionRealVariable("m", PDF_H, Params, -0.001, 0.001); 
  std::vector<RooRealVar*> s_vars = ProtectionRealVariable("s", PDF_H, Params, 0.0001, 0.001); 
 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  
  if (Params["fft_cache"].size() != 0)
  {
    x -> setBins(Params["fft_cache"][0], "fft"); 
    x -> setBins(Params["fft_cache"][0], "cache"); 
  }
  
  // Create the resolution model: Gaussian 
  std::vector<TString> g_N = NameGenerator(PDF_H, "_Gx");
  std::vector<RooGaussian*> g_vars = RooGaussianVariable(g_N, x, m_vars, s_vars); 
  
  // Convert the PDFs to RooPDFs
  std::vector<TString> pdf_N_D = NameGenerator(PDF_H, "_D"); 
  std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, PDF_H); 
  std::vector<TString> pdf_N_P = NameGenerator(PDF_H, "_P"); 
  std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 
  
  // Convolve the Gaussian and the PDFs
  std::vector<TString> pxg_N = NameGenerator(PDF_H, "_PxG"); 
  std::vector<RooFFTConvPdf*> PxG_vars = RooFFTVariable(pxg_N, x, g_vars, pdf_P); 
  
  RooArgList L; 
  RooArgList PxG; 
  for (int i(0); i < PxG_vars.size(); i++)
  {
    PxG_vars[i] -> setBufferFraction(1); 
    PxG.add(*PxG_vars[i]); 
    L.add(*l_vars[i]); 
  }
  
  // Fitting the PDFs to the Data 
  RooAddPdf model("model", "model", PxG, L); 
  RooDataHist D("Data", "Data", *x, Data); 
  
  RooFitResult* re; 
  std::vector<float> State; 
  std::vector<float> Bools; 
  
  // Make it such that only the current ntrk-ntru is allowed to float the fit variables 
  // First fit ==== Release m and not s
  Bools = std::vector<float>(PDF_H.size(), 1); 
  VariableConstant(Bools, s_vars); 
  Bools = std::vector<float>(PDF_H.size(), 0); 
  VariableConstant(Bools, m_vars); 
 
  re = MinimizationCustom(model, &D, Params, x); 
  State.push_back(re -> status()); 
  delete re; 
  
  // Second fit ===== Fix m and let s float 
  Bools = std::vector<float>(PDF_H.size(), 0); 
  VariableConstant(Bools, s_vars); 
  Bools = std::vector<float>(PDF_H.size(), 1); 
  VariableConstant(Bools, m_vars); 
  
  re = MinimizationCustom(model, &D, Params, x); 
  State.push_back(re -> status()); 
  delete re;
  
  // Third fit ==== Fix s and m and let only L float
  Bools = std::vector<float>(PDF_H.size(), 1); 
  VariableConstant(Bools, m_vars); 
  Bools = std::vector<float>(PDF_H.size(), 1); 
  VariableConstant(Bools, s_vars); 
  
  re = MinimizationCustom(model, &D, Params, x); 
  State.push_back(re -> status()); 

  std::map<TString, std::vector<float>> Output;
  CaptureResults(re, &Output);

  int status = 0; 
  for (float x : State){ if (x != 0){ status = x; } }
  Output["fit_status"][0] = status; 
  
  TString plot_t = Name; plot_t += "PULL.pdf"; 
  RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  IntoOutput(&Output, l_vars, m_vars, s_vars, PDF_H, PxG_vars, x, Data); 
  
  BulkDelete(l_vars); 
  BulkDelete(m_vars); 
  BulkDelete(s_vars); 
  BulkDelete(g_vars); 
  BulkDelete(pdf_D); 
  BulkDelete(pdf_P); 
  BulkDelete(PxG_vars); 
  delete x; 

  return Output; 
}

std::map<TString, std::vector<float>> SimultaneousFFT(std::vector<TH1F*> Data, std::vector<std::vector<TH1F*>> PDF_H_V, std::map<TString, std::vector<float>> Params, TString Name)
{
  float x_min = Data[0] -> GetXaxis() -> GetXmin(); 
  float x_max = Data[0] -> GetXaxis() -> GetXmax(); 
  int bins = Data[0] -> GetNbinsX(); 
  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 

  if (Params["fft_cache"].size() != 0)
  {
    x -> setBins(Params["fft_cache"][0], "fft"); 
    x -> setBins(Params["fft_cache"][0], "cache");  
  }

  std::vector<TH1F*> PDF_H; 
  for (int i(0); i < Data.size(); i++){PDF_H.push_back(PDF_H_V[i][i]);}
  
  std::vector<std::vector<RooRealVar*>> Lumi; 
  for (int i(0); i < PDF_H_V.size(); i++)
  {
    std::vector<RooRealVar*> l_var = ProtectionRealVariable("l", PDF_H_V[i], Params, 1e-9, Data[i] -> Integral() +1); 
    Lumi.push_back(l_var); 
    Average(PDF_H_V[i]); 
  } 

  std::vector<RooRealVar*> s_var = ProtectionRealVariable("s", PDF_H, Params, 0.001, 0.002); 
  std::vector<RooRealVar*> m_var = ProtectionRealVariable("m", PDF_H, Params, -0.0001, 0.0001); 

  // Create the resolution model: Gaussian 
  std::vector<TString> g_N = NameGenerator(PDF_H, "_Gx");
  std::vector<RooGaussian*> g_vars = RooGaussianVariable(g_N, x, m_var, s_var); 
  
  // Create the PDFs for the model 
  std::vector<TString> pdf_N_D = NameGenerator(PDF_H, "_D"); 
  std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, PDF_H); 
  std::vector<TString> pdf_N_P = NameGenerator(PDF_H, "_P"); 
  std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 

  // Convolve the Gaussian and the PDFs
  std::vector<TString> pxg_N = NameGenerator(PDF_H, "_PxG"); 
  std::vector<RooFFTConvPdf*> PxG_vars = RooFFTVariable(pxg_N, x, g_vars, pdf_P); 

  RooArgList PxG; 
  std::vector<RooArgList> L_Args; 
  for (int i(0); i < PxG_vars.size(); i++)
  {
    RooArgList L; 
    for (int j(0); j < Lumi[i].size(); j++){L.add(*Lumi[i][j]);}
    L_Args.push_back(L); 

    PxG_vars[i] -> setBufferFraction(1); 
    PxG.add(*PxG_vars[i]); 
  }
  
  RooCategory sample("Sample", "Sample"); 
  for (int i(0); i < Data.size(); i++)
  {
    TString name = "trk"; name += (i+1); 
    sample.defineType(name); 
  }
  
  RooSimultaneous simPdf("simPdf", "simPdf", sample);
  std::vector<RooAddPdf*> models; 
  for (int i(0); i < Data.size(); i++)
  {
    TString name = "trk"; name += (i+1); name += ("_F"); 
    RooAddPdf* trk_M = new RooAddPdf(name, name, PxG, L_Args[i]); 

    name = "trk"; name += (i+1);  
    simPdf.addPdf(*trk_M, name); 
    models.push_back(trk_M); 
  }

  std::vector<TString> name_D = NameGenerator(Data, "_F"); 
  std::vector<RooDataHist*> Data_D  = RooDataVariable(name_D, x, Data);
  std::map<std::string, RooDataHist*> Data_V; 
  RooArgList x_var; 
  for (int i(0); i < Data.size(); i++)
  {
    std::string name = "trk"; name += std::to_string(i+1); 
    x_var.add(*x); 
    
    Data_V.insert(std::pair<std::string, RooDataHist*>(name, Data_D[i])); 
  }
  
  RooDataHist* ComData = new RooDataHist("ComData", "ComData", x_var, sample, Data_V, 1.0); 
  RooFitResult* re = MinimizationCustom(simPdf, ComData, Params, x); 
  
  int stat = re -> status(); 
  delete re; 

  for (int i(0); i < PDF_H_V.size(); i++)
  {
    for (int j(0); j < PDF_H_V[i].size(); j++)
    {
      CopyPDFToTH1F(PxG_vars[j], x, PDF_H_V[i][j], Data[0]); 
      float e = Lumi[i][j] -> getVal(); 
      Normalize(PDF_H_V[i][j]); 
      PDF_H_V[i][j] -> Scale(e);
    } 
  }
  
  std::map<TString, std::vector<float>> Output; 
  for (int i(0); i < Lumi.size(); i++)
  {
    for (int j(0); j < Lumi[i].size(); j++)
    {
      TString name_t = "trk"; name_t += (i+1); name_t += ("_"); 
      Output[name_t + "Normalization"].push_back( Lumi[i][j] -> getVal() ); 
      Output[name_t + "Normalization_Error"].push_back( Lumi[i][j] -> getError() ); 
    }
    
    Output["Shift"].push_back(m_var[i] -> getVal()); 
    Output["Shift_Error"].push_back(m_var[i] -> getError()); 
    Output["Stdev"].push_back(s_var[i] -> getVal()); 
    Output["Stdev_Error"].push_back(s_var[i] -> getError()); 
  }
  Output["fit_status"].push_back(stat); 

  delete x; 
  BulkDelete(s_var); 
  BulkDelete(m_var); 
  BulkDelete(g_vars); 
  BulkDelete(pdf_D); 
  BulkDelete(pdf_P); 
  BulkDelete(PxG_vars); 
  for (int i(0); i < Lumi.size(); i++){ BulkDelete(Lumi[i]); }
  for (int i(0); i < models.size(); i++){ delete models[i]; }
  BulkDelete(Data_D); 
  delete ComData;

  return Output; 
}


std::map<TString, std::vector<float>> FractionFitter(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{
  std::map<TString, std::vector<float>> Output; 
  TObjArray* MC = new TObjArray(PDF_H.size()); 
  for (int i(0); i < PDF_H.size(); i++){MC -> Add(PDF_H[i]);}
  TFractionFitter* fit = new TFractionFitter(Data, MC);

  for (int i(0); i < PDF_H.size(); i++)
  {
    fit -> Constrain(i, 1e-9, Data -> Integral() +1);
  }
  int s = fit -> Fit(); 
  Double_t f[PDF_H.size()]; 
  Double_t e[PDF_H.size()]; 

  for (int i(0); i < PDF_H.size(); i++)
  {
    fit -> GetResult(i+1, f[i], e[i]); 
    PDF_H[i] -> Scale(f[i]); 
  }

  std::cout << "status: " << s << std::endl;
  


  return Output;


}

