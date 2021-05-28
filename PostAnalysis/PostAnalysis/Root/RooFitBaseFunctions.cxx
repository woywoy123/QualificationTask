#include<PostAnalysis/RooFitBaseFunctions.h>

std::map<TString, std::vector<float>> Normalization(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{
  Average(PDF_H);
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 

  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  if (Params["Range"].size() != 0){x -> setRange("fit", Params["Range"][0], Params["Range"][1]);}
  
  float r = 1; 
  if (Params["r_value"].size() != 0){ r = Params["r_value"][0]; }

  // Normalization variables
  std::vector<TString> l_N = NameGenerator(PDF_H, "_L"); 
  std::vector<float> l_s(l_N.size(), 0.); 
  std::vector<float> l_e(l_N.size(), r*Data -> Integral()); 
  std::vector<RooRealVar*> l_vars = RooRealVariable(l_N, l_s, l_e); 

  // PDF Data variables 
  std::vector<TString> pdf_N_D = NameGenerator(PDF_H, "_D"); 
  std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, PDF_H); 
  
  // PDF PDF variables
  std::vector<TString> pdf_N_P = NameGenerator(PDF_H, "_P"); 
  std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 

  RooArgList N; 
  RooArgList PDFs; 
  for (int i(0); i < l_N.size(); i++)
  {
    PDFs.add(*pdf_P[i]); 
    N.add(*l_vars[i]); 
  }
  RooAddPdf model("model", "model", PDFs, N); 

  // Fitting the PDFs to the Data 
  RooDataHist D("data", "data", *x, Import(*Data)); 
  RooFitResult* re; 
  if (Params["Minimizer"].size() == 0)
  {
    re = model.fitTo(D, Range("fit"), SumW2Error(true), Extended(true), Save()); 
  }
  else
  {
    RooAbsReal* nll = model.createNLL(D, Extended(true)); 
    re = MinimizationCustom(nll, Params);  
  }
  std::map<TString, std::vector<float>> Output;  
  CaptureResults(re, &Output);  
 
  if (Name != "")
  {
    TString plot_t = Name + "PULL.pdf"; 
    RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  }

  Normalize(PDF_H); 
  for (int i(0); i < l_N.size(); i++)
  {
    float n = l_vars[i] -> getVal(); 
    float n_e = l_vars[i] -> getError(); 
    Output["Normalization"].push_back(n); 
    Output["Normalization_Error"].push_back(n_e); 

    PDF_H[i] -> Scale(n); 

    for (int p(0); p < Data -> GetNbinsX(); p++){PDF_H[i] -> SetBinError(p+1, 0.);}
  }

  BulkDelete(l_vars); 
  BulkDelete(pdf_D); 
  BulkDelete(pdf_P); 
  delete x; 
  
  return Output; 
}


std::map<TString, std::vector<float>> NormalizationShift(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{

  float Lum = Data -> Integral(); 
  TH1F* Copy_D = (TH1F*)Data -> Clone("Temp");  

  Normalize(PDF_H); 
  //Average(PDF_H);
  //Normalize(Copy_D); 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 

  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  //if (Params["Range"].size() != 0){x -> setRange("fit", Params["Range"][0], Params["Range"][1]);}
  if (Params["fft_cache"].size() != 0){x -> setBins(Params["fft_cache"][0], "cache"); } 
  
  // Do a preliminary normalization fit:
  std::map<TString, std::vector<float>> Pre = Normalization(Data, PDF_H, Params);

  // Normalization variables
  std::vector<TString> l_N = NameGenerator(PDF_H, "_L"); 
  std::vector<float> l_s = MultiplyByConstant(Pre["Normalization"], 0.001); 
  std::vector<float> l_e = MultiplyByConstant(Pre["Normalization"], 100); 
  std::vector<RooRealVar*> l_vars = RooRealVariable(l_N, l_s, l_e); 

  // Shift variables
  std::vector<TString> d_N = NameGenerator(PDF_H, "_Dx"); 
  std::vector<float> d_s(PDF_H.size(), -0.1); 
  std::vector<float> d_e(PDF_H.size(), 0.1); 
  for (int i(0); i < Params["dx"].size(); i++)
  {
    if (i >= PDF_H.size()){continue;}
    d_s[i] = -Params["dx"][i]; 
    d_e[i] = Params["dx"][i];
  }
  std::vector<RooRealVar*> d_vars; 
  if (Params["dx_G"].size() != 0){d_vars = RooRealVariable(d_N, Params["dx_G"], d_s, d_e);}
  else {d_vars = RooRealVariable(d_N, d_s, d_e);}
  
  RooArgList N; 
  RooArgList PDFs; 
  std::vector<RooDataHist*> pdf_D; 
  std::vector<RooHistPdf*> pdf_P; 
  std::vector<RooFormulaVar*> pdf_F; 

  std::vector<TString> data_n = NameGenerator(PDF_H, "_D"); 
  std::vector<TString> pdf_n = NameGenerator(PDF_H, "_P"); 
  for (int i(0); i < PDF_H.size(); i++)
  {
    RooFormulaVar* dx = new RooFormulaVar(d_N[i] + "dx", "@0 - @1", RooArgList(*x, *d_vars[i])); 
    RooDataHist* DH = new RooDataHist(data_n[i], data_n[i], *x, PDF_H[i]); 
    RooHistPdf* PH = new RooHistPdf(pdf_n[i], pdf_n[i], *dx, *x, *DH); 

    pdf_D.push_back(DH); 
    pdf_P.push_back(PH); 
    pdf_F.push_back(dx); 

    PDFs.add(*PH); 
    N.add(*l_vars[i]);
  }



  // Fitting the PDFs to the Data 
  RooDataHist D("data", "data", *x, Copy_D); 
  RooAddPdf model("model", "model", PDFs, N); 
  std::map<TString, std::vector<float>> Output;  
  RooFitResult* re; 
  if (Params["Minimizer"].size() == 0)
  {
    re = model.fitTo(D, Range("fit"), SumW2Error(true), Extended(true), Save()); 
  }
  else
  {
    RooAbsReal* nll = model.createNLL(D, Extended(true)); 
    re = MinimizationCustom(nll, Params);  
  }
  CaptureResults(re, &Output);
  
  if (Name != "")
  {
    TString plot_t = Name + "PULL.pdf"; 
    RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  }

  for (int i(0); i < l_N.size(); i++)
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
    PDF_H[i] -> Scale(n); //*Lum); 
    
    for (int p(0); p < Data -> GetNbinsX(); p++){PDF_H[i] -> SetBinError(p+1, 0.);}
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
  TH1F* Data = (TH1F*)Data_Org -> Clone("temp"); 
  float Lum = Data -> Integral(); 
  Data -> Sumw2(); 

  //Normalize(Data); 
  Normalize(PDF_H); 
    
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 

  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  if (Params["Range"].size() != 0){x -> setRange("fit", Params["Range"][0], Params["Range"][1]);}
  if (Params["fft_cache"].size() != 0)
  {
    x -> setBins(Params["fft_cache"][0], "fft"); 
    x -> setBins(Params["fft_cache"][0], "cache");  
  }

  // Do a preliminary normalization fit:
  std::map<TString, std::vector<float>> Pre = Normalization(Data, PDF_H, Params);


  // Normalization variables
  std::vector<TString> l_N = NameGenerator(PDF_H, "_L"); 
  std::vector<float> l_s = MultiplyByConstant(Pre["Normalization"], 0.01); 
  std::vector<float> l_e = MultiplyByConstant(Pre["Normalization"], 100); 
  std::vector<RooRealVar*> l_vars; 
  if (Params["l_G"].size() != 0){l_vars = RooRealVariable(l_N, Params["l_G"], l_s, l_e);}
  else {l_vars = RooRealVariable(l_N, l_s, l_e);}

  // Mean variables
  std::vector<TString> m_N = NameGenerator(PDF_H, "_M"); 
  std::vector<float> m_s(PDF_H.size(), -0.1); 
  std::vector<float> m_e(PDF_H.size(), 0.1); 
  bool setconst = false; 
  for (int i(0); i < Params["m"].size(); i++)
  {
    if (i >= PDF_H.size()){continue;}
    if (Params["m"][i] == 0){setconst = true;}
    m_s[i] = -Params["m"][i]; 
    m_e[i] = Params["m"][i];
  }
  std::vector<RooRealVar*> m_vars; 
  if (Params["m_G"].size() != 0){m_vars = RooRealVariable(m_N, Params["m_G"], m_s, m_e);}
  else {m_vars = RooRealVariable(m_N, m_s, m_e);}
  if (setconst){ for (int i(0); i < m_vars.size(); i++){ m_vars[i] -> setConstant(true); } }

  // Standard Deviation variables
  std::vector<TString> s_N = NameGenerator(PDF_H, "_S"); 
  std::vector<float> s_s(PDF_H.size(), 0.001); 
  std::vector<float> s_e(PDF_H.size(), 0.0012); 
  setconst = false;
  
  if (Params["s_s"].size() == 0){setconst = true;}
  for (int i(0); i < Params["s_s"].size(); i++)
  {
    if (i >= PDF_H.size()){continue;}
    if (Params["s_s"][i] == 0)
    {
      setconst = true;
      continue;
    }
    
    s_s[i] = Params["s_s"][i]; 
    s_e[i] = Params["s_e"][i];
  }
  std::vector<RooRealVar*> s_vars; 
  if (Params["s_G"].size() != 0){s_vars = RooRealVariable(s_N, Params["s_G"], s_s, s_e);}
  else {s_vars = RooRealVariable(s_N, s_s, s_e);}
  if (setconst){ for (int i(0); i < s_vars.size(); i++){ s_vars[i] -> setConstant(true); } }

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
  RooDataHist D("data", "data", *x, Data); 
  RooAddPdf model("model", "model", PxG, L); 
  
  std::map<TString, std::vector<float>> Output;  
  RooFitResult* re; 
  if (Params["Minimizer"].size() == 0)
  {
    re = model.fitTo(D, Range("fit"), SumW2Error(true), Extended(true), Save()); 
  }
  else
  {
    RooAbsReal* nll = model.createNLL(D, Range("fit"), Extended(true)); 
    re = MinimizationCustom(nll, Params);  
  }
  CaptureResults(re, &Output); 

  if (Name != "")
  {
    TString plot_t = Name + "PULL.pdf"; 
    RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
  }

  for (int i(0); i < l_N.size(); i++)
  {
    float n = l_vars[i] -> getVal(); 
    float n_er = l_vars[i] -> getError(); 
    float m = m_vars[i] -> getVal(); 
    float m_er = m_vars[i] -> getError(); 
    float s = s_vars[i] -> getVal(); 
    float s_er = s_vars[i] -> getError(); 

    
    Output["Normalization"].push_back(n); 
    Output["Normalization_Error"].push_back(n_er); 
    Output["Mean"].push_back(m); 
    Output["Mean_Error"].push_back(m_er); 
    Output["Stdev"].push_back(s); 
    Output["Stdev_Error"].push_back(s_er); 

    CopyPDFToTH1F(PxG_vars[i], x, PDF_H[i], Data); 
    Normalize(PDF_H[i]); 
    PDF_H[i] -> Scale(n); 
    for (int p(0); p < Data -> GetNbinsX(); p++){PDF_H[i] -> SetBinError(p+1, 0.);}
  }

  delete Data; 
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


std::map<TString, std::vector<float>> DeConvolutionFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name)
{

  std::vector<TString> Names_Dec; 
  float r = 0.1; 
  int bins = Data -> GetNbinsX(); 
  float min = Data -> GetXaxis() -> GetXmin(); 
  float max = Data -> GetXaxis() -> GetXmax(); 
  float w = (max - min) / float(bins); 
  float new_min = min - w*bins*r; 
  float new_max = max + w*bins*r; 

  std::vector<TH1F*> PDF_D;
  std::vector<TH1F*> PSF;  
  std::vector<TH1F*> PDF; 
  for (int i(0); i < PDF_H.size(); i++)
  {
    TString nameG = "Gx_"; nameG+=(i+1); 
    TH1F* Gaus = Gaussian(Params["G_Mean"][i], Params["G_Stdev"][i], bins+2*bins*r, new_min, new_max, nameG); 
    PSF.push_back(Gaus);  

    TString name = "Dec_"; name += (PDF_H[i] -> GetTitle()); 
    TH1F* H = new TH1F(name, name, bins+2*bins*r, new_min, new_max); 
    PDF_D.push_back(H); 
  
    TString name_L = "L_"; name_L += (PDF_H[i] -> GetTitle()); 
    TH1F* X = new TH1F(name_L, name_L, bins+2*bins*r, new_min, new_max); 
    PDF.push_back(X); 

    for (int j(0); j < bins; j++)
    {
      X -> SetBinContent(j+1+r*bins, PDF_H[i] -> GetBinContent(j+1)); 
    }
  }
 
  Normalize(PDF); 
  Normalize(PSF); 
  MultiThreadingDeconvolution(PDF, PSF, PDF_D, Params["LR"][0]); 
  
  for (int i(0); i < PDF_H.size(); i++)
  {
    std::vector<float> vec; 
    for (int j(0); j < bins; j++)
    {
      float e = PDF_D[i] -> GetBinContent(j+1+r*bins); 
      PDF_H[i] -> SetBinContent(j+1, e);    
    }

    delete PDF[i]; 
    delete PSF[i]; 
    delete PDF_D[i]; 
  }
  std::map<TString, std::vector<float>> Output = ConvolutionFFT(Data, PDF_H, Params);
  
  return Output;
}


std::map<TString, std::vector<float>> SimultaneousFFT(std::vector<TH1F*> Data, std::vector<std::vector<TH1F*>> PDF_H_V, std::map<TString, std::vector<float>> Params, TString Name)
{
  float x_min = Data[0] -> GetXaxis() -> GetXmin(); 
  float x_max = Data[0] -> GetXaxis() -> GetXmax(); 
  int bins = Data[0] -> GetNbinsX(); 

  RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
  if (Params["Range"].size() != 0){x -> setRange("fit", Params["Range"][0], Params["Range"][1]);}
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
    // Do a preliminary normalization fit:
    std::map<TString, std::vector<float>> Pre = Normalization(Data[i], PDF_H_V[i], Params);
    Params["l_s"] = MultiplyByConstant(Pre["Normalization"], 0.5); 
    Params["l_e"] = MultiplyByConstant(Pre["Normalization"], 2); 
    std::vector<RooRealVar*> l_var = ProtectionRealVariable("l", PDF_H_V[i], Params, 0, Data[i] -> Integral()); 
    Lumi.push_back(l_var); 
  } 

  std::vector<RooRealVar*> s_var = ProtectionRealVariable("s", PDF_H, Params, 0.0001, 0.0002); 
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
  RooAbsReal* nll = simPdf.createNLL(*ComData, Range("fit")); 
  RooFitResult* re = MinimizationCustom(nll, Params); 
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
  
  for (int i(0); i < PDF_H_V.size(); i++)
  {
    for (int j(0); j < PDF_H_V[i].size(); j++)
    {
      for (int k(0); k < PDF_H_V[i][j] -> GetNbinsX(); k++){PDF_H_V[i][j] -> SetBinError(k+1, 0.);}
    }
  }


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
