#include<PostAnalysis/AlgorithmFunctions.h>

std::vector<TH1F*> BuildNtrkNtru(int n, TH1F* trk1_start, TString extension)
{
  std::vector<TString> ntrk_ntru_names; 
  for (int i(0); i < n; i++)
  {
    TString base = "dEdx_ntrk_"; base += (i+1); base += ("_ntru_"); base += (i+1); 
    ntrk_ntru_names.push_back(base); 
  }
  std::vector<TH1F*> ntrk_ntru_H = ConvolveNTimes(trk1_start, ntrk_ntru_names.size(), ntrk_ntru_names, extension); 
  return ntrk_ntru_H; 
}

// basic n-track, n-truth fits 
std::vector<TH1F*> Normalization_Fit(std::vector<TH1F*> Data, TH1F* trk1_Start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_Normal"; 
  std::vector<TH1F*> ntrk_ntru_H = BuildNtrkNtru(Data.size(), trk1_Start, ext); 
  for (int i(0); i < ntrk_ntru_H.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 

    TH1F* ntrk = ntrk_ntru_H[i];
    TH1F* ntrk_ntru_T = Data[i]; 

    TString base = "Fit_"; base += (i+1); base += (ext); 
    Normalization(ntrk_ntru_T, {ntrk}, Params, base); 
  }
 
  WriteHistsToFile(ntrk_ntru_H, JE + "/Normal"); 
  return ntrk_ntru_H; 
}

std::vector<TH1F*> NormalizationShift_Fit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  
  TString ext = "_" + JE + "_NormalShift"; 
  std::vector<TH1F*> ntrk_ntru_H = BuildNtrkNtru(Data.size(), trk1_start, ext); 
  
  for (int i(0); i < ntrk_ntru_H.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    TH1F* ntrk = ntrk_ntru_H[i]; 
    TH1F* ntrk_ntru_T = Data[i]; 

    TString base = "Fit_"; base += (i+1); base += (ext); 
    NormalizationShift(ntrk_ntru_T, {ntrk}, Params, base); 
  }
  WriteHistsToFile(ntrk_ntru_H, JE + "/ShiftNormal"); 

  return ntrk_ntru_H; 
}

std::vector<TH1F*> NormalizationShiftFFT_Fit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  
  TString ext = "_" + JE + "_NormalShiftFFT"; 
  std::vector<TH1F*> ntrk_ntru_H = BuildNtrkNtru(Data.size(), trk1_start, ext); 

  for (int i(0); i < ntrk_ntru_H.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    TH1F* ntrk = ntrk_ntru_H[i]; 
    TH1F* ntrk_ntru_T = Data[i]; 

    TString base = "Fit_"; base += (i+1); base += (ext); 
    ConvolutionFFT(ntrk_ntru_T, {ntrk}, Params, base); 
  }
  
  WriteHistsToFile(ntrk_ntru_H, JE + "/ShiftNormalFFT"); 

  return ntrk_ntru_H; 
}


std::vector<TH1F*> NormalizationShiftWidthFFT_Fit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  
  TString ext = "_" + JE + "_NormalShiftWidthFFT"; 
  std::vector<TH1F*> ntrk_ntru_H = BuildNtrkNtru(Data.size(), trk1_start, ext); 

  for (int i(0); i < ntrk_ntru_H.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    TH1F* ntrk = ntrk_ntru_H[i]; 
    TH1F* ntrk_ntru_T = Data[i]; 

    TString base = "Fit_"; base += (i+1); base += (ext); 
    ConvolutionFFT(ntrk_ntru_T, {ntrk}, Params, base); 
  }
  
  WriteHistsToFile(ntrk_ntru_H, JE + "/ShiftNormalWidthFFT"); 

  return ntrk_ntru_H; 
}

std::vector<std::vector<TH1F*>> BuildNtrkMtru(int n, TH1F* trk1_start, TString extension, int tru)
{
  std::vector<std::vector<TString>> ntrk_ntru_names; 
  for (int i(0); i < n; i++)
  {
    std::vector<TString> ntrk_ntru_n; 
    for (int j(0); j < tru; j++)
    {
      TString base = "dEdx_ntrk_"; base += (i+1); base += ("_ntru_"); base += (j+1); 
      ntrk_ntru_n.push_back(base); 
    }
    ntrk_ntru_names.push_back(ntrk_ntru_n); 
  }
  
  std::vector<std::vector<TH1F*>> ntrk_ntru_templates; 
  for (int i(0); i < n; i++)
  {
    std::vector<TH1F*> ntrk_ntru_H = ConvolveNTimes(trk1_start, ntrk_ntru_names[i].size(), ntrk_ntru_names[i], extension); 
    ntrk_ntru_templates.push_back(ntrk_ntru_H); 
  }
  return ntrk_ntru_templates; 
}

std::vector<std::vector<TH1F*>> Normalization_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_Normal_NtrkMtru"; 
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext);
  
  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/Normal"); 

  for (int i(0); i < Data.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i]; 
    TH1F* ntrk_Measure = Data[i]; 
    
    TString base = "Fit_"; base += (i+1); base += (ext); 
    std::map<TString, std::vector<float>> Map = Normalization(ntrk_Measure, ntrk_Template, Params, base); 
    WriteHistsToFile(ntrk_Template, JE + "/Normal");
    
    TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
    WriteOutputMapToFile(Map, JE + "/Normal", trk_n); 
  }
  
  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> NormalizationShift_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_ShiftNormal_NtrkMtru"; 
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext);

  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/ShiftNormal"); 

  for (int i(0); i < Data.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i]; 
    TH1F* ntrk_Measure = Data[i]; 
    
    TString base = "Fit_"; base += (i+1); base += (ext); 
    std::map<TString, std::vector<float>> Map = NormalizationShift(ntrk_Measure, ntrk_Template, Params, base); 
    
    WriteHistsToFile(ntrk_Template, JE + "/ShiftNormal"); 

    TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
    WriteOutputMapToFile(Map, JE + "/ShiftNormal", trk_n); 
  }
  
  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> NormalizationShiftFFT_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_ShiftNormalFFT_NtrkMtru"; 
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext);

  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/ShiftNormalFFT"); 

  for (int i(0); i < Data.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i]; 
    TH1F* ntrk_Measure = Data[i]; 
    
    TString base = "Fit_"; base += (i+1); base += (ext); 
    std::map<TString, std::vector<float>> Map = ConvolutionFFT(ntrk_Measure, ntrk_Template, Params, base); 
    
    WriteHistsToFile(ntrk_Template, JE + "/ShiftNormalFFT"); 
 
    TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
    WriteOutputMapToFile(Map, JE + "/ShiftNormalFFT", trk_n);  
  }
  
  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> NormalizationShiftWidthFFT_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_ShiftNormalWidthFFT_NtrkMtru"; 
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext);

  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/ShiftNormalWidthFFT"); 

  for (int i(0); i < Data.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 
    std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i]; 
    TH1F* ntrk_Measure = Data[i]; 
    
    TString base = "Fit_"; base += (i+1); base += (ext); 
    std::map<TString, std::vector<float>> Map = ConvolutionFFT(ntrk_Measure, ntrk_Template, Params, base); 
    
    WriteHistsToFile(ntrk_Template, JE + "/ShiftNormalWidthFFT"); 
   
    TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
    WriteOutputMapToFile(Map, JE + "/ShiftNormalWidthFFT", trk_n);  
  
  }
  
  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> Experimental_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_Experimental_NtrkMtru"; 
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext);
  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/Experimental"); 
  
  int iter = 3; 
  for (int x(0); x < iter; x++)
  {
    if (x > iter)
    {
      for (int t(0); t < Data.size(); t++)
      {
        TString name_temp = Data[t] -> GetTitle(); 
        TH1F* ntrk_Data = (TH1F*)Data[t] -> Clone(name_temp + "_C"); 
        for (int j(0); j < ntrk_mtru_H[t].size(); j++)
        {
          if (j != t){ ntrk_Data -> Add(ntrk_mtru_H[t][j], -1); }
        }
        
        Average(ntrk_Data); 
        Flush({ntrk_Data}, {ntrk_mtru_H[t][t]}); 
        delete ntrk_Data;  
      }

      for (int t(0); t < Data.size(); t++)
      {
        for (int j(0); j < Data.size(); j++)
        {
          Flush({ntrk_mtru_H[t][t]}, {ntrk_mtru_H[j][t]}); 
        }
      }
    }   
     
    for (int i(0); i < Data.size(); i++)
    {
      TString r_name = "Range_ntrk_"; r_name += (i+1); 
      Params["Range"] = Params[r_name]; 
      std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i]; 
      TH1F* ntrk_Measure = (TH1F*)Data[i] -> Clone("trk_cop"); 
      
      TString base = "Fit_"; base += (i+1); base += (ext); 
      
      Normalize(ntrk_Measure); 
      std::map<TString, std::vector<float>> prefit = Normalization(ntrk_Measure, ntrk_Template, Params, base); 
      std::vector<float> Norm = prefit["Normalization"]; 
      Params["l_G"] = Norm;
      std::map<TString, std::vector<float>> Map = ConvolutionFFT(ntrk_Measure, ntrk_Template, Params, base); 
      float L = Data[i] -> Integral(); 
      for (TH1F* H : ntrk_Template){ H -> Scale(L); }
      delete ntrk_Measure;
      if (x == iter -1)
      {
        TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
        WriteOutputMapToFile(Map, JE + "/Experimental", trk_n);  
      }
    }
  }

  for (int i(0); i < ntrk_mtru_H.size(); i++){ WriteHistsToFile(ntrk_mtru_H[i], JE + "/Experimental"); }


  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> Simultaneous_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_Simultaneous";
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext, Data.size()); 
  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/Simultaneous"); 
  
  std::map<TString, std::vector<float>> Fit_Res = SimultaneousFFT(Data, ntrk_mtru_H, Params, ext);  

  
  for (int i(0); i < Data.size(); i++)
  {
    std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i];  
    WriteHistsToFile(ntrk_Template, JE + "/Simultaneous"); 
   
    TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 

    std::map<TString, std::vector<float>> Map; 
    TString trk_e = "trk"; trk_e += (i+1); trk_e += ("_"); 
    Map.insert(std::pair<TString, std::vector<float>>("Normalization", Fit_Res[trk_e + "Normalization"])); 
    Map.insert(std::pair<TString, std::vector<float>>("Normalization_Error", Fit_Res[trk_e + "Normalization_Error"])); 
    Map.insert(std::pair<TString, std::vector<float>>("Shift", Fit_Res["Shift"])); 
    Map.insert(std::pair<TString, std::vector<float>>("Shift_Error", Fit_Res["Shift_Error"])); 
    Map.insert(std::pair<TString, std::vector<float>>("Stdev", Fit_Res["Stdev"])); 
    Map.insert(std::pair<TString, std::vector<float>>("Stdev_Error", Fit_Res["Stdev_Error"])); 
    Map.insert(std::pair<TString, std::vector<float>>("fit_status", Fit_Res["fit_status"])); 

    WriteOutputMapToFile(Map, JE + "/Simultaneous", trk_n);
  }

  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> BruceMethod(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  auto Minimization =[&] (RooAbsReal* nll, std::map<TString, std::vector<float>> Params)
  {
    RooMinimizer* pg = new RooMinimizer(*nll); 
    pg -> setMaxIterations(Params["Minimizer"][0]); 
    pg -> setMaxFunctionCalls(Params["Minimizer"][0]); 
    pg -> migrad(); 
    pg -> minos();
    pg -> hesse();
    pg -> improve(); 
    pg -> optimizeConst(true); 
    //pg -> setEvalErrorWall(true); 
    pg -> setEps(1e-6); 
    pg -> setPrintLevel(0); 
    RooFitResult* re = pg -> fit("r"); 
    pg -> cleanup(); 
    delete pg; 
    delete nll;
    
    int p = re -> status(); 
    delete re; 
    return p; 
  }; 

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
      for (int p(0); p < Data -> GetNbinsX(); p++){PDF_H[i] -> SetBinError(p+1, 0.);}
    }
  }; 

  auto GenerateVariableVectors =[&] (std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Ext)
  {
    std::vector<TString> m_N = NameGenerator(PDF_H, "_" + Ext); 
    std::vector<RooRealVar*> m_vars = RooRealVariable(m_N, Params[Ext + "_G"], Params[Ext + "_s"], Params[Ext + "_e"]);
    return m_vars; 
  }; 

  auto SetConstantVar =[&] (std::vector<RooRealVar*> Vars, std::vector<bool> OFF)
  {
    for (int i(0); i < Vars.size(); i++)
    { 
      Vars[i] -> setConstant(OFF[i]); 
      
     }

  };

  std::vector<bool> OFF; 
  std::vector<std::vector<float>> Ranges = {{0.5, 1.5}, {1.6, 2.5}, {3, 4.5}, {4.5, 8}};
  
  // Normalization Shift Width FFT parameters
  float m = 0.5;
  float l_s = 0; 
  float l_e = Data[0] -> Integral(); 
  float l_g = l_e * 0.5;
  
  std::map<TString, std::vector<float>> Params_Bruce; 
  Params_Bruce["r_value"] = {1};
  Params_Bruce["Range_ntrk_1"] = Ranges[0]; 
  Params_Bruce["Range_ntrk_2"] = Ranges[1];   
  Params_Bruce["Range_ntrk_3"] = Ranges[2];  
  Params_Bruce["Range_ntrk_4"] = Ranges[3]; 
  Params_Bruce["fft_cache"] = {10000}; 
  Params_Bruce["Minimizer"] = {10000}; 
  
  Params_Bruce["m_e"] = {m, m, m, m};
  Params_Bruce["m_G"] = {0, 0, 0, 0}; 
  Params_Bruce["m_s"] = {-m, -m, -m, -m}; 

  Params_Bruce["s_s"] = {0.0001, 0.0001, 0.0001, 0.0001};
  Params_Bruce["s_G"] = {0.02, 0.02, 0.02, 0.02};
  Params_Bruce["s_e"] = {0.1, 0.1, 0.1, 0.1};
  
  Params_Bruce["l_s"] = {l_s, l_s, l_s, l_s}; 
  Params_Bruce["l_G"] = {l_g, l_g, l_g, l_g}; 
  Params_Bruce["l_e"] = {l_e, l_e, l_e, l_e}; 

  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(4, trk1_start, JE); 

  float x_min = Data[0] -> GetXaxis() -> GetXmin(); 
  float x_max = Data[0] -> GetXaxis() -> GetXmax(); 
 
  for (int t(0); t < Data.size(); t++)
  {
    std::vector<TH1F*> PDF_H_trk = ntrk_mtru_H[t]; 
    
    RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
    x -> setBins(Params_Bruce["fft_cache"][0], "fft"); 
    x -> setBins(Params_Bruce["fft_cache"][0], "cache"); 
    x -> setRange("trk_1", Params_Bruce["Range_ntrk_1"][0], Params_Bruce["Range_ntrk_1"][1]); 
    x -> setRange("trk_2", Params_Bruce["Range_ntrk_2"][0], Params_Bruce["Range_ntrk_2"][1]); 
    x -> setRange("trk_3", Params_Bruce["Range_ntrk_3"][0], Params_Bruce["Range_ntrk_3"][1]); 
    x -> setRange("trk_4", Params_Bruce["Range_ntrk_4"][0], Params_Bruce["Range_ntrk_4"][1]); 

    // Normalization variables
    std::vector<RooRealVar*> l_vars = GenerateVariableVectors(PDF_H_trk, Params_Bruce, "l"); 
    std::vector<RooRealVar*> m_vars = GenerateVariableVectors(PDF_H_trk, Params_Bruce, "m"); 
    std::vector<RooRealVar*> s_vars = GenerateVariableVectors(PDF_H_trk, Params_Bruce, "s");  
    
    // Create the resolution model: Gaussian 
    std::vector<TString> g_N = NameGenerator(PDF_H_trk, "_Gx");
    std::vector<RooGaussian*> g_vars = RooGaussianVariable(g_N, x, m_vars, s_vars); 

    // Convert the PDFs to RooPDFs
    std::vector<TString> pdf_N_D = NameGenerator(PDF_H_trk, "_D"); 
    std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, PDF_H_trk); 
    std::vector<TString> pdf_N_P = NameGenerator(PDF_H_trk, "_P"); 
    std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 
    
    // Convolve the Gaussian and the PDFs
    std::vector<TString> pxg_N = NameGenerator(PDF_H_trk, "_PxG"); 
    std::vector<RooFFTConvPdf*> PxG_vars = RooFFTVariable(pxg_N, x, g_vars, pdf_P); 

    RooArgList L; 
    RooArgList PxG; 
    for (int i(0); i < PxG_vars.size(); i++)
    {
      PxG.add(*PxG_vars[i]); 
      L.add(*l_vars[i]); 
    }
    // Fitting the PDFs to the Data 
    RooDataHist D("Data", "Data", *x, Data[t]); 
    RooAddPdf model("model", "model", PxG, L); 
    RooAbsReal* nll; // = model.createNLL(D, Extended(true)); 
    int re; // = Minimization(nll, Params_Bruce); 

    //OFF = {false, false, false, false}; 
    //SetConstantVar(l_vars, OFF); 
    //SetConstantVar(m_vars, OFF); 
    //SetConstantVar(s_vars, OFF); 
    nll = model.createNLL(D, Range("trk_1, trk_2, trk_3, trk_4"), Extended(true)); 
    re = Minimization(nll, Params_Bruce); 


    // Fit only the 1-Track, 1-Truth
    OFF = {false, true, true, true}; 
    SetConstantVar(l_vars, OFF); 
    SetConstantVar(m_vars, OFF); 
    SetConstantVar(s_vars, OFF); 
    model.fitTo(D, Range("trk_1"));

    // Fit only the 1-Track, 1-Truth
    OFF = {true, false, true, true}; 
    SetConstantVar(l_vars, OFF); 
    SetConstantVar(m_vars, OFF); 
    SetConstantVar(s_vars, OFF); 
    model.fitTo(D, Range("trk_2")); 

    // Fit only the 1-Track, 1-Truth
    OFF = {true, true, false, true}; 
    SetConstantVar(l_vars, OFF); 
    SetConstantVar(m_vars, OFF); 
    SetConstantVar(s_vars, OFF); 
    model.fitTo(D, Range("trk_3")); 

    // Fit only the 1-Track, 1-Truth
    OFF = {true, true, true, false}; 
    SetConstantVar(l_vars, OFF); 
    SetConstantVar(m_vars, OFF); 
    SetConstantVar(s_vars, OFF); 
    model.fitTo(D, Range("trk_4"));


    TString plot_t = JE + "PULL.pdf"; 
    RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
    
    std::map<TString, std::vector<float>> Output;  
    IntoOutput(&Output, l_vars, m_vars, s_vars, PDF_H_trk, PxG_vars, x, Data[t]);  

    BulkDelete(l_vars); 
    BulkDelete(m_vars); 
    BulkDelete(s_vars); 
    BulkDelete(g_vars); 
    BulkDelete(pdf_D); 
    BulkDelete(pdf_P); 
    BulkDelete(PxG_vars); 
    delete x; 
  }


  return ntrk_mtru_H; 

}



std::vector<std::vector<TH1F*>> IncrementalFit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  auto Minimization =[&] (RooAbsReal* nll, std::map<TString, std::vector<float>> Params)
  {
    RooMinimizer* pg = new RooMinimizer(*nll); 
    pg -> setMaxIterations(Params["Minimizer"][0]); 
    pg -> setMaxFunctionCalls(Params["Minimizer"][0]); 
    pg -> migrad(); 
    pg -> minos();
    pg -> hesse();
    pg -> improve(); 
    pg -> optimizeConst(true); 
    pg -> setEps(1e-6); 
    pg -> setPrintLevel(0); 
    RooFitResult* re = pg -> fit("r"); 
    pg -> cleanup(); 
    delete pg; 
    delete nll;
    
    return re; 
  }; 

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
      for (int p(0); p < Data -> GetNbinsX(); p++){PDF_H[i] -> SetBinError(p+1, 0.);}
    }
  }; 

  auto Algorithm =[&] (std::map<TString, std::vector<float>> Params, 
                       std::vector<std::vector<TH1F*>> ntrk_mtru_H, 
                       std::vector<TH1F*> Data)
  {
    
    float r = 1; 
    if (Params["r_value"].size() != 0){ r = Params["r_value"][0]; } 
    
    float x_min = Data[0] -> GetXaxis() -> GetXmin(); 
    float x_max = Data[0] -> GetXaxis() -> GetXmax(); 

    std::vector<std::map<TString, std::vector<float>>> Output_List; 
    for (int i(0); i < Data.size(); i++)
    {
      std::map<TString, std::vector<float>> Map; 
      Output_List.push_back(Map); 
    }

    for (int i(0); i < Data.size(); i++)
    {
      std::vector<RooRealVar*> l_vars = ProtectionRealVariable("l", ntrk_mtru_H[i], Params, 0, r*(Data[i] -> Integral())); 
      std::vector<RooRealVar*> m_vars = ProtectionRealVariable("m", ntrk_mtru_H[i], Params, -0.001, 0.001); 
      std::vector<RooRealVar*> s_vars = ProtectionRealVariable("s", ntrk_mtru_H[i], Params, 0.0001, 0.001); 

      RooRealVar* x = new RooRealVar("x", "x", x_min, x_max); 
      TString rang;  
      if (Params["fft_cache"].size() != 0)
      {
        x -> setBins(Params["fft_cache"][0], "fft"); 
        x -> setBins(Params["fft_cache"][0], "cache"); 
        rang = "Range_ntrk_"; rang += (i+1); 
        //x -> setRange(rang, Params[rang][0], Params[rang][1]); 
      }

      // Create the resolution model: Gaussian 
      std::vector<TString> g_N = NameGenerator(ntrk_mtru_H[i], "_Gx");
      std::vector<RooGaussian*> g_vars = RooGaussianVariable(g_N, x, m_vars, s_vars); 

      // Convert the PDFs to RooPDFs
      std::vector<TString> pdf_N_D = NameGenerator(ntrk_mtru_H[i], "_D"); 
      std::vector<RooDataHist*> pdf_D = RooDataVariable(pdf_N_D, x, ntrk_mtru_H[i]); 
      std::vector<TString> pdf_N_P = NameGenerator(ntrk_mtru_H[i], "_P"); 
      std::vector<RooHistPdf*> pdf_P = RooPdfVariable(pdf_N_P, x, pdf_D); 
      
      // Convolve the Gaussian and the PDFs
      std::vector<TString> pxg_N = NameGenerator(ntrk_mtru_H[i], "_PxG"); 
      std::vector<RooFFTConvPdf*> PxG_vars = RooFFTVariable(pxg_N, x, g_vars, pdf_P); 

      RooArgList L; 
      RooArgList PxG; 
      for (int i(0); i < PxG_vars.size(); i++)
      {
        PxG.add(*PxG_vars[i]); 
        L.add(*l_vars[i]); 
      }

      // Fitting the PDFs to the Data 
      RooAddPdf model("model", "model", PxG, L); 
      RooDataHist D("Data", "Data", *x, Data[i]); 
      
      std::vector<float> Bools; 
      RooAbsReal* nll;
      RooFitResult* re; 
      // Make it such that only the current ntrk-ntru is allowed to float the fit variables 
      // First fit ==== Release m and not s
      Bools = std::vector<float>(ntrk_mtru_H[i].size(), 1); 
      VariableConstant(Bools, s_vars); 
      Bools = std::vector<float>(ntrk_mtru_H[i].size(), 0); 
      VariableConstant(Bools, m_vars); 

      nll = model.createNLL(D); 
      re = Minimization(nll, Params); 
      delete re; 

      // Second fit ===== Fix m and let s float 
      Bools = std::vector<float>(ntrk_mtru_H[i].size(), 0); 
      VariableConstant(Bools, s_vars); 
      Bools = std::vector<float>(ntrk_mtru_H[i].size(), 1); 
      VariableConstant(Bools, m_vars); 
      
      nll = model.createNLL(D); 
      re = Minimization(nll, Params); 
      delete re;

      // Third fit ==== Fix s and m and let only L float
      Bools = std::vector<float>(ntrk_mtru_H[i].size(), 1); 
      VariableConstant(Bools, m_vars); 
      Bools = std::vector<float>(ntrk_mtru_H[i].size(), 1); 
      VariableConstant(Bools, s_vars); 

      nll = model.createNLL(D); 
      re = Minimization(nll, Params); 
      CaptureResults(re, &(Output_List[i]));

      TString plot_t = JE + "_trk_"; plot_t += (i+1); plot_t += "PULL.pdf"; 
      RooFitPullPlot(model, x, pdf_P, &D, plot_t); 
      
      IntoOutput(&(Output_List[i]), l_vars, m_vars, s_vars, ntrk_mtru_H[i], PxG_vars, x, Data[i]); 

      BulkDelete(l_vars); 
      BulkDelete(m_vars); 
      BulkDelete(s_vars); 
      BulkDelete(g_vars); 
      BulkDelete(pdf_D); 
      BulkDelete(pdf_P); 
      BulkDelete(PxG_vars); 
      delete x; 
    }

    return Output_List;
  };

  
  auto FitBins =[&] (std::vector<TH1F*> ntrk, TH1F* Data)
  {

    for (int i(0); i < Data -> GetNbinsX(); i++)
    {
      float e = Data -> GetBinContent(i+1); 
      float sum = 0; 
      for (int x(0); x < ntrk.size(); x++)
      {
        sum += ntrk[x] -> GetBinContent(x+1); 
      } 
      
      if (sum == 0){continue;}
      float r = e/sum; 
      for (int x(0); x < ntrk.size(); x++)
      {
        float p = ntrk[x] -> GetBinContent(x+1); 
        ntrk[x] -> SetBinContent(x+1, r*p); 
      } 
    }
  }; 
  


  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, JE); 
  std::vector<std::map<TString, std::vector<float>>> Output_List;   
  for (int i(0); i < 2; i++)
  {
    Output_List = Algorithm(Params, ntrk_mtru_H, Data);
  }

  TString ext = "_" + JE + "_Incremental";
  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/Incremental"); 
  
  for (int i(0); i < Data.size(); i++)
  {
    TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
    WriteHistsToFile(ntrk_mtru_H[i], JE + "/Incremental"); 
    WriteOutputMapToFile(Output_List[i], JE + "/Incremental", trk_n); 
  }

  return ntrk_mtru_H; 
}


