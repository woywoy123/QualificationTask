#include<PostAnalysis/AlgorithmFunctions.h>
#include<TGraphSmooth.h>

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
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext, 4);

  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/ShiftNormal"); 

  for (int i(0); i < Data.size(); i++)
  {
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

  auto Smooth =[&] (std::vector<TH1F*> Hists)
  {
    for (int i(0); i < Hists.size(); i++)
    {
      TGraph* gr = new TGraph(Hists[i]); 
      TGraphSmooth* g = new TGraphSmooth("Smoother"); 
      TGraph* h = g -> SmoothKern(gr, "normal", 0.05); 

      int n = h -> GetN(); 
      for (int p(0); p < n; p++)
      {
        double x, y; 
        h -> GetPoint(p, x, y); 
        Hists[i] -> SetBinContent(p+1, y);
      }
      Average(Hists[i]);

      delete gr; 
      delete h;
    }
  }; 

  TString ext = "_" + JE + "_Experimental_NtrkMtru"; 
  gDirectory -> cd("/"); 
  gDirectory -> mkdir(JE + "/Experimental"); 
  std::vector<MVF> Params_V(Data.size(), Params);
  
  TH1F* trk1 = (TH1F*)trk1_start -> Clone("x"); 
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1, ext, Data.size());
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy();
  int iter = 3;
  for (int x(0); x < iter; x++)
  {
   
    for (int i(0); i < Data.size(); i++)
    {
      TH1F* D_t = (TH1F*)Data[i] -> Clone("D"); 
      MVF Map = NormalizationShift(D_t, ntrk_mtru_H[i], Params, JE); 
      PlotHists(D_t, ntrk_mtru_H[i], can); 
      can -> Print("Fit.pdf"); 
      if (x == iter -1)
      {
        TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
        WriteOutputMapToFile(Map, JE + "/Experimental", trk_n);  
      }
      delete D_t;
    }
    
    for (int i(0); i < ntrk_mtru_H.size(); i++)
    {
      for (int j(0); j < ntrk_mtru_H[i].size(); j++)
      {
        if (j > Data.size()-1){ continue; }
        TH1F* Update = (TH1F*)ntrk_mtru_H[i][i] -> Clone("CX"); 

        float L = ntrk_mtru_H[j][i] -> Integral(); 
        if (L <= 100){ continue; }
        ntrk_mtru_H[j][i] -> Reset(); 
        ntrk_mtru_H[j][i] -> FillRandom(Update, L); 
        delete Update;
      }
    }  
    
    for (int i(0); i < Data.size(); i++)
    {
      TH1F* D_t = (TH1F*)Data[i] -> Clone("D"); 
      SubtractData(ntrk_mtru_H[i], D_t, i, false); 
      Average(D_t);
      ntrk_mtru_H[i][i] -> Reset(); 
      ntrk_mtru_H[i][i] -> FillRandom(D_t, D_t -> Integral()); 
      delete D_t; 
    }
  }

  for (int i(0); i < ntrk_mtru_H.size(); i++){ WriteHistsToFile(ntrk_mtru_H[i], JE + "/Experimental"); }
  
  delete trk1; 

  return ntrk_mtru_H; 
}

std::vector<std::vector<TH1F*>> Simultaneous_Fit_NtrkMtru(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  TString ext = "_" + JE + "_Simultaneous_NtrkMtru";
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

std::vector<std::vector<TH1F*>> IncrementalFit(std::vector<TH1F*> Data, TH1F* trk1_start, std::map<TString, std::vector<float>> Params, TString JE)
{
  auto Algorithm =[&] (std::map<TString, std::vector<float>> Params, 
                       std::vector<std::vector<TH1F*>> ntrk_mtru_H, 
                       std::vector<TH1F*> Data)
  {
    std::vector<std::map<TString, std::vector<float>>> Output_List; 
    for (int i(0); i < Data.size(); i++)
    {
      TString rang = "Range_ntrk_"; rang += (i+1); 
      Params["Range"] = Params[rang]; 
      std::map<TString, std::vector<float>> Map = IncrementalFFT(Data[i], ntrk_mtru_H[i], Params); 
      Output_List.push_back(Map); 
    }
    return Output_List;
  };

  TString ext = "_" + JE + "_Incremental_NtrkMtru";
  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(Data.size(), trk1_start, ext); 
  std::vector<std::map<TString, std::vector<float>>> Output_List;   
  Output_List = Algorithm(Params, ntrk_mtru_H, Data);

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
