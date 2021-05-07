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

std::vector<std::vector<TH1F*>> BuildNtrkMtru(int n, TH1F* trk1_start, TString extension)
{
  std::vector<std::vector<TString>> ntrk_ntru_names; 
  for (int i(0); i < n; i++)
  {
    std::vector<TString> ntrk_ntru_n; 
    for (int j(0); j < 4; j++)
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

  for (int x(0); x < 3; x++)
  {
    //if (x == 1)
    //{
    //  for (int t(0); t < Data.size(); t++)
    //  {
    //    TString name_temp = Data[t] -> GetTitle(); 
    //    TH1F* ntrk_Data = (TH1F*)Data[t] -> Clone(name_temp + "_C"); 
    //    for (int j(0); j < ntrk_mtru_H[t].size(); j++)
    //    {
    //      if (j != t){ ntrk_Data -> Add(ntrk_mtru_H[t][j], -1); }
    //    }
    //    
    //    Flush({ntrk_Data}, {ntrk_mtru_H[t][t]}); 
    //    delete ntrk_Data;  
    //  }

    //  for (int t(0); t < Data.size(); t++)
    //  {
    //    for (int j(0); j < Data.size(); j++)
    //    {
    //      Flush({ntrk_mtru_H[t][t]}, {ntrk_mtru_H[j][t]}); 
    //    }
    //    Normalize(ntrk_mtru_H[t]);
    //  }
    //}   
    
    for (int i(0); i < Data.size(); i++)
    {
      TString r_name = "Range_ntrk_"; r_name += (i+1); 
      Params["Range"] = Params[r_name]; 
      std::vector<TH1F*> ntrk_Template = ntrk_mtru_H[i]; 
      TH1F* ntrk_Measure = Data[i]; 
      
      TString base = "Fit_"; base += (i+1); base += (ext); 
      
      std::map<TString, std::vector<float>> Map = ConvolutionFFT(ntrk_Measure, ntrk_Template, Params, base); 

      if (x == 2)
      {
        TString trk_n = "ntrk_"; trk_n += (i+1); trk_n += ("_error"); 
        WriteOutputMapToFile(Map, JE + "/Experimental", trk_n);  
      }
    }
  }

  for (int i(0); i < ntrk_mtru_H.size(); i++){ WriteHistsToFile(ntrk_mtru_H[i], JE + "/Experimental"); }


  return ntrk_mtru_H; 
}
