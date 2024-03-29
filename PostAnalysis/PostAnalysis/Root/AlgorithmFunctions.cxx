#include<PostAnalysis/AlgorithmFunctions.h>
#include<TGraphSmooth.h>

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

void Experimental(std::vector<TH1F*> Data, std::vector<std::vector<TH1F*>> ntrk_mtru, std::map<TString, std::vector<float>> Params)
{

  for (int it(0); it < 3; it++)
  {
    for (int i(0); i < Data.size(); i++)
    {
      TH1F* data = (TH1F*)Data[i] -> Clone("Data_TMP"); 

      ConvolutionFFT(data, ntrk_mtru[i], Params); 

      // Check if the bins are messed 
      Average(ntrk_mtru[i]);
      
      SubtractData(ntrk_mtru[i], data, i, false); 
      Average(data); 
    
      float L = ntrk_mtru[i][i] -> Integral(); 
      ntrk_mtru[i][i] -> Reset(); 
      ntrk_mtru[i][i] -> FillRandom(data, L); 

      delete data;
    }
    
    for (int i(0); i < ntrk_mtru.size(); i++)
    {
      
      TH1F* ntrkmtru = (TH1F*)ntrk_mtru[i][i] -> Clone("TMP");
      for (int x(0); x < ntrk_mtru.size(); x++)
      {
        float L = ntrk_mtru[x][i] -> Integral(); 
        ntrk_mtru[x][i] -> Reset(); 
        ntrk_mtru[x][i] -> FillRandom(ntrkmtru, L); 
      }
      delete ntrkmtru; 
    }
  }
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
