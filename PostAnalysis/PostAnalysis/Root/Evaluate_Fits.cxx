#include<PostAnalysis/Evaluate_Fits.h>
#include<iostream>
#include<fstream>

bool TruthLinker(TH1F* H, TH1F* V)
{
  bool Link = false; 
  TString Truth = H -> GetTitle(); 
  TString Recon = V -> GetTitle(); 
  for (int i(0); i < 4; i++)
  {
    for (int j(0); j < 4; j++)
    {
      TString Temp = "NTRU_"; Temp += (i+1); Temp += "_NTRK"; Temp += (j+1); 
      if (Recon.Contains(Temp))
      {
        TString Temp2 = "ntrk_"; Temp2 += (j+1); 
        TString Temp3 = "ntru_"; Temp3 += (i+1); 
        if (Truth.Contains(Temp2) && Truth.Contains(Temp3)){Link = true;}
      }
    }
  }
  return Link; 
}

void CoutText(TString *Input, int v, TString Text = "")
{
  
  for (int c(0); c < v; c++)
  {
    *Input += (Text); 
  }
}

TString PrecisionString(float number, int precision)
{
  TString out; 
  std::ostringstream o; 
  o.precision(precision); 
  o << std::scientific << number; 
  out += (o.str()); 
  return out; 
}

void Plotting_Fits(std::map<TString, std::vector<TH1F*>> Maps, std::map<TString, std::vector<TH1F*>> truth, TString name)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print(name + ".pdf["); 
  for (iz x = Maps.begin(); x != Maps.end(); x++)
  {
    TString na = x -> first; 
    std::vector<TH1F*> Hist = Maps[na]; 
    std::vector<TH1F*> tru = truth[na]; 

    PlotHists(Hist, tru, can); 
    can -> Print(name + ".pdf"); 

  }
  can -> Print(name + ".pdf]"); 
  delete can; 
}

void Evaluate_Fits(TString Filename)
{

  std::vector<TString> Layer_H = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy_H = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                  "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                  "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                  "higher_GeV"};  
  
  std::vector<TString> Algorithms = {"Normal", "ShiftNormal", "ShiftNormalFFT", "NormalShiftWidthFFT"};  


  // Define the different Algorithm maps
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Normal; 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> ShiftNormal; 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> ShiftNormalFFT; 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> NormalShiftWidthFFT; 

  std::map<TString, std::map<TString, std::vector<float>>> Normal_Error; 
  std::map<TString, std::map<TString, std::vector<float>>> ShiftNormal_Error; 
  std::map<TString, std::map<TString, std::vector<float>>> ShiftNormalFFT_Error; 
  std::map<TString, std::map<TString, std::vector<float>>> NormalShiftWidthFFT_Error; 

  std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> Results = ReadResults(Filename);
  std::map<TString, std::vector<TH1F*>> Truth = MC_Reader_All("Merged_MC.root"); 
  for (it i = Results.begin(); i != Results.end(); i++)
  {
    TString key_L_JE_AL = i -> first; 

    TString layer; 
    TString jetenergy; 
    TString algorithm;  

    // Get the attributes of the string 
    for (TString L : Layer_H)
    {
      for (TString JE : JetEnergy_H)
      {
        for (TString Alg : Algorithms)
        {
          TString Name = L + "_" + JE + "_" + Alg; 
          if(key_L_JE_AL == Name)
          {
            layer = L; 
            jetenergy = JE; 
            algorithm = Alg; 
          }
        }
      }
    }
   
    std::vector<TH1F*> v_T = Truth[layer + "_" + jetenergy + "_radius_Less_Truth"]; 
    std::vector<std::pair<TString, std::vector<TH1F*>>> v_ = Results[key_L_JE_AL]; 
  
    // This section links the truth histograms to the reconstructed histograms 
    std::map<TString, std::vector<TH1F*>> Temp_map_R; 
    std::map<TString, std::vector<TH1F*>> Temp_map_T; 

    std::map<TString, std::vector<float>> Temp_Flost_R; 
    std::map<TString, std::vector<float>> Temp_Flost_T; 
    std::vector<std::vector<TH1F*>> ntrk_ntru_R; 
    std::vector<std::vector<TH1F*>> ntrk_ntru_T; 
    std::vector<TH1F*> Err;   
    for (int p(0); p < v_.size(); p++)
    {
      TString ntrk = "dEdx_ntrk_"; ntrk += (p+1); 
      TString ntrk_T = "dEdx_ntrk_"; ntrk_T += (p+1); ntrk_T += ("_T"); 

      std::pair<TString, std::vector<TH1F*>> pairs_ = v_[p]; 
     
      std::vector<TH1F*> Temp_V; 
      std::vector<TH1F*> Temp_T; 
      for (TH1F* H_R : pairs_.second)
      {
        for (TH1F* H_T : v_T)
        {
          bool cap = TruthLinker(H_T, H_R); 
          if (cap == true)
          {
            Temp_V.push_back(H_R); 
            Temp_T.push_back(H_T); 
          }
        }
        TString n_H = H_R -> GetTitle(); 
        if (n_H.Contains("Error"))
        {
          Err.push_back(H_R); 
        }
      }
      Temp_map_R[ntrk] = Temp_V; 
      Temp_map_T[ntrk] = Temp_T; 

      ntrk_ntru_R.push_back(Temp_V); 
      ntrk_ntru_T.push_back(Temp_T); 
    }
   
    // Evaluate the Flost 2 and Flost 3
    Temp_Flost_R["FLost2"] = Flost2(ntrk_ntru_R, Err); 
    Temp_Flost_T["FLost2"] = Flost2(ntrk_ntru_T, {}); 
    Temp_Flost_R["FLost3"] = Flost3(ntrk_ntru_R, Err); 
    Temp_Flost_T["FLost3"] = Flost3(ntrk_ntru_T, {}); 

    // Populate the maps 
    if (algorithm == "Normal")
    {
      Normal[layer + "_" + jetenergy] = Temp_map_R; 
      Normal[layer + "_" + jetenergy + "_T"] = Temp_map_T; 
      
      Normal_Error[layer + "_" + jetenergy] = Temp_Flost_R; 
      Normal_Error[layer + "_" + jetenergy + "_T"] = Temp_Flost_T; 
    }
    if (algorithm == "ShiftNormal")
    {
      ShiftNormal[layer + "_" + jetenergy] = Temp_map_R; 
      ShiftNormal[layer + "_" + jetenergy + "_T"] = Temp_map_T; 

      ShiftNormal_Error[layer + "_" + jetenergy] = Temp_Flost_R; 
      ShiftNormal_Error[layer + "_" + jetenergy + "_T"] = Temp_Flost_T; 
    }
    if (algorithm == "ShiftNormalFFT")
    {
      ShiftNormalFFT[layer + "_" + jetenergy] = Temp_map_R; 
      ShiftNormalFFT[layer + "_" + jetenergy + "_T"] = Temp_map_T; 
 
      ShiftNormalFFT_Error[layer + "_" + jetenergy] = Temp_Flost_R; 
      ShiftNormalFFT_Error[layer + "_" + jetenergy + "_T"] = Temp_Flost_T; 
    }
    if (algorithm == "NormalShiftWidthFFT")
    {
      NormalShiftWidthFFT[layer + "_" + jetenergy] = Temp_map_R; 
      NormalShiftWidthFFT[layer + "_" + jetenergy + "_T"] = Temp_map_T; 
  
      NormalShiftWidthFFT_Error[layer + "_" + jetenergy] = Temp_Flost_R; 
      NormalShiftWidthFFT_Error[layer + "_" + jetenergy + "_T"] = Temp_Flost_T; 
    }
  }
 
  std::vector<TString> COUT_String; 
  std::vector<int> bounds; 
  for (ip z = Normal.begin(); z != Normal.end(); z++)
  {
    TString name = z -> first; 
    if (name.Contains("_T")){continue;}
    
    // Reconstructed Maps 
    std::map<TString, std::vector<TH1F*>> Normal_Map = Normal[name]; 
    std::map<TString, std::vector<TH1F*>> ShiftNormal_Map = ShiftNormal[name]; 
    std::map<TString, std::vector<TH1F*>> ShiftNormalFFT_Map = ShiftNormalFFT[name]; 
    std::map<TString, std::vector<TH1F*>> NormalShiftWidthFFT_Map = NormalShiftWidthFFT[name]; 

    std::map<TString, std::vector<float>> Normal_Map_E = Normal_Error[name]; 
    std::map<TString, std::vector<float>> ShiftNormal_Map_E = ShiftNormal_Error[name]; 
    std::map<TString, std::vector<float>> ShiftNormalFFT_Map_E = ShiftNormalFFT_Error[name]; 
    std::map<TString, std::vector<float>> NormalShiftWidthFFT_Map_E = NormalShiftWidthFFT_Error[name]; 

   
    // Truth Maps 
    std::map<TString, std::vector<TH1F*>> Normal_Map_T = Normal[name + "_T"]; 
    std::map<TString, std::vector<TH1F*>> ShiftNormal_Map_T = ShiftNormal[name + "_T"]; 
    std::map<TString, std::vector<TH1F*>> ShiftNormalFFT_Map_T = ShiftNormalFFT[name + "_T"]; 
    std::map<TString, std::vector<TH1F*>> NormalShiftWidthFFT_Map_T = NormalShiftWidthFFT[name + "_T"]; 
    
    std::map<TString, std::vector<float>> Normal_Map_E_T = Normal_Error[name + "_T"]; 
    std::map<TString, std::vector<float>> ShiftNormal_Map_E_T = ShiftNormal_Error[name + "_T"]; 
    std::map<TString, std::vector<float>> ShiftNormalFFT_Map_E_T = ShiftNormalFFT_Error[name + "_T"]; 
    std::map<TString, std::vector<float>> NormalShiftWidthFFT_Map_E_T = NormalShiftWidthFFT_Error[name + "_T"]; 

    // Plot the hists
    Plotting_Fits(Normal_Map, Normal_Map_T, name+"_Normal"); 
    Plotting_Fits(ShiftNormal_Map, ShiftNormal_Map_T, name+"_ShiftNormal"); 
    Plotting_Fits(ShiftNormalFFT_Map, ShiftNormalFFT_Map_T, name+"_ShiftNormalFFT"); 
    Plotting_Fits(NormalShiftWidthFFT_Map, NormalShiftWidthFFT_Map_T, name+"_NormalShiftWidthFFT"); 

    // Return the FLost values
    // ==== Measured
    float FLost2_Normal_R = Normal_Map_E["FLost2"][0]; // FLost2 value
    float FLost2_Normal_R_E = Normal_Map_E["FLost2"][1]; // FLost2 Error value 

    float FLost2_ShiftNormal_R = ShiftNormal_Map_E["FLost2"][0]; 
    float FLost2_ShiftNormal_R_E = ShiftNormal_Map_E["FLost2"][1]; 

    float FLost2_ShiftNormalFFT_R = ShiftNormalFFT_Map_E["FLost2"][0]; 
    float FLost2_ShiftNormalFFT_R_E = ShiftNormalFFT_Map_E["FLost2"][1]; 

    float FLost2_NormalShiftWidthFFT_R = NormalShiftWidthFFT_Map_E["FLost2"][0]; 
    float FLost2_NormalShiftWidthFFT_R_E = NormalShiftWidthFFT_Map_E["FLost2"][1]; 


    float FLost3_Normal_R = Normal_Map_E["FLost3"][0]; // FLost3 value
    float FLost3_Normal_R_E = Normal_Map_E["FLost3"][1]; // FLost3 Error value 

    float FLost3_ShiftNormal_R = ShiftNormal_Map_E["FLost3"][0]; 
    float FLost3_ShiftNormal_R_E = ShiftNormal_Map_E["FLost3"][1]; 

    float FLost3_ShiftNormalFFT_R = ShiftNormalFFT_Map_E["FLost3"][0]; 
    float FLost3_ShiftNormalFFT_R_E = ShiftNormalFFT_Map_E["FLost3"][1]; 

    float FLost3_NormalShiftWidthFFT_R = NormalShiftWidthFFT_Map_E["FLost3"][0]; 
    float FLost3_NormalShiftWidthFFT_R_E = NormalShiftWidthFFT_Map_E["FLost3"][1]; 


    // ==== Truth
    float FLost2_Normal_T = Normal_Map_E_T["FLost2"][0]; 
    float FLost2_ShiftNormal_T = ShiftNormal_Map_E_T["FLost2"][0]; 
    float FLost2_ShiftNormalFFT_T = ShiftNormalFFT_Map_E_T["FLost2"][0]; 
    float FLost2_NormalShiftWidthFFT_T = NormalShiftWidthFFT_Map_E_T["FLost2"][0]; 

    float FLost3_Normal_T = Normal_Map_E_T["FLost3"][0]; 
    float FLost3_ShiftNormal_T = ShiftNormal_Map_E_T["FLost3"][0]; 
    float FLost3_ShiftNormalFFT_T = ShiftNormalFFT_Map_E_T["FLost3"][0]; 
    float FLost3_NormalShiftWidthFFT_T = NormalShiftWidthFFT_Map_E_T["FLost3"][0]; 

    TString f = ("|"); 
    TString out = "";  
    if (COUT_String.size() < 1)
    { 
      CoutText(&out, 25, " "); 
      out += (f);      
      int Len = 20;
      for (TString al : Algorithms)
      {
        int i_s = out.Sizeof(); 
        int l_s = Len - al.Sizeof() / 2;  
        CoutText(&out, l_s, " "); 
        out += (al);  

        int l_e = Len - al.Sizeof() / 2; 
        CoutText(&out, l_e, " "); 
        
        int i_e = out.Sizeof(); 
        int del = i_e - i_s; 
        bounds.push_back(del);
       
        out += (f); 
      }
      COUT_String.push_back(out); 
      out.Clear(); 

      TString out = ""; 
      CoutText(&out, 25, " "); 
      out += (f); 
      for (int a(0); a < Algorithms.size(); a++)
      {
        int d1 = out.Sizeof(); 
        int b = bounds[a];   
        
        std::vector<TString> stats_test = {"KS", "CS", "IE"};
        for (TString t : stats_test)
        {
          int buf = b / 6; 
          CoutText(&out, buf-1, " "); 
          out += (t); 
          CoutText(&out, buf-1, " "); 
        }
        int dif = std::abs(out.Sizeof() - d1 -b); 
        CoutText(&out, dif, " "); 
        out+=(f);
      }
      COUT_String.push_back(out); 
      out.Clear();
    } 
    
    int Len = 25 - name.Sizeof(); 
    CoutText(&out, Len, " "); 
    out += (name); 
    out += (" ");
    out += (f); 
 
    for (int al(0); al < Algorithms.size(); al++)
    {
      int b = bounds[al]; 
      CoutText(&out, b, " ");
      out += (f);
    }
    COUT_String.push_back(out); 

    out.Clear(); 
    for (iz v = Normal_Map.begin(); v != Normal_Map.end(); v++)
    {
      TString ntrk = v -> first; 
      std::vector<TH1F*> ntrks_Normal = Normal_Map[ntrk]; 
      std::vector<TH1F*> ntrks_ShiftNormal = ShiftNormal_Map[ntrk]; 
      std::vector<TH1F*> ntrks_ShiftNormalFFT = ShiftNormalFFT_Map[ntrk]; 
      std::vector<TH1F*> ntrks_NormalShiftWidthFFT = NormalShiftWidthFFT_Map[ntrk]; 

      std::vector<TH1F*> ntrk_ntru = Normal_Map_T[ntrk]; 
      
      for (int H_i(0); H_i < ntrk_ntru.size(); H_i++)
      {
        TString ntrk_ntru_s = ntrk + "_tru_"; ntrk_ntru_s += (H_i +1); 
        Len = 25 - ntrk_ntru_s.Sizeof(); 
        CoutText(&out, Len, " ");  
        out += (ntrk_ntru_s);  
        out += (" ");
        out += (f); 
        
        TH1F* truth = ntrk_ntru[H_i]; 
        TH1F* normal_R = ntrks_Normal[H_i];  
        TH1F* shiftnormal_R = ntrks_ShiftNormal[H_i];  
        TH1F* shiftnormalfft_R = ntrks_ShiftNormalFFT[H_i];  
        TH1F* normalshiftwidthfft_R = ntrks_NormalShiftWidthFFT[H_i];  
        
        float ks_n = truth -> KolmogorovTest(normal_R); 
        float ks_sn = truth -> KolmogorovTest(shiftnormal_R); 
        float ks_snf = truth -> KolmogorovTest(shiftnormalfft_R); 
        float ks_nswf = truth -> KolmogorovTest(normalshiftwidthfft_R); 

        float chi_n = truth -> Chi2Test(normal_R); 
        float chi_sn = truth -> Chi2Test(shiftnormal_R); 
        float chi_snf = truth -> Chi2Test(shiftnormalfft_R); 
        float chi_nswf = truth -> Chi2Test(normalshiftwidthfft_R); 

        float ie_n = ErrorByIntegral(normal_R, truth); 
        float ie_sn = ErrorByIntegral(shiftnormal_R, truth); 
        float ie_snf = ErrorByIntegral(shiftnormalfft_R, truth); 
        float ie_nswf = ErrorByIntegral(normalshiftwidthfft_R, truth); 
        
        std::vector<float> Alg_ks = {ks_n, ks_sn, ks_snf, ks_nswf};  
        std::vector<float> Alg_chi = {chi_n, chi_sn, chi_snf, chi_nswf};  
        std::vector<float> Alg_ie = {ie_n, ie_sn, ie_snf, ie_nswf};  
        
        for (int c(0); c < Alg_ks.size(); c++)
        {
          TString ks; 
          TString chi; 
          TString ie; 

          std::ostringstream ssks; 
          std::ostringstream sschi; 
          std::ostringstream ssie; 
         
          ssks.precision(4); 
          sschi.precision(4); 
          ssie.precision(4); 
         
          ssks << std::scientific << Alg_ks[c];
          sschi << std::scientific << Alg_chi[c];
          ssie << std::scientific << Alg_ie[c];
          
          TString temp_s = ""; 
          CoutText(&temp_s, 2, " ");          
          ks += (ssks.str()); 
          temp_s += (ks); 
          CoutText(&temp_s, 2, " "); 
          chi += (sschi.str()); 
          temp_s += (chi); 
          CoutText(&temp_s, 2, " "); 
          ie += (ssie.str()); 
          temp_s += (ie);
          CoutText(&temp_s, 2, " "); 
          
          int del = bounds[c] - temp_s.Sizeof(); 
          out += (temp_s); 
          CoutText(&out, del+1, " "); 
          out += ("|");
        }
        COUT_String.push_back(out); 
        out.Clear();
      }
    }
    out += ("FLost 2 Truth: "); 
    out += (PrecisionString(FLost2_Normal_T, 4));  
    Len = 26 - out.Sizeof(); 
    CoutText(&out, Len, " ");  
    out += (f); 
   
    std::vector<float> FLost2_V = {FLost2_Normal_R, FLost2_ShiftNormal_R, FLost2_ShiftNormalFFT_R, FLost2_NormalShiftWidthFFT_R}; 
    std::vector<float> FLost2_V_E = {FLost2_Normal_R_E, FLost2_ShiftNormal_R_E, FLost2_ShiftNormalFFT_R_E, FLost2_NormalShiftWidthFFT_R_E}; 
    
    bool skip = false; 
    for (int x(0); x < FLost2_V.size(); x++)
    {
      if (std::isnan(FLost2_V[x])){skip = true;} 
    }
    if (skip) {continue;}

    for (int x(0); x < FLost2_V.size(); x++)
    {
      TString s_ = ""; 
      s_ += ("  ");  
      s_ += ("FLost2: "); 
  
      TString fl2 = PrecisionString(FLost2_V[x], 4); 
      TString fl2_e = PrecisionString(FLost2_V_E[x], 4); 
      s_ += (fl2); 
      s_ += (" +- "); 
      s_ += (fl2_e); 
      
      int del = bounds[x] - s_.Sizeof(); 
      out += (s_); 
      CoutText(&out, del, " "); 
      s_.Clear(); 
      
      out += (" |"); 
    }
    COUT_String.push_back(out); 
    out.Clear(); 

    out += ("FLost 3 Truth: "); 
    out += (PrecisionString(FLost3_Normal_T, 4));  
    Len = 26 - out.Sizeof(); 
    CoutText(&out, Len, " ");  
    out += (f); 
   
    std::vector<float> FLost3_V = {FLost3_Normal_R, FLost3_ShiftNormal_R, FLost3_ShiftNormalFFT_R, FLost3_NormalShiftWidthFFT_R}; 
    std::vector<float> FLost3_V_E = {FLost3_Normal_R_E, FLost3_ShiftNormal_R_E, FLost3_ShiftNormalFFT_R_E, FLost3_NormalShiftWidthFFT_R_E}; 
  
    skip = false; 
    for (int x(0); x < FLost3_V.size(); x++)
    {
      if (std::isnan(FLost3_V[x])){skip = true;} 
    }
    if (skip) {continue;}

  
    for (int x(0); x < FLost3_V.size(); x++)
    {
      TString s_ = ""; 
      s_ += ("  ");  
      s_ += ("FLost3: "); 
  
      TString fl3 = PrecisionString(FLost3_V[x], 4); 
      TString fl3_e = PrecisionString(FLost3_V_E[x], 4); 
      s_ += (fl3); 
      s_ += (" +- "); 
      s_ += (fl3_e); 
      
      int del = bounds[x] - s_.Sizeof(); 
      out += (s_); 
      CoutText(&out, del, " "); 
      s_.Clear(); 
      
      out += (" |"); 
    }
    COUT_String.push_back(out); 
    out.Clear();
    
    COUT_String.push_back(out); 
  }

  for (int x(0); x < COUT_String.size(); x++)
  {
    std::cout << COUT_String[x] << std::endl;
  }


}
