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
  o << std::fixed << number; 
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

  std::ofstream myfile; 
  myfile.open("results.txt");  
  for (int x(0); x < COUT_String.size(); x++)
  {
    std::cout << COUT_String[x] << std::endl;
    
    myfile << COUT_String[x] << "\n"; 
  }
  myfile.close();

}

void EvaluateErrorImpact(TString Filename)
{
  std::vector<TString> Layer_H = {"IBL", "Blayer", "layer1", "layer2"}; 
  std::vector<TString> JetEnergy_H = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                  "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                  "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                  "higher_GeV"};  
  
  std::vector<TString> Algorithms = {"Normal", "ShiftNormal", "ShiftNormalFFT", "NormalShiftWidthFFT", "ExperimentalShift"};  
  std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> Results = ReadResults(Filename);
  std::map<TString, std::vector<TH1F*>> Truth = MC_Reader_All("Merged_MC.root"); 
  
  std::map<TString, std::vector<std::pair<TString, std::map<TString, std::vector<float>>>>> Comparison; 
  for (it i = Results.begin(); i != Results.end(); i++)
  {
    TString LJEA_Name = i -> first; 
    
    // Check the jet energy, Algorithm, layer 
    TString Alg; 
    TString Je; 
    TString L; 
    for (int x(0); x < Layer_H.size(); x++)
    {
      if (!LJEA_Name.Contains(Layer_H.at(x))){ continue; }
      for (int j(0); j < JetEnergy_H.size(); j++)
      {
        if (!LJEA_Name.Contains(JetEnergy_H.at(j))){ continue; }
        for (int p(0); p < Algorithms.size(); p++)
        {
          if (!LJEA_Name.Contains(Algorithms.at(p))){ continue; }
          
          Alg = Algorithms.at(p); 
          Je = JetEnergy_H.at(j); 
          L = Layer_H.at(x); 

        }
      }
    }
    //std::cout << LJEA_Name << " ::: "  << L << " " << Je << " " << Alg << std::endl;
    
    // Get the truth histograms 
    std::vector<TH1F*> Truth_H = Truth[L + "_" + Je + "_radius_Less_Truth"]; 
    std::vector<std::pair<TString, std::vector<TH1F*>>>  Result_H = Results[LJEA_Name]; 
 
    std::map<TString, std::vector<float>> Fits; 
    for (int x(0); x < Result_H.size(); x++)
    {
      TString ntrk = "dEdx_ntrk_"; ntrk += (x+1); 
      std::pair<TString, std::vector<TH1F*>> Result_trk = Result_H[x]; 
      
      std::vector<TH1F*> Temp_R; 
      std::vector<TH1F*> Temp_T; 
      for (TH1F* H_R : Result_trk.second)
      {
        for (TH1F* H_T : Truth_H)
        {
          bool cap = TruthLinker(H_T, H_R); 
          if (cap == true)
          {
            Temp_R.push_back(H_R); 
            Temp_T.push_back(H_T); 
          }
        }
      }

      // Fit comparison 
      TH1F* Truth_Sum = SumHists(Temp_T, "Truth_Sum"); 
      TH1F* Recon_Sum = SumHists(Temp_R, "Recon_Sum"); 
      float inte_T = Truth_Sum -> Integral(); 
      float inte_R = Recon_Sum -> Integral(); 
      float r = inte_R / inte_T; 
      float diff = std::abs(inte_T - inte_R) / inte_T;  

      Fits[ntrk + "_Integral"].push_back(inte_T);
      Fits[ntrk + "_Integral"].push_back(inte_R);
      Fits[ntrk + "_Integral"].push_back(r);
      Fits[ntrk + "_Integral"].push_back(diff); 
      Fits[ntrk + "_Integral"].push_back(ErrorByIntegral(Recon_Sum, Truth_Sum)); 
      delete Truth_Sum; 
      delete Recon_Sum; 

      for (int j(0); j < Temp_T.size(); j++)
      {
        TString trk_tru = ntrk + "_ntru_"; trk_tru += (j+1); 
        float integ_R = Temp_R[j] -> Integral(); 
        float integ_T = Temp_T[j] -> Integral(); 
        float diff = std::abs(integ_T - integ_R) / inte_T; 
        
        Fits[trk_tru + "_Integral"].push_back(integ_R); 
        Fits[trk_tru + "_Integral"].push_back(integ_T); 
        Fits[trk_tru + "_Integral"].push_back(diff); 
      }
    }

    std::pair<TString, std::map<TString, std::vector<float>>> Alg_Comp_Pair; 
    Alg_Comp_Pair.first = Alg; 
    Alg_Comp_Pair.second = Fits; 
    Comparison[L + "_" + Je].push_back(Alg_Comp_Pair); 
  }


  std::map<TString, std::vector<int>> Perform_Matrix;  
  for (TString alg : Algorithms)
  {
    std::vector<int> trk_S; 
    for (int t(0); t < 4; t++)
    {
      trk_S.push_back(0); 
    }
    Perform_Matrix[alg] = trk_S; 
  }


  for (xi z = Comparison.begin(); z != Comparison.end(); z++)
  {
    std::vector<std::pair<TString, std::map<TString, std::vector<float>>>> Algs_Result = z -> second; 

    std::vector<float> error_V = {100, 100, 100, 100, 100}; 
    std::vector<TString> Alg = {"", "", "", "", ""}; 
    for (std::pair<TString, std::map<TString, std::vector<float>>> u : Algs_Result)
    {
      std::map<TString, std::vector<float>> Error_Map = u.second; 
     
      int p = 0; 
      for (ix h = Error_Map.begin(); h != Error_Map.end(); h++)
      {
        float error = error_V[p]; 
        TString ntrk_Int = h -> first; 
        std::vector<float> g = Error_Map[ntrk_Int]; 
        
        if (ntrk_Int.Contains("tru")){ continue; }
        //std::cout << g[4] << "  " << g[3] << "  " << ntrk_Int << std::endl;
        
        if (error > g[4])
        { 
          error_V[p] = g[4]; 
          Alg[p] = u.first;  
        }
        p++; 
      }
    }
    
    // Update the Performance matrix
    for (int o(0); o < Alg.size(); o++)
    {
      TString alg = Alg[o]; 
      float er = error_V[o]; 
      if (er == 100){continue;} 

      std::vector<int> vec = Perform_Matrix[alg]; 
      vec[o]++; 
      Perform_Matrix[alg] = vec; 
    }
  }

  for (ui u = Perform_Matrix.begin(); u != Perform_Matrix.end(); u++)
  {
    std::cout << u -> first << "    "; 
    for (int t : u -> second)
    {
      std::cout << t << "   "; 
    }
    
    std::cout << std::endl;
  }
}


void Evaluate_nTrackFits(TString Filename)
{
  auto ReturnCurrentDirs =[] ()
  {
    std::vector<TString> Output; 
    TDirectory* dir = gDirectory; 
    for (TObject* key : *dir -> GetListOfKeys())
    {
      auto k = dynamic_cast<TKey*>(key);
      TString dir = (TString)k -> GetName(); 
      Output.push_back(dir); 
    }
    return Output; 
  };
  
  std::map<TString, std::vector<std::pair<TString, std::vector<TH1F*>>>> Map; 
  std::map<TString, std::vector<TH1F*>> Truth_Map; 

  TFile* F = new TFile(Filename); 
  std::vector<TString> Dir_Roots = ReturnCurrentDirs(); 
  for (TString H : Dir_Roots)
  {
    F -> cd(H); 
    std::vector<TString> Algos = ReturnCurrentDirs(); 
    std::vector<std::pair<TString, std::vector<TH1F*>>> Temp2; 
    for (TString H2 : Algos)
    {
      F -> cd(H + "/" + H2); 

      std::vector<TString> Dir_Next_Next = ReturnCurrentDirs(); 
      for (TString H3 : Dir_Next_Next)
      {
        F -> cd(H + "/" + H2 + "/" + H3); 
       

        std::vector<TString> Hists = ReturnCurrentDirs(); 
        std::vector<TH1F*> Temp; 
        for (TString Hist : Hists)
        {
          TH1F* Hp = (TH1F*)gDirectory -> Get(Hist); 
          Temp.push_back(Hp);           
        }

        if (H3.Contains("Truth"))
        {
          Truth_Map[H] = Temp; 
          continue;
        }

        Temp2.push_back(std::pair<TString, std::vector<TH1F*>>(H2, Temp)); 
      }
    }
    Map[H] = Temp2;  
    
    F -> cd(); 
  }

  std::vector<TString> Layer_Energy; 
  std::vector<std::map<TString, std::vector<float>>> Stat_Results; 
  std::vector<std::map<TString, std::vector<float>>> KS_Results;
  std::vector<TString> Algos_S; 
  for (it i = Map.begin(); i != Map.end(); i++)
  {
    TString layer_energy = i -> first; 

    std::vector<TH1F*> Truth = Truth_Map[layer_energy]; 
    std::vector<std::pair<TString, std::vector<TH1F*>>> Algos = i -> second; 
    
    std::map<TString, std::vector<float>> Stats_Map; 
    std::map<TString, std::vector<float>> KS; 
    
    //TCanvas* can = new TCanvas(); 
    //can -> Print("./Temp/" + layer_energy + ".pdf["); 
    for (std::pair<TString, std::vector<TH1F*>> Alg_Res : Algos)
    {
      std::vector<TH1F*> Algs_V = Alg_Res.second; 
      for (int x(0); x < Algs_V.size(); x++)
      {
        TH1F* H_R = Algs_V[x]; 
        TH1F* H_T = Truth[x]; 

        Stats_Map[Alg_Res.first].push_back(ErrorByIntegral(H_R, H_T) * 100); 
        KS[Alg_Res.first].push_back(H_T -> KolmogorovTest(H_R));
        
        //RatioPlot(H_R, H_T, can); 
        //if (layer_energy.Contains("Blayer_600_800_GeV") || layer_energy.Contains("IBL_2200_2400_GeV")){continue;}
        //can -> Print("./Temp/" + layer_energy + ".pdf"); 
        
        if (i == Map.begin() && x == 0){Algos_S.push_back(Alg_Res.first);}
      }
    }
    //can -> Print("./Temp/" + layer_energy + ".pdf]"); 
    //delete can;
    
    Stat_Results.push_back(Stats_Map); 
    KS_Results.push_back(KS); 
    Layer_Energy.push_back(layer_energy); 
  }
  std::map<TString, std::vector<int>> BestFit; 

  std::vector<TString> couting; 
  std::vector<int> col; 
  int Margin = 25; 
  int Col_Margin = 25; 
  TString sym = " | "; 
  TString out = ""; 
  CoutText(&out, Margin -1, " "); 
  out += (sym); 
  for (int i(0); i < Algos_S.size(); i++)
  {
    col.push_back(out.Sizeof());  
    TString Al = Algos_S[i]; 
    TString out2 = ""; 
    out2 += Al;  
    CoutText(&out2, Col_Margin - Al.Sizeof(), " "); 
    out += out2; 
    out += (sym); 
    BestFit[Al] = {0, 0, 0, 0};
  }
  couting.push_back(out); 
  
  for (int i(0); i < Stat_Results.size(); i++)
  {
    TString LJE = Layer_Energy[i]; 
    std::map<TString, std::vector<float>> Stats = Stat_Results[i]; 
    std::map<TString, std::vector<float>> KS = KS_Results[i]; 
    
    out = "";
    CoutText(&out, Margin - LJE.Sizeof(), " ");
    out += (LJE);
    out += (sym);
    couting.push_back(out);  
    out = ""; 
    std::vector<TString> Temp; 
    std::vector<float> err; 
    std::vector<TString> Algo_err;  
    for (int p(0); p < Algos_S.size(); p++)
    {
      TString Alg = Algos_S[p]; 
      std::vector<float> IE = Stats[Alg]; 
      std::vector<float> ks = KS[Alg]; 
     
      if(p == 0)
      {
        for (int k(0); k < ks.size(); k++)
        {
          TString dk = " -> dEdx_ntrk_"; dk+=(k+1); dk += ("_ntru_"); dk += (k+1); 
          TString dk_s = ""; 
          CoutText(&dk_s, Margin - dk.Sizeof(), " "); 
          dk_s += dk; 
          dk_s += (sym); 
          Temp.push_back(dk_s); 
          err.push_back(100); 
          Algo_err.push_back("NAN"); 
        }
      }

      for (int trk(0); trk < ks.size(); trk++)
      {
        int l = col[p];
        TString T = ""; 
        TString g = PrecisionString(IE[trk], 4); 
        CoutText(&T, Margin - g.Sizeof(), " "); 
        T += (g); 
        
        Temp[trk] += T; 
        Temp[trk] += (sym);

        if(err[trk] > IE[trk])
        {
          Algo_err[trk] = Alg; 
          err[trk] = IE[trk]; 
        }
      }
    }

    for (int e(0); e < err.size(); e++)
    {
      BestFit[Algo_err[e]][e]++; 
    }

    for (TString k : Temp)
    {
      couting.push_back(k); 
    }
    
    couting.push_back(out); 

  }

  couting.push_back(""); 
  couting.push_back(""); 
  couting.push_back(""); 
  couting.push_back(""); 

  std::vector<TString> Temp; 
  for (ui ki = BestFit.begin(); ki != BestFit.end(); ki++)
  {
    TString alg = ""; 
    TString g = ki -> first; 
    CoutText(&alg, Margin -g.Sizeof(), " "); 
    alg +=(g); 
    std::vector<int> res = ki -> second; 
    for (int f(0); f < res.size(); f++)
    {
      alg += (sym); alg+=(res[f]);  
    }
    couting.push_back(alg); 
  
  }
 
  std::ofstream myfile; 
  myfile.open("results.txt");  
  for (int x(0); x < couting.size(); x++)
  {
    std::cout << couting[x] << std::endl;
    
    myfile << couting[x] << "\n"; 
  }
  myfile.close();




}
