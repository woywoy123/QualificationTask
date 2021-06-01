#include<PostAnalysis/SimpleExampleFits.h>

void FitTemplateToTruth( std::vector<std::vector<TH1F*>> Truth, TH1F* trk1_Start, std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, TString Mode, TString JE)
{
  auto Compare =[&] (std::map<TString, std::vector<float>> Fit, std::map<TString, std::vector<float>> Fit_Truth, TString K1, TString K2, std::vector<TString>* out)
  {
    std::vector<float> V1 = Fit[K1]; 
    std::vector<float> V2 = Fit_Truth[K2]; 
    for (int i(0); i < V2.size(); i++)
    {
      float e1 = V1[i]; 
      float e2 = V2[i]; 
      float del = e1 - e2; 

      float state = Fit["fit_status"][0]; 

      TString x = "tru_"; x += (i+1); 
      x += (" :: Key -> "); x+= (K1); 
      x+=(" :: Truth -> "); x+=(PrecisionString(e2, 4, true)); 
      x+=(" :: Fit -> "); x+=(PrecisionString(e1, 4, true)); 
      x+= (" :: Dif (%) -> "); 
      x+= (PrecisionString((del / e2)*100, 2, false)); 
      x += (" :: Fit Status -> "); 
      x += (state); 
      (*out).push_back(x); 
    }
  }; 



  std::vector<std::vector<TH1F*>> ntrk_mtru = BuildNtrkMtru(4, trk1_Start, "_Test_" + Mode, 4); 
  std::vector<std::vector<TH1F*>> ntrk_mtru_template = BuildNtrkMtru(4, trk1_Start, "_Template_" + Mode, 4); 

  std::vector<TString> Out; 
  for (int i(0); i < Truth.size(); i++)
  {
    Params["dx_G"] = {}; 
    Params["l_G"]  = {}; 
    Params["m_G"]  = {}; 
    Params["s_G"]  = {}; 

    std::vector<TH1F*> ntru_P = ntrk_mtru[i];
    std::vector<TH1F*> ntru_T = Truth[i]; 
    
    // Fit to individual Truth
    for (int j(0); j < ntru_T.size(); j++)
    {
      
      TH1F* trk_tru = ntru_P[j]; 
      TH1F* trk_tru_T = ntru_T[j]; 
      
      // Fit to the truth and get the parameters
      std::map<TString, std::vector<float>> Pred_Tru;
      
      if ( Mode == "NormalShift"){    Pred_Tru = NormalizationShift(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "ConvolutionFFT"){ Pred_Tru = ConvolutionFFT(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "IncrementalFFT"){ Pred_Tru = IncrementalFFT(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "Normalization"){  Pred_Tru = Normalization(trk_tru_T, {trk_tru}, Params, "_Test"); }

      // Update the parameters of the PARAMS as guess 
      if (Pred_Tru["Shift"].size() != 0)        { Params["dx_G"].push_back(Pred_Tru["Shift"][0]); }
      if (Pred_Tru["Normalization"].size() != 0){ Params["l_G"].push_back(Pred_Tru["Normalization"][0]); }
      if (Pred_Tru["Mean"].size() != 0)         { Params["m_G"].push_back(Pred_Tru["Mean"][0]); }
      if (Pred_Tru["Stdev"].size() != 0)        { Params["s_G"].push_back(Pred_Tru["Stdev"][0]); }
    }
    
    // Perform the full fit on the data
    std::vector<TH1F*> ntru_temp = ntrk_mtru_template[i]; 
    TH1F* ntrk_D = Data[i]; 

    std::map<TString, std::vector<float>> Pred_Tru;   
    if ( Mode == "NormalShift"){ Pred_Tru = NormalizationShift(ntrk_D, ntru_temp, Params, "_Test"); }
    if ( Mode == "ConvolutionFFT"){ Pred_Tru = ConvolutionFFT(ntrk_D, ntru_temp, Params, "_Test"); }
    if ( Mode == "IncrementalFFT"){ Pred_Tru = IncrementalFFT(ntrk_D, ntru_temp, Params, "_Test"); }
    if ( Mode == "Normalization"){ Pred_Tru = Normalization(ntrk_D, ntru_temp, Params, "_Test"); }
  
    // Get the values from the fits and find the delta of the fit. 
    Out.push_back(""); 
    TString Head = " ------- Performing fit under mode: "; Head += (Mode); Head+= (" :: Layer + Jet Energy "); Head += (JE); Head+= (" :: trk_"); Head += (i+1); Head += (" -------");
    Out.push_back(Head); 
    if (Pred_Tru["Normalization"].size() != 0){Compare(Pred_Tru, Params, "Normalization", "l_G", &Out);}
    if (Pred_Tru["Shift"].size() != 0){Compare(Pred_Tru, Params, "Shift", "dx_G", &Out);}
    if (Pred_Tru["Mean"].size() != 0){Compare(Pred_Tru, Params, "Mean", "m_G", &Out);}
    if (Pred_Tru["Stdev"].size() != 0){Compare(Pred_Tru, Params, "Stdev", "s_G", &Out);}
  }
 
  std::vector<std::vector<float>> f; 
  float F2_P = Flost2(ntrk_mtru_template, f)[0]; 
  float F2_T = Flost2(Truth, f)[0]; 
  TString Flost2P = "====== Flost2 Performance; Fit -> "; Flost2P += (PrecisionString(F2_P, 4, false)); 
  Flost2P += (" Truth -> "); Flost2P += (PrecisionString(F2_T, 4, false)); 
  Flost2P += (" Ratio % -> "); Flost2P += (PrecisionString( (( F2_P - F2_T ) / F2_T)*100, 4, false)); 

  float F3_P = Flost3(ntrk_mtru_template, f)[0]; 
  float F3_T = Flost3(Truth, f)[0]; 
  TString Flost3P = "====== Flost3 Performance; Fit -> "; Flost3P += (PrecisionString(F3_P, 4, false)); 
  Flost3P += (" Truth -> "); Flost3P += (PrecisionString(F3_T, 4, false)); 
  Flost3P += (" Ratio % -> "); Flost3P += (PrecisionString( (( F3_P - F3_T ) / F3_T)*100, 4, false)); 

  Out.push_back(Flost2P); 
  Out.push_back(Flost3P); 

  std::ofstream myfile; 
  myfile.open("Results.txt", std::ios_base::app); 
  for (TString S : Out)
  {
    std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 



}
