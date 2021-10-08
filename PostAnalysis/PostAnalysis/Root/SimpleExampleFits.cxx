#include<PostAnalysis/SimpleExampleFits.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/Plotting.h>

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

  //TCanvas* can = new TCanvas(); 
  //can -> SetLogy(); 
  //can -> Print(Mode + ".pdf["); 

  std::vector<std::vector<TH1F*>> ntrk_mtru; 
  std::vector<std::vector<TH1F*>> ntrk_mtru_template; 
  if (Mode.Contains("_TRUTH"))
  {
    for (int trk(0); trk < 4; trk++)
    {
      std::vector<TString> names_Test; 
      std::vector<TString> names_Templ;
      for (int tru(0); tru < 4; tru++)
      {
        TString base = "dEdx_ntrk_"; base += (trk+1); base += ("_ntru_"); base += (tru+1); 
        names_Test.push_back(base + "_Test_" + Mode); 
        names_Templ.push_back(base + "_Template_" + Mode); 
      }
      ntrk_mtru.push_back(BulkClone(Truth[trk], names_Test)); 
      ntrk_mtru_template.push_back(BulkClone(Truth[trk], names_Templ));
    }
    Mode = Mode.ReplaceAll("_TRUTH", "");
  }
  else
  {
    ntrk_mtru = BuildNtrkMtru(4, trk1_Start, "_Test_" + Mode, 4); 
    ntrk_mtru_template = BuildNtrkMtru(4, trk1_Start, "_Template_" + Mode, 4); 
  }


  std::vector<TString> Out; 
  for (int i(0); i < Data.size(); i++)
  {
    MVF Params_TFit = Params; 

    std::vector<TH1F*> ntru_P = ntrk_mtru[i];
    std::vector<TH1F*> ntru_T = Truth[i]; 
    
    // Fit to individual Truth
    //TString na = Mode + "_"; na += (i+1); na += (".pdf");
    //Can -> Print(na + "[");  
    for (int j(0); j < ntru_T.size(); j++)
    {
      TH1F* trk_tru = ntru_P[j];  
      TH1F* trk_tru_T = ntru_T[j]; 

      if (trk_tru_T -> GetEntries() < 1000){continue;}

      // Fit to the truth and get the parameters
      std::map<TString, std::vector<float>> Pred_Tru;
      
      if ( Mode == "NormalShift"){    Pred_Tru = NormalizationShift(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "ShiftNormalFFT"){ Pred_Tru = ConvolutionFFT(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "ConvolutionFFT"){ Pred_Tru = ConvolutionFFT(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "IncrementalFFT"){ Pred_Tru = IncrementalFFT(trk_tru_T, {trk_tru}, Params, "_Test"); }
      if ( Mode == "Normalization"){  Pred_Tru = Normalization(trk_tru_T, {trk_tru}, Params, "_Test"); }

      //PlotHists({trk_tru}, {trk_tru_T}, can); 
      //can -> Print(Mode + ".pdf");
      //can -> Print(na);
      //can -> Print("Debug.pdf");

    }
   
    // Perform the full fit on the data
    std::vector<TH1F*> ntru_temp = ntrk_mtru_template[i]; 
    TH1F* ntrk_D = Data[i];

    std::map<TString, std::vector<float>> Pred;   
    if ( Mode == "NormalShift"){ Pred = NormalizationShift(ntrk_D, ntru_temp, Params_TFit, "_Template"); }
    if ( Mode == "ShiftNormalFFT"){ Pred = ConvolutionFFT(ntrk_D, ntru_temp, Params_TFit, "_Template"); }
    if ( Mode == "ConvolutionFFT"){ Pred = ConvolutionFFT(ntrk_D, ntru_temp, Params_TFit, "_Template"); }
    if ( Mode == "IncrementalFFT"){ Pred = IncrementalFFT(ntrk_D, ntru_temp, Params_TFit, "_Template"); }
    if ( Mode == "Normalization"){ Pred = Normalization(ntrk_D, ntru_temp, Params_TFit, "_Template"); }
  
    // Get the values from the fits and find the delta of the fit. 
    Out.push_back(""); 
    TString Head = " ------- Performing fit under mode: "; Head += (Mode); Head+= (" :: Layer + Jet Energy "); Head += (JE); Head+= (" :: trk_"); Head += (i+1); Head += (" -------");
    Out.push_back(Head); 
    if (Pred["Normalization"].size() != 0){Compare(Pred, Params_TFit, "Normalization", "l_G", &Out);}
    if (Pred["Shift"].size() != 0){Compare(Pred, Params_TFit, "Shift", "dx_G", &Out);}
    if (Pred["Mean"].size() != 0){Compare(Pred, Params_TFit, "Mean", "m_G", &Out);}
    if (Pred["Stdev"].size() != 0){Compare(Pred, Params_TFit, "Stdev", "s_G", &Out);}

    //PlotHists(Truth[i], ntrk_mtru_template[i], can); 
    //can -> Print(Mode + ".pdf"); 
    //can -> Print(na);
    //can -> Clear();
    //can -> Print(na + "]");
  }
  
  if ( Mode == "Experimental" ){ Experimental(Data, ntrk_mtru_template, Params); }

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
  myfile.open(Mode + ".txt", std::ios_base::app); 
  for (TString S : Out)
  {
    std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 
  
  //can -> Print(Mode + ".pdf]"); 
}
