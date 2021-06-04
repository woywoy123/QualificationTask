#include<PostAnalysis/Evaluate.h>


void CompareToTruth(TString dir)
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Fits = ReadAlgorithmResults(dir); 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Truth = ReadCTIDE("Merged_MC.root");
  std::vector<std::map<TString, std::map<TString, std::vector<float>>>> Collector; 
  std::vector<TString> JetEnergy; 
  std::vector<TString> Algo_Strings; 

  int index = 0; 
  for (MMVi r = Fits.begin(); r != Fits.end(); r++)
  {
    TString current = r -> first; 

    std::map<TString, std::vector<TH1F*>> Post_Fits_M = r -> second; 
    std::map<TString, std::vector<TH1F*>> Truth_M = Truth[current]; 
  
    for (MVi z = Post_Fits_M.begin(); z != Post_Fits_M.end(); z++)
    {
      TString s = z -> first; 
      int split = s.Index("_"); 
      s.Resize(split); 
     
      bool skip = false; 
      for (int i(0); i < Algo_Strings.size(); i++){ if ( s == Algo_Strings[i] ) { skip = true; }}
      if (skip) {continue;} 

      Algo_Strings.push_back(s); 
    }
    
    // Truth inside the jet core
    std::vector<TH1F*> ntrk_1_T = Truth_M["ntrk_1_T_I"]; 
    std::vector<TH1F*> ntrk_2_T = Truth_M["ntrk_2_T_I"]; 
    std::vector<TH1F*> ntrk_3_T = Truth_M["ntrk_3_T_I"]; 
    std::vector<TH1F*> ntrk_4_T = Truth_M["ntrk_4_T_I"]; 
    std::vector<std::vector<TH1F*>> All_T = {ntrk_1_T, ntrk_2_T, ntrk_3_T, ntrk_4_T};
  
    // Get the fit resutls 
    std::map<TString, std::map<TString, std::vector<float>>> Fit_Stats; 
    for (int i(0); i < Algo_Strings.size(); i++)
    {
      TString Algo = Algo_Strings[i]; 
      std::vector<TH1F*> ntrk_1_R = Post_Fits_M[Algo + "_ntrk_1"]; 
      std::vector<TH1F*> ntrk_2_R = Post_Fits_M[Algo + "_ntrk_2"]; 
      std::vector<TH1F*> ntrk_3_R = Post_Fits_M[Algo + "_ntrk_3"]; 
      std::vector<TH1F*> ntrk_4_R = Post_Fits_M[Algo + "_ntrk_4"]; 
      std::vector<std::vector<TH1F*>> All_R = {ntrk_1_R, ntrk_2_R, ntrk_3_R, ntrk_4_R};
     
      //TCanvas* can = new TCanvas(); 
      //TString Fname = current + "_" + Algo + ".pdf"; 
      //can -> Print(Fname + "["); 
      
      std::map<TString, std::vector<float>> Stats_i; 
      for (int trk(0); trk < All_R.size(); trk++)
      {
        //can -> SetLogy();
        //PlotHists(All_R[trk], All_T[trk], current + "_" + Algo, can); 
        //can -> Print(Fname); 
        
        for (int tru(0); tru < All_R[trk].size(); tru++)
        {
          TString name = All_R[trk][tru] -> GetTitle(); 
          
          int n_tru;  
          for (int x(0); x < All_T[trk].size(); x++)
          {
            TString dEdx_S = "dEdx_ntrk_"; dEdx_S += (trk+1); dEdx_S += ("_ntru_"); dEdx_S += (x+1);  
            if (name.Contains(dEdx_S)){ n_tru = x; }
          }

          if (trk == n_tru)
          {
            //RatioPlot(All_R[trk][tru], All_T[trk][n_tru], can); 
            //can -> Print(Fname); 
            
            // Do the stats: 
            TString dEdx_S = "dEdx_ntrk_"; dEdx_S += (trk+1); dEdx_S += ("_ntru_"); dEdx_S += (n_tru+1);  
            Stats_i[dEdx_S].push_back(All_T[trk][n_tru] -> Chi2Test( All_R[trk][tru], "WW" )); 
            Stats_i[dEdx_S].push_back(All_T[trk][n_tru] -> Chi2Test( All_R[trk][tru], "WW CHI2/NDF" ));
           
            Stats_i[dEdx_S].push_back(All_T[trk][n_tru] -> KolmogorovTest( All_R[trk][tru], "" )); 
            Stats_i[dEdx_S].push_back(All_T[trk][n_tru] -> KolmogorovTest( All_R[trk][tru], "M" ));
            
            //Stats_i[dEdx_S].push_back(All_T[trk][n_tru] -> AndersonDarlingTest( All_R[trk][tru] ));
            Stats_i[dEdx_S].push_back(ErrorByIntegral(All_R[trk][tru], All_T[trk][n_tru])); 
          }
        }
      }
      Fit_Stats[Algo] = Stats_i; 
      //can -> Print(Fname + "]"); 
      //delete can; 
    }
    JetEnergy.push_back(current); 
    Collector.push_back(Fit_Stats); 
  
  }

  CompileCout(JetEnergy, Algo_Strings, Collector); 
}

void MultiTrackTruthComparison(TString dir)
{

  auto Plotting =[&] (std::vector<std::vector<TH1F*>> Truth, std::vector<std::vector<TH1F*>> Prediction, TString name, bool skip)
  {
    if (!skip){return;}
    TCanvas* can = new TCanvas(); 
    can -> SetLogy(); 
    can -> Print(name + ".pdf["); 
    for (int i(0); i < Prediction.size(); i++)
    {
      std::vector<TH1F*> ntrk_T = Truth[i]; 
      std::vector<TH1F*> ntrk_P = Prediction[i]; 
    
      PlotHists(ntrk_P, ntrk_T, name, can); 
      can -> Print(name + ".pdf"); 
      for (int p(0); p < Prediction[i].size(); p++)
      {
        RatioPlot(ntrk_T[p], ntrk_P[p], can);
        can -> Print(name + ".pdf"); 
      }
    }
    can -> Print(name + ".pdf]"); 
    delete can;
  }; 

  auto Statistics =[&] (std::vector<std::vector<TH1F*>> Truth, std::vector<std::vector<TH1F*>> Prediction)
  {
    std::map<TString, std::vector<float>> Map; 
    for (int i(0); i < Prediction.size(); i++)
    {
      TH1F* temp = SumHists(Truth[i], "temp"); 
      float L = temp -> Integral(); 
      delete temp; 

      for (int p(0); p < Prediction[i].size(); p++)
      {
        TH1F* trk_P = Prediction[i][p]; 
        TH1F* trk_T = Truth[i][p]; 
        
        TString dEdx_S = "dEdx_ntrk_"; dEdx_S += (i+1); dEdx_S += ("_ntru_"); dEdx_S += (p+1);  
        Map[dEdx_S].push_back(0); //trk_P -> Chi2Test( trk_T, "WW" )); 
        Map[dEdx_S].push_back(0); //trk_P -> Chi2Test( trk_T, "WW CHI2/NDF" ));
       
        Map[dEdx_S].push_back(0); //trk_P -> KolmogorovTest( trk_T, "" )); 
        Map[dEdx_S].push_back(0); //trk_P -> KolmogorovTest( trk_T, "M" ));
        
        float r = trk_T -> Integral() / L;
        Map[dEdx_S].push_back(r*ErrorByIntegral(trk_P, trk_T)); 
      }
    }
    return Map; 
  }; 

  auto CreateErrorBoard =[&] (std::map<TString, std::vector<float>>* Scores, std::vector<std::vector<TH1F*>> Truth, std::vector<std::vector<TH1F*>> Prediction)
  {
    for (int i(0); i < Prediction.size(); i++)
    {
      TH1F* temp = SumHists(Truth[i], "temp"); 
      float L = temp -> Integral(); 
      delete temp; 

      for (int p(0); p < Prediction[i].size(); p++)
      {
        TH1F* trk_P = Prediction[i][p]; 
        TH1F* trk_T = Truth[i][p]; 
        
        TString dEdx_S = "dEdx_ntrk_"; dEdx_S += (i+1); dEdx_S += ("_ntru_"); dEdx_S += (p+1);  
        float r = trk_T -> Integral() / L; 
       
        // Simple test cases to validate result table
        //TString name = trk_P -> GetTitle(); 
        //if (name.Contains("Experimental")){ r = 0.1; }
        //if (name.Contains("Normal")){ r = 0.01; }
        //if (name.Contains("ShiftNormal")){ r = 0.1; }
        //if (name.Contains("ShiftNormalFFT")){ r = 0.1; }
        //if (name.Contains("ShiftNormalWidthFFT")){ r = 0.1; }
        
        //if (name.Contains("ShiftNormalWidthFFT_")) { r = 0; }
        float err = r*ErrorByIntegral(trk_P, trk_T); 
        (*Scores)[dEdx_S].push_back(err);
      }
    }
  }; 

  auto GetScores =[&] (std::map<TString, std::vector<float>> Map, std::map<TString, std::vector<float>>* Scores)
  {
    
    std::vector<std::pair<TString, float>> Best; 
    for (MVFi x = Map.begin(); x != Map.end(); x++)
    {
      float err = 1000; 
      int index = -1; 
      for (int i(0); i < Map[x -> first].size(); i++)
      {
        float e = Map[x -> first][i]; 
        
        if (e < err)
        {
          err = e;
          index = i; 
        }
      }
      Best.push_back(std::pair<TString, float>(x -> first, index)); 
    }
    
    for (int i(0); i < Best.size(); i++)
    {
      if ((*Scores)[Best[i].first].size() == 0){(*Scores)[Best[i].first] = std::vector<float>(8, 0);}
      (*Scores)[Best[i].first][Best[i].second]++; 
    }
  };
  
  auto GetNormalizationError =[&] (std::map<TString, std::vector<float>> Err)
  {
    std::vector<std::vector<float>> Output; 
    for (int i(0); i < 4; i++)
    {
      for (MVFi p = Err.begin(); p != Err.end(); p++)
      {
        TString trk = "ntrk_"; trk += (i+1); 
        TString metric = p -> first; 
          
        if (!metric.Contains(trk) || !metric.Contains("Normalization_Error")){continue;}
        Output.push_back(p -> second);    
      } 
    }
    return Output; 
  };
  

  int n_energy = 99;
  bool PlotOn = false; 
  MMVT Results = ReadAlgorithmResults(dir); 
  MMMVF Errors = ReadOutputFileToMap(dir); 
  MMVF Stats;
  MVF Scores; 
  MVF Flost2_Map;
  MVF Flost3_Map;
 // VS Algo_Strings = {"Normal", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT", "Incremental", "Experimental", "Simultaneous"};
  VS Algo_Strings = {"Normal", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT", "Experimental"};
 
  int ip = 0; 
  for (MMVi x = Results.begin(); x != Results.end(); x++)
  {
    if (ip > n_energy && n_energy != -1){continue;}
    ip++;  
 
    std::map<TString, std::vector<float>> Error_Map; 
    TString LE = x -> first; 

    // Return the truth 
    VT ntrk_1_T = Results[LE]["ntrk_1_Truth"]; 
    VT ntrk_2_T = Results[LE]["ntrk_2_Truth"]; 
    VT ntrk_3_T = Results[LE]["ntrk_3_Truth"]; 
    VT ntrk_4_T = Results[LE]["ntrk_4_Truth"]; 
    VVT All_Truth = {ntrk_1_T, ntrk_2_T, ntrk_3_T, ntrk_4_T}; 
    VVF N_Err = GetNormalizationError(Errors[LE]["Truth"]); 
    
    Flost2_Map[LE + "_Truth"] = Flost2(All_Truth, N_Err);  
    Flost3_Map[LE + "_Truth"] = Flost3(All_Truth, N_Err); 
   
    // Go over the algorithms 
    for (int i(0); i < Algo_Strings.size(); i++)
    {
      TString Alg = Algo_Strings[i];  
      VT ntrk_1_Algo = Results[LE][Alg + "_ntrk_1"];
      VT ntrk_2_Algo = Results[LE][Alg + "_ntrk_2"];
      VT ntrk_3_Algo = Results[LE][Alg + "_ntrk_3"];
      VT ntrk_4_Algo = Results[LE][Alg + "_ntrk_4"];

      VVT All_Algo = {ntrk_1_Algo, ntrk_2_Algo, ntrk_3_Algo, ntrk_4_Algo}; 
      Plotting(All_Truth, All_Algo, LE + Alg, PlotOn); 
      Stats[LE + "_" + Alg] = Statistics(All_Truth, All_Algo); 
      CreateErrorBoard(&Error_Map, All_Truth, All_Algo); 
      VVF Algo_Error = GetNormalizationError(Errors[LE][Alg]); 
      
      // Debugging line 
      if (Alg != "x")
      {
        Flost2_Map[LE + "_" + Alg] = Flost2(All_Algo, Algo_Error);  
        Flost3_Map[LE + "_" + Alg] = Flost3(All_Algo, Algo_Error); 
      }
      else
      {
        Flost2_Map[LE + "_" + Alg] = Flost2(All_Truth, N_Err);  
        Flost3_Map[LE + "_" + Alg] = Flost3(All_Truth, N_Err); 
      }
      GetScores(Error_Map, &Scores); 
    }
  }

  std::vector<TString> cout_V; 
  
  // ============= COUT STUFF ===================== //
  CompileTruthTrackEvaluation(Results, &cout_V);
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 

  // ===== Status codes from individual fits ===== //
  CompileStatusCodeTable(Errors, &cout_V);
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 

  // ========= FLOST comparisons to truth ============= //
  CompileFLostTable(Results, Errors, &cout_V, "FL2"); 
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 
  CompileFLostTable(Results, Errors, &cout_V, "FL3"); 

  std::ofstream myfile; 
  myfile.open("Results.txt"); 
  for (TString S : cout_V)
  {
    std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 
}

void CompileCout(std::vector<TString> JetEnergy, std::vector<TString> Algo_Strings, std::vector<std::map<TString, std::map<TString, std::vector<float>>>> Collector)
{
  
  std::vector<TString> cout_V; 
  TString out; 
  int margin = 49; 
  TString sep = " | "; 
  
  CoutText(&out, margin/2, " "); 
  for (int i(0); i < Algo_Strings.size(); i++)
  {
    if (Algo_Strings[i].Contains("ntrk")){continue;}
    out += (sep); 
    out += Algo_Strings[i]; 
    CoutText(&out, margin - Algo_Strings[i].Sizeof()+9, " "); 
  }
  out += (sep); 
  cout_V.push_back(out); 
  out.Clear(); 

  CoutText(&out, margin/2, " "); 
  std::vector<TString> Stats_N = {"Chi2-WW", "Chi2-NDF", "KS", "KS-M", "IE"}; 
  
  int sub_M = 10;
  for (int i(0); i < Algo_Strings.size(); i++)
  {
    for (TString st : Stats_N)
    {
      out += (sep); 
      out += (st); 
      CoutText(&out, sub_M - st.Sizeof(), " "); 
    }
  }
  out += (sep);
  cout_V.push_back(out); 
  
  for (int j(0); j < JetEnergy.size(); j++)
  {
    out.Clear(); 
    CoutText(&out, margin / 2 - JetEnergy[j].Sizeof() +1, " "); 
    out += (JetEnergy[j]); 
    out += (sep); 

    cout_V.push_back(out); 
    out.Clear(); 

    std::map<TString, std::map<TString, std::vector<float>>> E = Collector[j]; 
   
    std::map<TString, TString> Cout_ntrk; 
    for (TString Alg : Algo_Strings)
    {
      std::map<TString, std::vector<float>> ntrk_ntru = E[Alg]; 
      for (MVFi tt = E[Alg].begin(); tt != E[Alg].end(); tt++)
      {
        for (float x : tt -> second)
        {
          TString y = tt -> first; 
          TString num = PrecisionString(x, 3, true); 
          out += (num); 
          CoutText(&out, sub_M - num.Sizeof(), " "); 
          out += (sep); 
          
          if (Cout_ntrk[y].Sizeof() == 1)
          { 
            TString out2; 
            CoutText(&out2, margin / 2 - y.Sizeof() -1, " "); 
            out2 += ("->");  
            out2 += (y); 
            out2 += (sep);  
            out2 += (out); 

            Cout_ntrk[y] += (out2); 
            out.Clear(); 
          }
          else 
          {
            Cout_ntrk[y] += (out); 
            out.Clear(); 
          }
        }
      }
    }

    for (std::map<TString, TString>::iterator xp = Cout_ntrk.begin(); xp != Cout_ntrk.end(); xp++){cout_V.push_back(xp -> second);}
    cout_V.push_back(" "); 
  }

  std::vector<int> ScoreOverall(Algo_Strings.size(), 0); 
  std::map<TString, std::vector<float>> ScorePerTrack; 
  for (int j(0); j < JetEnergy.size(); j++)
  {
    std::map<TString, std::map<TString, std::vector<float>>> E = Collector[j]; 
    
    std::vector<float> temp_E(4, 100); 
    std::vector<int> temp_id(4, -1);  
    std::vector<float> temp(Algo_Strings.size(), -1);  
    for (int k(0); k < Algo_Strings.size(); k++)
    {
      std::map<TString, std::vector<float>> re = E[Algo_Strings[k]]; 
      
      int id = 0; 
      for (MVFi x = re.begin(); x != re.end(); x++)
      {
        if (ScorePerTrack[x -> first].size() == 0)
        {
          ScorePerTrack[x -> first] = std::vector<float>(Algo_Strings.size(), 0);  
        }
        
        float IE = (x -> second)[4]; // Return Integral Error
        temp[k] += IE; 
        
        if (temp_E[id] > IE)
        { 
          temp_E[id] = IE; 
          temp_id[id] = k;
        }
        id++; 
      }
    }
    
    // Update the score 
    int id = 0; 
    for (MVFi sc = ScorePerTrack.begin(); sc != ScorePerTrack.end(); sc++)
    {
      int id_s = temp_id[id]; 
      if (id_s == -1){ continue; }
      ScorePerTrack[sc -> first][id_s]++; 
      id++;
    }

    int ind = -1; 
    float err_t = 1000; 
    for (int k(0); k < temp.size(); k++)
    {
      if (temp[k] < err_t)
      {
        err_t = temp[k]; 
        ind = k; 
      }
    }
    ScoreOverall[ind]++; 
  }
  
  
  cout_V.push_back(" "); 
  cout_V.push_back(" "); 

  out.Clear(); 
  CoutText(&out, 18, " "); 
  out += (sep);  

  for (TString Alg : Algo_Strings)
  {
    out += (Alg); 
    out += (sep); 
  }
  cout_V.push_back(out); 

  for (MVFi sc = ScorePerTrack.begin(); sc != ScorePerTrack.end(); sc++)
  {
    out.Clear(); 
    out += (sc -> first); 
    out += (sep); 
    for (int res(0); res < (sc -> second).size(); res++)
    {
      TString out2 = ""; 
      out2 += (sc -> second)[res]; 
      CoutText(&out2, Algo_Strings[res].Sizeof() - out2.Sizeof(), " "); 
      out += ( out2 ); 
      out += (sep); 
    }
    cout_V.push_back(out); 
  }
  out.Clear(); 
  
  CoutText(&out, 18, " "); 
  out += (sep); 
  for (int sc(0); sc < ScoreOverall.size(); sc++)
  {
    TString out2 = ""; 
    TString xout = ""; 
    xout += (ScoreOverall[sc]);  
    CoutText(&out2, Algo_Strings[sc].Sizeof()-xout.Sizeof(), " "); 
    
    out += (xout); 
    out += (out2);  
    out += (sep); 
  }
  cout_V.push_back(out); 

  std::ofstream myfile; 
  myfile.open("Results.txt"); 
  for (TString S : cout_V)
  {
    std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 

}

