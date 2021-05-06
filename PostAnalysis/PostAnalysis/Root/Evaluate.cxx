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
        Map[dEdx_S].push_back(trk_P -> Chi2Test( trk_T, "WW" )); 
        Map[dEdx_S].push_back(trk_P -> Chi2Test( trk_T, "WW CHI2/NDF" ));
       
        Map[dEdx_S].push_back(trk_P -> KolmogorovTest( trk_T, "" )); 
        Map[dEdx_S].push_back(trk_P -> KolmogorovTest( trk_T, "M" ));
        
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
        float r = 1; //trk_T -> Integral() / L; 
        
        TString name = trk_P -> GetTitle(); 
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
        
        if (e <= err)
        {
          err = e;
          index = i; 
        }
      }
      //std::cout << x -> first << std::endl;
      Best.push_back(std::pair<TString, float>(x -> first, index)); 
    }
    
    for (int i(0); i < Best.size(); i++)
    {
      if ((*Scores)[Best[i].first].size() == 0){(*Scores)[Best[i].first] = std::vector<float>(5, 0);}
      (*Scores)[Best[i].first][Best[i].second]++; 
    }
  };

  bool PlotOn = false; 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> Results = ReadAlgorithmResults(dir); 
  std::map<TString, std::map<TString, std::map<TString, std::vector<float>>>> Errors = ReadOutputFileToMap(dir); 
  std::map<TString, std::map<TString, std::vector<float>>> Stats;
  std::map<TString, std::vector<float>> Scores; 
  for (MMVi x = Results.begin(); x != Results.end(); x++)
  {
  
    std::map<TString, std::vector<float>> Error_Map; 
    TString LE = x -> first; 

    // Return the truth 
    std::vector<TH1F*> ntrk_1_T = Results[x -> first]["ntrk_1_Truth"]; 
    std::vector<TH1F*> ntrk_2_T = Results[x -> first]["ntrk_2_Truth"]; 
    std::vector<TH1F*> ntrk_3_T = Results[x -> first]["ntrk_3_Truth"]; 
    std::vector<TH1F*> ntrk_4_T = Results[x -> first]["ntrk_4_Truth"]; 
    std::vector<std::vector<TH1F*>> All_Truth = {ntrk_1_T, ntrk_2_T, ntrk_3_T, ntrk_4_T};  
    
    // Return the Algos:
    // Normal
    std::vector<TH1F*> ntrk_1_Normal = Results[x -> first]["Normal_ntrk_1"]; 
    std::vector<TH1F*> ntrk_2_Normal = Results[x -> first]["Normal_ntrk_2"]; 
    std::vector<TH1F*> ntrk_3_Normal = Results[x -> first]["Normal_ntrk_3"]; 
    std::vector<TH1F*> ntrk_4_Normal = Results[x -> first]["Normal_ntrk_4"]; 
    std::vector<std::vector<TH1F*>> All_Normal = {ntrk_1_Normal, ntrk_2_Normal, ntrk_3_Normal, ntrk_4_Normal}; 
    Plotting(All_Truth, All_Normal, LE + "_Normal", PlotOn); 
    Stats[LE + "_Normal"] = Statistics(All_Truth, All_Normal); 
    CreateErrorBoard(&Error_Map, All_Truth, All_Normal);  

    // ShiftNormal
    std::vector<TH1F*> ntrk_1_ShiftNormal = Results[x -> first]["ShiftNormal_ntrk_1"]; 
    std::vector<TH1F*> ntrk_2_ShiftNormal = Results[x -> first]["ShiftNormal_ntrk_2"]; 
    std::vector<TH1F*> ntrk_3_ShiftNormal = Results[x -> first]["ShiftNormal_ntrk_3"]; 
    std::vector<TH1F*> ntrk_4_ShiftNormal = Results[x -> first]["ShiftNormal_ntrk_4"]; 
    std::vector<std::vector<TH1F*>> All_ShiftNormal = {ntrk_1_ShiftNormal, ntrk_2_ShiftNormal, ntrk_3_ShiftNormal, ntrk_4_ShiftNormal}; 
    Plotting(All_Truth, All_ShiftNormal, LE + "_ShiftNormal", PlotOn); 
    Stats[LE + "_ShiftNormal"] = Statistics(All_Truth, All_ShiftNormal);
    CreateErrorBoard(&Error_Map, All_Truth, All_ShiftNormal);  

    // ShiftNormalFFT
    std::vector<TH1F*> ntrk_1_ShiftNormalFFT = Results[x -> first]["ShiftNormalFFT_ntrk_1"]; 
    std::vector<TH1F*> ntrk_2_ShiftNormalFFT = Results[x -> first]["ShiftNormalFFT_ntrk_2"]; 
    std::vector<TH1F*> ntrk_3_ShiftNormalFFT = Results[x -> first]["ShiftNormalFFT_ntrk_3"]; 
    std::vector<TH1F*> ntrk_4_ShiftNormalFFT = Results[x -> first]["ShiftNormalFFT_ntrk_4"]; 
    std::vector<std::vector<TH1F*>> All_ShiftNormalFFT = {ntrk_1_ShiftNormalFFT, ntrk_2_ShiftNormalFFT, ntrk_3_ShiftNormalFFT, ntrk_4_ShiftNormalFFT}; 
    Plotting(All_Truth, All_ShiftNormalFFT, LE + "_ShiftNormalFFT", PlotOn); 
    Stats[LE + "_ShiftNormalFFT"] = Statistics(All_Truth, All_ShiftNormalFFT);
    CreateErrorBoard(&Error_Map, All_Truth, All_ShiftNormalFFT);  
   
    // ShiftNormalWidthFFT
    std::vector<TH1F*> ntrk_1_ShiftNormalWidthFFT = Results[x -> first]["ShiftNormalWidthFFT_ntrk_1"]; 
    std::vector<TH1F*> ntrk_2_ShiftNormalWidthFFT = Results[x -> first]["ShiftNormalWidthFFT_ntrk_2"]; 
    std::vector<TH1F*> ntrk_3_ShiftNormalWidthFFT = Results[x -> first]["ShiftNormalWidthFFT_ntrk_3"]; 
    std::vector<TH1F*> ntrk_4_ShiftNormalWidthFFT = Results[x -> first]["ShiftNormalWidthFFT_ntrk_4"]; 
    std::vector<std::vector<TH1F*>> All_ShiftNormalWidthFFT = {ntrk_1_ShiftNormalWidthFFT, ntrk_2_ShiftNormalWidthFFT, ntrk_3_ShiftNormalWidthFFT, ntrk_4_ShiftNormalWidthFFT}; 
    Plotting(All_Truth, All_ShiftNormalWidthFFT, LE + "_ShiftNormalWidthFFT", PlotOn); 
    Stats[LE + "_ShiftNormalWidthFFT"] = Statistics(All_Truth, All_ShiftNormalWidthFFT);
    CreateErrorBoard(&Error_Map, All_Truth, All_ShiftNormalWidthFFT);  
    
    GetScores(Error_Map, &Scores); 
    
  }


  std::vector<TString> cout_V; 
  TString out; 
  int margin = 20; 
  TString sep = " | "; 
  
  CoutText(&out, margin-1, " "); 
  
  std::vector<TString> Algo_Strings = {"Normal", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT"}; 
  for (int i(0); i < Algo_Strings.size(); i++)
  {
    out += (sep); 
    out += Algo_Strings[i]; 
    CoutText(&out, margin - Algo_Strings[i].Sizeof(), " "); 
  }
  out += (sep); 
  cout_V.push_back(out); 
  out.Clear();  

  for (MVFi p = Scores.begin(); p != Scores.end(); p++)
  {
    TString name = p -> first; 
    out += (name); 
    CoutText(&out, margin - name.Sizeof(), " "); 
    out += (sep);
    
    int total = 0; 
    for (int i(0); i < Algo_Strings.size(); i++)
    {
      int x = (p -> second)[i];  
      TString sub = "   "; 
      sub += (x); 
      CoutText(&sub, margin - sub.Sizeof(), " "); 
      out += (sub); 
      out += (sep); 
      total += x; 
    }
    out += (total); 

    cout_V.push_back(out); 
    out.Clear();
  }

  std::map<TString, std::vector<float>> verdict; 
  for (MVFi p = Scores.begin(); p != Scores.end(); p++)
  {
    auto SumVect =[] (std::vector<float> *A, std::vector<float> B)
    {
      if (A -> size() == 0){*A = std::vector<float>(5, 0);}
      for (int i(0); i < A -> size(); i++)
      {
        if (i == B.size()){continue;}
        (*A)[i] += B[i]; 
      }
    }; 

    TString name = p -> first; 
    if (name.Contains("ntrk_1")) { SumVect(&verdict["dEdx_ntrk_1"], (p -> second)); }
    if (name.Contains("ntrk_2")) { SumVect(&verdict["dEdx_ntrk_2"], (p -> second)); }
    if (name.Contains("ntrk_3")) { SumVect(&verdict["dEdx_ntrk_3"], (p -> second)); }
    if (name.Contains("ntrk_4")) { SumVect(&verdict["dEdx_ntrk_4"], (p -> second)); }
  }


  CoutText(&out, 120, ":"); 
  cout_V.push_back(out); 
  out.Clear();

  for (MVFi p = verdict.begin(); p != verdict.end(); p++)
  {
    TString name = p -> first; 
    out += (name); 
    CoutText(&out, margin - name.Sizeof(), " "); 
    out += (sep);
    
    int total = 0; 
    for (int i(0); i < Algo_Strings.size(); i++)
    {
      int x = (p -> second)[i];  
      TString sub = "   "; 
      sub += (x); 
      CoutText(&sub, margin - sub.Sizeof(), " "); 
      out += (sub); 
      out += (sep); 
      total += x; 
    }
    out += (total); 

    cout_V.push_back(out); 
    out.Clear();
  }

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

