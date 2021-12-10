#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/SimpleExampleFits.h>

void TestFits_AllTruth_ToTrack(TString JE, TString Mode, TString MCFile)
{
  if (MCFile == "x"){ MCFile = "Merged_MC.root"; }
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(MCFile); 
  
  std::vector<float> k1 = {0.0, 13.6}; 
  std::vector<std::vector<float>> Ranges = {k1}; 

  float m = 0.3; 
  float s_e = 0.005; 
  // Normalization parameters
  std::map<TString, std::vector<float>> Params_N; 
  Params_N["Minimizer"] = {10000};

  // Normalization Shift parameters
  int Minim = 10000; //50000; 
  std::map<TString, std::vector<float>> Params_NS; 
  Params_NS["dx_s"] = {-m, -m, -m, -m}; 
  Params_NS["dx_G"] = {0, 0, 0, 0}; 
  Params_NS["dx_e"] = {m, m, m, m}; 
  Params_NS["Seek"] = {1};
  
  // Normalization Shift FFT parameters
  std::map<TString, std::vector<float>> Params_FFT; 
  Params_FFT["m_s"] = {-m, -m, -m, -m};
  Params_FFT["m_e"] = {m, m, m, m};
  Params_FFT["s_C"] = {1, 1, 1, 1};
  Params_FFT["fft_cache"] = {Minim}; 

  // Normalization Shift Width FFT parameters
  std::map<TString, std::vector<float>> Params_WidthFFT; 
  Params_WidthFFT["m_s"] = {-m, -m, -m, -m};
  Params_WidthFFT["m_e"] = {m, m, m, m};
  Params_WidthFFT["s_s"] = {0.001, 0.001, 0.001, 0.001};
  Params_WidthFFT["s_e"] = {s_e, s_e, s_e, s_e};
  Params_WidthFFT["fft_cache"] = {Minim}; 
  
  // Simultaneous Fitting method 
  std::map<TString, std::vector<float>> Params_Sim; 
  Params_Sim["m_e"] = {m, m, m, m};
  Params_Sim["m_s"] = {-m, -m, -m, -m};
  Params_Sim["s_s"] = {0.0005, 0.0005, 0.0005, 0.0005};
  Params_Sim["s_e"] = {s_e, s_e, s_e, s_e};
  Params_Sim["fft_cache"] = {Minim}; 

  // Experimental Fitting method 
  std::map<TString, std::vector<float>> Params_Exp; 
  Params_Exp["m_e"] = {m, m, m, m};
  Params_Exp["m_G"] = {0, 0, 0, 0};
  Params_Exp["m_s"] = {-m, -m, -m, -m};
  Params_Exp["s_e"] = {s_e, s_e, s_e, s_e};
  Params_Exp["s_s"] = {0.0005, 0.0005, 0.0005, 0.0005};
  Params_Exp["Minimizer"] = {Minim}; 

  TFile* X = new TFile("Fit_Tracks.root", "RECREATE"); 
  int p = 0; 
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    auto Plotter =[&] (std::vector<std::vector<TH1F*>> Fits, std::vector<std::vector<TH1F*>> Truth, std::vector<TH1F*> Data, TString current, TString Mode)
    {
      TString name = current + "_" + Mode + ".pdf"; 
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(name + "["); 
      for (int i(0); i < Fits.size(); i++)
      {
        PlotHists(Data[i], Fits[i], Truth[i], can); 
        can -> Print(name); 
      }
      can -> Print(name + "]"); 
      for (int i(0); i < Fits.size(); i++){BulkDelete(Fits[i]);}
    };

    TString current = x -> first; 
    if (!JE.Contains("x")){if (JE != current){continue;}} // This is for parallel computing capabilities on a cluster
    std::cout << "++++++++" << current << " " << Mode << std::endl;
    
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 
    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    std::vector<TH1F*> ntrk_2_T = M["ntrk_2_T_I"]; 
    std::vector<TH1F*> ntrk_3_T = M["ntrk_3_T_I"]; 
    std::vector<TH1F*> ntrk_4_T = M["ntrk_4_T_I"]; 
    std::vector<std::vector<TH1F*>> TruthVector = { ntrk_1_T,  ntrk_2_T, ntrk_3_T, ntrk_4_T };
    
    TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
    TH1F* ntrk_2_M = M["ntrk_2_M_I"][0]; 
    TH1F* ntrk_3_M = M["ntrk_3_M_I"][0]; 
    TH1F* ntrk_4_M = M["ntrk_4_M_I"][0]; 
    std::vector<TH1F*> Proposed = {ntrk_1_M, ntrk_2_M, ntrk_3_M, ntrk_4_M};  
   
    TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
    TH1F* trk2_start = M["ntrk_2_M_O"][0];
    TH1F* trk3_start = M["ntrk_3_M_O"][0];
    TH1F* trk4_start = M["ntrk_4_M_O"][0];
    std::vector<TH1F*> starter = {trk1_start, trk2_start, trk3_start, trk4_start};  

    std::vector<TH1F*> ToBeUsed; 
    for ( int i(0); i < Proposed.size(); i++)
    {
      //if (Proposed[i] -> GetEntries() < 100){continue;}
      ToBeUsed.push_back(Proposed[i]); 
    }
  
    //if (ToBeUsed.size() < 2){continue;}
    
    X -> mkdir(current); 
    X -> cd(current); 

    std::vector<TString> ntrk_String_T = {"ntrk_1_T", "ntrk_2_T", "ntrk_3_T", "ntrk_4_T"};  
    if (Mode == "Truth")
    { 
      for (int p(0); p < TruthVector.size(); p++)
      {
        gDirectory -> mkdir(ntrk_String_T[p]); 
        gDirectory -> cd(ntrk_String_T[p]); 
        BulkWrite(TruthVector[p]); 
        gDirectory -> cd("/");
        gDirectory -> cd(current); 
      }
    }
    

    std::vector<std::vector<TH1F*>> Fits; 
    bool All = false; 
    if (Mode == "All"){All = true;}

    if (All){Mode = "Normal"; }
    if (Mode == "Normal")
    {
      Fits = Normalization_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_N, current);
      Plotter(Fits, TruthVector, ToBeUsed, current, Mode); 
    } 
    
    if (All){Mode = "ShiftNormal"; }
    if (Mode == "ShiftNormal")
    {
      Fits = NormalizationShift_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_NS, current);
      Plotter(Fits, TruthVector, ToBeUsed, current, Mode); 
    }

    if (All){Mode = "ShiftNormalFFT"; }
    if (Mode == "ShiftNormalFFT")
    {
      Fits = NormalizationShiftFFT_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_FFT, current);
      Plotter(Fits, TruthVector, ToBeUsed, current, Mode); 
    }   
    
    if (All){Mode = "ShiftNormalWidthFFT"; }
    if (Mode == "ShiftNormalWidthFFT")
    {
      Fits = NormalizationShiftWidthFFT_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_WidthFFT, current);
      Plotter(Fits, TruthVector, ToBeUsed, current, Mode); 
    } 
    
    if (All){Mode = "Incremental"; }
    if (Mode == "Incremental")
    {
      Fits = IncrementalFit(ToBeUsed, trk1_start, Params_Sim, current); 
      Plotter(Fits, TruthVector, ToBeUsed, current, Mode); 
    }

    if (All){Mode = "Simultaneous"; }
    if (Mode == "Simultaneous")
    {
      Fits = Simultaneous_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_Sim, current); 
      Plotter(Fits, TruthVector, ToBeUsed, current, Mode);
    }
    

    // ====== Simple tester Fits ====== //
    // String settings 
    if (Mode.Contains("Smooth")){Smooth(trk1_start, 0.1); }
    if (Mode.Contains("FitTo"))
    {
      Params_N["Minimizer"] = {}; 
      Params_NS["Minimizer"] = {}; 
      Params_FFT["Minimizer"] = {}; 
      Params_WidthFFT["Minimizer"] = {}; 
    }
    if (Mode.Contains("Minimizer"))
    {
      Params_N["Minimizer"] = {Minim}; 
      Params_NS["Minimizer"] = {Minim}; 
      Params_FFT["Minimizer"] = {Minim}; 
      Params_WidthFFT["Minimizer"] = {Minim}; 
    }
    if (Mode.Contains("Range"))
    {
      Params_N["Range_ntrk_1"] = Ranges[0];
      Params_NS["Range_ntrk_1"] = Ranges[0];
      Params_FFT["Range_ntrk_1"] = Ranges[0];
      Params_WidthFFT["Range_ntrk_1"] = Ranges[0];
    }

    if (Mode.Contains("Subtract")){ SubtractData(starter, trk1_start, 0, false); }
    
    TString Perfect = ""; 
    if (Mode.Contains("TRUTH")){ Perfect = "_TRUTH"; }

    // ------- Fit calls
    if (Mode.Contains("FitT_Normal_"))
    {
      gDirectory -> mkdir("Normalization"+Perfect);
      gDirectory -> cd("Normalization"+Perfect);
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_N, "Normalization"+Perfect, current);
    }

    else if (Mode.Contains("FitT_ShiftNormalFFT_"))
    {
      gDirectory -> mkdir("ShiftNormalFFT"+Perfect);
      gDirectory -> cd("ShiftNormalFFT"+Perfect);
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_FFT, "ShiftNormalFFT"+Perfect, (current + " + ShiftOnly")); 
    }

    else if (Mode.Contains("FitT_ShiftNormal_"))
    {
      gDirectory -> mkdir("ShiftNormal"+Perfect);
      gDirectory -> cd("ShiftNormal"+Perfect);
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_NS, "NormalShift"+Perfect, current); 
    }

    else if (Mode.Contains("FitT_ShiftNormalWidthFFT_"))
    {
      gDirectory -> mkdir("ShiftNormalWidthFFT"+Perfect);
      gDirectory -> cd("ShiftNormalWidthFFT"+Perfect);
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_WidthFFT, "ConvolutionFFT"+Perfect, (current + "+ ShiftWidth")); 
    }
    
    else if (Mode.Contains("FitT_Incremental_"))
    {
      gDirectory -> mkdir("Incremental"+Perfect);
      gDirectory -> cd("Incremental"+Perfect);
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_WidthFFT, "IncrementalFFT"+Perfect, current); 
    }
    
    else if (Mode.Contains("FitT_Experimental_"))
    {
      gDirectory -> mkdir("Experimental" + Perfect); 
      gDirectory -> cd("Experimental" + Perfect); 
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_WidthFFT, "Experimental" + Perfect, current); 
    }

    p++;
    X -> Write();
  }

  X -> Close(); 
}


void CompareToTruth(TString Mode, TString Energy)
{

  auto Plotter =[&] (std::vector<std::vector<TH1F*>> Fits, std::vector<std::vector<TH1F*>> Truth, std::vector<TH1F*> Data, TString current, TString Mode)
  {
    TString name = current + "_" + Mode + ".pdf"; 
    TCanvas* can = new TCanvas(); 
    can -> SetLogy(); 
    can -> Print(name + "["); 
    for (int i(0); i < Fits.size(); i++)
    {
      PlotHists(Data[i], Fits[i], Truth[i], can); 
      can -> Print(name); 
    }
    can -> Print(name + "]"); 
  };
  std::vector<std::vector<TH1F*>> Fits; 
 
  float m = 0.4; 
  float s_e = 0.01; 
  std::vector<float> k = {0.5, 9.5}; 
  std::vector<std::vector<float>> Ranges = {k, k, k, k}; 


  // Experimental Fitting method 
  std::map<TString, std::vector<float>> Params_Exp; 
  Params_Exp["Range_ntrk_1"] = Ranges[0];
  Params_Exp["m_e"] = {m, m, m, m};
  Params_Exp["m_G"] = {0, 0, 0, 0};
  Params_Exp["m_s"] = {-m, -m, -m, -m};
  Params_Exp["s_e"] = {s_e, s_e, s_e, s_e};
  Params_Exp["s_s"] = {0.0005, 0.0005, 0.0005, 0.0005};
  //Params_Exp["s_G"] = {0.05, 0.05, 0.05, 0.05};
  Params_Exp["Minimizer"] = {10000}; 
  Params_Exp["fft_cache"] = {10000};
  //Params_Exp["Seek"] = {1};
  Params_Exp["Print"] = {1}; 


  TString MCFile = "Merged_MC.root";
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(MCFile); 
  std::map<TString, std::vector<TH1F*>> M = F[Energy]; 
  
  std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
  std::vector<TH1F*> ntrk_2_T = M["ntrk_2_T_I"]; 
  std::vector<TH1F*> ntrk_3_T = M["ntrk_3_T_I"]; 
  std::vector<TH1F*> ntrk_4_T = M["ntrk_4_T_I"]; 
  std::vector<std::vector<TH1F*>> Truth = { ntrk_1_T,  ntrk_2_T, ntrk_3_T, ntrk_4_T };

  TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
  TH1F* ntrk_2_M = M["ntrk_2_M_I"][0]; 
  TH1F* ntrk_3_M = M["ntrk_3_M_I"][0]; 
  TH1F* ntrk_4_M = M["ntrk_4_M_I"][0]; 
  std::vector<TH1F*> Proposed = {ntrk_1_M, ntrk_2_M, ntrk_3_M, ntrk_4_M};  
  
  TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
  TH1F* trk2_start = M["ntrk_2_M_O"][0];
  TH1F* trk3_start = M["ntrk_3_M_O"][0];
  TH1F* trk4_start = M["ntrk_4_M_O"][0];
  std::vector<TH1F*> starter = {trk1_start, trk2_start, trk3_start, trk4_start};  
  //SubtractData(starter, trk1_start, 0, false); 

  std::vector<TH1F*> ToBeUsed; 
  for ( int i(0); i < Proposed.size(); i++)
  {
    if (Proposed[i] -> GetEntries() < 200){continue;}
    ToBeUsed.push_back(Proposed[i]); 
  }

  //auto Reb =[&] (std::vector<TH1F*> D)
  //{
  //  for (int i(0); i < D.size(); i++){ D[i] -> Rebin(2); }
  //};
  //
  //Reb( Proposed ); 
  //Reb( starter ); 
  //for (int x(0); x < Truth.size(); x++){ Reb( Truth[x] ); }

  std::vector<std::vector<float>> Err;  


  Plotter(Fits, Truth, ToBeUsed, Energy, Mode); 
  

  std::cout << "Flost 3 Predicted: " << Flost3(Fits, Err)[0] << " Flost 3 Truth: " << Flost3(Truth, Err)[0] << std::endl; 
  std::cout << "Flost 2 Predicted: " << Flost2(Fits, Err)[0] << " Flost 2 Truth: " << Flost2(Truth, Err)[0] << std::endl;
}

void FastFits(TString JE, TString Mode, TString MCFile)
{
  if (MCFile == "x"){ MCFile = "Merged_MC.root"; }
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(MCFile); 
  
  std::vector<float> k1 = {0.0, 13.6}; 
  std::vector<std::vector<float>> Ranges = {k1}; 

  float m = 0.4; 
  float s_e = 0.005; 
  // Normalization parameters
  std::map<TString, std::vector<float>> Params_N; 
  Params_N["Minimizer"] = {10000};

  // Normalization Shift parameters
  int Minim = 10000; //50000; 
  std::map<TString, std::vector<float>> Params_NS; 
  Params_NS["dx_s"] = {-m, -m, -m, -m}; 
  Params_NS["dx_G"] = {0, 0, 0, 0}; 
  Params_NS["dx_e"] = {m, m, m, m}; 
  Params_NS["Seek"] = {1};
  
  // Normalization Shift FFT parameters
  std::map<TString, std::vector<float>> Params_FFT; 
  Params_FFT["m_s"] = {-m, -m, -m, -m};
  Params_FFT["m_e"] = {m, m, m, m};
  Params_FFT["s_C"] = {1, 1, 1, 1};
  Params_FFT["fft_cache"] = {Minim}; 

  // Normalization Shift Width FFT parameters
  std::map<TString, std::vector<float>> Params_WidthFFT; 
  Params_WidthFFT["m_s"] = {-m, -m, -m, -m};
  Params_WidthFFT["m_e"] = {m, m, m, m};
  Params_WidthFFT["s_s"] = {0.001, 0.001, 0.001, 0.001};
  Params_WidthFFT["s_e"] = {s_e, s_e, s_e, s_e};
  Params_WidthFFT["fft_cache"] = {Minim}; 
  

  TFile* X = new TFile("Fit_Tracks.root", "RECREATE"); 
  int p = 0; 
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    
    TString current = x -> first; 
    if (!JE.Contains("x")){if (JE != current){continue;}} // This is for parallel computing capabilities on a cluster
    std::cout << "++++++++" << current << " " << Mode << std::endl;
    
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 
    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    std::vector<TH1F*> ntrk_2_T = M["ntrk_2_T_I"]; 
    std::vector<TH1F*> ntrk_3_T = M["ntrk_3_T_I"]; 
    std::vector<TH1F*> ntrk_4_T = M["ntrk_4_T_I"]; 
    std::vector<std::vector<TH1F*>> TruthVector = { ntrk_1_T,  ntrk_2_T, ntrk_3_T, ntrk_4_T };
    
    TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
    TH1F* ntrk_2_M = M["ntrk_2_M_I"][0]; 
    TH1F* ntrk_3_M = M["ntrk_3_M_I"][0]; 
    TH1F* ntrk_4_M = M["ntrk_4_M_I"][0]; 
    std::vector<TH1F*> DataVector = {ntrk_1_M, ntrk_2_M, ntrk_3_M, ntrk_4_M};  
   
    TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
    TH1F* trk2_start = M["ntrk_2_M_O"][0];
    TH1F* trk3_start = M["ntrk_3_M_O"][0];
    TH1F* trk4_start = M["ntrk_4_M_O"][0];
    std::vector<TH1F*> starter = {trk1_start, trk2_start, trk3_start, trk4_start};  

    X -> mkdir(current); 
    X -> cd(current); 

    std::vector<TString> ntrk_String_T = {"ntrk_1_T", "ntrk_2_T", "ntrk_3_T", "ntrk_4_T"};  
    if (Mode == "Truth")
    { 
      for (int p(0); p < TruthVector.size(); p++)
      {
        gDirectory -> mkdir(ntrk_String_T[p]); 
        gDirectory -> cd(ntrk_String_T[p]); 
        BulkWrite(TruthVector[p]); 
        gDirectory -> cd("/");
        gDirectory -> cd(current); 
        X -> Write();
      }
      continue;
    }
    
    TString Settings = "";
    int Tracks = -1; 
    if (Mode.Contains("_ntrk1")){Tracks = 0;}
    if (Mode.Contains("_ntrk2")){Tracks = 1;}
    if (Mode.Contains("_ntrk3")){Tracks = 2;}
    if (Mode.Contains("_ntrk4")){Tracks = 3;}
    if (Mode.Contains("_ntrkmtru")){Tracks = -1;}

    // Get the correct Algo 
    TString alg = ""; 
    if (Mode.Contains("Experimental")){ alg = "Experimental"; }
    else if (Mode.Contains("ShiftNormalWidthFFT")){ alg = "ShiftNormalWidthFFT"; }
    else if (Mode.Contains("ShiftNormalFFT")){ alg = "ShiftNormalFFT"; }
    else if (Mode.Contains("ShiftNormal")){ alg = "ShiftNormal"; }
    else if (Mode.Contains("Incremental")){ alg = "Incremental"; }
    else if (Mode.Contains("Normal")){ alg = "Normalization"; }
    Settings += (alg + "_"); 

    // Minimizer and Range
    if (Mode.Contains("FitTo"))
    {
      Params_N["Minimizer"] = {}; 
      Params_NS["Minimizer"] = {}; 
      Params_FFT["Minimizer"] = {}; 
      Params_WidthFFT["Minimizer"] = {}; 
      Settings += ("FitTo"); 
    }
    if (Mode.Contains("Minimizer"))
    {
      Params_N["Minimizer"] = {Minim}; 
      Params_NS["Minimizer"] = {Minim}; 
      Params_FFT["Minimizer"] = {Minim}; 
      Params_WidthFFT["Minimizer"] = {Minim}; 
      Settings += ("Minimizer"); 
    }

    if (Mode.Contains("Range"))
    {
      Params_N["Range_ntrk_1"] = Ranges[0];
      Params_NS["Range_ntrk_1"] = Ranges[0];
      Params_FFT["Range_ntrk_1"] = Ranges[0];
      Params_WidthFFT["Range_ntrk_1"] = Ranges[0];
      Settings += ("_Range");
    }
    
    // Create the templates 
    if (Mode.Contains("Subtract")){ SubtractData(starter, trk1_start, 0, false); Settings += ("_Subtract"); }
    

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
          names_Test.push_back(base + "_Test_" + Settings + "_TRUTH"); 
          names_Templ.push_back(base + "_Template_" + Settings + "_TRUTH"); 
        }
        ntrk_mtru.push_back(BulkClone(TruthVector[trk], names_Test)); 
        ntrk_mtru_template.push_back(BulkClone(TruthVector[trk], names_Templ));
      }
    }
    else
    {
      ntrk_mtru = BuildNtrkMtru(4, trk1_start, "_Test_" + Settings, 4); 
      ntrk_mtru_template = BuildNtrkMtru(4, trk1_start, "_Template_" + Settings, 4); 
    }

    if (Mode.Contains("Smooth"))
    {
      for (int trk(0); trk < ntrk_mtru_template.size(); trk++)
      { 
        Smooth(ntrk_mtru_template[trk][0], 0.1);  
        Smooth(ntrk_mtru[trk][0], 0.1);  
      }
      Settings += ("_Smooth");
    }
    
    std::cout << "here" << std::endl;
    if (Tracks > -1)
    {
      for (int trk(0); trk < ntrk_mtru_template.size(); trk++)
      {
        if (trk != Tracks)
        { 
          BulkDelete(ntrk_mtru_template[trk]); 
          BulkDelete(ntrk_mtru[trk]); 
        }
      }
      
      for (int trk(0); trk < ntrk_mtru[Tracks].size(); trk++)
      {
        if ( alg == "ShiftNormal"){    NormalizationShift(TruthVector[Tracks][trk], {ntrk_mtru[Tracks][trk]}, Params_NS, "_Test"); }
        if ( alg == "ShiftNormalFFT"){ ConvolutionFFT(TruthVector[Tracks][trk], {ntrk_mtru[Tracks][trk]}, Params_FFT, "_Test"); }
        if ( alg == "ShiftNormalWidthFFT"){ ConvolutionFFT(TruthVector[Tracks][trk], {ntrk_mtru[Tracks][trk]}, Params_WidthFFT, "_Test"); }
        if ( alg == "Incremental"){ IncrementalFFT(TruthVector[Tracks][trk], {ntrk_mtru[Tracks][trk]}, Params_WidthFFT, "_Test"); }
        if ( alg == "Normalization"){  Normalization(TruthVector[Tracks][trk], {ntrk_mtru[Tracks][trk]}, Params_N, "_Test"); }
      } 
      if ( alg == "ShiftNormal"){    NormalizationShift(DataVector[Tracks], ntrk_mtru_template[Tracks], Params_NS, "_Template"); }
      if ( alg == "ShiftNormalFFT"){ ConvolutionFFT(DataVector[Tracks], ntrk_mtru_template[Tracks], Params_FFT, "_Template"); }
      if ( alg == "ShiftNormalWidthFFT"){ ConvolutionFFT(DataVector[Tracks], ntrk_mtru_template[Tracks], Params_WidthFFT, "_Template"); }
      if ( alg == "Incremental"){ IncrementalFFT(DataVector[Tracks], ntrk_mtru_template[Tracks], Params_WidthFFT, "_Template"); }
      if ( alg == "Normalization"){  Normalization(DataVector[Tracks], ntrk_mtru_template[Tracks], Params_N, "_Template"); }
      WriteHistsToFile(ntrk_mtru_template[Tracks], current + "/" + alg); 
      WriteHistsToFile(ntrk_mtru[Tracks], current + "/" + alg); 
    }
    else
    {
      if (alg == "Experimental")
      {
        Experimental(DataVector, ntrk_mtru_template, Params_WidthFFT);
        
        for (int c(0); c < ntrk_mtru_template.size(); c++)
        {
          WriteHistsToFile(ntrk_mtru_template[c], current + "/" + alg); 
          WriteHistsToFile(ntrk_mtru[c], current + "/" + alg); 
        }
      }
    }


    p++;
  }

  X -> Close(); 
}


