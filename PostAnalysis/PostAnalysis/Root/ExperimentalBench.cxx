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
  
  std::vector<float> k1 = {0.5, 9.5}; 
  std::vector<std::vector<float>> Ranges = {k1}; 

  float m = 0.45; 
  float s_e = 0.01; 
  // Normalization parameters
  std::map<TString, std::vector<float>> Params_N; 
  Params_N["Minimizer"] = {10000};

  // Normalization Shift parameters
  std::map<TString, std::vector<float>> Params_NS; 
  Params_NS["Range_ntrk_1"] = Ranges[0]; 
  Params_NS["dx_s"] = {-m, -m, -m, -m}; 
  //Params_NS["dx_G"] = {0, 0, 0, 0}; 
  Params_NS["dx_e"] = {m, m, m, m}; 
  //Params_NS["Seek"] = {1};
  Params_NS["Minimizer"] = {10000};
  //Params_NS["GSL"] = {1};

  // Normalization Shift FFT parameters
  std::map<TString, std::vector<float>> Params_FFT; 
  Params_FFT["Range_ntrk_1"] = Ranges[0]; 
  Params_FFT["m_s"] = {-m, -m, -m, -m};
  //Params_FFT["m_G"] = {0, 0, 0, 0};
  Params_FFT["m_e"] = {m, m, m, m};
  Params_FFT["s_C"] = {1, 1, 1, 1};
  Params_FFT["fft_cache"] = {10000}; 
  Params_FFT["Minimizer"] = {10000}; 

  // Normalization Shift Width FFT parameters
  std::map<TString, std::vector<float>> Params_WidthFFT; 
  Params_WidthFFT["Range_ntrk_1"] = Ranges[0]; 
  Params_WidthFFT["m_s"] = {-m, -m, -m, -m};
  Params_WidthFFT["m_G"] = {0, 0, 0, 0};
  Params_WidthFFT["m_e"] = {m, m, m, m};
  Params_WidthFFT["s_s"] = {0.001, 0.001, 0.001, 0.001};
  Params_WidthFFT["s_e"] = {s_e, s_e, s_e, s_e};
  Params_WidthFFT["fft_cache"] = {10000}; 
  Params_WidthFFT["Minimizer"] = {10000};

  // Simultaneous Fitting method 
  std::map<TString, std::vector<float>> Params_Sim; 
  Params_Sim["Range_ntrk_1"] = Ranges[0]; 
  Params_Sim["m_e"] = {m, m, m, m};
  Params_Sim["m_G"] = {0, 0, 0, 0};
  Params_Sim["m_s"] = {-m, -m, -m, -m};
  Params_Sim["s_s"] = {0.0005, 0.0005, 0.0005, 0.0005};
  Params_Sim["s_e"] = {s_e, s_e, s_e, s_e};
  Params_Sim["fft_cache"] = {10000}; 
  Params_Sim["Minimizer"] = {10000}; 

  // Experimental Fitting method 
  std::map<TString, std::vector<float>> Params_Exp; 
  //Params_Exp["Range_ntrk_1"] = Ranges[0];
  Params_Exp["m_e"] = {m, m, m, m};
  Params_Exp["m_G"] = {0, 0, 0, 0};
  Params_Exp["m_s"] = {-m, -m, -m, -m};
  Params_Exp["s_e"] = {s_e, s_e, s_e, s_e};
  Params_Exp["s_s"] = {0.0005, 0.0005, 0.0005, 0.0005};
  Params_Exp["Minimizer"] = {10000}; 
  Params_Exp["fft_cache"] = {10000};

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
    SubtractData(starter, trk1_start, 0, false); 

    std::vector<TH1F*> ToBeUsed; 
    for ( int i(0); i < Proposed.size(); i++)
    {
      if (Proposed[i] -> GetEntries() < 100){continue;}
      ToBeUsed.push_back(Proposed[i]); 
    }
    if (ToBeUsed.size() < 2){continue;}
    
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
    
    if (All){Mode = "Experimental"; }
    if (Mode == "Experimental")
    {
      Fits = Experimental_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_Exp, current);
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
    if (Mode == "FitT_Normal")
    {
      gDirectory -> mkdir("Normalization");
      gDirectory -> cd("Normalization");
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_N, "Normalization", current);
    }
    
    if (Mode == "FitT_ShiftNormal")
    {
      gDirectory -> mkdir("ShiftNormal");
      gDirectory -> cd("ShiftNormal");
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_NS, "NormalShift", current); 
    }

    if (Mode == "FitT_ShiftNormalFFT")
    {
      gDirectory -> mkdir("ShiftNormalFFT");
      gDirectory -> cd("ShiftNormalFFT");
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_FFT, "ShiftNormalFFT", (current + " + ShiftOnly")); 
    }

    if (Mode == "FitT_ShiftNormalWidthFFT")
    {
      gDirectory -> mkdir("ShiftNormalWidthFFT");
      gDirectory -> cd("ShiftNormalWidthFFT");
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_WidthFFT, "ConvolutionFFT", (current + "+ ShiftWidth")); 
    }
    
    if (Mode == "FitT_Incremental")
    {
      gDirectory -> mkdir("Incremental");
      gDirectory -> cd("Incremental");
      FitTemplateToTruth( TruthVector, trk1_start, ToBeUsed, Params_WidthFFT, "IncrementalFFT", current); 
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



  Fits = Experimental_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_Exp, Energy);
  Plotter(Fits, Truth, ToBeUsed, Energy, Mode); 
  

  std::cout << "Flost 3 Predicted: " << Flost3(Fits, Err)[0] << " Flost 3 Truth: " << Flost3(Truth, Err)[0] << std::endl; 
  std::cout << "Flost 2 Predicted: " << Flost2(Fits, Err)[0] << " Flost 2 Truth: " << Flost2(Truth, Err)[0] << std::endl;
    






}

