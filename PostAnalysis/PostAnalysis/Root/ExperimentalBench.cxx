#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/Plotting.h>

void IOTest(){TestReadCTIDE();}
void TestRead(){TestReadAlgorithm();}
void RooFitBaseFunctionTest()
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE("Merged_MC.root");
  
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 

    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 

    //// Test the Normalization function out
    std::vector<TString> ntrk1_Name = NameGenerator(ntrk_1_T.size(), "dEdx_ntrk_1_tru_"); 
    std::vector<TH1F*> PDF_H = BulkClone(ntrk_1_T, ntrk1_Name); 
    
    std::map<TString, std::vector<float>> Params_N; 
    //Params_N["Range"] = {0, 10}; 
    //Params_N["r_value"] = {1.1};
    
    Normalization(ntrk_1_M, PDF_H, Params_N, "Normal"); 
  
    TCanvas* can = new TCanvas(); 
    can -> SetLogy();
    PlotHists(ntrk_1_T, PDF_H, can); 
    can -> Print("Example.pdf"); 
    delete can; 

    //std::map<TString, std::vector<float>> Params_NS; 
    //Params_NS["Range"] = {1, 8}; 
    //Params_NS["r_value"] = {1.1};
    //Params_NS["dx"] = {0.5, 0.5, 0.5, 0.5}; 
    //Params_NS["dx_G"] = {0, 0, 0, 0};
    //NormalizationShift(ntrk_1_M, PDF_H, Params_NS, "NormalShift"); 
    //
    //TCanvas* can1 = new TCanvas(); 
    //can1 -> SetLogy();
    //PlotHists(ntrk_1_T, PDF_H, can1); 
    //can1 -> Print("ExampleS.pdf"); 
    //delete can1; 


    ntrk1_Name = NameGenerator(ntrk_1_T.size(), "FFTdEdx_ntrk_1_tru_"); 
    std::vector<TH1F*> PDF_H_FFT = BulkClone(ntrk_1_T, ntrk1_Name); 
    std::map<TString, std::vector<float>> Params_FFT; 
    Params_FFT["LR"] = {20};  
    Params_FFT["G_Mean"] = {0, 0, 0, 0}; 
    Params_FFT["G_Stdev"] = {0.01, 0.01, 0.01, 0.01}; 
    Params_FFT["Range"] = {0.2, 8}; 
    Params_FFT["r_value"] = {1.1};
    Params_FFT["m"] = {0.25, 0.25, 0.25, 0.25};
    Params_FFT["m_G"] = {0, 0, 0, 0}; 
    Params_FFT["s_s"] = {0.0001, 0.0001, 0.0001, 0.0001};
    Params_FFT["s_e"] = {0.05, 0.05, 0.05, 0.05};
    Params_FFT["fft_cache"] = {10000}; 
    Params_FFT["Minimizer"] = {100000}; 

    //ConvolutionFFT(ntrk_1_M, PDF_H_FFT, Params_FFT, "NormalShiftFFT"); 
    //TCanvas* can2 = new TCanvas(); 
    //can2 -> SetLogy();
    //PlotHists(ntrk_1_T, PDF_H_FFT, can2); 
    //can2 -> Print("ExampleSFFT.pdf"); 
    //delete can2; 

    //TCanvas* can3 = new TCanvas(); 
    //can3-> SetLogy();
    //for (int i(0); i < 10; i++)
    //{
    //  DeConvolutionFFT(ntrk_1_M, PDF_H_FFT, Params_FFT, "DeConvolutionFFT"); 
    //  PlotHists(ntrk_1_T, PDF_H_FFT, can3); 
    //  can3 -> Print("ExampleSFFT.pdf"); 
    //}


    break;
  }
}

void TestFits_NTruth_NTrack()
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE("Merged_MC.root"); 
  
  // Normalization parameters
  std::map<TString, std::vector<float>> Params_N; 
  Params_N["Range"] = {0, 10}; 
  Params_N["r_value"] = {1.1};
  Params_N["Range_ntrk_1"] = {0, 6}; 
  Params_N["Range_ntrk_2"] = {1, 6}; 
  Params_N["Range_ntrk_3"] = {2.5, 7}; 
  Params_N["Range_ntrk_4"] = {4, 8}; 

  // Normalization Shift parameters
  std::map<TString, std::vector<float>> Params_NS; 
  Params_NS["Range"] = {0.2, 8}; 
  Params_NS["r_value"] = {1.1};
  Params_NS["Range_ntrk_1"] = {0, 6}; 
  Params_NS["Range_ntrk_2"] = {1, 6}; 
  Params_NS["Range_ntrk_3"] = {2.5, 7}; 
  Params_NS["Range_ntrk_4"] = {4, 8}; 
  Params_NS["dx"] = {0.5, 0.5, 0.5, 0.5}; 
  Params_NS["dx_G"] = {0, 0, 0, 0};
  //Params_NS["Minimizer"] = {100000}; 
  
  // Normalization Shift FFT parameters
  std::map<TString, std::vector<float>> Params_FFT; 
  Params_FFT["Range"] = {0.2, 8}; 
  Params_FFT["r_value"] = {1.1};
  Params_FFT["Range_ntrk_1"] = {0, 6}; 
  Params_FFT["Range_ntrk_2"] = {1, 6}; 
  Params_FFT["Range_ntrk_3"] = {2.5, 7}; 
  Params_FFT["Range_ntrk_4"] = {4, 8}; 
  Params_FFT["m"] = {0.5, 0.5, 0.5, 0.5};
  Params_FFT["m_G"] = {0, 0, 0, 0}; 
  Params_FFT["s_s"] = {0, 0, 0, 0};
  Params_FFT["s_e"] = {0.1, 0.1, 0.1, 0.1};
  Params_FFT["fft_cache"] = {10000}; 
  //Params_FFT["Minimizer"] = {100000}; 

  // Normalization Shift Width FFT parameters
  std::map<TString, std::vector<float>> Params_WidthFFT; 
  Params_WidthFFT["Range"] = {0.2, 8}; 
  Params_WidthFFT["r_value"] = {1.1};
  Params_WidthFFT["Range_ntrk_1"] = {0, 6}; 
  Params_WidthFFT["Range_ntrk_2"] = {1, 6}; 
  Params_WidthFFT["Range_ntrk_3"] = {2.5, 7}; 
  Params_WidthFFT["Range_ntrk_4"] = {4, 8}; 
  Params_WidthFFT["m"] = {0.5, 0.5, 0.5, 0.5};
  Params_WidthFFT["m_G"] = {0, 0, 0, 0}; 
  Params_WidthFFT["s_s"] = {0.0001, 0.0001, 0.0001, 0.0001};
  Params_WidthFFT["s_e"] = {0.1, 0.1, 0.1, 0.1};
  Params_WidthFFT["fft_cache"] = {10000}; 
  Params_WidthFFT["Minimizer"] = {100000}; 

  TFile* X = new TFile("ntrk_ntru.root", "RECREATE"); 
  int p = 0;
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    TString current = x -> first;  
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 
    
    TH1F* ntrk_1_T = M["ntrk_1_T_I"][0]; 
    TH1F* ntrk_2_T = M["ntrk_2_T_I"][1]; 
    TH1F* ntrk_3_T = M["ntrk_3_T_I"][2]; 
    TH1F* ntrk_4_T = M["ntrk_4_T_I"][3]; 
    TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
    
    std::vector<TH1F*> Proposed = { ntrk_1_T,  ntrk_2_T, ntrk_3_T, ntrk_4_T };
    std::vector<TH1F*> ToBeUsed; 
    
    for ( int i(0); i < Proposed.size(); i++)
    {
      if (Proposed[i] -> GetEntries() < 5000){continue;}
      ToBeUsed.push_back(Proposed[i]); 
    }
    if (ToBeUsed.size() == 0){continue;}
    X -> mkdir(current); 
    X -> cd(current); 

    std::vector<TH1F*> Normal_Fits = Normalization_Fit(ToBeUsed, trk1_start, Params_N, current); 
    std::vector<TH1F*> ShiftNormal_Fits = NormalizationShift_Fit(ToBeUsed, trk1_start, Params_NS, current); 
    std::vector<TH1F*> ShiftNormalFFT_Fits = NormalizationShiftFFT_Fit(ToBeUsed, trk1_start, Params_FFT, current); 
    std::vector<TH1F*> ShiftNormalWidthFFT_Fits = NormalizationShiftWidthFFT_Fit(ToBeUsed, trk1_start, Params_WidthFFT, current); 
   
    TCanvas* can = new TCanvas(); 
    can -> SetLogy(); 
    can -> Print(current + ".pdf["); 
    PlotHists(Normal_Fits, ToBeUsed, can); 
    can -> Print(current + ".pdf"); 

    PlotHists(ShiftNormal_Fits, ToBeUsed, can); 
    can -> Print(current + ".pdf"); 
    
    PlotHists(ShiftNormalFFT_Fits, ToBeUsed, can); 
    can -> Print(current + ".pdf"); 

    PlotHists(ShiftNormalWidthFFT_Fits, ToBeUsed, can); 
    can -> Print(current + ".pdf"); 
    can -> Print(current + ".pdf]"); 


    BulkDelete(Normal_Fits); 
    BulkDelete(ShiftNormal_Fits); 
    BulkDelete(ShiftNormalFFT_Fits); 
    BulkDelete(ShiftNormalWidthFFT_Fits); 
    
    X -> cd(); 
    p++;
  }
}

void TestFits_AllTruth_ToTrack(TString JE, TString Mode, TString MCFile)
{
  if (MCFile == ""){ MCFile = "Merged_MC.root"; }
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(MCFile); 
 
  std::vector<std::vector<float>> Ranges = {{0.2, 6}, {0.5, 8}, {2, 8}, {2, 8}}; 

  float m = 0.2; 
  // Normalization parameters
  std::map<TString, std::vector<float>> Params_N; 
  Params_N["Range"] = {0, 10}; 
  Params_N["r_value"] = {1.5};
  Params_N["Range_ntrk_1"] = Ranges[0]; 
  Params_N["Range_ntrk_2"] = Ranges[1]; 
  Params_N["Range_ntrk_3"] = Ranges[2]; 
  Params_N["Range_ntrk_4"] = Ranges[3]; 

  // Normalization Shift parameters
  std::map<TString, std::vector<float>> Params_NS; 
  Params_NS["Range"] = {0.2, 8}; 
  Params_NS["r_value"] = {1.5};
  Params_NS["Range_ntrk_1"] = Ranges[0]; 
  Params_NS["Range_ntrk_2"] = Ranges[1]; 
  Params_NS["Range_ntrk_3"] = Ranges[2];   
  Params_NS["Range_ntrk_4"] = Ranges[3]; 
  Params_NS["dx"] = {m, m, m, m}; 
  Params_NS["dx_G"] = {0, 0, 0, 0};
  Params_NS["Minimizer"] = {100000}; 
  
  // Normalization Shift FFT parameters
  std::map<TString, std::vector<float>> Params_FFT; 
  Params_FFT["Range"] = {0.2, 8}; 
  Params_FFT["r_value"] = {1.5};
  Params_FFT["Range_ntrk_1"] = Ranges[0]; 
  Params_FFT["Range_ntrk_2"] = Ranges[1]; 
  Params_FFT["Range_ntrk_3"] = Ranges[2];   
  Params_FFT["Range_ntrk_4"] = Ranges[3]; 
  Params_FFT["m"] = {m, m, m, m};
  Params_FFT["m_G"] = {0, 0, 0, 0}; 
  Params_FFT["s_s"] = {0, 0, 0, 0};
  Params_FFT["s_e"] = {0.1, 0.1, 0.1, 0.1};
  Params_FFT["fft_cache"] = {50000}; 
  Params_FFT["Minimizer"] = {100000}; 

  // Normalization Shift Width FFT parameters
  std::map<TString, std::vector<float>> Params_WidthFFT; 
  Params_WidthFFT["Range"] = {0.2, 8}; 
  Params_WidthFFT["r_value"] = {1.5};
  Params_WidthFFT["Range_ntrk_1"] = Ranges[0]; 
  Params_WidthFFT["Range_ntrk_2"] = Ranges[1];   
  Params_WidthFFT["Range_ntrk_3"] = Ranges[2];  
  Params_WidthFFT["Range_ntrk_4"] = Ranges[3]; 
  Params_WidthFFT["m"] = {m, m, m, m};
  Params_WidthFFT["m_G"] = {0, 0, 0, 0}; 
  Params_WidthFFT["s_s"] = {0.0001, 0.0001, 0.0001, 0.0001};
  Params_WidthFFT["s_e"] = {0.1, 0.1, 0.1, 0.1};
  Params_WidthFFT["s_G"] = {0.001, 0.001, 0.001, 0.001};
  Params_WidthFFT["fft_cache"] = {50000}; 
  Params_WidthFFT["Minimizer"] = {100000}; 

  TFile* X = new TFile("Fit_Tracks.root", "RECREATE"); 
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    TString current = x -> first;  
    if (JE != ""){if (JE != current){continue;}} // This is for parallel computing capabilities on a cluster
 
    
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 
    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    std::vector<TH1F*> ntrk_2_T = M["ntrk_2_T_I"]; 
    std::vector<TH1F*> ntrk_3_T = M["ntrk_3_T_I"]; 
    std::vector<TH1F*> ntrk_4_T = M["ntrk_4_T_I"]; 
    TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
    std::vector<std::vector<TH1F*>> TruthVector = { ntrk_1_T,  ntrk_2_T, ntrk_3_T, ntrk_4_T };
   
    TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
    TH1F* ntrk_2_M = M["ntrk_2_M_I"][0]; 
    TH1F* ntrk_3_M = M["ntrk_3_M_I"][0]; 
    TH1F* ntrk_4_M = M["ntrk_4_M_I"][0]; 
    std::vector<TH1F*> Proposed = {ntrk_1_M, ntrk_2_M, ntrk_3_M, ntrk_4_M};  

    std::vector<TH1F*> ToBeUsed; 
    for ( int i(0); i < Proposed.size(); i++)
    {
      if (Proposed[i] -> GetEntries() < 5000){continue;}
      ToBeUsed.push_back(Proposed[i]); 
    }
    
    if (ToBeUsed.size() == 0){continue;}
    X -> mkdir(current); 
    X -> cd(current); 
    std::vector<TString> ntrk_String_T = {"ntrk_1_T", "ntrk_2_T", "ntrk_3_T", "ntrk_4_T"};  
    for (int p(0); p < TruthVector.size(); p++)
    {
      gDirectory -> mkdir(ntrk_String_T[p]); 
      gDirectory -> cd(ntrk_String_T[p]); 
      BulkWrite(TruthVector[p]); 
      gDirectory -> cd("../");
    }

    std::cout << "++++++++" << current << " " << Mode << std::endl;


    if (Mode == "")
    {
      std::vector<std::vector<TH1F*>> Normal_Fits = Normalization_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_N, current);
      std::vector<std::vector<TH1F*>> ShiftNormal_Fits = NormalizationShift_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_NS, current);   
      std::vector<std::vector<TH1F*>> ShiftNormalFFT_Fits = NormalizationShiftFFT_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_FFT, current);
      std::vector<std::vector<TH1F*>> ShiftNormalWidthFFT_Fits = NormalizationShiftWidthFFT_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_WidthFFT, current);
      std::vector<std::vector<TH1F*>> Experimental_Fits = Experimental_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_WidthFFT, current);

      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(current + ".pdf["); 
      
      for (int i(0); i < ToBeUsed.size(); i++)
      {
        PlotHists(Normal_Fits[i], TruthVector[i], can); 
        can -> Print(current + ".pdf"); 

        PlotHists(ShiftNormal_Fits[i], TruthVector[i], can); 
        can -> Print(current + ".pdf"); 
        
        PlotHists(ShiftNormalFFT_Fits[i], TruthVector[i], can); 
        can -> Print(current + ".pdf"); 

        PlotHists(ShiftNormalWidthFFT_Fits[i], TruthVector[i], can); 
        can -> Print(current + ".pdf"); 

        PlotHists(Experimental_Fits[i], TruthVector[i], can); 
        can -> Print(current + ".pdf"); 

        BulkDelete(Normal_Fits[i]); 
        BulkDelete(ShiftNormal_Fits[i]); 
        BulkDelete(ShiftNormalFFT_Fits[i]); 
        BulkDelete(ShiftNormalWidthFFT_Fits[i]); 
        BulkDelete(Experimental_Fits[i]); 
      }
      can -> Print(current + ".pdf]"); 
    }
    

    if (Mode == "Normal")
    {
      TString name = current + "_" + Mode + ".pdf"; 
      std::vector<std::vector<TH1F*>> Fits = Normalization_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_N, current);
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(name + "["); 
      
      for (int i(0); i < ToBeUsed.size(); i++)
      {
        PlotHists(Fits[i], TruthVector[i], can); 
        can -> Print(name); 
      }
      can -> Print(name + "]"); 
      for (int i(0); i < Fits.size(); i++){BulkDelete(Fits[i]);}
    }

    if (Mode == "ShiftNormal")
    {
      TString name = current + "_" + Mode + ".pdf"; 
      std::vector<std::vector<TH1F*>> Fits = NormalizationShift_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_NS, current);
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(name + "["); 
      
      for (int i(0); i < ToBeUsed.size(); i++)
      {
        PlotHists(Fits[i], TruthVector[i], can); 
        can -> Print(name); 
      }
      can -> Print(name + "]"); 
      for (int i(0); i < Fits.size(); i++){BulkDelete(Fits[i]);}
    }


    if (Mode == "ShiftNormalFFT")
    {

      TString name = current + "_" + Mode + ".pdf"; 
      std::vector<std::vector<TH1F*>> Fits = NormalizationShiftFFT_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_FFT, current);
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(name + "["); 
      
      for (int i(0); i < ToBeUsed.size(); i++)
      {
        PlotHists(Fits[i], TruthVector[i], can); 
        can -> Print(name); 
      }
      can -> Print(name + "]"); 
      for (int i(0); i < Fits.size(); i++){BulkDelete(Fits[i]);}
    }



    if (Mode == "ShiftNormalWidthFFT")
    {
      TString name = current + "_" + Mode + ".pdf"; 
      std::vector<std::vector<TH1F*>> Fits = NormalizationShiftWidthFFT_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_WidthFFT, current);
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(name + "["); 
      
      for (int i(0); i < ToBeUsed.size(); i++)
      {
        PlotHists(Fits[i], TruthVector[i], can); 
        can -> Print(name); 
      }
      can -> Print(name + "]"); 
      for (int i(0); i < Fits.size(); i++){BulkDelete(Fits[i]);}
    }

    if (Mode == "Experimental")
    {
      TString name = current + "_" + Mode + ".pdf"; 
      std::vector<std::vector<TH1F*>> Fits = Experimental_Fit_NtrkMtru(ToBeUsed, trk1_start, Params_WidthFFT, current);
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      can -> Print(name + "["); 
      
      for (int i(0); i < ToBeUsed.size(); i++)
      {
        PlotHists(Fits[i], TruthVector[i], can); 
        can -> Print(name); 
      }
      can -> Print(name + "]"); 
      for (int i(0); i < Fits.size(); i++){BulkDelete(Fits[i]);}
    }
    X -> cd(current); 
    X -> Write();
    X -> cd(); 
  }

  X -> Close(); 













}
