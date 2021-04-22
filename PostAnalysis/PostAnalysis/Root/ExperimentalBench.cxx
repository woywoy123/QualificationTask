#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/Plotting.h>

void IOTest(){TestReadCTIDE();}
void RooFitBaseFunctionTest()
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE("Merged_MC.root");
  
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 

    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 

    // Test the Normalization function out
    std::vector<TString> ntrk1_Name = NameGenerator(ntrk_1_T.size(), "dEdx_ntrk_1_tru_"); 
    std::vector<TH1F*> PDF_H = BulkClone(ntrk_1_T, ntrk1_Name); 
    
    std::map<TString, std::vector<float>> Params_N; 
    Params_N["Range"] = {1, 8}; 
    Params_N["r_value"] = {1.1};
    
    Normalization(ntrk_1_M, PDF_H, Params_N, "Normal"); 
  
    TCanvas* can = new TCanvas(); 
    can -> SetLogy();
    PlotHists(ntrk_1_T, PDF_H, can); 
    can -> Print("Example.pdf"); 
    delete can; 

    std::map<TString, std::vector<float>> Params_NS; 
    Params_NS["Range"] = {1, 8}; 
    Params_NS["r_value"] = {1.1};
    Params_NS["dx"] = {0.5, 0.5, 0.5, 0.5}; 
    Params_NS["dx_G"] = {0, 0, 0, 0};
    NormalizationShift(ntrk_1_M, PDF_H, Params_NS, "NormalShift"); 
    
    TCanvas* can1 = new TCanvas(); 
    can1 -> SetLogy();
    PlotHists(ntrk_1_T, PDF_H, can1); 
    can1 -> Print("ExampleS.pdf"); 
    delete can1; 


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

    ConvolutionFFT(ntrk_1_M, PDF_H_FFT, Params_FFT, "NormalShiftFFT"); 
    TCanvas* can2 = new TCanvas(); 
    can2 -> SetLogy();
    PlotHists(ntrk_1_T, PDF_H_FFT, can2); 
    can2 -> Print("ExampleSFFT.pdf"); 
    delete can2; 


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

  TFile* F = new TFile("ntrk_ntru.root", "RECREATE"); 
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    TString current = x -> first;  
    std::cout << current << std::endl;
    F -> mkdir(current); 
    F -> cd(current); 

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
      if (Proposed[i] -> GetEntries() < 1000){continue;}
      ToBeUsed.push_back(Proposed[i]); 
    }
    if (ToBeUsed.size() == 0){continue;}

    std::vector<TH1F*> Normal_Fits = Normalization_Fit(ToBeUsed, trk1_start, Params_N, current); 
   







    F -> cd(); 

  }

}


