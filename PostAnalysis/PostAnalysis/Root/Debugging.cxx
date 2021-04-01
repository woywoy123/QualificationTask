#include<PostAnalysis/Debugging.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/Plotting.h>

void Entry()
{
  int mode = 4; 

  if (mode == 4){StepByStepFit();}
}

void StepByStepFit()
{
  std::map<TString, std::vector<TH1F*>> MC = MC_Reader("Merged_MC.root"); 
  std::vector<TH1F*> All = MC["All"]; 

  // Ideal PDF being used 
  std::vector<TH1F*> Truth = {All[0], All[1], All[2], All[3]}; 
  TH1F* Summed_MC = SumHists(Truth, "Data_track_1"); 
  std::vector<TH1F*> ntrk_C = ConvolveNTimes(Summed_MC, 4, "_C"); 
  std::vector<TH1F*> trk1 = {ntrk_C[0], ntrk_C[1], ntrk_C[2], ntrk_C[3]}; 

  // Creating the PSF
  int bins = trk1[0] -> GetNbinsX(); 
  float min = trk1[0] -> GetXaxis() -> GetXmin(); 
  float max = trk1[0] -> GetXaxis() -> GetXmax(); 
  
  float m = 0.5; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.01, 0.01, 0.01, 0.01}; 
  Params["cache"] = {100000}; 
  Params["x_range"] = {0.1, max-1}; 

  // ====== WARNING! REPLACED CONVOLVED HISTS WITH TRUTH!!!! ============= //


  // Empty TH1F used for the deconvolution using the PSF
  std::vector<TString> names = NameGenerator(Truth.size(), "_trk_D"); 
  std::vector<TH1F*> trk_D = CloneTH1F(Truth[0], names); 

  TH1F* Gaus = Gaussian(Params["G_Stdev"][0], Params["G_Stdev"][1], bins, min, max, "psf"); 
  std::vector<TH1F*> psf_vector(Truth.size(), Gaus); 

  MultiThreadingDeconvolutionExperimental(Truth, psf_vector, trk_D, Params["LR_iterations"][0]); 
  //TH1F* trk = (TH1F*)Summed_MC -> Clone("trk"); 
  //std::vector<TH1F*> trk_1 = {trk}; 


  // Deconvolve the PDF with the PSF 
  //std::vector<TH1F*> Output = IterativeFitting(Summed_MC, trk_D, Params, Params["LR_iterations"][0], Params["LR_iterations"][0]); 
  
  std::vector<std::pair<TH1F*, std::vector<float>>> Tracks_Result = FitDeconvolutionPerformance(Summed_MC, trk_D, Params, Params["cache"][0], Params["cache"][0]);

  std::vector<TH1F*> Output; 
  for (int p(0); p < Tracks_Result.size(); p++)
  {
    TH1F* H = Tracks_Result[p].first;     
    Output.push_back(H); 
  }


  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  PlotHists(Summed_MC, Truth, Output, can); 
  can -> Print("debugging.pdf"); 

}



