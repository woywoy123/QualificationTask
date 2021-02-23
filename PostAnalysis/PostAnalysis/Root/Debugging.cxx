#include<PostAnalysis/Debugging.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/Plotting.h>

void Entry()
{
  int mode = 4; 

  if (mode == 1){FitWithConstraint();}
  if (mode == 2){ExplicitConstrain();}
  if (mode == 3){ExplicitConstrainExternal();}
  if (mode == 4){StepByStepFit();}
}

void FitWithConstraint()
{
  auto Loops =[] (std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params)
  {
    int bins = Data -> GetNbinsX(); 
    float min = Data -> GetXaxis() -> GetXmin(); 
    float max = Data -> GetXaxis() -> GetXmax(); 
    
    std::vector<TH1F*> PSF; 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      TString name = "Gaussian_"; name += (x+1); 
      float m = Params["G_Mean"][x]; 
      float s = Params["G_Stdev"][x];
      TH1F* Gaus = Gaussian(m ,s, bins, min, max, name); 
      PSF.push_back(Gaus);
    }
    
    std::vector<TString> Names_Dec; 
    for (int i(0); i < ntrk_Conv.size(); i++)
    {
      TString name = "Dec_"; name += (ntrk_Conv[i] -> GetTitle()); 
      Names_Dec.push_back(name);
    }
    std::vector<TH1F*> PDF_D = CloneTH1F(Data, Names_Dec); 
    MultiThreadingDeconvolutionExperimental(ntrk_Conv, PSF, PDF_D, Params["LR_iterations"][0]); 
   std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitWithConstraint(Data, PDF_D, Params, Params["cache"][0], Params["cache"][0]);
    BulkDelete(PSF); 
    BulkDelete(PDF_D); 

    std::vector<TH1F*> F_C;
    for (int i(0); i < trk_Fit.size(); i++){F_C.push_back(trk_Fit[i].first); }
    return F_C; 
  };


  std::map<TString, std::vector<TH1F*>> MC = Experimental_MC_Reader("Merged_MC.root"); 
  std::vector<TH1F*> All = MC["All"]; 
  
  TH1F* Trk1_ideal = All[0]; 
  
  std::vector<TH1F*> Trk2 = {All[3], All[4], All[5], All[6]}; 
  TH1F* Data = SumHists(Trk2, "Data_2Trk"); 
  
  float m = 1; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params["x_range"] = {0.01, 11.5}; 
  Params["iterations"] = {30}; 
  Params["LR_iterations"] = {150}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 
  
  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(Data, 4, "C"); 
  std::vector<TH1F*> Result = Loops(ntrk_Conv, Data, Params);
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  PlotHists(Data, Trk2, Result, can); 
  can -> Print("Constraint_Debugging.pdf"); 
}

void ExplicitConstrain()
{
  float m = 0.001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params["x_range"] = {0.01, 11.5}; 
  Params["iterations"] = {30}; 
  Params["LR_iterations"] = {150}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 

  std::cout << "Here" << std::endl;
  std::map<TString, std::vector<TH1F*>> MC = Experimental_MC_Reader("Merged_MC.root"); 
  std::vector<TH1F*> All = MC["All"]; 

  // Ideal PDF being used 
  TH1F* trk1 = All[0]; 

  // Empty TH1F used for the deconvolution using the PSF
  TH1F* trk_D = CloneTH1F(trk1, {"trk_D"})[0]; 

  // Creating the PSF
  int bins = trk1 -> GetNbinsX(); 
  float min = trk1 -> GetXaxis() -> GetXmin(); 
  float max = trk1 -> GetXaxis() -> GetXmax(); 
  TH1F* Gaus = Gaussian(Params["G_Stdev"][0], Params["G_Stdev"][1], bins, min, max, "psf"); 
  
  // Deconvolve the PDF with the PSF 
  DeconvolutionExperimental(trk1, Gaus, trk_D, Params["LR_iterations"][0]); 
  
  TH1F* Output = ExplicitConstraining(trk1, trk_D, Params); 

}

void ExplicitConstrainExternal()
{
  float m = 0.001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params["x_range"] = {0.01, 11.5}; 
  Params["iterations"] = {30}; 
  Params["LR_iterations"] = {150}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 

  std::cout << "Here" << std::endl;
  std::map<TString, std::vector<TH1F*>> MC = Experimental_MC_Reader("Merged_MC.root"); 
  std::vector<TH1F*> All = MC["All"]; 

  // Ideal PDF being used 
  TH1F* trk1 = All[0]; 

  // Empty TH1F used for the deconvolution using the PSF
  TH1F* trk_D = CloneTH1F(trk1, {"trk_D"})[0]; 

  // Creating the PSF
  int bins = trk1 -> GetNbinsX(); 
  float min = trk1 -> GetXaxis() -> GetXmin(); 
  float max = trk1 -> GetXaxis() -> GetXmax(); 
  TH1F* Gaus = Gaussian(Params["G_Stdev"][0], Params["G_Stdev"][1], bins, min, max, "psf"); 
  
  // Deconvolve the PDF with the PSF 
  DeconvolutionExperimental(trk1, Gaus, trk_D, Params["LR_iterations"][0]); 
  
  TH1F* Output = ExplicitConstrainingExternal(trk1, trk_D, Params); 

}

void StepByStepFit()
{
  std::map<TString, std::vector<TH1F*>> MC = Experimental_MC_Reader("Merged_MC.root"); 
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
  
  float m = 0.1; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params["s_e"] = {0.05, 0.05, 0.1, 0.1};  
  Params["iterations"] = {30}; 
  Params["LR_iterations"] = {400}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.025, 0.025, 0.025, 0.025}; 
  Params["cache"] = {100000}; 
  Params["x_range"] = {0.1, max-1}; 

  // Empty TH1F used for the deconvolution using the PSF
  std::vector<TString> names = NameGenerator(trk1.size(), "_trk_D"); 
  std::vector<TH1F*> trk_D = CloneTH1F(trk1[0], names); 

  TH1F* Gaus = Gaussian(Params["G_Stdev"][0], Params["G_Stdev"][1], bins, min, max, "psf"); 
  std::vector<TH1F*> psf_vector(trk1.size(), Gaus); 

  MultiThreadingDeconvolutionExperimental(trk1, psf_vector, trk_D, Params["LR_iterations"][0]); 
  //TH1F* trk = (TH1F*)Summed_MC -> Clone("trk"); 
  //std::vector<TH1F*> trk_1 = {trk}; 


  // Deconvolve the PDF with the PSF 
  std::vector<TH1F*> Output = IterativeFitting(Summed_MC, trk_D, Params, Params["LR_iterations"][0], Params["LR_iterations"][0]); 

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  PlotHists(Summed_MC, Truth, Output, can); 
  can -> Print("debugging.pdf"); 

}



