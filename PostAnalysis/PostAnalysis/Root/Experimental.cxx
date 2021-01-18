#include<PostAnalysis/Experimental.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth, int trk_Data)
{
  auto ShapeSigmoid =[] (TH1F* trk_Fit, TH1F* ntrk_Conv)
  {
    for (int z(0); z < trk_Fit -> GetNbinsX(); z++)
    {
      float e = trk_Fit -> GetBinContent(z+1); 
      float f = ntrk_Conv -> GetBinContent(z+1); 
      
      if ( e == 0 && f == 0) {continue;}
      float dif = std::pow(-1 * std::abs(e-f)/(e+f), 1); 
      float sig = 1./(1+std::exp(dif)); 
      ntrk_Conv -> SetBinContent(z+1, e*(1-sig)+f*sig); 
    } 
  };

  auto LoopGen =[&](std::vector<TH1F*> ntrk_Conv, std::vector<TH1F*> PSF, TH1F* Data, std::map<TString, std::vector<float>> Params)
  {
    std::vector<TString> Names_Dec; 
    int index = 0;  
    for (int i(0); i < ntrk_Conv.size(); i++)
    {
      TString name = "Dec_"; name += (ntrk_Conv[i] -> GetTitle()); 
      Names_Dec.push_back(name);
    }
    std::vector<TH1F*> PDF_D = CloneTH1F(Data, Names_Dec);
    
    MultiThreadingDeconvolutionExperimental(ntrk_Conv, PSF, PDF_D, Params["LR_iterations"][0]); 
    std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data, PDF_D, Params, 10000, 10000);
    
    std::vector<TH1F*> F_C;
    std::vector<float> Lumi; 
    for (int i(0); i < trk_Fit.size(); i++)
    {
      F_C.push_back(trk_Fit[i].first); 
      Lumi.push_back(trk_Fit[i].first -> Integral()); 
    }
    
    if (ntrk_Conv.size() != 1)
    {
      //Scale(Data, F_C); 
      ScaleShape(Data, F_C, 1); 
    }
    for (int i(0); i < ntrk_Conv.size(); i++){ShapeSigmoid(F_C[i], ntrk_Conv[i]); } 

    for (int i(0); i < PDF_D.size(); i++){delete PDF_D[i];}
    for (int i(0); i < F_C.size(); i++){delete F_C[i];}
    return Lumi; 
  };

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int iterations = Params["iterations"][0]; 
  int LucyRichardson = Params["LR_iterations"][0]; 
  int bins = Data[0] -> GetNbinsX(); 
  float min = Data[0] -> GetXaxis() -> GetXmin(); 
  float max = Data[0] -> GetXaxis() -> GetXmax(); 
  float width = (max - min) / float(bins); 
  min += width / 2.; 
  max += width / 2.; 
  
  std::vector<TH1F*> PSF; 
  for (int x(0); x < Data.size(); x++)
  {
    TString name = "Gaussian_"; name += (x+1); 
    float m = Params["G_Mean"][x]; 
    float s = Params["G_Stdev"][x];
    TH1F* Gaus = Gaussian(m ,s, bins, min, max, name); 
    PSF.push_back(Gaus);
  }
  PlotHists(PSF, can);  
  can -> Print("Gaus.pdf"); 

  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(Data[0], Data.size(), "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
  std::vector<float> L = LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
  
  for (int i(0); i < L.size(); i++){ ntrk_Conv[i] -> Scale(L[i]); }

  for (int i(0); i < iterations; i++)
  {
    LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print("Debug.pdf"); 

    std::vector<TH1F*> Delta;
    std::vector<TH1F*> Delta_PSF; 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x != trk_Data){Data_Copy -> Add(ntrk_Conv[x], -1);}
      else 
      { 
        Delta.push_back(ntrk_Conv[x]); 
        Delta_PSF.push_back(PSF[x]); 
      }
    }
    
    LoopGen(Delta, Delta_PSF, Data_Copy, Params); 
    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print("Debug.pdf"); 

    delete Data_Copy; 
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
  }
  std::map<TString, std::vector<TH1F*>> Out; 
  return Out; 
}

void Shifting(TH1F* H)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int LucyRichardson = 100; 
  int bins = H -> GetNbinsX(); 
  float min = H -> GetXaxis() -> GetXmin(); 
  float max = H -> GetXaxis() -> GetXmax(); 
  float width = (max - min) / float(bins); 
  min += width / 2.; 
  max += width / 2.; 
  
  std::vector<TH1F*> PSF; 
  float m = 0; 
  float s = 0.05;
  TH1F* Gaus = Gaussian(m ,s, bins, min, max, "Gaus"); 
  TH1F* Dec = (TH1F*)H -> Clone("DecP"); 
  TH1F* Temp = (TH1F*)H -> Clone("Temp"); 

  Gaus -> Draw("HIST"); 
  can -> Print("Debug.pdf"); 

  H -> GetYaxis() -> SetRangeUser(1e-6, 1);
  for (int i(0); i < 100; i++)
  {
    DeconvolutionExperimental(Temp, Gaus, Dec, LucyRichardson); 
    Convolution(Dec, Gaus, Temp); 
    Normalize(Temp); 
    Normalize(H); 
    H -> SetLineStyle(kSolid); 
    H -> Draw("HIST");
    Temp -> SetLineColor(kRed);
    Temp -> SetLineStyle(kDashed); 
    Temp -> Draw("SAMEHIST"); 
    can -> Print("Debug.pdf"); 
  }

}

void AlgorithmMonteCarlo()
{
  auto Sum_Hist =[] (std::vector<TH1F*> Hists, TString Name)
  {
    TH1F* Hist = (TH1F*)Hists[0] -> Clone(Name); 
    Hist -> Reset();  
    Hist -> SetTitle(Name); 
    for (int i(0); i < Hists.size(); i++)
    {
      Hist -> Add(Hists[i]); 
    }
    return Hist; 
  };

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.001, -0.001, -0.001, -0.001}; 
  Params["m_e"] = {0.001, 0.001, 0.001, 0.001}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.15, 0.15, 0.15, 0.15};  
  Params["x_range"] = {0.1, 10}; 
  Params["iterations"] = {150}; 
  Params["LR_iterations"] = {150}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 

  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarloLayerEnergy(Dir); 
  std::vector<TH1F*> Track1 = MC["Track_1_All"];
  std::vector<TH1F*> Track2 = MC["Track_2_All"];
  std::vector<TH1F*> Track3 = MC["Track_3_All"];
  std::vector<TH1F*> Track4 = MC["Track_4_All"];
  
  TH1F* Trk1 = Sum_Hist(Track1, "trk1_data"); 
  TH1F* Trk2 = Sum_Hist(Track2, "trk2_data"); 
  TH1F* Trk3 = Sum_Hist(Track3, "trk3_data"); 
  TH1F* Trk4 = Sum_Hist(Track4, "trk4_data"); 
  std::vector<TH1F*> Data = {Trk1, Trk2, Trk3, Trk4}; 
  MainAlgorithm(Data, Params, Track2, 1); 
  //Shifting(Data[0]); 
  
}
