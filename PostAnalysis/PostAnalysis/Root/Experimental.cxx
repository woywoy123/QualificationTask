#include<PostAnalysis/Experimental.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth)
{
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

  std::vector<TString> Names_Dec = NameGenerator(Data.size(), "Dec_Trk_");  
  std::vector<TH1F*> PDF_D = CloneTH1F(Data[0], Names_Dec);  

  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(Data[0], Data.size(), "C"); 
  
  for (int i(0); i < iterations; i++)
  {
    MultiThreadingDeconvolution(ntrk_Conv, PSF, PDF_D, Params["LR_iterations"][0]); 

    std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data[0], PDF_D, Params, 10000, 10000);
    
    if (i == 0)
    {
      for (int p(0); p < trk_Fit.size(); p++)
      {
        TH1F* New = trk_Fit[p].first;
        ntrk_Conv[p] -> Reset(); 
        ntrk_Conv[p] -> Add(New); 
      }
    }

    //ScaleShape(Data[0], ntrk_Conv); 
    for (int p(0); p < trk_Fit.size(); p++)
    {
      TH1F* New = trk_Fit[p].first; 
      for (int z(0); z < bins; z++)
      {
        float e = New -> GetBinContent(z+1); 
        float f = ntrk_Conv[p] -> GetBinContent(z+1); 
        if (i > 0)
        {
          ntrk_Conv[p] -> SetBinContent(z+1, (e+f)/2.); 
        }
      } 
      ntrk_Conv[p] -> Reset(); 
      ntrk_Conv[p] -> Add(New); 
      delete New;
    }


    PlotHists(Data[0], Truth, ntrk_Conv, can); 
    can -> Print("Debug.pdf"); 

  }
  std::map<TString, std::vector<TH1F*>> Out; 
  return Out; 
}




void Shifting(TH1F* H)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int LucyRichardson = 500; 
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
  TH1F* Dec = (TH1F*)H -> Clone("Dec"); 
  Dec -> Reset(); 
  
  TH1F* Temp = (TH1F*)H -> Clone("Temp"); 

  Gaus -> Draw("HIST"); 
  can -> Print("Debug.pdf"); 

  H -> GetYaxis() -> SetRangeUser(1e-6, 1);
  for (int i(0); i < 100; i++)
  {
    MultiThreadingDeconvolution({Temp}, {Gaus}, {Dec}, LucyRichardson); 
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
  Params["m_s"] = {-0.1, -0.1, -0.1, -0.1}; 
  Params["m_e"] = {0.1, 0.1, 0.1, 0.1}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params["x_range"] = {0, 10}; 
  Params["iterations"] = {50}; 
  Params["LR_iterations"] = {200}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.1, 0.1, 0.1, 0.1}; 

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
  //MainAlgorithm(Data, Params, Track1); 
  Shifting(Data[0]); 
}
