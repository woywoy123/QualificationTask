#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/AlgorithmFunctions.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth, int trk_Data)
{
  auto Smear =[](TH1F* Data, float stdev)
  {
    float lumi = Data -> Integral(); 
    int bins = Data -> GetNbinsX(); 
    float min = Data -> GetXaxis() -> GetXmin(); 
    float max = Data -> GetXaxis() -> GetXmax(); 
    float width = (max - min) / float(bins); 
    min += width / 2.; 
    max += width / 2.; 
    TH1F* Gaus = Gaussian(0 ,stdev , bins, min, max, "Tempo"); 
    Convolution(Data, Gaus, Data); 
    Normalize(Data);
    Data -> Scale(lumi); 
    delete Gaus;  
  }; 

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int iterations = Params["iterations"][0]; 
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

  std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
  Flush(F_C, ntrk_Conv, false); 

  TString name = Data[trk_Data] -> GetTitle(); name += ("_Experimental.pdf"); 
  std::vector<TString> Names_Dec; 
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString name = "Temp_"; name += (ntrk_Conv[i] -> GetTitle()); 
    Names_Dec.push_back(name);
  }
  
  for (int i(0); i < iterations; i++)
  {
    float s = 0.02*(0.95 - float(i+1) / float (iterations)); 
    
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");

    std::vector<TH1F*> Not_trk; 
    std::vector<TH1F*> Not_PSF;
    if (s > 0){Smear(Data_Copy, s);} 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x == trk_Data){ Data_Copy -> Add(ntrk_Conv[x], -1); }
      else
      {
        Not_trk.push_back(ntrk_Conv[x]); 
        Not_PSF.push_back(PSF[x]);
      }
    }
    //Average(Data_Copy);  
    F_C = LoopGen(Not_trk, Not_PSF, Data_Copy, Params);  
    //ScaleShape(Data_Copy, F_C);  
    Flush(F_C, Not_trk, true);

    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print(name); 
    delete Data_Copy; 

    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
    
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
   
    if (s > 0){Smear(Data_Copy, s);} 
    //Average(Data_Copy); 
    F_C = LoopGen(Delta, Delta_PSF, Data_Copy, Params); 
    ScaleShape(Data_Copy, F_C); 
    Flush(F_C, Delta, true); 

    delete Data_Copy; 

    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
    if (s > 0){Smear(Data_Copy, s);} 
    
    F_C = CloneTH1F(Data_Copy, Names_Dec);
    for (int x(0); x < F_C.size(); x++){F_C[x] -> Add(ntrk_Conv[x], 1);}
    //Average(Data_Copy); 
    
    ScaleShape(Data_Copy, F_C); 
    Flush(F_C, ntrk_Conv, true); 
    delete Data_Copy;
    
  }
  std::map<TString, std::vector<TH1F*>> Out; 
  return Out; 
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

  float m = 0.001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.04, 0.04, 0.04, 0.04};  
  Params["x_range"] = {0.05, 11.5}; 
  Params["iterations"] = {50}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.02, 0.02, 0.02, 0.02}; 
  Params["cache"] = {10000}; 

  TString Dir = "Merged_MC.root"; 
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
  MainAlgorithm(Data, Params, Track1, 0); 
  MainAlgorithm(Data, Params, Track2, 1); 
  MainAlgorithm(Data, Params, Track3, 2); 
  MainAlgorithm(Data, Params, Track3, 3);  
  
}


void Shifting(TH1F* H)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int LucyRichardson = 300; 
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
  for (int i(0); i < 50; i++)
  {
    Convolution(Dec, Gaus, Temp); 
    DeconvolutionExperimental(Temp, Gaus, Dec, LucyRichardson); 
     
    
    
    //Convolution(Dec, Gaus, Temp); 
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


