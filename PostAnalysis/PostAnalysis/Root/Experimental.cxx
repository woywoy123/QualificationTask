#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/AlgorithmFunctions.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth, int trk_Data)
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

  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(Data[0], Data.size(), "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 

  std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
  Flush(F_C, ntrk_Conv, false); 

  TString name = Data[trk_Data] -> GetTitle(); name += (".pdf"); 

  std::vector<TString> Names_Dec; 
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString name = "Temp_"; name += (ntrk_Conv[i] -> GetTitle()); 
    Names_Dec.push_back(name);
  }
 
  float Error = 1e10; 
  for (int i(0); i < iterations; i++)
  {

    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
    F_C = CloneTH1F(Data_Copy, Names_Dec);
    for (int x(0); x < F_C.size(); x++){F_C[x] -> Add(ntrk_Conv[x], 1);}
    ScaleShape(Data_Copy, F_C); 

    //float E = FitError(Data_Copy, ntrk_Conv); 
    //if (E < Error){Flush(F_C, ntrk_Conv, true); Error = E;} 
    //else{BulkDelete(F_C);}  
    Flush(F_C, ntrk_Conv, true); 

    std::vector<TH1F*> Not_trk; 
    std::vector<TH1F*> Not_PSF; 
    std::vector<TH1F*> Truth_Not;  
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x == trk_Data){ Data_Copy -> Add(ntrk_Conv[x], -1); }
      else
      {
        Not_trk.push_back(ntrk_Conv[x]); 
        Not_PSF.push_back(PSF[x]);
        Truth_Not.push_back(Truth[x]); 
      }
    }
    //Average(Data_Copy); 

    F_C = LoopGen(Not_trk, Not_PSF, Data_Copy, Params);  
    Flush(F_C, Not_trk, true);
    
    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print(name); 

    delete Data_Copy; 
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 

    std::vector<TH1F*> Delta;
    std::vector<TH1F*> Delta_PSF; 

    float heat = 0; //1 - float(1.)/float(iterations-i);
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      if (x != trk_Data){Data_Copy -> Add(ntrk_Conv[x], (-1+0.01*heat));}
     
      if ( x == trk_Data )
      {
        Delta.push_back(ntrk_Conv[x]); 
        Delta_PSF.push_back(PSF[x]); 
      }
    }
    //Average(Data_Copy); 
   
    F_C = LoopGen(Delta, Delta_PSF, Data_Copy, Params);     
    //ScaleShape(Data_Copy, F_C); 
    //E = FitError(Data_Copy, Delta); 
    //if (E < Error){Flush(F_C, Delta, true); Error = E; std::cout << "LOWER!!!" << std::endl;} 
    //else{BulkDelete(F_C); std::cout << "Higher " << E << std::endl;}
    Flush(F_C, Delta, true); 


    delete Data_Copy; 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      for (int t(0); t < ntrk_Conv.size(); t++)
      {

        if (x == t){continue;}
        Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
        
        std::vector<TH1F*> Temp; 
        std::vector<TH1F*> Temp_PSF; 
        for (int p(0); p < ntrk_Conv.size(); p++)
        {
          if (p == x || p == t) { Data_Copy -> Add(ntrk_Conv[p], (-1+0.01*heat)); } 
          else
          {
            Temp.push_back(ntrk_Conv[p]); 
            Temp_PSF.push_back(PSF[p]); 
          }
        }
        //Average(Data_Copy); 
        F_C = LoopGen(Temp, Temp_PSF, Data_Copy, Params); 
        //ScaleShape(Data_Copy, F_C); 
        float E = FitError(Data_Copy, F_C); 
        if (E < Error){Flush(F_C, Temp, true); Error = E; std::cout << "LOWER!!!" << std::endl;} 
        else{BulkDelete(F_C); std::cout << "Higher " << E << std::endl;}
        //Flush(F_C, Temp, true); 
    
        PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
        can -> Print(name); 

        delete Data_Copy; 
      }
    }
    
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

  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-0.001, -0.001, -0.001, -0.001}; 
  Params["m_e"] = {0.01, 0.01, 0.01, 0.01}; 
  Params["s_s"] = {0.04, 0.04, 0.04, 0.04};
  Params["s_e"] = {0.06, 0.06, 0.06, 0.06};  
  Params["x_range"] = {0.01, 9.8}; 
  Params["iterations"] = {100}; 
  Params["LR_iterations"] = {100}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 

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
  MainAlgorithm(Data, Params, Track1, 0); 
  MainAlgorithm(Data, Params, Track2, 1); 
  MainAlgorithm(Data, Params, Track3, 2); 
  MainAlgorithm(Data, Params, Track3, 3);  
  
  Shifting(Data[0]); 
  
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
    Deconvolution(Temp, Gaus, Dec, LucyRichardson); 
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


