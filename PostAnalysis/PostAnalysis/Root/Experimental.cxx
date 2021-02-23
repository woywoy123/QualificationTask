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

  std::vector<TH1F*> PSF; 
  for (int x(0); x < Data.size(); x++)
  {
    TString name = "Gaussian_"; name += (x+1); 
    float m = Params["G_Mean"][x]; 
    float s = Params["G_Stdev"][x];
    TH1F* Gaus = Gaussian(m ,s, bins, min, max, name); 
    PSF.push_back(Gaus);
  }

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

  PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
  can -> Print(name); 

  for (int i(0); i < iterations; i++)
  {
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> not_trk; 
    std::vector<TH1F*> not_psf;
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p != trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else
      { 
        not_trk.push_back(ntrk_Conv[p]); 
        not_psf.push_back(PSF[p]);  
      }
    }
    F_C = LoopGen(not_trk, not_psf, Data_Copy, Params); 
    ScaleShape(Data_Copy, not_trk);  
    Flush(F_C, not_trk, true); 
    PlotHists(Data_Copy, {Truth[1], Truth[2],Truth[3]}, not_trk, can); 
    can -> Print(name); 
    delete Data_Copy;    
   
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    Smear(Data_Copy, 0.05); 
    F_C = LoopGen(ntrk_Conv, PSF, Data_Copy, Params); 
    Flush(F_C, ntrk_Conv, true); 

    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print(name); 
 
    delete Data_Copy; 


  }
  BulkDelete(PSF);
  std::map<TString, std::vector<TH1F*>> Out; 
  return Out; 
}

void AlgorithmMonteCarlo()
{
  int nsdo = 4;
  float m = 0.1; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params["s_e"] = {0.05, 0.05, 0.05, 0.05};  
  Params["iterations"] = {10}; 
  Params["LR_iterations"] = {400}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 
  Params["x_range"] = {0.1, 11.5}; 

  TString Dir = "Merged_MC.root"; 
 
  std::map<TString, std::vector<TH1F*>> Output = Experimental_MC_Reader(Dir);

  typedef std::map<TString, std::vector<TH1F*>>::iterator M;
  for (M x = Output.begin(); x != Output.end(); x++)
  {
    std::vector<TH1F*> Trks = x -> second; 
    std::cout << x -> first << std::endl; 
    std::vector<TH1F*> T1;  
    std::vector<TH1F*> T2;  
    std::vector<TH1F*> T3;  
    std::vector<TH1F*> T4;  
    for (int i(0); i < nsdo; i++){T1.push_back(Trks[i]);}
    for (int i(4); i < nsdo+4; i++){T2.push_back(Trks[i]);}
    for (int i(8); i < nsdo+8; i++){T3.push_back(Trks[i]);}
    for (int i(12); i < nsdo+12; i++){T4.push_back(Trks[i]);}
    TH1F* Trk1 = SumHists(T1, x -> first + "Track1");
    TH1F* Trk2 = SumHists(T2, x -> first + "Track2");
    TH1F* Trk3 = SumHists(T3, x -> first + "Track3");
    TH1F* Trk4 = SumHists(T4, x -> first + "Track4");
    std::vector<TH1F*> Data = {Trk1, Trk2, Trk3, Trk4}; 
  
    bool kill = false; 
    for (TH1F* O : Data)
    {
      std::cout << O -> GetEntries() << std::endl;
      if (O -> GetEntries() < 1000){kill = true; }
    }
    if (kill == true)
    {
      BulkDelete(Trks); 
      continue;
    }
    MainAlgorithm(Data, Params, T1, 0); 
    break;
    MainAlgorithm(Data, Params, T2, 1); 
    MainAlgorithm(Data, Params, T3, 2); 
    MainAlgorithm(Data, Params, T3, 3);  

  }
}

void PlotInsideOutsideJet()
{
  TString Dir = "Merged_MC.root";
  std::map<TString, std::vector<TH1F*>> Output = Experimental_MC_Reader(Dir);
  
  std::vector<TH1F*> R_out = Output["All_Greater"]; 
  std::vector<TH1F*> R_in = Output["All_Less"]; 
  
  std::cout << R_out.size() << std::endl;



  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  PlotHists(R_out, can); 
  can -> Print("Out_jetcore.pdf");
  delete can; 

  can = new TCanvas(); 
  can -> SetLogy(); 
  PlotHists(R_in, can); 
  can -> Print("In_jetcore.pdf");

  delete can;
}

void Shifting(TH1F* H)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int LucyRichardson = 300; 
  int bins = H -> GetNbinsX(); 
  float min = H -> GetXaxis() -> GetXmin(); 
  float max = H -> GetXaxis() -> GetXmax(); 
  
  std::vector<TH1F*> PSF; 
  float m = 0; 
  float s = 0.03;
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


