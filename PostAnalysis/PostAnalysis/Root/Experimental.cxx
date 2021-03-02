#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/AlgorithmFunctions.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(TH1F* InputTrk1, std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth, int trk_Data)
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

  std::vector<TH1F*> ntrk_Conv = ConvolveNTimes(InputTrk1, 4, "C"); 
  TH1F* Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy"); 
  std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, Data_Copy, Params); 
  Flush(F_C, ntrk_Conv, false); 

  TString name = Data[trk_Data] -> GetTitle(); name += ("_Experimental.pdf"); 
  PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
  can -> Print(name); 

  for (int i(0); i < iterations; i++)
  {
    
    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> trk; 
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p == trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else{ trk.push_back(ntrk_Conv[p]); }
    }
    Average(Data_Copy); 
    Smear(Data_Copy, 0.05); 
    F_C = LoopGen(trk, Data_Copy, Params); 
    Flush(F_C, trk, true);
    
    PlotHists(Data_Copy, Truth, trk, can); 
    can -> Print(name); 
    delete Data_Copy; 


    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> not_trk; 
    for (int p(0); p < ntrk_Conv.size(); p++)
    {
      if (p != trk_Data){ Data_Copy -> Add(ntrk_Conv[p], -1); }
      else{ not_trk.push_back(ntrk_Conv[p]); }
    }
    F_C = LoopGen(not_trk, Data_Copy, Params); 
    ScaleShape(Data_Copy, not_trk);  
    Flush(F_C, not_trk, true); 

    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print(name); 
    delete Data_Copy;    

    Data_Copy = (TH1F*)Data[trk_Data] -> Clone("Data_Copy");    
    std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, Data_Copy, Params); 
    ScaleShape(Data_Copy, F_C);  
    Flush(F_C, ntrk_Conv, true); 

    PlotHists(Data_Copy, Truth, ntrk_Conv, can); 
    can -> Print(name); 
    delete Data_Copy;    
  
  }
  BulkDelete(ntrk_Conv);

  std::map<TString, std::vector<TH1F*>> Out; 
  return Out; 
}

void AlgorithmMonteCarlo()
{
  int nsdo = 4;
  float m = 0.0001; 
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {0.04, 0.04, 0.04, 0.04};  
  Params["iterations"] = {2}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.02, 0.02, 0.02, 0.02}; 
  Params["cache"] = {10000}; 
  Params["x_range"] = {0., 10.}; 

  TString Dir = "Merged_MC.root"; 
 
  std::map<TString, std::vector<TH1F*>> Output = MC_Reader(Dir);

  typedef std::map<TString, std::vector<TH1F*>>::iterator M;
  for (M x = Output.begin(); x != Output.end(); x++)
  {
    TString name = x -> first; 
    std::vector<TH1F*> Trks = x -> second; 
   
    if (name.Contains("_Less") || name.Contains("Greater") || name.Contains("_Data")){continue;}
    std::vector<TH1F*> Outside_Core = Output[name + "_radius_Greater"]; 
    
    std::vector<TH1F*> T1;  
    std::vector<TH1F*> T2;  
    std::vector<TH1F*> T3;  
    std::vector<TH1F*> T4;  
    for (int i(0); i < nsdo; i++){T1.push_back(Trks[i]);}
    for (int i(4); i < nsdo+4; i++){T2.push_back(Trks[i]);}
    for (int i(8); i < nsdo+8; i++){T3.push_back(Trks[i]);}
    for (int i(12); i < nsdo+12; i++){T4.push_back(Trks[i]);} 
    
    TH1F* Outside_JetCore = SumHists(Outside_Core, name + "_Outside_JetCore");  
    TH1F* Trk1 = Output[name + "_Data"][0]; 
    TH1F* Trk2 = Output[name + "_Data"][1]; 
    TH1F* Trk3 = Output[name + "_Data"][2]; 
    TH1F* Trk4 = Output[name + "_Data"][3]; 
    std::vector<TH1F*> Data = {Trk1, Trk2, Trk3, Trk4}; 
  
    std::cout << " :::: " << name << std::endl;  
    bool kill = false; 
    for (TH1F* O : Data)
    {
      std::cout << O -> GetEntries() << std::endl;
      if (O -> GetEntries() < 1000){kill = true; }
    }
    if (kill == true)
    {
      BulkDelete(Trks); 
      BulkDelete(Outside_Core); 
      delete Outside_JetCore;
      continue;
    }
    
    //Shifting(Outside_JetCore); 
    MainAlgorithm(Outside_JetCore, Data, Params, T1, 0); 
    MainAlgorithm(Outside_JetCore, Data, Params, T2, 1); 
    MainAlgorithm(Outside_JetCore, Data, Params, T3, 2); 
    MainAlgorithm(Outside_JetCore, Data, Params, T3, 3);  
  }
}

void PlotInsideOutsideJet()
{
  TString Dir = "Merged_MC.root";
  std::map<TString, std::vector<TH1F*>> Output = MC_Reader(Dir);
  
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

  int LucyRichardson = 20; 
  int bins = H -> GetNbinsX(); 
  float min = H -> GetXaxis() -> GetXmin(); 
  float max = H -> GetXaxis() -> GetXmax(); 
  
  std::vector<TH1F*> PSF; 
  float m = 0; 
  float s = 0.01;
  float w = (max - min)/float(bins); 
  
  float new_min = min - w*120; 
  float new_max = max + w*120; 
  
  TH1F* PDF = new TH1F("PDF", "PDF", bins+2*120, new_min, new_max); 
  TH1F* Gaus = Gaussian(m ,s, bins+2*120, new_min, new_max, "Gaus"); 
  
  for (int i(0); i < bins; i++)
  {
    float e = H -> GetBinContent(i+1); 
    PDF -> SetBinContent(120+i+1, e); 
  }
  
  TH1F* Dec = (TH1F*)PDF -> Clone("DecP"); 
  TH1F* Temp = (TH1F*)PDF -> Clone("Temp"); 

  Gaus -> Draw("HIST"); 
  can -> Print("Debug.pdf"); 

  for (int i(0); i < 50; i++)
  {
    DeconvolutionExperimental(Temp, Gaus, Dec, LucyRichardson); 
    Convolution(Dec, Gaus, Temp); 
    
    //Convolution(Dec, Gaus, Temp); 
    Normalize(Temp); 
    Normalize(PDF); 
    PDF -> SetLineStyle(kSolid); 
    PDF -> Draw("HIST");
    Temp -> SetLineColor(kRed);
    Temp -> SetLineStyle(kDashed); 
    Temp -> Draw("SAMEHIST"); 
    can -> Print("Debug.pdf"); 
  }

}


