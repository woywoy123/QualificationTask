#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/AlgorithmFunctions.h>

std::map<TString, std::vector<TH1F*>> MainAlgorithm(TH1F* InputTrk1, std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, std::vector<TH1F*> Truth, int trk_Data)
{
  std::vector<TH1F*> T1;  
  std::vector<TH1F*> T2;  
  std::vector<TH1F*> T3;  
  std::vector<TH1F*> T4;  
  for (int i(0); i < trk_Data; i++){T1.push_back(Truth[i]);}
  for (int i(4); i < trk_Data+4; i++){T2.push_back(Truth[i]);}
  for (int i(8); i < trk_Data+8; i++){T3.push_back(Truth[i]);}
  for (int i(12); i < trk_Data+12; i++){T4.push_back(Truth[i]);} 

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int iterations = Params["iterations"][0]; 
  std::vector<TH1F*> ntrk_Conv1 = ConvolveNTimes(InputTrk1, 4, "C1"); 
  std::vector<TH1F*> ntrk_Conv2 = ConvolveNTimes(InputTrk1, 4, "C2"); 
  std::vector<TH1F*> ntrk_Conv3 = ConvolveNTimes(InputTrk1, 4, "C3"); 
  std::vector<TH1F*> ntrk_Conv4 = ConvolveNTimes(InputTrk1, 4, "C4"); 
  TString name = Data[0] -> GetTitle(); name += ("_Experimental.pdf"); 

  for (int i(0); i < iterations; i++)
  {
    auto Clone =[&] (TH1F* H, TH1F* O)
    {
      float HL = H -> Integral(); 
      float OL = O -> Integral(); 
    
      H -> Scale(1. / HL); 
      O -> Scale(1. / OL); 
      for (int x(0); x < H -> GetNbinsX(); x++)
      {
        float e = H -> GetBinContent(x+1); 
        float f = O -> GetBinContent(x+1); 
        O -> SetBinContent(x+1, 0.5*f+0.5*e);
      }
      float T = O -> Integral(); 
      O -> Scale(1. / T); 
      H -> Scale(HL); 
      O -> Scale(OL); 
    };

    auto Subtract =[&] (std::vector<TH1F*> ntrk_Conv, TH1F* Data, int trk)
    {
      for (int i(0); i < ntrk_Conv.size(); i++){if (i != trk-1){ Data -> Add(ntrk_Conv[i], -1); }}
      Average(Data); 
      Flush({Data}, {ntrk_Conv[trk-1]}); 
    }; 

    auto FitTrack =[&] (std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params, int trk)
    {
      std::vector<TH1F*> F_C = LoopGen(ntrk_Conv, Data, Params); 
      ScaleShape(Data, ntrk_Conv); 
      Subtract(F_C, Data, trk);
      Flush(F_C, ntrk_Conv); 
    }; 


    TH1F* trk1 = (TH1F*)Data[0] -> Clone("Trk1_"); 
    TH1F* trk2 = (TH1F*)Data[1] -> Clone("Trk2_"); 
    TH1F* trk3 = (TH1F*)Data[2] -> Clone("Trk3_"); 
    TH1F* trk4 = (TH1F*)Data[3] -> Clone("Trk4_"); 

    //ScalingFit(trk1, ntrk_Conv1); 
    //ScalingFit(trk2, ntrk_Conv2); 
    //ScalingFit(trk3, ntrk_Conv3); 
    //ScalingFit(trk4, ntrk_Conv4); 
    
    //Subtract(ntrk_Conv1, trk1, 1); 
    //Subtract(ntrk_Conv2, trk2, 2); 
    //Subtract(ntrk_Conv3, trk3, 3); 
    //Subtract(ntrk_Conv4, trk4, 4); 

    //trk1 = (TH1F*)Data[0] -> Clone("Trk1_"); 
    //trk2 = (TH1F*)Data[1] -> Clone("Trk2_"); 
    //trk3 = (TH1F*)Data[2] -> Clone("Trk3_"); 
    //trk4 = (TH1F*)Data[3] -> Clone("Trk4_"); 

    //std::vector<TH1F*> ntrk_Conv_temp1 = ConvolveNTimes(ntrk_Conv1[0], 4, "Temp1"); 
    //std::vector<TH1F*> ntrk_Conv_temp2 = ConvolveNTimes(ntrk_Conv1[0], 4, "Temp2"); 
    //std::vector<TH1F*> ntrk_Conv_temp3 = ConvolveNTimes(ntrk_Conv1[0], 4, "Temp3"); 
    //std::vector<TH1F*> ntrk_Conv_temp4 = ConvolveNTimes(ntrk_Conv1[0], 4, "Temp4"); 
   
    //for (int x(0); x < ntrk_Conv_temp1.size(); x++)
    //{
    //  Clone(ntrk_Conv_temp1[x], ntrk_Conv1[x]); 
    //  Clone(ntrk_Conv_temp2[x], ntrk_Conv2[x]); 
    //  Clone(ntrk_Conv_temp3[x], ntrk_Conv3[x]); 
    //  Clone(ntrk_Conv_temp4[x], ntrk_Conv4[x]); 
    //} 
    //
    //BulkDelete(ntrk_Conv_temp1); 
    //BulkDelete(ntrk_Conv_temp2); 
    //BulkDelete(ntrk_Conv_temp3); 
    //BulkDelete(ntrk_Conv_temp4); 
    
    ScalingShift(trk1, ntrk_Conv1); 
    ScalingShift(trk2, ntrk_Conv2); 
    ScalingShift(trk3, ntrk_Conv3); 
    ScalingShift(trk4, ntrk_Conv4); 

    Subtract(ntrk_Conv1, trk1, 1); 
    Subtract(ntrk_Conv2, trk2, 2); 
    Subtract(ntrk_Conv3, trk3, 3); 
    Subtract(ntrk_Conv4, trk4, 4);

    if (i < 1)
    {
      Clone(ntrk_Conv2[1], ntrk_Conv1[1]); 
      Clone(ntrk_Conv3[2], ntrk_Conv1[2]); 
      Clone(ntrk_Conv4[3], ntrk_Conv1[3]); 

      Clone(ntrk_Conv1[0], ntrk_Conv2[0]); 
      Clone(ntrk_Conv3[2], ntrk_Conv2[2]); 
      Clone(ntrk_Conv4[3], ntrk_Conv2[3]); 

      Clone(ntrk_Conv1[0], ntrk_Conv3[0]); 
      Clone(ntrk_Conv2[1], ntrk_Conv3[1]); 
      Clone(ntrk_Conv4[3], ntrk_Conv3[3]); 

      Clone(ntrk_Conv1[0], ntrk_Conv4[0]); 
      Clone(ntrk_Conv2[1], ntrk_Conv4[1]); 
      Clone(ntrk_Conv3[2], ntrk_Conv4[2]); 
    }


    //trk1 = (TH1F*)Data[0] -> Clone("Trk1_"); 
    //trk2 = (TH1F*)Data[1] -> Clone("Trk2_"); 
    //trk3 = (TH1F*)Data[2] -> Clone("Trk3_"); 
    //trk4 = (TH1F*)Data[3] -> Clone("Trk4_");     


    //FitTrack(ntrk_Conv1, trk1, Params, 1); 
    //FitTrack(ntrk_Conv2, trk2, Params, 2); 
    //FitTrack(ntrk_Conv3, trk3, Params, 3); 
    //FitTrack(ntrk_Conv4, trk4, Params, 4); 

    trk1 = (TH1F*)Data[0] -> Clone("Trk1_"); 
    trk2 = (TH1F*)Data[1] -> Clone("Trk2_"); 
    trk3 = (TH1F*)Data[2] -> Clone("Trk3_"); 
    trk4 = (TH1F*)Data[3] -> Clone("Trk4_");     

    ScalingFit(trk1, ntrk_Conv1); 
    ScalingFit(trk2, ntrk_Conv2); 
    ScalingFit(trk3, ntrk_Conv3); 
    ScalingFit(trk4, ntrk_Conv4);   

    ScaleShape(trk1, ntrk_Conv1);
    ScaleShape(trk2, ntrk_Conv2);
    ScaleShape(trk3, ntrk_Conv3);
    ScaleShape(trk4, ntrk_Conv4);

    delete trk1; 
    delete trk2; 
    delete trk3; 
    delete trk4; 

    can -> Print(name + "["); 
    PlotHists(T1, ntrk_Conv1, can);
    can -> Print(name); 
    PlotHists(T2, ntrk_Conv2, can);
    can -> Print(name); 
    PlotHists(T3, ntrk_Conv3, can);
    can -> Print(name); 
    PlotHists(T4, ntrk_Conv4, can);
    can -> Print(name); 
    can -> Print(name + "]"); 

  }
   
  delete can; 
  std::map<TString, std::vector<TH1F*>> Out; 
  BulkDelete(ntrk_Conv1); 
  BulkDelete(ntrk_Conv2); 
  BulkDelete(ntrk_Conv3); 
  BulkDelete(ntrk_Conv4);  
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
  Params["s_e"] = {0.2, 0.2, 0.2, 0.2};  
  Params["iterations"] = {1}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.05, 0.05, 0.05, 0.05}; 
  Params["cache"] = {10000}; 
  Params["x_range"] = {0.1, 9.0}; 

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
    //std::vector<TH1F*> T2;  
    //std::vector<TH1F*> T3;  
    //std::vector<TH1F*> T4;  
    for (int i(0); i < nsdo+12; i++){T1.push_back(Trks[i]);}
    //for (int i(4); i < nsdo+4; i++){T2.push_back(Trks[i]);}
    //for (int i(8); i < nsdo+8; i++){T3.push_back(Trks[i]);}
    //for (int i(12); i < nsdo+12; i++){T4.push_back(Trks[i]);} 
    
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
    MainAlgorithm(Outside_JetCore, Data, Params, T1, nsdo); 
    //MainAlgorithm(Outside_JetCore, Data, Params, T2, 1); 
    //MainAlgorithm(Outside_JetCore, Data, Params, T3, 2); 
    //MainAlgorithm(Outside_JetCore, Data, Params, T3, 3);  
    //break;
  }
}

void Shifting(TH1F* H)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  float m = 0.5;
  std::map<TString, std::vector<float>> Params; 
  Params["m_s"] = {-m, -m, -m, -m}; 
  Params["m_e"] = {m, m, m, m}; 
  Params["s_s"] = {0.005, 0.005, 0.005, 0.005};
  Params["s_e"] = {0.1, 0.1, 0.1, 0.1};  
  Params["iterations"] = {15}; 
  Params["LR_iterations"] = {50}; 
  Params["G_Mean"] = {0, 0, 0, 0}; 
  Params["G_Stdev"] = {0.025, 0.025, 0.025, 0.025}; 
  Params["cache"] = {10000}; 
  Params["x_range"] = {0.5, 9.0}; 


  int LucyRichardson = 20; 
  int bins = H -> GetNbinsX(); 
  float min = H -> GetXaxis() -> GetXmin(); 
  float max = H -> GetXaxis() -> GetXmax(); 
  
  std::vector<TH1F*> PSF; 
  //float m = 0; 
  float s = 0.1;
  float w = (max - min)/float(bins); 
 
  int r = 0; 
  float new_min = min - w*r; 
  float new_max = max + w*r; 
  
  TH1F* PDF = new TH1F("PDF", "PDF", bins+2*r, new_min, new_max); 
  //TH1F* Gaus = Gaussian(m ,s, bins+2*r, new_min, new_max, "Gaus"); 
  
  for (int i(0); i < bins; i++)
  {
    float e = H -> GetBinContent(i+1); 
    PDF -> SetBinContent(r+i+1, e); 
  }
  
  TH1F* Dec = (TH1F*)PDF -> Clone("DecP"); 
  TH1F* Temp = (TH1F*)PDF -> Clone("Temp"); 
  TH1F* Residual = (TH1F*)PDF -> Clone("Residual"); 

  //Gaus -> Draw("HIST"); 
  can -> Print("Debug.pdf"); 

  for (int i(0); i < 50; i++)
  {
    std::vector<TH1F*> t = LoopGen({Temp}, PDF, Params); 
    Temp -> Reset(); 
    Temp -> Add(t[0]); 
    delete t[0]; 

    //DeconvolutionExperimental(Temp, Gaus, Dec, LucyRichardson); 
    //Convolution(Dec, Gaus, Temp); 

    Normalize(Temp); 
    Normalize(PDF); 

    Residual -> Reset(); 
    Residual -> Add(PDF); 
    Residual -> Add(Temp, -1); 
    
    PDF -> SetLineStyle(kSolid); 
    //PDF -> Draw("HIST");
    Temp -> SetLineColor(kRed);
    Temp -> SetLineStyle(kDashed); 
    //Temp -> Draw("SAMEHIST"); 
    
    Residual -> Draw("HIST"); 
    
    
    can -> Print("Debug.pdf"); 
  }

}


