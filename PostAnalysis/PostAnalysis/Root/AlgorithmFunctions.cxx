#include<PostAnalysis/AlgorithmFunctions.h>

void ScaleShape(TH1F* Data, std::vector<TH1F*> ntrk)
{
  for (int i(0); i < Data ->GetNbinsX(); i++)
  {
    float e = Data -> GetBinContent(i+1); 
    float sum(0);  
    for (int c(0); c < ntrk.size(); c++)
    {
      sum += ntrk[c] -> GetBinContent(i+1);  
    }
    for (int z(0); z < ntrk.size(); z++)
    {
      float point = ntrk[z] -> GetBinContent(i+1); 
      float ed = point*(e/sum);  
      float log = ed;  
      if (std::isnan(log)){log = 0;}
      if (std::isinf(log)){log = 0;} 
      ntrk[z] -> SetBinContent(i+1, log); 
    }     
  }
}

float FitError(TH1F* Data, std::vector<TH1F*> ntrk)
{
  int bins = Data -> GetNbinsX(); 
  float sum_D = 0; 
  float sum_E = 0;  
  for (int i(0); i < bins; i++)
  {
    float e = Data -> GetBinContent(i+1); 
    float sum_H = 0; 
    for (int p(0); p < ntrk.size(); p++){sum_H += std::abs(ntrk[p] -> GetBinContent(i+1));}
    sum_E += std::abs(sum_H - e);  
  }
  sum_E = sum_E/Data -> Integral();
  return sum_E;
}



void Flush(std::vector<TH1F*> F_C, std::vector<TH1F*> ntrk_Conv, bool sig)
{
  if (F_C.size() != ntrk_Conv.size()) {std::cout << "!!!!!!!!!!!!!!!!!" << std::endl; }
  for (int i(0); i < F_C.size(); i++)
  {
    float HL = F_C[i] -> Integral(); 
    float OL = ntrk_Conv[i] -> Integral(); 
    
    F_C[i] -> Scale(1. / HL); 
    ntrk_Conv[i] -> Scale(1. / OL); 
    for (int x(0); x < F_C[i] -> GetNbinsX(); x++)
    {
      float e = F_C[i] -> GetBinContent(x+1); 
      float f = ntrk_Conv[i] -> GetBinContent(x+1); 
      ntrk_Conv[i] -> SetBinContent(x+1, (e+f)/2.);
    }
    
    F_C[i] -> Scale(HL); 
    ntrk_Conv[i] -> Scale(HL); 
    delete F_C[i];  
  }
}

void Average(TH1F* Data)
{
  std::vector<float> D; 
  for (int i(0); i < Data -> GetNbinsX(); i++)
  {
    float y = Data -> GetBinContent(i+1); 
    float sum = 0; 
    int p = 0; 
    for (int x(0); x < 1; x++)
    {
      float e = Data -> GetBinContent(i + x + 1); 
      if (e < 0){ e = 0; }
      sum += e;
      p++; 
    }
    if ( sum < 1e-10 ){sum = 1e-10;}
    sum = sum / float(p); 
    D.push_back(sum);  
  }
  for (int i(0); i < Data -> GetNbinsX(); i++)
  {
    Data -> SetBinContent(i+1, D[i]); 
  }
}


std::vector<TH1F*> LoopGen(std::vector<TH1F*> ntrk_Conv, TH1F* Data, std::map<TString, std::vector<float>> Params)
{
  std::vector<TString> Names_Dec; 
  float r = 0.2; 
  int bins = Data -> GetNbinsX(); 
  float min = Data -> GetXaxis() -> GetXmin(); 
  float max = Data -> GetXaxis() -> GetXmax(); 
  float w = (max - min) / float(bins); 
  float new_min = min - w*bins*r; 
  float new_max = max + w*bins*r; 

  std::vector<TH1F*> PDF_D;
  std::vector<TH1F*> PSF;  
  std::vector<TH1F*> PDF; 
  for (int i(0); i < ntrk_Conv.size(); i++)
  {
    TString nameG = "Gx_"; nameG+=(i+1); 
    TH1F* Gaus = Gaussian(Params["G_Mean"][i], Params["G_Stdev"][i], bins+2*bins*r, new_min, new_max, nameG); 
    PSF.push_back(Gaus);  

    TString name = "Dec_"; name += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* H = new TH1F(name, name, bins+2*bins*r, new_min, new_max); 
    PDF_D.push_back(H); 
  
    TString name_L = "L_"; name_L += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* X = new TH1F(name_L, name_L, bins+2*bins*r, new_min, new_max); 
    PDF.push_back(X); 

    for (int j(0); j < bins; j++)
    {
      X -> SetBinContent(j+1+r*bins, ntrk_Conv[i] -> GetBinContent(j+1)); 
    }
  }
  Normalize(PDF); 
  Normalize(PSF); 
  MultiThreadingDeconvolutionExperimental(PDF, PSF, PDF_D, Params["LR_iterations"][0]); 

  TH1F* Data_D = new TH1F("DATA_D", "DATA_D", bins+2*bins*r, new_min, new_max); 
  for (int j(0); j < bins; j++)
  {
    Data_D -> SetBinContent(j+1+r*bins, Data -> GetBinContent(j+1)); 
  }

  std::vector<std::pair<TH1F*, std::vector<float>>> trk_Fit = FitDeconvolutionPerformance(Data_D, PDF_D, Params, Params["cache"][0], Params["cache"][0]);
  delete Data_D; 
  std::vector<TH1F*> F_C;
  float bin_min = PDF[0] -> GetXaxis() -> FindBin(min) -1; 
  for (int i(0); i < trk_Fit.size(); i++)
  {
    TString name_N = "N_"; name_N += (ntrk_Conv[i] -> GetTitle()); 
    TH1F* H = new TH1F(name_N, name_N, bins, min, max); 
    F_C.push_back(H);  
   
    for (int j(0); j < bins; j++)
    {
      H -> SetBinContent(j+1, trk_Fit[i].first -> GetBinContent(j+1 + bin_min)); 
    }
    delete trk_Fit[i].first;
  } 
  
  BulkDelete(PDF_D); 
  BulkDelete(PSF); 
  BulkDelete(PDF); 
  return F_C; 
}


std::map<TString, std::vector<TH1F*>> MainAlgorithm(TH1F* InputTrk1, std::vector<TH1F*> Data, std::map<TString, std::vector<float>> Params, int trk_Data)
{
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 

  int iterations = Params["iterations"][0]; 
  std::vector<TH1F*> ntrk_Conv1 = ConvolveNTimes(InputTrk1, trk_Data, "C1"); 
  std::vector<TH1F*> ntrk_Conv2 = ConvolveNTimes(InputTrk1, trk_Data, "C2"); 
  std::vector<TH1F*> ntrk_Conv3 = ConvolveNTimes(InputTrk1, trk_Data, "C3"); 
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
        O -> SetBinContent(x+1, (e+f)/2.);
      }
      
      H -> Scale(HL); 
      O -> Scale(OL); 
    };


    Clone(ntrk_Conv2[1], ntrk_Conv1[1]); 
    Clone(ntrk_Conv3[2], ntrk_Conv1[2]); 

    Clone(ntrk_Conv1[0], ntrk_Conv2[0]); 
    Clone(ntrk_Conv3[2], ntrk_Conv2[2]); 

    Clone(ntrk_Conv1[0], ntrk_Conv3[0]); 
    Clone(ntrk_Conv2[1], ntrk_Conv3[1]); 

    for (int c(0); c < ntrk_Conv1.size(); c++)
    {
      Average(ntrk_Conv1[c]); 
      Average(ntrk_Conv2[c]); 
      Average(ntrk_Conv3[c]); 
    }
 
    TH1F* trk1_D = (TH1F*)Data[0] -> Clone("Trk1_temp"); 
    ScalingFit(trk1_D, ntrk_Conv1); 
    trk1_D -> Add(ntrk_Conv1[1], -1); 
    trk1_D -> Add(ntrk_Conv1[2], -1); 
    std::vector<TH1F*> ntrk_Conv_temp = ConvolveNTimes(trk1_D, 3, "Temp1"); 
    delete trk1_D; 

    Clone(ntrk_Conv_temp[0], ntrk_Conv1[0]); 
    Clone(ntrk_Conv_temp[1], ntrk_Conv2[1]); 
    Clone(ntrk_Conv_temp[2], ntrk_Conv3[2]); 
    BulkDelete(ntrk_Conv_temp); 
    
    Normalize(ntrk_Conv1); 
    Normalize(ntrk_Conv2);
    Normalize(ntrk_Conv3); 

    TH1F* trk1 = (TH1F*)Data[0] -> Clone("Trk1_D"); 
    TH1F* trk2 = (TH1F*)Data[1] -> Clone("Trk2_D"); 
    TH1F* trk3 = (TH1F*)Data[2] -> Clone("Trk3_D"); 
    
    std::vector<TH1F*> F_C1 = LoopGen(ntrk_Conv1, trk1, Params); 
    ScaleShape(trk1, F_C1); 
    trk1 -> Add(F_C1[1], -1);  
    trk1 -> Add(F_C1[2], -1); 
    Average(trk1); 
    Flush(F_C1, ntrk_Conv1, true); 
    Flush({trk1}, {ntrk_Conv1[0]}); 
    
    std::vector<TH1F*> F_C2 = LoopGen(ntrk_Conv2, trk2, Params); 
    ScaleShape(trk2, F_C2); 
    trk2 -> Add(F_C2[0], -1);  
    trk2 -> Add(F_C2[2], -1); 
    Average(trk2); 
    Flush(F_C2, ntrk_Conv2, true); 
    Flush({trk2}, {ntrk_Conv2[1]}); 
    
    std::vector<TH1F*> F_C3 = LoopGen(ntrk_Conv3, trk3, Params); 
    ScaleShape(trk3, F_C3); 
    trk3 -> Add(F_C3[0], -1);  
    trk3 -> Add(F_C3[1], -1); 
    Average(trk3); 
    Flush(F_C3, ntrk_Conv3, true); 
    Flush({trk3}, {ntrk_Conv3[2]}); 

    Clone(ntrk_Conv2[1], ntrk_Conv1[1]); 
    Clone(ntrk_Conv3[2], ntrk_Conv1[2]); 

    Clone(ntrk_Conv1[0], ntrk_Conv2[0]); 
    Clone(ntrk_Conv3[2], ntrk_Conv2[2]); 

    Clone(ntrk_Conv1[0], ntrk_Conv3[0]); 
    Clone(ntrk_Conv2[1], ntrk_Conv3[1]); 

    for (int c(0); c < ntrk_Conv1.size(); c++)
    {
      Average(ntrk_Conv1[c]); 
      Average(ntrk_Conv2[c]); 
      Average(ntrk_Conv3[c]); 
    }
  
    Normalize(ntrk_Conv1); 
    Normalize(ntrk_Conv2);
    Normalize(ntrk_Conv3);  
 
    trk1 = (TH1F*)Data[0] -> Clone("Trk1_D"); 
    trk2 = (TH1F*)Data[1] -> Clone("Trk2_D"); 
    trk3 = (TH1F*)Data[2] -> Clone("Trk3_D");  
    ScalingFit(trk1, ntrk_Conv1); 
    ScalingFit(trk2, ntrk_Conv2); 
    ScalingFit(trk3, ntrk_Conv3); 
 
    can -> Print("t.pdf["); 
    PlotHists(trk1, ntrk_Conv1, can);
    can -> Print("t.pdf"); 
    PlotHists(trk2, ntrk_Conv2, can);
    can -> Print("t.pdf"); 
    PlotHists(trk3, ntrk_Conv3, can);
    can -> Print("t.pdf"); 
    can -> Print("t.pdf]"); 

    delete trk1; 
    delete trk2; 
    delete trk3; 
  }

  std::map<TString, std::vector<TH1F*>> Out; 
  std::vector<std::vector<TH1F*>> All = {ntrk_Conv1, ntrk_Conv2, ntrk_Conv3}; 
  for (int d(0); d < trk_Data; d++)
  {
    std::vector<TH1F*> ntrk_Conv = All[d]; 
    for (int x(0); x < ntrk_Conv.size(); x++)
    {
      TString name2 = Data[d]->GetTitle();   
      name2.Remove(name2.Length() - 5, name2.Length()); 
      name2.Remove(0, 12); 
      name = "dEdx_ntrk_"; name += (d+1); name += ("_ntru_"); name += (x+1); name += ("_"); 
      name2 += ("_R");
      name = name + name2;
      TH1F* H = (TH1F*)ntrk_Conv[x] -> Clone(name); 
      H -> SetTitle(name); 
       
      TString OutName = "Result"; OutName += (d+1); 
      Out[OutName].push_back(H); 
    }
  }

  Out["Data"] = Data;
  
  for (int i(0); i < All.size(); i++){BulkDelete(All[i]);}
  delete can;

  return Out; 
}


