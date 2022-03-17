#include "BaseFunctions.h"

VTH1F GetSubVector(std::vector<TH1F*> Input, int start, int end)
{
  if (Input.size() == 0){ return {}; }
  VTH1F Out = {Input.begin() + start, Input.end() - Input.size() + end}; 
  
  for (TH1F* H : Out)
  {
    for (int i = 0; i < H -> GetNbinsX(); i++)
    {
      float k = H -> GetBinContent(i+1);
      if (std::isnan(k)) { H -> SetBinContent(i+1, 1e-9); }
    }
  }
  return Out; 
}

VF ShapeComparison(VTH1F Hists, VTH1F Truth)
{
  VF Output; 
  int n = 0;
  if (Hists.size() < 4){ n = Hists.size(); }
  else {n = 4;}
  for (int i(0); i < n; i++)
  {
    TH1F* H = (TH1F*)Hists[i] -> Clone("H"); 
    TH1F* H_T = (TH1F*)Truth[i] -> Clone("T"); 
    
    float l_H = H -> Integral();
    H -> Scale(1./l_H); 

    float l_T = H_T -> Integral();
    H_T -> Scale(1./l_T); 

    float chi2 = ChiSquareDistanceSym(H, H_T);  
    Output.push_back(chi2); 
    
    delete H, H_T; 
  }
  return Output;
}

float ChiSquareDistanceSym(TH1F* H1, TH1F* H2)
{
  float chi2 = 0; 
  for (int i = 0; i < H1 -> GetNbinsX(); i++)
  {
    float x = H1 -> GetBinContent(i+1); 
    float y = H2 -> GetBinContent(i+1); 
    float sum = x + y; 
    float dif = x-y; 
    if (sum == 0){continue;}
    chi2 += std::pow(dif, 2)/sum; 
  }
  chi2 = 0.5*chi2; 
  return chi2;
}

float WeightedComparisonToTruth(VTH1F Hists, VTH1F Truth)
{
  VF Output = ShapeComparison(Hists, Truth);
  
  if (Output.size() == 0){ return -1; }
  // Get the sum of the data we are testing on 
  float sum = 0; 
  for (int i(0); i < Truth.size(); i++){sum += Truth[i] -> Integral();}
  
  // Check the impact the shape difference made to the overall data size
  float err = 0;  
  for (int i(0); i < Truth.size(); i++)
  {
    err += Output[i]*((Truth[i] -> Integral()) / sum);
  }
  return err;
}

float Flost2(std::vector<std::vector<TH1F*>> ntrk)
{
  std::vector<float> fi = {0, 0}; 
  for (int i(0); i < ntrk.size(); i++)
  {
    if (ntrk[i].size() > 1 && i < fi.size())
    {  
      fi[i] = ntrk[i][1] -> Integral(); 
    }
  }

  // Convention: n_trk_tru
  float n_1_2 = fi[0];
  float n_2_2 = fi[1];
  
  if (n_1_2 == 0 && n_2_2 == 0){ return -1.0; }

  float FL2 = float(n_1_2) / float(2.*(n_1_2 + n_2_2)); 
  return FL2;
}

float Flost3(std::vector<std::vector<TH1F*>> ntrk)
{
  std::vector<float> fi = {0, 0, 0, 0}; 
  for (int i(0); i < ntrk.size(); i++)
  {
    if (ntrk[i].size() > 2 && i < fi.size())
    {  
      fi[i] = ntrk[i][2] -> Integral(); 
    }
  }

  // Convention: n_trk_tru
  float n_1_3 = fi[0];
  float n_2_3 = fi[1];
  float n_3_3 = fi[2];
  float n_4_3 = fi[3]; 

  if (n_1_3 == 0 && n_2_3 == 0){ return -1.0; }

  float FL3 = float(2*n_1_3 + n_2_3) / float(3.*(n_1_3 + n_2_3 + n_3_3 + n_4_3)); 
  
  return FL3;
}

MVTH1F PrintFlost(MMMF FLostMap, TString Lay, TString dir)
{
  TString FL; 

  if (dir.Contains("FLost2")){ FL = "FLost2"; }
  if (dir.Contains("FLost3")){ FL = "FLost3"; }
  MVTH1F Output; 

  MMF Layer = FLostMap[Lay]; 
  std::vector<float> Errors; 
  std::vector<float> Perc; 
  std::vector<TString> Best_Alg; 
  std::vector<float> Flost; 

  std::vector<TString> out; 
  out.push_back(dir);
  TString tmp = "Algorithm"; 
  CoutText(&tmp, 12, " "); 
  tmp += " | "; 
  for (TString JE : JetEnergy)
  { 
    tmp += JE.ReplaceAll("_GeV", "").ReplaceAll("_", "-"); 
    tmp +=(" | "); 
    Errors.push_back(-2); 
    Perc.push_back(0); 
    Best_Alg.push_back("x"); 
    Flost.push_back(0); 
  }

  out.push_back(tmp); 
  
  for (MMFi alg = Layer.begin(); alg != Layer.end(); alg++)
  {
    tmp.Clear();  
    tmp += (alg -> first);
    CoutText(&tmp, 22 - tmp.Sizeof(), " "); 
    tmp += (" | "); 
    
    TString Title_E = "Error_" + FL + "_"; Title_E += (alg->first); Title_E += ("_" + Lay); 
    TH1F* Hist_E = new TH1F(Title_E, Title_E, JetEnergy.size(), 0, JetEnergy.size()); 

    TString Title = "Original_" + FL + "_"; Title += (alg->first); Title += ("_" + Lay); 
    TH1F* Hist = new TH1F(Title, Title, JetEnergy.size(), 0, JetEnergy.size());

    int x = 0; 
    for (TString JE : JetEnergy)
    {
      float FL = Layer[alg -> first][JE]; 
      TString fl = PrecisionString(FL, 4, false); 
      
      float Tru = Layer["Truth"][JE]; 
      TString l = JE.ReplaceAll("_GeV", "").ReplaceAll("_", "-");
      
      CoutText(&tmp, l.Sizeof() - fl.Sizeof(), " "); 
      tmp += fl; 
      tmp += (" | ");
      
      Hist -> SetBinContent(x+1, FL);
    
      if (Tru != 0){ Hist_E -> SetBinContent(x+1, float(std::abs(Tru - FL)/Tru)*100); }
      else{ Hist_E -> SetBinContent(x+1, 0); }
      
      if (Errors[x] == -2 || Errors[x] > std::abs(Tru - FL))
      {
        if (alg -> first == "Truth" || FL <= 0){ }
        else
        {
          Errors[x] = std::abs(Tru - FL); 
          Best_Alg[x] = alg -> first; 
          Perc[x] = (Errors[x] / Tru)*100; 
          Flost[x] = FL; 
        }
      }
      x++; 
    }

    Output["Original_" + FL].push_back(Hist); 
    Output["Error_" + FL].push_back(Hist_E); 

    out.push_back(tmp); 
    tmp.Clear(); 
  }

  CoutText(&tmp, 100, "_"); 
  out.push_back(tmp); 

  tmp.Clear();
  for (int a(0); a < Best_Alg.size(); a++)
  { 
    tmp += JetEnergy[a]; 
    tmp += (" -> "); 
    tmp += Best_Alg[a]; 
    tmp += (" | ");
    tmp += (PrecisionString(Flost[a], 3, false)); 
    tmp += (" +- "); 
    tmp += (PrecisionString(Perc[a], 2, false)); 
    tmp += ("% \n"); 
  }
  out.push_back(tmp); 
  
  std::ofstream print; 
  print.open(dir); 
  for (TString a : out)
  {
    print << a << "\n"; 
  }
  print.close(); 

  return Output; 
}

TString PrecisionString(float number, int precision, bool sci)
{
  TString out; 
  std::ostringstream p;
  p.precision(precision);
  if (sci){p << std::scientific << number;}
  else {p << std::fixed << number;}
  out += (p.str()); 
  return out;
}

void CoutText(TString *Input, int v, TString Text)
{
  for (int c(0); c < v; c++){ *Input += (Text); }
}

void BulkWrite(std::vector<TH1F*> H)
{
  for (TH1F* x : H){ x -> Write(); }  
}

MF ReadTH1F(TH1F* H)
{
  MF O; 
  for (int i(0); i < H -> GetNbinsX(); i++)
  {
    float e = H -> GetBinContent(i+1);
    TString je = JetEnergy[i]; 
    O[je] = e; 
  }
  return O; 
}

TH1F* SumHist(std::vector<TH1F*> H, TString Title)
{
  TH1F* Data = (TH1F*)H[0] -> Clone(Title); 
  Data -> SetTitle(Title); 
  
  for (int i(1); i < H.size(); i++){Data -> Add(H[i]);}

  return Data; 
}


