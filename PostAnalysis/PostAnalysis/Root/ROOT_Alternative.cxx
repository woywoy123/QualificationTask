#include<PostAnalysis/ROOT_Alternative.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/AlgorithmFunctions.h>

#include<TH1F.h>
#include<TH1.h>
#include<TF1.h>
#include<TF1Convolution.h>
#include<TRandom.h>
#include<TGraph.h>
#include<TMath.h>


void FitWithRoot(std::vector<TH1F*> Data, TH1F* trk1_start, std::vector<std::vector<TH1F*>> Truth)
{
 

  std::vector<std::vector<TH1F*>> ntrk_mtru_H = BuildNtrkMtru(1, trk1_start, "Test"); 
  
  std::vector<TH1F*> trk1_tpts = ntrk_mtru_H[0]; 
  TH1F* trk1 = Data[0]; 
  std::vector<TH1F*> trk1_truth = Truth[0];  

  Normalize(trk1_start); 
  
  
  int bins = trk1 -> GetNbinsX(); 
  float x_min = trk1 -> GetXaxis() -> GetXmin(); 
  float x_max = trk1 -> GetXaxis() -> GetXmax(); 



  TGraph* trk1_gr = new TGraph(); 
  for (int i(0); i < bins; i++){trk1_gr -> SetPoint(i, trk1_start -> GetBinCenter(i+1), trk1_start -> GetBinContent(i+1));}
  
  TF1* trk1_F1 = new TF1("trk1_F1", [&](double*x, double*){ return trk1_gr -> Eval(x[0]); }, x_min, x_max, 1); 

  
  TCanvas* can = new TCanvas(); 
  //can -> SetLogy(); 
  can -> Print("Test.pdf["); 
  trk1_F1 -> Draw();
  can -> Print("Test.pdf"); 
  can -> Clear(); 
 
  TF1* Gaus = new TF1("G_TX", "TMath::Gaus(x, [0], [1])", x_min, x_max); 
  Gaus -> SetParameters(1, 0.1);  
  //Gaus -> SetParameter(0, 1);
  //Gaus -> SetParameter(1, 1); 
  //Gaus -> SetParameter(2, 10); 
  //Gaus -> SetParName(0, "Mean"); 
  //Gaus -> SetParName(1, "sdev"); 
  //Gaus -> SetParName(2, "Norm"); 
  
  TF1Convolution* C_Gaus = new TF1Convolution(Gaus, Gaus, x_min, x_max, true); 
  C_Gaus -> SetNofPointsFFT(10000); 
  
  TF1* GxT = new TF1("GxT", *C_Gaus, x_min, x_max, C_Gaus -> GetNpar());
  //trk1 -> Fit("GxT");  
  
  GxT -> Draw();
  //Gaus -> Draw();  
  //C_Gaus -> Draw(); 
  can -> Print("Test.pdf"); 

  can -> Print("Test.pdf]"); 
  


}
