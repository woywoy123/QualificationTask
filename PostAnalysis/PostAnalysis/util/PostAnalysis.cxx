// Including commonly reused functions
#include<PostAnalysis/Functions.h>
#include<PostAnalysis/Verification.h>
#include<PostAnalysis/UnitClosures.h>
#include<TF1.h>

// Including standard C++ libraries 
#include<iostream>

using namespace RooFit;

void PostAnalysis()
{
  auto FillHist = [](TH1F* Hist, std::vector<float> Comps, std::vector<float> Params)
  {
    TF1 Lan("Lan", "landau", 0, 20);
    for (int i(0); i < Params.size(); i++)
    {
      Lan.SetParameter(i, Params.at(i));
    }
  
    for ( int i(0); i < 500000; i++)
    {
      double r1 = Lan.GetRandom();
      double r2 = Lan.GetRandom();
      double r3 = Lan.GetRandom();
      double r4 = Lan.GetRandom();

      Hist -> Fill(r1, Comps[0]);
      Hist -> Fill(r1+r2, Comps[1]);
      Hist -> Fill(r1 + r2 + r3, Comps[2]);
      Hist -> Fill(r1 + r2 + r3 + r4, Comps[3]);
    }
  };

  auto LandauGenerator = [](std::vector<TH1F*> Hists, std::vector<float> Params)
  {
    TF1 Lan("Lan", "landau", 0, 20);
    for (int i(0); i < Params.size(); i++)
    {
      Lan.SetParameter(i, Params.at(i));
    }
  
    for ( int i(0); i < 500000; i++)
    {
      double r1 = Lan.GetRandom();
      double r2 = Lan.GetRandom();
      double r3 = Lan.GetRandom();
      double r4 = Lan.GetRandom();
      Hists.at(0) -> Fill(r1);  
      Hists.at(1) -> Fill(r1 + r2); 
      Hists.at(2) -> Fill(r1 + r2 + r3); 
      Hists.at(3) -> Fill(r1 + r2 + r3 + r4); 
    }
  };

  Verification V; 
  Fit_Functions f;
  Functions F;
  UnitClosures U;
  Presentation P;

  // ====================== Generate the data =================== //  
  // Histograms: Pure with no cross contamination
  std::vector<TString> Hist_Names = {"trk1_P", "trk2_P", "trk3_P", "trk4_P"};
  std::vector<TH1F*> Hists = F.MakeTH1F(Hist_Names, 500, 0, 20);
  LandauGenerator(Hists, {1, 0.9, 0.1}); 
 
  // Create some Datasets which have some contamination
  // === COMPN is the composition of cross contamination
  std::vector<float> COMP1 = {1.  , 0. , 0.  , 0.  };  
  std::vector<float> COMP2 = {0.  , 1. , 0.  , 0.  };  
  std::vector<float> COMP3 = {0.01, 0.2, 0.59, 0.2 };  
  std::vector<float> COMP4 = {0.02, 0.2, 0.2 , 0.58};
  std::vector<std::vector<float>> Closure = {COMP1, COMP2, COMP3, COMP4};
 
  // === Create the data sets
  // State the names explicitly 
  std::vector<TString> Data_Names = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> Data_ntrk = F.MakeTH1F(Data_Names, 500, 0, 20); 
  TH1F* trk1 = Data_ntrk.at(0);
  TH1F* trk2 = Data_ntrk.at(1);
  TH1F* trk3 = Data_ntrk.at(2);
  TH1F* trk4 = Data_ntrk.at(3);
  std::vector<TH1F*> Data = {trk1, trk2, trk3, trk4};
 
  // Fill the histograms with the composition fractions  
  FillHist(trk1, COMP1, {1, 0.9, 0.1});
  FillHist(trk2, COMP2, {1, 0.9, 0.1});
  FillHist(trk3, COMP3, {1, 0.9, 0.1});
  FillHist(trk4, COMP4, {1, 0.9, 0.1});

  // Add some styles 
  trk1 -> SetLineColor(kRed);
  trk2 -> SetLineColor(kBlue);
  trk3 -> SetLineColor(kOrange);
  trk4 -> SetLineColor(kGreen); 
 
  // ======================= End of Data Generation ================ //  

//  TCanvas* can = new TCanvas("Cant", "Cant", 800, 800);
//  can -> SetLogy();
//  trk1 ->Draw("SAMEHIST"); 
//  trk2 -> Draw("SAMEHIST");
//  trk3 -> Draw("SAMEHIST");
//  trk4 -> Draw("SAMEHIST");
//  trk2 -> Draw("SAMEHIST*");
//  can -> Update();

  
  //U.TestFit(Hists, Data_ntrk, 0, 20, Closure);
  //U.TestTailAndDeconv(Hists[0], Hists[1], 300, 0, 20);
  //U.TestDeconvolution(Hists[0], Hists[1], 300);
  //U.TestSubtraction(trk4, 4, Hists, 0, 20, COMP4);
  P.Threshold("/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/Merger/Merger.root");  
  //P.TestMinimalAlgorithm(Data, 0, 20, 0.1, Hists, Closure);
  //P.TestGaussianAlgorithm(Data, 0, 20, 0.1, Hists, Closure);
  
 }

void StandaloneApplications(int argc, char**argv)
{
  PostAnalysis();
}

int main(int argc, char** argv)
{
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplications(app.Argc(), app.Argv());
  app.Run();

  return 0;
}
