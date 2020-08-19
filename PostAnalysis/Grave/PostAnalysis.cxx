// Including commonly reused functions
#include<PostAnalysis/Functions.h>
#include<PostAnalysis/Verification.h>
#include<PostAnalysis/UnitClosures.h>
#include<iostream>

void PostAnalysis()
{

  Benchmark B;
  Verification V; 
  Fit_Functions f;
  Functions F;
  UnitClosures U;
  Presentation P;
  DataGeneration D; 

  // ======================== Base Variables ======================== // 
  TString Mode = "MC";
  int bins = 500;
  int min = 0;
  int max = 20;
  float offset = 0.1;
  int iter = 1000;
  
  // ===== Data vectors
  std::vector<TString> Data_Names = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> Data_ntrk = F.MakeTH1F(Data_Names, bins, min, max); 

  // ===== Monte Carlo dir
  TString dir = "/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/MonteCarlo/Merged.root"; 
   

  // ===== Landau Parameter
  std::vector<float> LP = {1, 0.9, 0.1};

  // Contamination Composition 
  std::vector<float> COMP1;
  std::vector<float> COMP2; 
  std::vector<float> COMP3; 
  std::vector<float> COMP4; 
  std::vector<std::vector<float>> Closure;
  
  // ====================== Get the Monte Carlo data ================ //
  if ( Mode == "MC")
  { 
    // ===== Define the names  
    std::vector<TString> trk_P_n = Constants::Pure_Names;
    std::vector<TString> trk_1_n = Constants::trk_1;
    std::vector<TString> trk_2_n = Constants::trk_2;
    std::vector<TString> trk_3_n = Constants::trk_3;
    std::vector<TString> trk_4_n = Constants::trk_4;

    // ===== Create the histograms
    std::vector<TH1F*> Pure_trks = F.MakeTH1F(trk_P_n, bins, min, max, "_P");
    std::vector<TH1F*> trk_1_V = F.MakeTH1F(trk_1_n, bins, min, max);
    std::vector<TH1F*> trk_2_V = F.MakeTH1F(trk_2_n, bins, min, max);
    std::vector<TH1F*> trk_3_V = F.MakeTH1F(trk_3_n, bins, min, max);
    std::vector<TH1F*> trk_4_V = F.MakeTH1F(trk_4_n, bins, min, max);
    std::vector<std::vector<TH1F*>> trk_P = {trk_1_V, trk_2_V, trk_3_V, trk_4_V};

    // ===== Fill the histograms through the MC 
    D.MonteCarlo(Pure_trks, dir, "_P");
    D.MonteCarlo(trk_1_V, dir);
    D.MonteCarlo(trk_2_V, dir);
    D.MonteCarlo(trk_3_V, dir);
    D.MonteCarlo(trk_4_V, dir);
  
    // ===== Merge the histograms and get the composition  
    COMP1 = D.MergeToys(trk_1_V, Data_ntrk.at(0));
    COMP2 = D.MergeToys(trk_2_V, Data_ntrk.at(1));   
    COMP3 = D.MergeToys(trk_3_V, Data_ntrk.at(2));   
    COMP4 = D.MergeToys(trk_4_V, Data_ntrk.at(3));   
    Closure = {COMP1, COMP2, COMP3, COMP4};
  
   
    //U.TestFit(trk_P, Data_ntrk, min, max, Closure);
    //U.TestTailAndDeconv(Pure_trks[0], Pure_trks[1], iter, min, max); 
    //U.TestDeconvolution(Pure_trks[0], Pure_trks[1], iter); 
    //U.TestSubtraction(Data_ntrk.at(2), 3, trk_3_V, min, max, COMP3); 
    //P.TestMinimalAlgorithm(Data_ntrk, min, max, offset, trk_P);
    //P.TestGaussianAlgorithm(Data_ntrk, min, max, offset); 
    TCanvas* Truth = B.ClosurePlot("Truth", Data_ntrk, trk_P);
    
    Debug D; 
    D.GaussianFlow(Data_ntrk, min, max, offset); 
  
  
  }
  
  // ======================= Toy Data ===================== //
  if ( Mode == "Landau" )
  {  

    // ===== Define the names  
    std::vector<TString> trk_P_n = Constants::Pure_Names;
    std::vector<TString> trk_1_n = Constants::trk_1;
    std::vector<TString> trk_2_n = Constants::trk_2;
    std::vector<TString> trk_3_n = Constants::trk_3;
    std::vector<TString> trk_4_n = Constants::trk_4;

    // ===== Create the histograms
    std::vector<TH1F*> trk_P = F.MakeTH1F(trk_P_n, bins, min, max, "_P");
    std::vector<TH1F*> trk_1_V = F.MakeTH1F(trk_1_n, bins, min, max);
    std::vector<TH1F*> trk_2_V = F.MakeTH1F(trk_2_n, bins, min, max);
    std::vector<TH1F*> trk_3_V = F.MakeTH1F(trk_3_n, bins, min, max);
    std::vector<TH1F*> trk_4_V = F.MakeTH1F(trk_4_n, bins, min, max);

    // ===== Histograms: Pure, 1trk, 2trk, 3trk, 4trk
    D.IdealLandau(trk_P, {1, 1, 1, 1}, LP);
    D.IdealLandau(trk_1_V, {0.8  , 0.2 , 0.  , 0.  }, LP);
    D.IdealLandau(trk_2_V, {0.  , 1. , 0.  , 0.  }, LP);
    D.IdealLandau(trk_3_V, {0.01, 0.2, 0.59, 0.2 }, LP);
    D.IdealLandau(trk_4_V, {0.02, 0.2, 0.2 , 0.58}, LP);
    
    //  ===== Merge the vectors into a single TH1F for data 
    COMP1 = D.MergeToys(trk_1_V, Data_ntrk[0]);
    COMP2 = D.MergeToys(trk_2_V, Data_ntrk[1]);
    COMP3 = D.MergeToys(trk_3_V, Data_ntrk[2]);
    COMP4 = D.MergeToys(trk_4_V, Data_ntrk[3]);
    Closure = {COMP1, COMP2, COMP3, COMP4};
  }

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
