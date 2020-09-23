#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DerivedFunctions.h>
#include<PostAnalysis/Constants.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/UnitTest.h>
#include<TApplication.h>

using namespace Constants;

void PostAnalysis()
{
  // ==== Classes being imported ==== //
  BaseFunctions B;
  DistributionGenerators D; 
  Presentation P; 
  BaseFunctionTest BFT; 
  DerivedFunctionTest DFT;
  DerivedFunctions DF;

  // ==== Constants used for the algorithm ==== //
  // Execution parameter 
  int Mode = 10;  // Change to 0 - MC, 1 - Toy, 2 - Data, 3 - Presentation
  bool Test = false; // Test Components 
  int Shift = 0;

  // Histogram parameters  
  float bins = 500; 
  float min = 0;
  float max = 20;

  // Gaussian Parameter
  float npts = 500000; 
  float mean = 0;
  float stdev = 0.01; 

  // Other parameters
  float offset = 0.1;
  int iter = 50;
  int cor_loop = 100; // Correction loop number 

  // ==== Forward declaration for Histograms ==== //
  std::vector<TH1F*> trk1_N;
  std::vector<TH1F*> trk2_N;
  std::vector<TH1F*> trk3_N;
  std::vector<TH1F*> trk4_N;
  std::vector<TH1F*> trk5_N;

  std::vector<TH1F*> ntrk_Data; 
  std::vector<std::vector<TH1F*>> Truth_Sets; 
  std::vector<float> CLS1;
  std::vector<float> CLS2;
  std::vector<float> CLS3;
  std::vector<float> CLS4; 
  std::vector<float> CLS5; 
  std::vector<std::vector<float>> Closure; 

  //==== Monte Carlo Parameters 
  std::map<TString, std::vector<float>> Params;
 
  // Monte Carlo Reading 
  if( Mode == 0 )
  {
    // Generate hists and fill with MC data
    trk1_N = D.FillTH1F(trk_1, MC_dir); 
    trk2_N = D.FillTH1F(trk_2, MC_dir); 
    trk3_N = D.FillTH1F(trk_3, MC_dir); 
    trk4_N = D.FillTH1F(trk_4, MC_dir);
    trk5_N = D.FillTH1F(trk_5, MC_dir);
    ntrk_Data = D.FillTH1F({"dEdx_ntrk_1", "dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4", "dEdx_ntrk_5"}, MC_dir);  
    Truth_Sets = {trk1_N, trk2_N, trk3_N, trk4_N, trk5_N}; 
   
    // Get Closure values and fill data 
    CLS1 = B.ClosureAndData(trk1_N, ntrk_Data[0]); 
    CLS2 = B.ClosureAndData(trk2_N, ntrk_Data[1]); 
    CLS3 = B.ClosureAndData(trk3_N, ntrk_Data[2]); 
    CLS4 = B.ClosureAndData(trk4_N, ntrk_Data[3]);   
    CLS5 = B.ClosureAndData(trk4_N, ntrk_Data[4]);   
    Closure = {CLS1, CLS2, CLS3, CLS4, CLS5};  
  }

  // Toy model 
  if( Mode == 1 )     
  {
    // Make Hists
    trk1_N = B.MakeTH1F(trk_1, bins, min, max);    
    trk2_N = B.MakeTH1F(trk_2, bins, min, max);    
    trk3_N = B.MakeTH1F(trk_3, bins, min, max);    
    trk4_N = B.MakeTH1F(trk_4, bins, min, max); 
    trk5_N = B.MakeTH1F(trk_5, bins, min, max); 
    ntrk_Data = B.MakeTH1F(Data_Names, bins, min, max);  
    Truth_Sets = {trk1_N, trk2_N, trk3_N, trk4_N, trk5_N};
   
    // Fill Hists 
    D.Landau(trk1_N, COMP1, LandauParameters, npts, min, max);
    D.Landau(trk2_N, COMP2, LandauParameters, npts, min, max);
    D.Landau(trk3_N, COMP3, LandauParameters, npts, min, max);
    D.Landau(trk4_N, COMP4, LandauParameters, npts, min, max); 
    D.Landau(trk5_N, COMP5, LandauParameters, npts, min, max); 
        
    // Get Closure values and fill data 
    CLS1 = B.ClosureAndData(trk1_N, ntrk_Data[0]); 
    CLS2 = B.ClosureAndData(trk2_N, ntrk_Data[1]); 
    CLS3 = B.ClosureAndData(trk3_N, ntrk_Data[2]); 
    CLS4 = B.ClosureAndData(trk4_N, ntrk_Data[3]);  
    CLS5 = B.ClosureAndData(trk5_N, ntrk_Data[4]); 
    Closure = {CLS1, CLS2, CLS3, CLS4, CLS5}; 
  }

  // Presentation Stuff
  if ( Mode == 3 )
  {
    Presentation P; 
    P.ThresholdEffects();  
    Test = false; 
  }

  // Component testing 
  if( Test == true)
  {
    //BFT.NormalFit(trk2_N, ntrk_Data[1], CLS2, 0.04, 20); 
    //DFT.NormalFit(trk1_N, ntrk_Data[0], CLS2, 0, 20); 
    //BFT.Convolve(trk1_N[1], trk1_N[1], trk1_N[3]); 
    //BFT.Deconvolve(trk2_N[1], trk1_N[0], offset, 25);
    //DFT.ShiftTest(trk1_N[0], Shift);
    //DFT.ReplaceShiftTail(trk1_N[0], trk1_N[1], Shift);
    //DFT.DeconvolveReconvolve(trk1_N, offset, iter);
    //DFT.DeconvolveGaussianFit(ntrk_Data[0], ntrk_Data[1], mean, stdev, offset, iter);
    //BFT.Constraint(); 
 
    //Gaussian Parameter used for deconvolution
//    Params["Gaussian"] = {0, 0.1};
//    Params["m_s"] = {-3, -10, -10, -10, -1};
//    Params["m_e"] = {3, 20, 20, 20, 15};        
//    Params["s_s"] = {0.01, 0.01, 0.01, 0.01, 0.01};
//    Params["s_e"] = {5, 20, 30, 40, 60};



 
    //P.MainAlgorithm(ntrk_Data, Params, offset, iter, cor_loop, Truth_Sets); 
  }
  else
  {
    ////Gaussian Parameter used for deconvolution
    //Params["Gaussian"] = {0, 0.1};
    //Params["m_s"] = {-5, -10, -20, -20, -20};
    //Params["m_e"] = {5, 30, 40, 40, 40};        
    //Params["s_s"] = {0.01, 0.01, 0.01, 0.01, 0.01};
    //Params["s_e"] = {5, 20, 40, 120, 120};

    // ===== Good parameters that have been tested (out.root)
    //Gaussian Parameter used for deconvolution
    Params["Gaussian"] = {0, 1};
    Params["m_e"] = {1, 1, 1, 1, 1};
    Params["m_s"] = {0, 0, 0, 0, 0};        
    Params["s_s"] = {0.6, 0.6, 0.6, 0.6, 0.6};
    Params["s_e"] = {1.5, 1.5, 1.5, 1.5, 1.5};

    //P.MainAlgorithm(ntrk_Data, Params, offset, iter, cor_loop, Truth_Sets);   
    P.DataAnalysis(Params, offset, iter, cor_loop, bins, min, max);   
    //P.ReconstructNTrack();

    //P.AlgorithmPlots("/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/AnalysisOutput/out.root", cor_loop); 
  }
 
  std::cout << "Fin" << std::endl; 
  
}






























void StandaloneApplications(int argc, char** argv){PostAnalysis();}
int main(int argc, char** argv)
{
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplications(app.Argc(), app.Argv());
  app.Run();
  return 0; 
}
