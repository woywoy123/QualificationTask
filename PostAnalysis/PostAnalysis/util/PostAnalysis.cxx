#include<PostAnalysis/BaseFunctions.h>
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
  Plotting P; 
  BaseFunctionTest BFT; 
  DerivedFunctionTest DFT;

  // ==== Constants used for the algorithm ==== //
  // Execution parameter 
  int Mode = 2;  // Change to 0 - MC, 1 - Toy, 2 - RealData
  bool Test = true; // Test Components 
  int Shift = 0;

  // Histogram parameters  
  float bins = 500; 
  float min = 0;
  float max = 20;

  // Gaussian Parameter
  float npts = 500000; 
  float mean = 0;
  float stdev = 0.01; 
  float m_s = -2;
  float m_e = 2;
  float s_s = 0.001;
  float s_e = 0.1;

  // Other parameters
  float offset = 0.1;
  float Gamma = 0.1;
  int iter = 250;
  int cor_loop = 3; // Correction loop number 
  std::vector<float> Params = {mean, stdev, m_s, m_e, s_s, s_e}; 

  // ==== Forward declaration for Histograms ==== //
  std::vector<TH1F*> Pure; 
  std::vector<TH1F*> trk1_N;
  std::vector<TH1F*> trk2_N;
  std::vector<TH1F*> trk3_N;
  std::vector<TH1F*> trk4_N;
  std::vector<TH1F*> ntrk_Data; 
  std::vector<std::vector<TH1F*>> Truth_Sets; 
  std::vector<float> CLS1;
  std::vector<float> CLS2;
  std::vector<float> CLS3;
  std::vector<float> CLS4; 
  std::vector<std::vector<float>> Closure; 

  // Monte Carlo Reading 
  if( Mode == 0 )
  {
    // Generate hists and fill with MC data
    trk1_N = D.FillTH1F(trk_1, MC_dir); 
    trk2_N = D.FillTH1F(trk_2, MC_dir); 
    trk3_N = D.FillTH1F(trk_3, MC_dir); 
    trk4_N = D.FillTH1F(trk_4, MC_dir);
    ntrk_Data = D.FillTH1F(Data_Names, MC_dir);  
    Truth_Sets = {trk1_N, trk2_N, trk3_N, trk4_N}; 
   
    // Get Closure values and fill data 
    CLS1 = B.ClosureAndData(trk1_N, ntrk_Data[0]); 
    CLS2 = B.ClosureAndData(trk2_N, ntrk_Data[1]); 
    CLS3 = B.ClosureAndData(trk3_N, ntrk_Data[2]); 
    CLS4 = B.ClosureAndData(trk4_N, ntrk_Data[3]);   
    Closure = {CLS1, CLS2, CLS3, CLS4};  
  }

  // Toy model 
  if( Mode == 1 )     
  {
    // Make Hists
    trk1_N = B.MakeTH1F(trk_1, bins, min, max);    
    trk2_N = B.MakeTH1F(trk_2, bins, min, max);    
    trk3_N = B.MakeTH1F(trk_3, bins, min, max);    
    trk4_N = B.MakeTH1F(trk_4, bins, min, max); 
    ntrk_Data = B.MakeTH1F(Data_Names, bins, min, max);  
    Truth_Sets = {trk1_N, trk2_N, trk3_N, trk4_N};
   
    // Fill Hists 
    D.Landau(trk1_N, COMP1, LandauParameters, npts, min, max);
    D.Landau(trk2_N, COMP2, LandauParameters, npts, min, max);
    D.Landau(trk3_N, COMP3, LandauParameters, npts, min, max);
    D.Landau(trk4_N, COMP4, LandauParameters, npts, min, max); 
    
    // Get Closure values and fill data 
    CLS1 = B.ClosureAndData(trk1_N, ntrk_Data[0]); 
    CLS2 = B.ClosureAndData(trk2_N, ntrk_Data[1]); 
    CLS3 = B.ClosureAndData(trk3_N, ntrk_Data[2]); 
    CLS4 = B.ClosureAndData(trk4_N, ntrk_Data[3]);  
    Closure = {CLS1, CLS2, CLS3, CLS4}; 
  }

  // Data
  if ( Mode == 2 )
  {
    TH1F* trk1 = D.FillTH1F("dEdx_out_ntrk1_calib", energies, MC_dir); 
    ntrk_Data.push_back(trk1);
    
    for (TH1F* H : D.FillTH1F({"dEdx_ntrk_2", "dEdx_ntrk_3", "dEdx_ntrk_4"} , MC_dir))
    {
      ntrk_Data.push_back(H);
    } 
 
    /// Delete below    
    trk1_N = D.FillTH1F(trk_1, MC_dir); 
    trk2_N = D.FillTH1F(trk_2, MC_dir); 
    trk3_N = D.FillTH1F(trk_3, MC_dir); 
    trk4_N = D.FillTH1F(trk_4, MC_dir);
    Truth_Sets = {trk1_N, trk2_N, trk3_N, trk4_N}; 
    // ==================  
     
    DFT.MainAlgorithm(ntrk_Data, Params, offset, iter, cor_loop, Gamma, Truth_Sets); // << Replace this later  
    
    Test = false;
  }

  // Component testing 
  if( Test == true)
  {
    // Test the subtraction 
    //BFT.Subtraction();  
    //BFT.NormalFit(trk2_N, ntrk_Data[1], CLS2, 0.04, 20); 
    //DFT.NormalFit(trk1_N, ntrk_Data[0], CLS2, 0, 20); 
    //BFT.Convolve(trk1_N[1], trk1_N[1], trk1_N[3]); 
    //BFT.Deconvolve(trk2_N[1], trk1_N[0], offset, 25);
    //DFT.ShiftTest(trk1_N[0], Shift);
    //DFT.ReplaceShiftTail(trk1_N[0], trk1_N[1], Shift);
    //DFT.DeconvolveReconvolve(trk1_N, offset, iter);
    //DFT.DeconvolveGaussianFit(ntrk_Data[0], ntrk_Data[1], mean, stdev, offset, iter);
  
    //P.PlotHists(Truth_Sets, ntrk_Data);
    //DFT.MainAlgorithm(ntrk_Data, Params, offset, iter, cor_loop, Gamma, Truth_Sets);   
  }
 
 
  std::cout << "Fin" << std::endl;
 
  
  // To do:
  // - Write the threshold function and do a fit to remove 1 trk  
 
 
 
 
 
 
 
 
 
 
 
  
}






























void StandaloneApplications(int argc, char** argv){PostAnalysis();}
int main(int argc, char** argv)
{
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplications(app.Argc(), app.Argv());
  app.Run();
  return 0; 
}
