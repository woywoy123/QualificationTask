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
  int Mode = 1; 
  bool Test = true;
  float bins = 500;
  float min = 0;
  float max = 20;
  float npts = 500000; 

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

  // Component testing 
  if( Test == true)
  {
    // Test the subtraction 
    //BFT.Subtraction();  
    //BFT.NormalFit(trk2_N, ntrk_Data[1], CLS2, 0.04, 20); 
    DFT.NormalFit(trk2_N, ntrk_Data[1], CLS2, 0, 20); 
  }
 
 
  std::cout << "Fin" << std::endl;
 
 
 
  
  // To do:
  // - write LR function 
  // - write RooFit class
  // - Define the ROOFIT tail replace class as experimental
  // - Write Gaussian function 
  // - write minimal where the tail replace happens after LR loop
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
