#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Evaluate.h>


int main(int argc, char* argv[])
{
  TString JE = argv[1]; 
  TString Mode = argv[2];  
  TString File = argv[3]; 
  std::cout << File << std::endl;
  //IOTest(); 
  //RooFitBaseFunctionTest();  
  //TestFits_NTruth_NTrack();  
  //TestRead(); 
  //CompareToTruth("ntrk_ntru.root");  
  TestFits_AllTruth_ToTrack(JE, Mode, File);  
  //ReadOutputFileToMap("Fit_Tracks.root"); 

  //MultiTrackTruthComparison("MultiTrackFit.root"); 
  //CompareToTruth("MultiTrackFit.root"); 





  //std::cout << "fin" << std::endl;
  //return 0; 
}
