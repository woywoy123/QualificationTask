#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Evaluate.h>


int main(int argc, char* argv[])
{
  TString JE = "Blayer"; //argv[1]; 
  TString Mode = "ShiftNormal"; // argv[2];  
  TString File = "Merged_MC.root"; //argv[3]; 
  //std::cout << File << std::endl;
  //IOTest(); 
  //RooFitBaseFunctionTest();  
  //TestFits_NTruth_NTrack();  
  //TestRead(); 
  //CompareToTruth("ntrk_ntru.root"); 
  std::cout << "+++" << JE << " " << Mode << std::endl;
  TestFits_AllTruth_ToTrack(JE, Mode, File);  
  //ReadOutputFileToMap("Fit_Tracks.root"); 

  //MultiTrackTruthComparison("MultiTrackFit.root"); 
  //CompareToTruth("MultiTrackFit.root"); 





  //std::cout << "fin" << std::endl;
  //return 0; 
}
