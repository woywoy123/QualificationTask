#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Evaluate.h>
#include<PostAnalysis/Debug.h>

int main(int argc, char* argv[])
{
  TString JE = argv[1]; 
  TString Mode = argv[2];  
  TString File = argv[3]; 

  //File = "./Merged_MC.root"; 
  //Mode = "Simultaneous"; 
  //JE = "x";
  //std::cout << File << std::endl;
  //IOTest(); 
  //RooFitBaseFunctionTest();  
  //TestFits_NTruth_NTrack();  
  //TestRead(); 
  //CompareToTruth("ntrk_ntru.root"); 
  //std::cout << "+++" << JE << " " << Mode << std::endl;
  //ReadOutputFileToMap("Fit_Tracks.root"); 

  TestFits_AllTruth_ToTrack(JE, Mode, File);  
  //MultiTrackTruthComparison("MultiTrackFit.root"); 
  
  //Proxy(); 
  //SmoothingTest();

  std::cout << "fin" << std::endl; 
}
