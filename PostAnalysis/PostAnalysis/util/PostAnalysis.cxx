#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Evaluate.h>


int main(int argc, char* argv[])
{
  TString JE = argv[1]; 
  //IOTest(); 
  //RooFitBaseFunctionTest();  
  //TestFits_NTruth_NTrack();  
  //TestRead(); 
  //CompareToTruth("ntrk_ntru.root");  
  TestFits_AllTruth_ToTrack(JE);  
 









  std::cout << "fin" << std::endl;
  return 0; 
}
