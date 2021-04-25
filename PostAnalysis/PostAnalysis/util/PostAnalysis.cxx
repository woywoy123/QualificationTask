#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Evaluate.h>


int main(int argc, char** argv)
{
  //IOTest(); 
  //RooFitBaseFunctionTest();  
  //TestFits_NTruth_NTrack();  
  //TestRead(); 
  CompareToTruth("ntrk_ntru.root");  


  std::cout << "fin" << std::endl;
  return 0; 
}
