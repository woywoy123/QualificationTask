#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Evaluate.h>
#include<PostAnalysis/Debug.h>
#include<PostAnalysis/ReportFigures.h>

int main(int argc, char* argv[])
{
  TString JE = argv[1]; 
  TString Mode = argv[2];  
  TString File = argv[3]; 

  //Figure_Proxy();

  //File = "./Merged_MC.root"; 
  //Mode = "Simultaneous"; 
  //JE = "x";
  //std::cout << File << std::endl;
  //IOTest(); 
  //RooFitBaseFunctionTest();  
  //TestFits_NTruth_NTrack();  
  //TestRead(); 
  //ReadOutputFileToMap("Fit_Tracks.root"); 
    
  //SplitTest(JE);
  
  
  if (Mode.Contains("Debug")){ Proxy(JE, File, Mode); } 
  else {FastFits(JE, Mode, File);} 
  
  //MultiTrackTruthComparison("MultiTrackFit.root"); 
  std::cout << "fin" << std::endl; 
}
