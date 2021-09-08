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
  //CompareToTruth("ntrk_ntru.root"); 
  //std::cout << "+++" << JE << " " << Mode << std::endl;
  //ReadOutputFileToMap("Fit_Tracks.root"); 
    
  //CompareToTruth(Mode, JE);  
  //SplitTest(JE);
  
  
  if (Mode.Contains("Debug")){ Proxy(JE, File, Mode); } 
  else {TestFits_AllTruth_ToTrack(JE, Mode, File); }

  //MultiTrackTruthComparison("MultiTrackFit.root"); 
  std::cout << "fin" << std::endl; 
}
