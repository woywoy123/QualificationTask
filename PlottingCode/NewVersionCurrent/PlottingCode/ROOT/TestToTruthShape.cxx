#include "../PlottingCode/TestToTruthShape.h"

int main(int argc, char* argv[])
{
  TString Fit_ROOT = argv[1]; 
  TString Merged_ROOT = argv[2]; 
  
  std::map<_TS, std::vector<TH1F*>> CTIDE = ReadCTIDE(Merged_ROOT);
  ReadPostAnalysis(Fit_ROOT); 





  return 0; 
}
