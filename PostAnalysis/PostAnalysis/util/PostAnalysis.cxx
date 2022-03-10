#include<PostAnalysis/ExperimentalBench.h>
#include<PostAnalysis/Debug.h>

int main(int argc, char* argv[])
{
  TString JE = argv[1]; 
  TString Mode = argv[2];  
  TString File = argv[3]; 

  if (Mode.Contains("Debug")){ Proxy(JE, File, Mode); } 
  else {FastFits(JE, Mode, File);} 
  
  std::cout << "fin" << std::endl; 
}
