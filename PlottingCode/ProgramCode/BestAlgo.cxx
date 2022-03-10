#include "BestAlgo.h"
#include "Evaluate.h"

int main(int argc, char* argv[])
{
  TString mini = argv[1];
  if (!mini.Contains(".root")){ std::cout << "#> " << mini << std::endl; WriteFLostFile(F, mini);}
  else if (mini.Contains("BestAlgos.root"))
  { 
    BestMinimizerAlgo(mini);
    ErrorFlostAlgo(mini);
  }
  //else if (mini.Contains("Minimizer"))
  //{
  //  TemplateVariants(mini);
  //}

  return 0;

}
