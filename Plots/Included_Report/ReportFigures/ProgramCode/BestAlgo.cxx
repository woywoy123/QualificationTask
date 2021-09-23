#include "BestAlgo.h"
#include "Evaluate.h"

int main(int argc, char* argv[])
{
  TString mini = argv[1];
  if (!mini.Contains(".root")){WriteFLostFile(F, mini);}
  else { ErrorFlostAlgo(mini); }

  return 0;

}
