#include "Evaluate.h"

int main(int argc, char* argv[])
{
  
  TString DataName = argv[1]; 
  std::cout << DataName << std::endl;

  Debug_Cases(DataName);
  return 0; 
}
