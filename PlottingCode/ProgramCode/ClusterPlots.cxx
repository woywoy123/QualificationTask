#include "ClusterPlots.h"

int main(int argc, char* argv[])
{
  
  TString Data = argv[1];
  if (Data.Contains(".root")){Clusters(Data);}

  return 0; 

}
