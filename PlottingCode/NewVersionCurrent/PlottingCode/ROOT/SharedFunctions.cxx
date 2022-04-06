#include "../PlottingCode/SharedFunctions.h"

std::vector<_TS> Split(_TS Input, _TS Sub)
{
  std::vector<_TS> Output;
  if (!Input.Contains(Sub)){ return Output; }
  while (true)
  {
    int l = Input.First(Sub); 
    Output.push_back(Input(0, l)); 
    Input.Remove(0, l+Sub.Length());
    if (!Input.Contains(Sub)){ Output.push_back(Input); break;}
  }
  return Output; 
}

void BulkDelete(std::vector<TGraph*> gr)
{
  for (TGraph* g : gr){ delete g; }
}
