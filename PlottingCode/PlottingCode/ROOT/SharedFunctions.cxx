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

float WeightedShapeError(std::vector<TH1F*> Pred, std::vector<TH1F*> Truth)
{
  float err = 0; 
  for (int i(0); i < Pred.size(); i++)
  {
    TH1F* P = Pred[i]; 
    TH1F* T = Truth[i]; 
    
    float err_t = 0; 
    for (int j(0); j < P -> GetNbinsX(); j++)
    {
      float sum = int(P -> GetBinContent(j+1)) + T -> GetBinContent(j+1);
      if (sum == 0){continue;}
      float diff = int(P -> GetBinContent(j+1)) - T -> GetBinContent(j+1); 

      err_t += std::pow(diff, 2) / sum;
    }
    err += 0.5*err_t*(T -> Integral());
  }
  
  float n_int = 0; 
  for (TH1F* H : Truth)
  {
    n_int += (H -> Integral());  
  }
  err = err/n_int;
  return err; 
}


