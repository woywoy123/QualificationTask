#include<PostAnalysis/RooFitFunctions.h>










void FitDeconvolution(TH1F* Data, std::vector<TH1F*> PDF_H)
{
  // First we get the domain of the Data histogram we are fitting 
  float x_min = Data -> GetXaxis() -> GetXmin(); 
  float x_max = Data -> GetXaxis() -> GetXmax(); 
  int bins = Data -> GetNbinsX(); 

  // Declare the domain using RooVariables 
  RooRealVar* x =  new RooRealVar("x", "x", x_min, x_max); 
  
























  delete x; 

}
