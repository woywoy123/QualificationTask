#include<PostAnalysis/DerivedFunctions.h>

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, 
                                                     TH1F* Data, 
                                                     RooRealVar* Domain, 
                                                     std::vector<float> Begin, 
                                                     std::vector<float> End, 
                                                     std::vector<TString> Names)
{
  BaseFunctions B;
  B.Normalize(Hists);
  std::vector<RooRealVar*> Var = B.RooVariables(Names, Begin, End);
  std::vector<RooHistPdf*> Hist = B.RooPDF(Hists, Domain);
  RooDataHist* data = B.RooData(Data, Domain);

  RooAddPdf model("Model", "Model", B.RooList(Hist), B.RooList(Var));
  model.fitTo(*data);
  
  for (int i(0); i < Var.size(); i++)
  {
    delete Hist[i];
  }
  delete data;
  
  return Var;
} 

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, TH1F* Data, float min, float max)
{
  RooRealVar* x = new RooRealVar("x", "x", min, max);
  std::vector<float> Begin(Hists.size(), 0);
  std::vector<float> End(Hists.size(), Data -> Integral());
  std::vector<RooRealVar*> Vars = FitToData(Hists, Data, x, Begin, End, Constants::FitNames);
  return Vars; 
}
















