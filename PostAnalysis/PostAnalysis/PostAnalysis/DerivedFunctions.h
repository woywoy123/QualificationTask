#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Constants.h>

#ifndef DERIVEDFUNCTIONS_H
#define DERIVEDFUNCTIONS_H

class DerivedFunctions
{
  public:
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data, RooRealVar* Domain, std::vector<float> Begin, std::vector<float> End, std::vector<TString> Names);
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data);
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data, float min, float max);
};

#endif
