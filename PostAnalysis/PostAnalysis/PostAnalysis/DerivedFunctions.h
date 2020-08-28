#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Constants.h>
#include<PostAnalysis/Plotting.h>

#ifndef DERIVEDFUNCTIONS_H
#define DERIVEDFUNCTIONS_H

class DerivedFunctions
{
  public:
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data, RooRealVar* Domain, std::vector<float> Begin, std::vector<float> End, std::vector<TString> Names);
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data, float min, float max);
    void RooShift(TH1F* H1, TH1F* H2);

    int NumericalShift(TH1F* H1, TH1F* H2); 
    void ReplaceShiftTail(TH1F* Source, TH1F* Target, float offset = 0); 
    std::vector<TH1F*> nTRKGenerator(TH1F* trk1, TH1F* trk2, float offset, int iter);    
    TH1F* GaussianConvolve(TH1F* Hist, float mean, float stdev, int Toys = Constants::GaussianToys);
    std::map<TString, float> FitGaussian(TH1F* GxTrk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> Params, float offset, int iter);
    void RemoveArtifact(TH1F* Conv); // Move to base later
    std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, float Gamma, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure); 

};


#endif
