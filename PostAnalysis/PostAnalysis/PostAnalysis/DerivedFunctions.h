#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Constants.h>
#include<PostAnalysis/Plotting.h>
#include<thread>
#include<future>
#include<unistd.h>
#include<Math/MinimizerOptions.h>
#include<TGraph.h>

#ifndef DERIVEDFUNCTIONS_H
#define DERIVEDFUNCTIONS_H

class DerivedFunctions
{
  public:
    
    // === Basic Fitting Functions
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data, RooRealVar* Domain, std::vector<float> Begin, std::vector<float> End, std::vector<TString> Names);
    std::vector<RooRealVar*> FitToData(std::vector<TH1F*> Hists, TH1F* Data, float min, float max);

    // === Tail Replace and scaling 
    int NumericalShift(TH1F* H1, TH1F* H2); 
    void ReplaceShiftTail(TH1F* Source, TH1F* Target, float offset = 0); 
    void SafeScale(std::vector<TH1F*> PDFs, TH1F* Data); 
    void SafeScaleNew(std::vector<TH1F*> PDFs, TH1F* Data); 

    // === PDF generators
    std::vector<TH1F*> nTRKGenerator(TH1F* trk1, TH1F* trk2, float offset, int iter);    
    std::vector<TH1F*> nTRKFrom1Trk(TH1F* trk1);    
    
    // === Gaussian Fitting stuff  
    TH1F* GaussianConvolve(TH1F* Hist, float mean, float stdev, int Toys = Constants::GaussianToys);
    std::vector<TH1F*> ConvolveFit(TH1F* GxTrk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> Params, float offset, int iter);  
   
    
    // === Main Algorithm  
    std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure); 

};


#endif
