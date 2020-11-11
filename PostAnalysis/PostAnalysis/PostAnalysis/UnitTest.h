#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DerivedFunctions.h>
#include<TH2F.h>


#ifndef UNITTEST_H
#define UNITTEST_H

class BaseFunctionTest
{
  public:
    void NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max);  
    void Convolve(TH1F* Hist1, TH1F* Hist2, TH1F* Expection);
    void Deconvolve(TH1F* Trk2, TH1F* Trk1, float offset, int iter); 
    void Constraint();  
};

class DerivedFunctionTest
{
  public:
    void NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max);
    void ShiftTest(TH1F* H1, int Shift);
    void ReplaceShiftTail(TH1F* Source, TH1F* Target, int Shift);
    void DeconvolveReconvolve(std::vector<TH1F*> trk1, float offset, int iter);
    void Deconvolve(TH1F* Hist);
    void AlgorithmTest();
};

class Presentation
{
  public:
    void ThresholdEffects(); 
    void MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure);
    void DataAnalysis(std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, int bins, float min, float max);
    void AlgorithmPlots(TString dir, int iter); 
    void ReconstructNTrack(std::vector<TH1F*> Hists);
    void ShowMonteCarloDistribution(); 
    void ShowMonteCarloClusterPosition(); 
};

#endif
