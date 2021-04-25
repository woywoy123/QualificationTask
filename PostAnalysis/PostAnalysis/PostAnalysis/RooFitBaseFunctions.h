#ifndef ROOFITBASEFUNCTIONS_H
#define ROOFITBASEFUNCTIONS_H
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooHistPdf.h>
#include<RooArgList.h>
#include<RooAddPdf.h>
#include<RooFormulaVar.h>
#include<TF1.h>
#include<RooGaussian.h>
#include<RooFitResult.h>
#include<RooFFTConvPdf.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Plotting.h>

#include<RooAbsReal.h>
#include<RooMinimizer.h>

using namespace RooFit;

// ==================  Basic RooFit Variables
static std::vector<RooRealVar*> RooRealVariable(std::vector<TString> Names, std::vector<float> Guess, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables; 
  for (int i(0); i < Names.size(); i++)
  {
    RooRealVar* r = new RooRealVar(Names[i], Names[i], Guess[i], Begin[i], End[i]); 
    Variables.push_back(r); 
  }
  return Variables; 
}

static std::vector<RooRealVar*> RooRealVariable(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables; 
  for (int i(0); i < Names.size(); i++)
  {
    RooRealVar* r = new RooRealVar(Names[i], Names[i], Begin[i], End[i]); 
    Variables.push_back(r); 
  }
  return Variables; 
}

// PDF + Data related Variables 
static std::vector<RooDataHist*> RooDataVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<TH1F*> Hist)
{
  std::vector<RooDataHist*> Out; 
  for (int i(0); i < Hist.size(); i++)
  {
    RooDataHist* H = new RooDataHist(Name[i], Name[i], *domain, RooFit::Import(*Hist[i])); 
    Out.push_back(H); 
  }
  return Out; 
}

static std::vector<RooHistPdf*> RooPdfVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<RooDataHist*> Hist)
{
  std::vector<RooHistPdf*> Out; 
  for (int i(0); i < Hist.size(); i++)
  {
    RooHistPdf* H = new RooHistPdf(Name[i], Name[i], *domain, *Hist[i]); 
    Out.push_back(H); 
  }
  return Out; 
}

static std::vector<RooGaussian*> RooGaussianVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<RooRealVar*> m, std::vector<RooRealVar*> s)
{
  std::vector<RooGaussian*> Out; 
  for (int i(0); i < Name.size(); i++)
  {
    RooGaussian* H = new RooGaussian(Name[i], Name[i], *domain, *m[i], *s[i]); 
    Out.push_back(H); 
  }
  return Out; 
}

static std::vector<RooFFTConvPdf*> RooFFTVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<RooGaussian*> g, std::vector<RooHistPdf*> p)
{
  std::vector<RooFFTConvPdf*> Out; 
  for (int i(0); i < Name.size(); i++)
  {
    RooFFTConvPdf* H = new RooFFTConvPdf(Name[i], Name[i], *domain, *p[i], *g[i]); 
    Out.push_back(H); 
  }
  return Out; 
}

static void CopyPDFToTH1F(RooHistPdf* pdf, RooRealVar* domain, TH1F* Out, TH1F* Data)
{
  auto T = pdf -> generateBinned(*domain, Data -> GetEntries(), true); 
  auto H = T -> createHistogram("temp", *domain, RooFit::Binning(Data -> GetNbinsX(), Data -> GetXaxis() -> GetXmin(), Data -> GetXaxis() -> GetXmax())); 

  Out -> Reset(); 
  Out -> Add(H, 1); 
  delete H;
}

static void CopyPDFToTH1F(RooFFTConvPdf* pdf, RooRealVar* domain, TH1F* Out, TH1F* Data)
{

  auto T = pdf -> generateBinned(*domain, Data -> GetEntries(), true); 
  auto H = T -> createHistogram("temp", *domain, RooFit::Binning(Data -> GetNbinsX(), Data -> GetXaxis() -> GetXmin(), Data -> GetXaxis() -> GetXmax())); 

  Out -> Reset(); 
  Out -> Add(H, 1);  
  delete H; 
}

static void BulkDelete(std::vector<RooDataHist*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooRealVar*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooHistPdf*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooFormulaVar*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooGaussian*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooFFTConvPdf*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}


// ==================== Base Functions 
std::map<TString, std::vector<float>> Normalization(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 
std::map<TString, std::vector<float>> NormalizationShift(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 
std::map<TString, std::vector<float>> ConvolutionFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = "");
std::map<TString, std::vector<float>> DeConvolutionFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = "");











#endif
