#ifndef ROOFITBASEFUNCTIONS_H
#define ROOFITBASEFUNCTIONS_H
#include<RooRealVar.h>
#include<RooDataHist.h>
#include<RooHistPdf.h>
#include<RooArgList.h>
#include<RooAddPdf.h>
#include<RooFormulaVar.h>
#include<TF1.h>
#include<RooGaussModel.h>
#include<RooLandau.h>
#include<RooFitResult.h>
#include<RooFFTConvPdf.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Plotting.h>
#include<TMatrixDSym.h>
#include<RooAbsReal.h>
#include<RooMinimizer.h>
#include<RooSimultaneous.h>
#include<RooDataSet.h>
#include<RooCategory.h>
#include<TFractionFitter.h>

using namespace RooFit;
const std::vector<TString> FitRanges_Names = {"Range_ntrk_1", "Range_ntrk_2", "Range_ntrk_3", "Range_ntrk_4"}; 
const int n_cpu = 1; 


// ==================  Basic RooFit Variables
static std::vector<RooRealVar*> RooRealVariable(std::vector<TString> Names, std::vector<float> Guess, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables; 
  for (int i(0); i < Names.size(); i++)
  {
    RooRealVar* r; 
    if (i < Guess.size()){r = new RooRealVar(Names[i], Names[i], Guess[i], Begin[i], End[i]);}
    else {r = new RooRealVar(Names[i], Names[i], Begin[i], End[i]);}
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

static void VariableConstant(std::vector<float> Bools, std::vector<RooRealVar*> Vars)
{
  for (int i(0); i < Vars.size(); i++)
  { 
    bool on; 
    if (Bools[i] == 1){ on = true; }
    else {on = false;}

    Vars[i] -> setConstant(on); 
  }
}


static std::vector<RooRealVar*> ProtectionRealVariable(TString ext, std::vector<TH1F*> PDF, std::map<TString, std::vector<float>> Params, float start, float end)
{
  std::vector<TString> var_N = NameGenerator(PDF, "_" + ext); 
  if (Params[ext + "_s"].size() == 0){Params[ext + "_s"] = std::vector<float>(PDF.size(), start);} 
  if (Params[ext + "_e"].size() == 0){Params[ext + "_e"] = std::vector<float>(PDF.size(), end);} 

  std::vector<RooRealVar*> v_vars;
  if (Params[ext + "_G"].size() != 0){v_vars = RooRealVariable(var_N, Params[ext + "_G"], Params[ext + "_s"], Params[ext + "_e"]);}
  else{v_vars = RooRealVariable(var_N, Params[ext + "_s"], Params[ext + "_e"]);}

  if (Params[ext + "_C"].size() != 0){ VariableConstant(Params[ext + "_C"], v_vars);}

  return v_vars; 
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

static std::vector<RooGaussModel*> RooGaussianVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<RooRealVar*> m, std::vector<RooRealVar*> s)
{
  std::vector<RooGaussModel*> Out; 
  for (int i(0); i < Name.size(); i++)
  {
    RooGaussModel* H = new RooGaussModel(Name[i], Name[i], *domain, *m[i], *s[i]); 
    Out.push_back(H); 
  }
  return Out; 
}

static std::vector<RooFFTConvPdf*> RooFFTVariable(std::vector<TString> Name, RooRealVar* domain, std::vector<RooGaussModel*> g, std::vector<RooHistPdf*> p)
{
  std::vector<RooFFTConvPdf*> Out; 
  for (int i(0); i < Name.size(); i++)
  {
    RooFFTConvPdf* H = new RooFFTConvPdf(Name[i], Name[i], *domain, *p[i], *g[i]); 
    //H -> setInterpolationOrder(2); 
    Out.push_back(H); 
  }
  return Out; 
}

static void CopyPDFToTH1F(RooHistPdf* pdf, RooRealVar* domain, TH1F* Out, TH1F* Data)
{
  TF1 D = TF1(*pdf -> asTF(RooArgList(*domain))); 
  D.SetNpx(Data -> GetNbinsX()); 
  TH1* H = D.CreateHistogram(); 

  Out -> Reset(); 
  Out -> Add(H, 1); 
  
}

static void CopyPDFToTH1F(RooFFTConvPdf* pdf, RooRealVar* domain, TH1F* Out, TH1F* Data)
{
  TF1 D = TF1(*pdf -> asTF(RooArgList(*domain))); 
  D.SetNpx(Data -> GetNbinsX()); 
  TH1* H = D.CreateHistogram(); 

  Out -> Reset(); 
  Out -> Add(H, 1);  
}

static void CaptureResults(RooFitResult* Re, std::map<TString, std::vector<float>>* Output)
{
  TMatrixDSym M = Re -> covarianceMatrix(); 
  
  for (int i(0); i < M.GetNcols(); i++)
  {
    std::vector<float> row;
    for (int j(0); j < M.GetNcols(); j++){row.push_back(M[j][i]);}
    TString R = "Covar_M_i"; R += (i+1); 
    Output -> insert({R, row}); 
  }

  float res = Re -> status(); 
  (*Output)["fit_status"].push_back(res); 
  delete Re; 
}

template<typename mod>
static RooFitResult* MinimizationCustom(mod model, RooDataHist* data, std::map<TString, std::vector<float>> Params, RooRealVar* x)
{

  RooFitResult* re; 
  TString Ranges; 
  for (TString F : FitRanges_Names){if (Params[F].size() != 0){x -> setRange(F, Params[F][0], Params[F][1]); Ranges += (F+ ",");}}
  if (Params["Minimizer"].size() == 0)
  {
    re = model.fitTo(*data, Range(Ranges), SumW2Error(true), NumCPU(n_cpu, 1), Save(), Extended(true));
  }
  else
  {
    RooAbsReal* nll = model.createNLL(*data, Range(Ranges), NumCPU(n_cpu, 1), Extended(true)); 

    RooMinimizer* pg = new RooMinimizer(*nll);
    int print = -1; 
    if (Params["Print"].size() != 0){print = Params["Print"][0]; }
    pg -> setPrintLevel(print); 
    pg -> setPrintEvalErrors(print);
    pg -> setMaxFunctionCalls(Params["Minimizer"][0]); 
    pg -> setMaxIterations(Params["Minimizer"][0]); 

    if (Params["Strategy"].size() != 0){pg -> setStrategy(Params["Strategy"][0]); }
    if (Params["GSL"].size() != 0){pg -> setMinimizerType("GSLMultiMin");}
    //pg -> setProfile(true);
    pg -> setOffsetting(true); 
    pg -> setEvalErrorWall(true); 
    pg -> optimizeConst(true); 
    
    //pg -> migrad();
    //pg -> simplex();
    //pg -> hesse();
    if (Params["Seek"].size() != 0){pg -> seek();}
    pg -> fit("mhr"); 
    for (int i(0); i < 20; i++)
    {
      //pg -> improve();
      //pg -> minos(); 
      pg -> migrad(); 
      re = pg -> save();
      //pg -> minos();
      int res = re -> status(); 
      if (res == 0){break;}
      //pg -> minos();
    }
    pg -> cleanup(); 
   
    delete pg; 
    delete nll;
  }

  return re; 
}

static void BulkDelete(std::vector<RooDataHist*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooRealVar*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooHistPdf*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooFormulaVar*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooGaussModel*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}
static void BulkDelete(std::vector<RooFFTConvPdf*> var){ for (int i(0); i < var.size(); i++){ delete var[i]; }}

static std::vector<float> MultiplyByConstant(std::vector<float> Vec, float c)
{
  std::vector<float> Out; 
  for (int i(0); i < Vec.size(); i++){ Out.push_back(Vec[i]*c); }
  return Out; 
}


// ==================== Base Functions 
std::map<TString, std::vector<float>> Normalization(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 
std::map<TString, std::vector<float>> NormalizationShift(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 
std::map<TString, std::vector<float>> ConvolutionFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = "");
std::map<TString, std::vector<float>> DeConvolutionFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = "");
std::map<TString, std::vector<float>> SimultaneousFFT(std::vector<TH1F*> Data, std::vector<std::vector<TH1F*>> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 
std::map<TString, std::vector<float>> IncrementalFFT(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 
std::map<TString, std::vector<float>> FractionFitter(TH1F* Data, std::vector<TH1F*> PDF_H, std::map<TString, std::vector<float>> Params, TString Name = ""); 


#endif
