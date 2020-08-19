#include<PostAnalysis/BaseFunctions.h>


std::vector<TH1F*> BaseFunctions::MakeTH1F(std::vector<TString> Names, int bins, float min, float max, TString Extension)
{
  std::vector<TH1F*> Histograms;
  for (TString name : Names)
  {
    if (Extension != ""){ name+=(Extension); }
    TH1F* hist = new TH1F(name, name, bins, min, max);
    Histograms.push_back(hist);
  }

  return Histograms;
}

std::vector<float> BaseFunctions::Ratio(std::vector<TH1F*> Hists, TH1F* Data)
{
  float Lumi = Data -> Integral(); 
  std::vector<float> Ratios;  
  for (TH1F* H : Hists)
  {
    float e = H -> Integral(); 
    Ratios.push_back(e/Lumi);  
  }
  return Ratios; 
}

std::vector<float> BaseFunctions::Ratio(std::vector<RooRealVar*> Vars, TH1F* Data)
{
  float Lumi = Data -> Integral(); 
  std::vector<float> Output; 
  for (RooRealVar* v : Vars)
  {
    Output.push_back((v -> getVal())/Lumi);   
    delete v; 
  }
  return Output;
}

std::vector<float> BaseFunctions::ClosureAndData(std::vector<TH1F*> Hists, TH1F* Data)
{
  for (TH1F* H : Hists){Data -> Add(H);}
  return Ratio(Hists, Data); 
}

void BaseFunctions::Normalize(TH1F* Hist)
{
  float e = Hist -> Integral(); 
  Hist -> Scale(1/e);
}

void BaseFunctions::Normalize(std::vector<TH1F*> Hist)
{
  for (TH1F* H : Hist)
  {
    Normalize(H); 
  }

}

void BaseFunctions::Subtraction(std::vector<TH1F*> ntrk, TH1F* Data, int Exclude, std::vector<float> Ratios)
{
  float lumi = Data -> Integral(); 
  for (int i(0); i < ntrk.size(); i++)
  {
    Normalize(ntrk[i]);
    ntrk[i] -> Scale(Ratios[i]*lumi);
    if ( i != Exclude -1 )
    {  
      Data -> Add(ntrk[i], -1); 
    } 
  }
}

std::vector<RooRealVar*> BaseFunctions::RooVariables(std::vector<TString> Names, std::vector<float> Begin, std::vector<float> End)
{
  std::vector<RooRealVar*> Variables(Names.size()); 
  for (int i(0); i < Names.size(); i++)
  {
    Variables[i] = new RooRealVar(Names[i]+"_Fit", Names[i], Begin[i], End[i]);  
  }
  return Variables;
}

std::vector<RooDataHist*> BaseFunctions::RooData(std::vector<TH1F*> Hist, RooRealVar* Domain)
{
  std::vector<RooDataHist*> DataHist(Hist.size());
  for (int i(0); i < Hist.size(); i++)
  {
    TString name = Hist[i] -> GetTitle(); 
    DataHist[i] = new RooDataHist(name, name, *Domain, Hist[i]);
  }
  return DataHist;
}

RooDataHist* BaseFunctions::RooData(TH1F* Hist, RooRealVar* Domain)
{
  std::vector<TH1F*> H = {Hist};
  std::vector<RooDataHist*> D = RooData(H, Domain);
  return D[0];
}

std::vector<RooHistPdf*> BaseFunctions::RooPDF(std::vector<TH1F*> Hist, RooRealVar* Domain)
{
  std::vector<RooDataHist*> Data = RooData(Hist, Domain);
  std::vector<RooHistPdf*> HistPdf(Hist.size());
  for (int i(0); i < Data.size(); i++)
  {
    TString name = Hist[i] -> GetTitle(); 
    HistPdf[i] = new RooHistPdf(name, name, *Domain, *Data[i]);
  }
  return HistPdf; 
}

RooArgList BaseFunctions::RooList(std::vector<RooRealVar*> Vector)
{
  RooArgList L;
  for (RooRealVar* a : Vector){L.add(*a);}
  return L; 
}

RooArgList BaseFunctions::RooList(std::vector<RooHistPdf*> Vector)
{
  RooArgList L;
  for (RooHistPdf* a : Vector){L.add(*a);}
  return L; 
}

float BaseFunctions::ChiSquare(std::vector<float> V1, std::vector<float> V2)
{
  float di(0);
  for (int i(0); i < V1.size(); i++)
  {
    di += pow(V1[i] - V2[i], 2);
  }
  return di;
}

void BaseFunctions::PredictionTruthPrint(std::vector<float> Truth, std::vector<float> Prediction)
{
  for (int i(0); i < Truth.size(); i++)
  {
    std::cout << "trk-" << i+1 << " Truth: " << Truth[i] << " Prediction: " << Prediction[i] << std::endl;
  }

}
