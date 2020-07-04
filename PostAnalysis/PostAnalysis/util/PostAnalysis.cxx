// Including commonly reused functions
#include<PostAnalysis/Functions.h>

// Including standard C++ libraries 
#include<iostream>

using namespace RooFit;

void PostAnalysis()
{
  // Init the function that hosts common tools
  Functions F;
  Fit_Functions Fit;
  Plot_Functions plot;

  // Check that the file is there  
  TString dir = "../PostAnalysisData/Merged.root";
  TFile *File = new TFile(dir);
  if (!File -> IsOpen()){ gSystem -> Exit(0); }

  // Here we create the histograms that are needed for the pure-templates 
  std::vector<TH1F*> Pure_Hists = F.MakeTH1F(Constants::Pure_Names, 50, 0, 20); 
  for (TString layer : Constants::Detector)
  {
    F.FillTH1F_From_File(Pure_Hists, File, layer); 
  }
 
  // Build PDFs and all the histograms needed for the fitting  
  RooRealVar* dEdx_range = new RooRealVar("dEdx_range", "dEdx_range", 0.4, 19.6);
  std::vector<RooHistPdf*> Pure_PDF = Fit.ConvertTH1FtoPDF(Pure_Hists, dEdx_range);

  std::vector<RooRealVar*> Variables = Fit.GenerateVariables(Constants::Variable_Names, Constants::Begin, Constants::End);
  RooArgList Vars = Fit.VectorToArgList(Variables);
  RooArgList PDFs = Fit.VectorToArgList(Pure_PDF);

  RooAddPdf model("model","model", PDFs, Vars);
  
  // 1. For Verification of scaling and correct implementation of RooFit functions
  // === 1.1: Create Scaled histograms and sum them to produce toy data sample 
  TH1F* Data_TH = (TH1F*)Pure_Hists.at(0) -> Clone("hnew"); 
  Data_TH -> Scale(0.001);
  Pure_Hists.at(1) -> Scale(0.1);
  Data_TH -> Add(Pure_Hists.at(1));

  // === 1.2: Convert to RooDataHist and fit 
  RooDataHist* Data = Fit.ConvertTH1FtoDataHist(Data_TH, dEdx_range);
  model.fitTo(*Data);
  
  // === 1.3: Plot the output
  TCanvas* Plot = plot.GeneratePlot("Addition of two scales histograms", dEdx_range, Data, model, Pure_PDF, Constants::Pure_Names);
  Plot -> Draw(); 

  // === 1.4: Calculate the Scalings
  float trk1_e = Pure_Hists.at(0) -> GetEntries(); 
  float trk2_e = Pure_Hists.at(1) -> GetEntries();
  float Pre1_e = Variables.at(0) -> getVal(); 
  float Pre2_e = Variables.at(1) -> getVal();

  std::cout << "Scale of the 1-Track: " << Pre1_e / trk1_e << std::endl;
  std::cout << "Scale of the 2-Track: " << Pre2_e / trk2_e << std::endl;
}

void StandaloneApplications(int argc, char**argv)
{
  PostAnalysis();
}

int main(int argc, char** argv)
{
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplications(app.Argc(), app.Argv());
  app.Run();

  return 0;
}
