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
  RooRealVar* dEdx_range = new RooRealVar("dEdx_range", "dEdx_range", 0, 20);
  std::vector<RooHistPdf*> Pure_PDF = Fit.ConvertTH1FtoPDF(Pure_Hists, dEdx_range);
  std::vector<RooRealVar*> Variables = Fit.GenerateVariables(Constants::Variable_Names, Constants::Begin, Constants::End);
  
  RooArgList Vars = Fit.VectorToArgList(Variables);
  RooArgList PDFs = Fit.VectorToArgList(Pure_PDF);

  RooAddPdf model("model","model", PDFs, Vars);

  // ================================================== Fit Verification steps ========================================================== // 
  // 1. For Verification of scaling and correct implementation of RooFit functions
  // === 1.1: Create Scaled histograms and sum them to produce toy data sample 
  TH1F* Data_TH = (TH1F*)Pure_Hists.at(0) -> Clone("hnew"); 
  Data_TH -> Scale(0.001);
  TH1F* Data_TH_2 = (TH1F*)Pure_Hists.at(1) -> Clone("hnew");
  Data_TH_2 -> Scale(0.1);
  Data_TH -> Add(Data_TH_2);

  // === 1.2: Convert to RooDataHist and fit 
  RooDataHist* Data = Fit.ConvertTH1toDataHist(Data_TH, dEdx_range);
//  model.fitTo(*Data);
  
  // === 1.3: Plot the output
//  TCanvas* Plot = plot.GeneratePlot("Addition of two scales histograms", dEdx_range, Data, model, Pure_PDF, Constants::Pure_Names);
//  Plot -> Draw(); 

  // === 1.4: Calculate the Scalings
  float trk1_e = Pure_Hists.at(0) -> GetEntries(); 
  float trk2_e = Pure_Hists.at(1) -> GetEntries();
  float Pre1_e = Variables.at(0) -> getVal(); 
  float Pre2_e = Variables.at(1) -> getVal();

  std::cout << "Scale of the 1-Track: " << Pre1_e / trk1_e << std::endl;
  std::cout << "Scale of the 2-Track: " << Pre2_e / trk2_e << std::endl;

  // 2. Stack the pure-ntrk templates and subtract at each stage
  // === 2.1: Create the stacked histograms for each n-track
  TH1F* Unscaled_Stack = new TH1F("Unscaled", "Unscaled", 50, 0, 20);
  Unscaled_Stack -> Add(Pure_Hists.at(0));
  Unscaled_Stack -> Add(Pure_Hists.at(1)); 
  Unscaled_Stack -> Add(Pure_Hists.at(2));
  Unscaled_Stack -> Add(Pure_Hists.at(3)); 
  RooDataHist* Data_Stack = Fit.ConvertTH1toDataHist(Unscaled_Stack, dEdx_range); 
  
  // Get the entries in the pure histograms  
  float trk1_mc = Pure_Hists.at(0) -> GetEntries();
  float trk2_mc = Pure_Hists.at(1) -> GetEntries();
  float trk3_mc = Pure_Hists.at(2) -> GetEntries();
  float trk4_mc = Pure_Hists.at(3) -> GetEntries();
  
  // === 2.2: Fit the model to the data. 
  model.fitTo(*Data_Stack);
  TCanvas* Plot_1 = plot.GeneratePlot("Histogram subtraction", dEdx_range, Data_Stack, model, Pure_PDF, Constants::Pure_Names);

  // Get the predicted values
  float trk1_fit_1 = Variables.at(0) -> getVal();
  float trk2_fit_1 = Variables.at(1) -> getVal();
  float trk3_fit_1 = Variables.at(2) -> getVal();
  float trk4_fit_1 = Variables.at(3) -> getVal();

  // === 2.3: Subtract the 4-trk histogram from data and fit  
  Unscaled_Stack -> Add(Pure_Hists.at(3), -1);
  RooDataHist* Data_123_H = Fit.ConvertTH1toDataHist(Unscaled_Stack, dEdx_range);
  model.fitTo(*Data_123_H);
  TCanvas* Plot_2 = plot.GeneratePlot("Histogram subtraction", dEdx_range, Data_123_H, model, Pure_PDF, Constants::Pure_Names);

  // Get the prediction values
  float trk1_fit_2 = Variables.at(0) -> getVal();
  float trk2_fit_2 = Variables.at(1) -> getVal();
  float trk3_fit_2 = Variables.at(2) -> getVal();
  float trk4_fit_2 = Variables.at(3) -> getVal();
 
  // === 2.5: Subtract the 3-trk histogram from data and fit
  Unscaled_Stack -> Add(Pure_Hists.at(2), -1);
  RooDataHist* Data_12_H = Fit.ConvertTH1toDataHist(Unscaled_Stack, dEdx_range);
  model.fitTo(*Data_12_H);
  TCanvas* Plot_3 = plot.GeneratePlot("Histogram subtraction", dEdx_range, Data_12_H, model, Pure_PDF, Constants::Pure_Names);

  // Get the prediction values
  float trk1_fit_3 = Variables.at(0) -> getVal();
  float trk2_fit_3 = Variables.at(1) -> getVal();
  float trk3_fit_3 = Variables.at(2) -> getVal();
  float trk4_fit_3 = Variables.at(3) -> getVal();
 

  std::cout << "========================== No subtraction results ===============================" << std::endl; 
  std::cout << "Fraction of 1-Track in Data from original: " << trk1_fit_1 / trk1_mc << std::endl;
  std::cout << "Fraction of 2-Track in Data from original: " << trk2_fit_1 / trk2_mc << std::endl;
  std::cout << "Fraction of 3-Track in Data from original: " << trk3_fit_1 / trk3_mc << std::endl;
  std::cout << "Fraction of 4-Track in Data from original: " << trk4_fit_1 / trk4_mc << std::endl;
 
  std::cout << "========================== 4 track subtraction  ===============================" << std::endl; 
  std::cout << "Fraction of 1-Track in Data from original: " << trk1_fit_2 / trk1_mc << std::endl;
  std::cout << "Fraction of 2-Track in Data from original: " << trk2_fit_2 / trk2_mc << std::endl;
  std::cout << "Fraction of 3-Track in Data from original: " << trk3_fit_2 / trk3_mc << std::endl;
  std::cout << "Fraction of 4-Track in Data from original: " << trk4_fit_2 / trk4_mc << std::endl;  

  std::cout << "========================== 3 track subtraction  ===============================" << std::endl; 
  std::cout << "Fraction of 1-Track in Data from original: " << trk1_fit_3 / trk1_mc << std::endl;
  std::cout << "Fraction of 2-Track in Data from original: " << trk2_fit_3 / trk2_mc << std::endl;
  std::cout << "Fraction of 3-Track in Data from original: " << trk3_fit_3 / trk3_mc << std::endl;
  std::cout << "Fraction of 4-Track in Data from original: " << trk4_fit_3 / trk4_mc << std::endl;
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
