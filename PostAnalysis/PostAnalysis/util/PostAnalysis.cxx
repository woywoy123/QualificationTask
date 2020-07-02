// Including commonly reused functions
#include<PostAnalysis/Functions.h>

// Including standard C++ libraries 
#include<iostream>

// Include ROOT Lib
#include<TCanvas.h>

void PostAnalysis()
{
  // Init the function that hosts common tools
  Functions F;
  Fitting Fit;

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

  RooRealVar* dEdx_range = new RooRealVar("dEdx_range", "dEdx_range", 0, 20);
  std::vector<RooHistPdf*> Pure_PDF = Fit.ConvertTH1FtoPDF(Pure_Hists, dEdx_range);
  RooDataHist* Data = Fit.ConvertTH1FtoDataHist(Pure_Hists.at(0), dEdx_range);
  std::vector<RooRealVar*> Variables = Fit.GenerateVariables(Constants::Variable_Names, Constants::Begin, Constants::End);
  RooArgList Vars = Fit.VectorToArgList(Variables);
  RooArgList PDFs = Fit.VectorToArgList(Pure_PDF);

  RooAddPdf model("model","model", PDFs, Vars);
  model.fitTo(*Data);

  
  
  
   
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  


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
