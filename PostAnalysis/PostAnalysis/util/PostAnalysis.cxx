// Including commonly reused functions
#include<PostAnalysis/Functions.h>
#include<PostAnalysis/Verification.h>

// Including standard C++ libraries 
#include<iostream>

using namespace RooFit;

void PostAnalysis()
{
  // Init the function that hosts common tools
  Functions F;
  Fit_Functions Fit;
  Plot_Functions plot;
  Verification V;

  // Check that the file is there  
  TString dir = "../PostAnalysisData/Merged.root";
  TFile *File = new TFile(dir);
  if (!File -> IsOpen()){ gSystem -> Exit(0); }

  // Create the required Histograms  
  std::vector<TH1F*> Pure_Hists = F.MakeTH1F(Constants::Pure_Names, 50, 0, 20); 
  std::vector<TH1F*> trk_1 = F.MakeTH1F(Constants::trk_1, 50, 0, 20); 
  std::vector<TH1F*> trk_2 = F.MakeTH1F(Constants::trk_2, 50, 0, 20); 
  std::vector<TH1F*> trk_3 = F.MakeTH1F(Constants::trk_3, 50, 0, 20); 
  std::vector<TH1F*> trk_4 = F.MakeTH1F(Constants::trk_4, 50, 0, 20); 

  for (TString layer : Constants::Detector)
  {
    F.FillTH1F_From_File(Pure_Hists, File, layer); 
    F.FillTH1F_From_File(trk_1, File, layer); 
    F.FillTH1F_From_File(trk_2, File, layer);
    F.FillTH1F_From_File(trk_3, File, layer);
    F.FillTH1F_From_File(trk_4, File, layer);
  }

  // Build PDFs and all the histograms needed for the fitting  
  RooRealVar* dEdx_range = new RooRealVar("dEdx_range", "dEdx_range", 0.4, 19.6);
  std::vector<RooHistPdf*> Pure_PDF = Fit.ConvertTH1FtoPDF(Pure_Hists, dEdx_range);
  std::vector<RooRealVar*> Variables = Fit.GenerateVariables(Constants::Variable_Names, Constants::Begin, Constants::End);
  
  RooArgList Vars = Fit.VectorToArgList(Variables);
  RooArgList PDFs = Fit.VectorToArgList(Pure_PDF);
  RooAddPdf model("model","model", PDFs, Vars);

  //V.RecoverScaling(model, Pure_Hists, Pure_PDF, dEdx_range, Variables, 0.1, 0.001);
  //V.Subtraction(Pure_Hists, dEdx_range, model, Pure_PDF, Variables); 
  V.Reconstruction(trk_1, trk_2, trk_3, trk_4, Variables, Pure_PDF, model, dEdx_range);

   
 
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
