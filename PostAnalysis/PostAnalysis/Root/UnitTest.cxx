#include<PostAnalysis/UnitTest.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Constants.h>

void BaseFunctionTest::Subtraction()
{
    DistributionGenerators D; 
    BaseFunctions B;
    Plotting P; 

    // Example parameters
    std::vector<float> LDau = {1, 0.9, 0.1};
    std::vector<float> ratio = {0.5, 0.5, 0.5, 0.5};
    float min = 0; 
    float max = 20; 
    int bins = 500; 
    int npts = 500000;
    std::vector<TString> trk_2 = {"dEdx_ntrk_2_ntru_1", "dEdx_ntrk_2_ntru_2", "dEdx_ntrk_2_ntru_3", "dEdx_ntrk_2_ntru_4"};

    // Make Hists
    std::vector<TH1F*> trk2_N = B.MakeTH1F(trk_2, bins, min, max); 
    TH1F* Data = new TH1F("Data", "Data", bins, min, max);      
    D.Landau(trk2_N, ratio, LDau, npts, min, max);
    
    // Get Closure values and fill data 
    std::vector<float> CLS2 = B.ClosureAndData(trk2_N, Data); 
    
    B.Subtraction(trk2_N, Data, 2, CLS2);
    
    TCanvas* can = P.PlotHists(trk2_N[1], Data); 
    can -> Draw(); 
}

void BaseFunctionTest::NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max )
{
  BaseFunctions B; 

  B.Normalize(Hists);
  RooRealVar* x = new RooRealVar("x", "x", min, max); 
 
  std::vector<float> Begin(Constants::trk_2.size(), 0); 
  std::vector<float> End(Constants::trk_2.size(), Data -> Integral()); 
  
  std::vector<RooRealVar*> Var = B.RooVariables(Constants::trk_2, Begin, End);
  std::vector<RooHistPdf*> Hist = B.RooPDF(Hists, x);
  RooDataHist* data = B.RooData(Data, x);

  RooAddPdf model("Model", "Model", B.RooList(Hist), B.RooList(Var));
  model.fitTo(*data);

  Plotting P;
  TCanvas* can = P.PlotHists(model, x, Hist, data);
  can -> Draw(); 
}

void BaseFunctionTest::Convolve(TH1F* Hist1, TH1F* Hist2, TH1F* Expectation)
{
  BaseFunctions B; 
  Plotting P;
   
  TH1F* Conv = (TH1F*)Expectation -> Clone("Conv"); 
  B.ConvolveHists(Hist1, Hist2, Conv);
  B.Normalize(Conv);
  Conv -> Scale(Expectation -> Integral());
  TCanvas* can = P.PlotHists(Conv, Expectation);
  can -> Draw();
}

void BaseFunctionTest::Deconvolve(TH1F* Trk2, TH1F* Trk1, float offset, int iter)
{
  BaseFunctions B;
  Plotting P;
   
  std::vector<float> Data_V = B.TH1FDataVector(Trk2, offset);
  std::vector<float> deconv(Data_V.size(), 0.5); 
  for (int i(0); i < iter; i++)
  {
    deconv = B.LucyRichardson(Data_V, deconv, deconv, 1);  
  }  
  
  TH1F* Hist = new TH1F("Hist", "Hist", 500, 0, 20); 
  B.ToTH1F(deconv, Hist);
 
  P.PlotHists(Hist, Trk1); 
}

void DerivedFunctionTest::NormalFit(std::vector<TH1F*> Hists, TH1F* Data, std::vector<float> CL, float min, float max)
{
  DerivedFunctions DF; 
  BaseFunctions B; 
  std::vector<RooRealVar*> Vars = DF.FitToData(Hists, Data, min, max);
  std::vector<float> pred = B.Ratio(Vars, Data);  
  float distance = B.ChiSquare(pred, CL);
  std::cout << distance << std::endl;
  B.PredictionTruthPrint(CL, pred); 
}


void DerivedFunctionTest::ShiftTest(TH1F* H1, int Shift)
{
  BaseFunctions B;
  DerivedFunctions D;
  Plotting P;

  int bins = H1 -> GetNbinsX();
  int Padding = bins/2;
 
  TH1F* Hist = new TH1F("Test", "Test", bins, 0, 20); //2*bins, -Padding, bins + Padding);
  B.ShiftExpandTH1F(H1, Hist, Shift);
  //D.RooShift(Hist, H1);
  int shift = D.NumericalShift(H1, Hist);
  std::cout << " The histograms are shifted by: " << shift << std::endl;
}

void DerivedFunctionTest::ReplaceShiftTail(TH1F* Source, TH1F* Target, int Shift)
{
  BaseFunctions B;
  DerivedFunctions D;
  Plotting P; 
 
  // Create a shifted version of the Target 
  TH1F* Copy = (TH1F*)Target -> Clone("Copy");  
  Copy -> Reset();
  Copy -> SetTitle("Copy");
  B.ShiftExpandTH1F(Target, Copy, Shift); 
  D.ReplaceShiftTail(Source, Copy);

  TCanvas* can = P.PlotHists({Source, Target}, Copy);
}

void DerivedFunctionTest::DeconvolveReconvolve(std::vector<TH1F*> ntrk, float offset, int iter)
{
  DerivedFunctions DF;
  Plotting P; 

  std::vector<TH1F*> PDFs = DF.nTRKGenerator(ntrk[0], ntrk[1], offset, iter);

  for (int i(0); i < ntrk.size(); i++)
  {
    PDFs[i] -> Scale(ntrk[i] -> Integral()); 
  }
  
  P.PlotHists(PDFs, ntrk); 
}

void DerivedFunctionTest::DeconvolveGaussianFit(TH1F* trk1, TH1F* trk2,  float mean, float stdev, float offset, int iter)
{
  DerivedFunctions DF; 
  BaseFunctions B; 
  Plotting P; 
 
  std::vector<TH1F*> PDFs = DF.nTRKGenerator(trk1, trk2, offset, iter);
  std::map<TString, float> Parameters = DF.FitGaussian(trk2, PDFs, mean, stdev, -1, 1, 0.01, 1, offset, iter);  
  std::vector<TString> Names = {"n_trk1", "n_trk2", "n_trk3", "n_trk4"};
  std::vector<TString> Stdev = {"s1", "s2", "s3", "s4"};
  std::vector<TString> Mean = {"m1", "m2", "m3", "m4"};

  for (int i(0); i < Names.size(); i++)
  {
    TH1F* H = DF.GaussianConvolve(PDFs[i], Parameters[Mean[i]], Parameters[Stdev[i]]);
    PDFs[i] -> Reset();
    B.ShiftExpandTH1F(H, PDFs[i]);
   
    float e = Parameters[Names[i]];  
    PDFs[i] -> Scale(e);   
  }
 
  P.PlotHists(PDFs, trk2);   
}

void DerivedFunctionTest::MainAlgorithm(std::vector<TH1F*> ntrk, std::vector<float> Params, float offset, int iter, int cor_loop, float Gamma, std::vector<std::vector<TH1F*>> Closure)
{
  DerivedFunctions DF; 
  BaseFunctions B;
   
  std::map<TH1F*, std::vector<TH1F*>> PDFs = DF.MainAlgorithm(ntrk, Params, offset, Gamma, iter, cor_loop);

  std::vector<TH1F*> truth;
  std::vector<std::vector<TH1F*>> result;

  int i = 0; 
  std::vector<TString> Names = {"trk1_pure", "trk2_pure", "trk3_pure", "trk4_pure"};
  std::vector<TString> PDF_Names = {"trk1_pdf", "trk2_pdf", "trk3_pdf", "trk4_pdf"};
  for (auto &x : PDFs)
  {
    std::vector<TH1F*> Hists = x.second;
    std::vector<TH1F*> Cl_Hists = Closure[i]; 
    
    TH1F* t = x.first;
    TH1F* H2 = Cl_Hists[i]; 

    TH1F* H1 = (TH1F*)H2 -> Clone(Names[i]);
    H1 -> Reset();
    H1 -> SetTitle(Names[i]);
    
    TH1F* H3 = (TH1F*)H2 -> Clone(PDF_Names[i]);
    H3 -> Reset();
    H3 -> SetTitle(PDF_Names[i]);   
    
    B.ShiftExpandTH1F(t, H1);
    B.ShiftExpandTH1F(Hists[i], H3);

    truth.push_back(H2);
    std::vector<TH1F*> Out;
    Out.push_back(H1);
    Out.push_back(H3);

    result.push_back(Out);
    i++;
  }

  Plotting P;
  TCanvas* can = P.PlotHists(result, truth);
}

























