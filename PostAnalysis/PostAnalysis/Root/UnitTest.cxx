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

void BaseFunctionTest::Constraint()
{
  BaseFunctions B; 

  std::vector<TString> Names = {"m", "s", "f"}; 
  std::vector<float> V1 = {0, 2, 0.5};
  std::vector<float> V2 = {-10, 0.1, 0.};
  std::vector<float> V3 = {10, 10, 1.};
  std::vector<RooRealVar*> Var = B.RooVariables(Names, V1, V2, V3); 
  RooRealVar* x = new RooRealVar("x", "x", -10, 10); 
  RooGaussian* Gaus = new RooGaussian("gaus", "gaus", *x, *Var[0], *Var[1]); 

  RooPolynomial poly("poly", "poly", *x); 
 
  RooAddPdf model("model", "model", RooArgSet(*Gaus, poly), *Var[2]);  
  
  RooDataSet *d = model.generate(*x, 50); 

  std::vector<Double_t> v1 = {0.8};
  std::vector<Double_t> v2 = {0.2}; 
  
  std::vector<RooGaussian*> G = B.RooVariables({"fconst"}, {Var[2]}, v1, v2); 

  RooProdPdf modelc("modelc", "modelc", RooArgSet(model, *G[0])); 

  RooFitResult *r2 = modelc.fitTo(*d, RooFit::Constrain(*Var[2]), RooFit::Save()); 


 

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
  //std::map<TString, float> Parameters = DF.FitGaussian(trk2, PDFs, mean, stdev, -1, 1, 0.01, 1, offset, iter);  
  std::vector<TString> Names = {"n_trk1", "n_trk2", "n_trk3", "n_trk4"};
  std::vector<TString> Stdev = {"s1", "s2", "s3", "s4"};
  std::vector<TString> Mean = {"m1", "m2", "m3", "m4"};

  //for (int i(0); i < Names.size(); i++)
  //{
  //  TH1F* H = DF.GaussianConvolve(PDFs[i], Parameters[Mean[i]], Parameters[Stdev[i]]);
  //  PDFs[i] -> Reset();
  //  B.ShiftExpandTH1F(H, PDFs[i]);
  // 
  //  float e = Parameters[Names[i]];  
  //  PDFs[i] -> Scale(e);   
  //}
 
  //P.PlotHists(PDFs, trk2);   
}

void DerivedFunctionTest::MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, float Gamma, std::vector<std::vector<TH1F*>> Closure)
{
  DerivedFunctions DF; 
  BaseFunctions B;
   
  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> PDFs = DF.MainAlgorithm(ntrk, Params, offset, Gamma, iter, cor_loop, Closure);

  std::vector<TH1F*> truth;
  std::vector<std::vector<TH1F*>> result;

  std::vector<TString> Names = {"trk1_pure", "trk2_pure", "trk3_pure", "trk4_pure"};
  std::vector<TString> PDF_Names = {"trk1_pdf", "trk2_pdf", "trk3_pdf", "trk4_pdf"};
  
   
  for (int i(0); i < Closure.size(); i++)
  {
    std::vector<TH1F*> Cl_Hists = Closure[i];  
    TH1F* H2 = Cl_Hists[i]; 

    std::pair<TH1F*, std::vector<TH1F*>> m = PDFs[i]; 
    TH1F* trkN = m.first; 
    std::vector<TH1F*> PDF_H = m.second;
    
    TH1F* H1 = (TH1F*)H2 -> Clone(Names[i]);
    H1 -> Reset();
    H1 -> SetTitle(Names[i]);
    
    TH1F* H3 = (TH1F*)H2 -> Clone(PDF_Names[i]);
    H3 -> Reset();
    H3 -> SetTitle(PDF_Names[i]);   
    
    B.ShiftExpandTH1F(trkN, H1);
    B.ShiftExpandTH1F(PDF_H[i], H3);
    
    std::vector<TH1F*> out;
    out.push_back(H1); 
    out.push_back(H3);
    result.push_back(out);
    truth.push_back(H2); 
  }

  Plotting P;
  TCanvas* can = P.PlotHists(result, truth);
}

void Presentation::ThresholdEffects()
{
  auto MakeHist = [](TString name, int bins, float min, float max, TString Detector, Color_t col)
  {
    TH1F* trk = new TH1F(name, name, bins, min, max); 
    trk -> SetTitle(name + Detector); 
    trk -> SetLineColor(col);
    trk -> GetXaxis() -> SetTitle("dEdx (MeV g^{-1} cm^{2})");
    return trk;
  };

  auto Fill =[](std::vector<std::vector<TString>> Batch, std::vector<TString> Data, std::vector<TString> Detector_Layer, TString Sample)
  {
    DistributionGenerators DG; 
    std::vector<std::vector<TH1F*>> trk; 
    for (TString Dat : Data)
    { 
      std::vector<TH1F*> temp;
      for (std::vector<TString> H : Batch)
      {
        temp.push_back(DG.FillTH1F(Sample, H, Dat, Detector_Layer)); 
      }
      trk.push_back(temp); 
    }  
    return trk; 
  };

  auto Combine =[](std::vector<std::vector<TH1F*>> trk, TH1F* All, std::vector<TH1F*> Sets)
  {
    for (int i(0); i < trk.size(); i++)
    {
       std::vector<TH1F*> En = trk[i];  
       for (int j(0); j < En.size(); j++)
       {
         All -> Add(En[j]); 
         Sets[i] -> Add(En[j]);
       }
    }
  };

  auto CompositionPlot = [](std::vector<TH1F*> Composition, TH1F* Data, TString Name)
  {
    Plotting P;
    std::vector<TString> Titles = { "<200", "200-600", "600-1200", "1200+"};
    for (int i(0); i < Titles.size(); i++)
    {
      Composition[i] -> SetTitle(Titles[i]);  
      
    }
   
    TCanvas* Comp = new TCanvas();
    Comp -> SetWindowSize(2400, 1200);  
    P.PlotHists(Composition, Data, Comp);
    Comp -> Print(Name);
    delete Comp;
  };

  auto Convolve = [](std::vector<TH1F*> trk1_Comps, std::vector<TH1F*> trk2_Comps, TString FileName, TString Name = "")
  {
    BaseFunctions BF; 
    Plotting P;
    std::vector<TString> Titles = { "<200", "200-600", "600-1200", "1200+"};
    
    std::vector<TH1F*> Convol;  
    for (int i(0); i < trk1_Comps.size(); i++)
    {
      TH1F* Conv = (TH1F*)trk2_Comps[i] -> Clone(Titles[i] + "-Convoluted");
      Conv -> SetTitle(Titles[i] + "-Convoluted");
      Conv -> Reset(); 
      BF.ConvolveHists(trk1_Comps[i], trk1_Comps[i], Conv); 
      BF.Normalize(Conv); 
      BF.ResidualRemove(Conv);  
      Conv -> SetAxisRange(1e-5, 0.4, "Y");
      Convol.push_back(Conv);
    }
  
    BF.Normalize(trk2_Comps);  
    TH1F* Temp = (TH1F*)trk2_Comps[0] -> Clone("Temp"); 
    Temp -> Reset(); 
    Temp -> SetTitle(Name);
    Temp -> SetAxisRange(1e-5, 0.4, "Y");
    TCanvas* can = new TCanvas();
    can -> SetWindowSize(2400, 1200); 
    can -> SetLogy(); 
    gStyle -> SetOptStat(0);
    Temp -> Draw("SAMEHIST"); 
    can -> Update();
    
    TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
    P.Populate(Convol, can, len);
    P.Populate(trk2_Comps, can, len, kSolid); 
     
    can -> Print(FileName);
    delete can;
  };


  // Import the relevant classes
  Plotting P; 
  DerivedFunctions DF; 

  // Constants
  int bins = 300; 
  float min = -0.5; 
  float max = 14.5; 
 
  std::vector<TString> Data = {Constants::Data2016, Constants::Data2017, Constants::Data2018}; 
  std::vector<TString> Detector_Layer = {"IBL"};
  std::vector<TString> E = Constants::energies;
  std::vector<std::vector<TString>> Batch = {{E[0]}, 
                                             {E[1], E[2]}, 
                                             {E[3], E[4], E[5]}, 
                                             {E[6], E[7], E[8], E[9], E[10], E[11], E[12], E[13], E[14], E[15]}};
 
  TH1F* trk1_2018 = MakeHist("Track1-2018 ",  bins, min, max, Detector_Layer[0], kRed);
  TH1F* trk1_2017 = MakeHist("Track1-2017 ",  bins, min, max, Detector_Layer[0], kBlue);
  TH1F* trk1_2016 = MakeHist("Track1-2016 ",  bins, min, max, Detector_Layer[0], kOrange);
  TH1F* trk1_All = MakeHist("Track1 ", bins, min, max, Detector_Layer[0], kBlack); 

  TH1F* trk2_2018 = MakeHist("Track2-2018 ",  bins, min, max, Detector_Layer[0], kRed);
  TH1F* trk2_2017 = MakeHist("Track2-2017 ",  bins, min, max, Detector_Layer[0], kBlue);
  TH1F* trk2_2016 = MakeHist("Track2-2016 ",  bins, min, max, Detector_Layer[0], kOrange);
  TH1F* trk2_All = MakeHist("Track2 ", bins, min, max, Detector_Layer[0], kBlack); 
 
  std::vector<TH1F*> Sets_1 = {trk1_2016, trk1_2017, trk1_2018};  
  std::vector<TH1F*> Sets_2 = {trk2_2016, trk2_2017, trk2_2018}; 

  std::vector<std::vector<TH1F*>> trk1 = Fill(Batch, Data, Detector_Layer, "dEdx_out_ntrk1_calib"); 
  std::vector<std::vector<TH1F*>> trk2 = Fill(Batch, Data, Detector_Layer, "dEdx_in_ntrk2_calib"); 

  Combine(trk1, trk1_All, Sets_1); 
  Combine(trk2, trk2_All, Sets_2); 


  // Declare the different TCanvas 
  TCanvas* Data_All1 = new TCanvas();  
  Data_All1 -> SetWindowSize(2400, 1200);  
  P.PlotHists(Sets_1, trk1_All, Data_All1);
  Data_All1 -> Print("Data_All_Trk1.pdf"); 

  TCanvas* Data_All2 = new TCanvas(); 
  Data_All2 -> SetWindowSize(2400, 1200);  
  P.PlotHists(Sets_2, trk2_All, Data_All2);
  Data_All2 -> Print("Data_All_Trk2.pdf"); 
 
  CompositionPlot(trk1[0], trk1_2016, "trk1_Composition_2016.pdf");  
  CompositionPlot(trk1[1], trk1_2017, "trk1_Composition_2017.pdf");  
  CompositionPlot(trk1[2], trk1_2018, "trk1_Composition_2018.pdf");  
  
  CompositionPlot(trk2[0], trk2_2016, "trk2_Composition_2016.pdf");  
  CompositionPlot(trk2[1], trk2_2017, "trk2_Composition_2017.pdf");  
  CompositionPlot(trk2[2], trk2_2018, "trk2_Composition_2018.pdf");  
 
  // Now we want to compare the the convolution of 1 track with the 2 track for different energies.
  Convolve(trk1[0], trk2[0], "Convolved_trk_2016.pdf", "Convolved 2016"); 
  Convolve(trk1[1], trk2[1], "Convolved_trk_2017.pdf", "Convolved 2017"); 
  Convolve(trk1[2], trk2[2], "Convolved_trk_2018.pdf", "Convolved 2018");

  Convolve({trk1_All}, {trk2_All}, "Convolved_trk_All.pdf", {"Convolved 1 track Data"});

  // Note to self: Fix the titles of the canvas it says <200-convoluted 

}























