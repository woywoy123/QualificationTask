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

  std::vector<float> v1 = {0.8};
  std::vector<float> v2 = {0.2}; 
  
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
  P.SimplePlot(Hist); 
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

  std::map<TString, std::vector<float>> Params; 
  Params["Gaussian"] = {0, 0.01};
  Params["m_s"] = {-3, -1, -1, -1, -1}; 
  Params["m_e"] = {3, 1, 1, 1, 1}; 
  Params["s_s"] = {0.01, 0.01, 0.01, 0.01, 0.01};
  Params["s_e"] = {1, 1, 1, 1, 1};
  std::map<TString, float> Parameters = DF.FitGaussian(trk2, PDFs, Params, offset, iter);  
  std::vector<TString> Names = {"n_trk1", "n_trk2", "n_trk3", "n_trk4"};
  std::vector<TString> Stdev = {"s1", "s2", "s3", "s4"};
  std::vector<TString> Mean = {"m1", "m2", "m3", "m4"};

  for (int i(0); i < Names.size(); i++)
  {
    TH1F* H = DF.GaussianConvolve(PDFs[i], Params["Gaussian"][0], Params["Gaussian"][1]); 
    PDFs[i] -> Reset();
    B.ShiftExpandTH1F(H, PDFs[i]);
   
    float e = Parameters[Names[i]];  
    PDFs[i] -> Scale(e);   
    delete H;
  }
  B.Normalize(trk2); 
  P.PlotHists(PDFs, trk2);   
}

void Presentation::ReconstructNTrack()
{
  auto GausConvolute = [](TH1F* PDF, float offset, int iter, TH1F* Hist, TH1F* Gaus)
  {
    BaseFunctions B;
    int bins = PDF -> GetNbinsX();
    std::vector<float> PDF_V = B.TH1FDataVector(PDF, offset);
    TString v = PDF -> GetTitle(); v += ("_Temp");
    TH1F* Temp = new TH1F(v, v, PDF_V.size(), 0, PDF_V.size()); 
    B.ToTH1F(PDF_V, Temp);
    
    B.ShiftExpandTH1F(Temp, Hist, bins*offset/2);
    Temp -> Reset();
    
    PDF_V = B.TH1FDataVector(Hist, 0);
    Hist -> Reset();

    std::vector<float> deconv(PDF_V.size(), 0.5);
    B.ShiftExpandTH1F(Gaus, Hist);

    std::vector<float> PSF_V = B.TH1FDataVector(Hist, 0); 
    Hist -> Reset(); 

    for (int x(0); x < iter; x++) 
    {
      deconv = B.LucyRichardson(PDF_V, PSF_V, deconv);
    }
    
    B.ToTH1F(deconv, Hist); 
       
    delete Temp;    
  }; 
  
  auto Fitting =[&GausConvolute](TH1F* Gaus, TH1F* PDF, TH1F* Hist, TH1F* FitTo, int iter)
  {
    BaseFunctions B; 
    DerivedFunctions DF; 

    GausConvolute(PDF, 0.2, iter, Hist, Gaus); 
    B.Normalize(Hist); 
    B.ResidualRemove(Hist);
    
    TH1F* Data = (TH1F*)Gaus -> Clone("trk_D"); 
    Data -> Reset(); 
    B.ShiftExpandTH1F(FitTo, Data); 
    B.Normalize(Data); 
 
    RooRealVar x("x", "x", 0, 500); 
    RooRealVar m("m", "m", -10, 10); 
    RooRealVar s("s", "s", -1, 1); 
    RooRealVar n("n", "n", 0, Data -> Integral()); 
 
    for (int i(0); i < Hist -> GetNbinsX(); i++)
    {
      std::cout << Hist -> GetBinContent(i+1) << std::endl;   
      
    }
  
    Hist -> Draw("SAMEHIST"); 
    RooDataHist pdf_H("pdf_H", "pdf_H", x, Hist); 
    RooDataHist data("data", "data", x, Data); 
  
    RooGaussian g("g", "g", x, m, s); 
    RooHistPdf pdf("pdf", "pdf", x, pdf_H);  
 
    std::vector<RooFFTConvPdf*> Conv_Vars = B.RooVariables({"gxtrk"}, {&pdf}, {&g}, &x); 
 
    RooAddPdf model("model", "model", RooArgSet(*Conv_Vars[0]), RooArgSet(n)); 
    model.fitTo(data, RooFit::SumW2Error(true)); 

    RooPlot* xframe = x.frame(); 
    model.plotOn(xframe, RooFit::LineStyle(kDotted)); 
    //gPad -> SetLogy(); 
    xframe -> SetMinimum(1e-9); 
    xframe -> Draw(); 

    float m_out = m.getVal(); 
    float s_out = s.getVal(); 
    float n_out = n.getVal(); 
  
    TH1F* Out_G = DF.GaussianConvolve(PDF, m_out, n_out); 
    TString name = "gX_"; name += (PDF -> GetTitle()); 
    TH1F* Out = (TH1F*)PDF -> Clone(name); 
    Out -> Reset(); 
    B.ShiftExpandTH1F(Out_G, Out); 
    Out -> Scale(n_out); 
    return Out;  
  };

  std::vector<TString> Detector_Layer = {"IBL", "Blayer", "layer1", "layer2"};
  std::vector<TString> E = Constants::energies;
  std::vector<std::vector<TString>> Batch = {{E[0]}, 
                                             {E[1], E[2]}, 
                                             {E[3], E[4], E[5]}, 
                                             {E[6], E[7], E[8], E[9], E[10], E[11], E[12], E[13], E[14], E[15]}};
  
  TH1F* tru1_trk1 = new TH1F("trk1", "trk1", 500, 0, 20); 
  TH1F* tru2_trk2 = new TH1F("trk2", "trk2", 500, 0, 20); 
  TH1F* tru3_trk3 = new TH1F("trk3", "trk3", 500, 0, 20); 
  TH1F* tru4_trk4 = new TH1F("trk4", "trk4", 500, 0, 20); 
  
  TFile* f = new TFile(Constants::MC_dir); 
  for (TString Lay : Detector_Layer)
  {
    for (std::vector<TString> Ba : Batch)
    {
      for ( TString B : Ba )
      {
        f -> cd(Lay+B); 
        TH1F* H1 = (TH1F*)gDirectory -> Get("dEdx_ntrk_1_ntru_1"); 
        TH1F* H2 = (TH1F*)gDirectory -> Get("dEdx_ntrk_2_ntru_2"); 
        TH1F* H3 = (TH1F*)gDirectory -> Get("dEdx_ntrk_3_ntru_3"); 
        TH1F* H4 = (TH1F*)gDirectory -> Get("dEdx_ntrk_4_ntru_4"); 

        tru1_trk1 -> Add(H1); 
        tru2_trk2 -> Add(H2);
        tru3_trk3 -> Add(H3);
        tru4_trk4 -> Add(H4);
      } 
    }
  }  

  //TCanvas* can = new TCanvas(); 
  //can -> SetLogy(); 
  //gStyle -> SetOptStat(0); 
  //tru1_trk1 -> Draw("SAMEHIST");
  //tru2_trk2 -> Draw("SAMEHIST");
  //tru3_trk3 -> Draw("SAMEHIST");
  //tru4_trk4 -> Draw("SAMEHIST");
  //can -> Update(); 

  DistributionGenerators DG; 
  BaseFunctions B; 
  DerivedFunctions DF; 
  Plotting P; 

  B.Normalize(tru1_trk1); 
  B.Normalize(tru2_trk2); 
  B.Normalize(tru3_trk3); 
  B.Normalize(tru4_trk4); 

  int iter = 50;
  std::vector<TH1F*> PDFs = DF.nTRKGenerator(tru1_trk1, tru2_trk2, 0.2, iter);    
   
  //P.PlotHists({PDFs[0], PDFs[1], PDFs[2], PDFs[3]} , {tru1_trk1, tru2_trk2, tru3_trk3, tru4_trk4}); 

  // ====== Doing the deconvolution and fit.
  std::vector<TString> names;  
  for (int i(0); i < PDFs.size(); i++)
  {
    TString name = PDFs[i] -> GetTitle(); name += ("Conv"); 
    names.push_back(name); 
  }
  std::vector<TH1F*> Hists = B.MakeTH1F(names, 500, 0, 500); 
  TH1F* Gaus = (TH1F*)Hists[0] -> Clone("Gaussian"); 
  Gaus -> Reset(); 
  DG.Gaussian(10, 1, Constants::GaussianToys, Gaus);  
  B.Normalize(Gaus); 


  TH1F* g1 = Fitting(Gaus, PDFs[1], Hists[1], tru2_trk2, iter); 
  
  P.PlotHists({g1}, {tru2_trk2});








}

void Presentation::MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure)
{
  DerivedFunctions DF;  
  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> PDFs = DF.MainAlgorithm(ntrk, Params, offset, iter, cor_loop, Closure);
}

void Presentation::DataAnalysis(std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, int bins, float min, float max)
{ 
  auto Generate =[](std::map<TString, TH1F*> &Dic, std::vector<TString> Detector, std::vector<std::vector<TString>> Batch, std::vector<TString> Names, TString Dir, std::vector<TString> ntrk, int bins, float min, float max)
  {
    std::map<TString, int> Test; 
    TFile* File = new TFile(Dir); 
    for (int i(0); i < ntrk.size(); i++)
    {
      TString trk = ntrk[i] + "_All";
      Dic[trk] = new TH1F(trk, trk, bins, min, max); 
      Test[trk] = 0; 
      for (TString v : Names)
      {
        TString trk_all = trk + "_" + v;
        Dic[trk_all] = new TH1F(trk_all, trk_all, bins, min, max); 
        Test[trk_all] = 0; 
      }

      for (int v(0); v < Detector.size(); v++)
      {
        TString Lay = Detector[v]; 
        TString trk_L = ntrk[i] + "_" + Lay; 
        Dic[trk_L] = new TH1F(trk_L, trk_L, bins, min, max);
        Test[trk_L] = 0;  
        for (int x(0); x < Batch.size(); x++)
        {
          std::vector<TString> Bat = Batch[x];
          TString BN = Names[x]; 
          TString trk_bn = trk_L + "_" + BN;
          TString trk_all = trk + "_" + BN;
          Dic[trk_bn] = new TH1F(trk_bn, trk_bn, bins, min, max);  
          Test[trk_bn] = 0; 
          for (int b(0); b < Bat.size(); b++)
          {
            TString dir = Lay + Bat[b]; 
            File -> cd(dir);  
            TH1F* H = (TH1F*)gDirectory -> Get(ntrk[i]);
            Dic[trk_all] -> Add(H); 
            Dic[trk] -> Add(H); 
            Dic[trk_L] -> Add(H); 
            Dic[trk_bn] -> Add(H); 

            Test[trk_all]++;
            Test[trk]++;
            Test[trk_L]++;
            Test[trk_bn]++;

            File -> cd(); 
          }
        }
      }
    }

    // Making sure that the counting makes sense for the histogram addition 
    //std::map<TString, int>::iterator it; 
    //for (it = Test.begin(); it != Test.end(); it++)
    //{
    //  std::cout << it -> first << " :: " << it -> second << std::endl; 
    //}
  };
 
  DerivedFunctions DF;  
  
  std::vector<TString> Detector_Layer = {"IBL", "Blayer", "layer1", "layer2"};
  std::vector<TString> E = Constants::energies;
  std::vector<std::vector<TString>> Batch = {{E[0]}, 
                                             {E[1], E[2]}, 
                                             {E[3], E[4], E[5]}, 
                                             {E[6], E[7], E[8], E[9], E[10], E[11], E[12], E[13], E[14], E[15]}};
  std::vector<TString> Batch_Names = {"200", "200-600", "600-1200", "1200+"};

  // Create the histograms for both the MC and the Data. The Data wont have truth but an empty set of histograms 
  std::map<TString, TH1F*> Data_MC;
  for (int i(0); i < Constants::All_MC.size(); i++)
  {
    std::vector<TString> ntrk = Constants::All_MC[i]; 
    Generate(Data_MC, Detector_Layer, Batch, Batch_Names, Constants::MC_dir, ntrk, bins, min, max);     
  }
  Generate(Data_MC, Detector_Layer, Batch, Batch_Names, Constants::MC_dir, Constants::Data_Names, bins, min, max);    

  Detector_Layer.push_back("All");
   
  // Iterate through the dictionary and build the data and truth sets used for the iteration 
  std::map<TString, std::vector<std::vector<TH1F*>>> Closure_Dic;
  std::map<TString, std::vector<TH1F*>> Data_Dic; 


  
  for (TString L : Detector_Layer)
  { 
    for (TString BN : Batch_Names)
    {
      for (std::vector<TString> tru : Constants::All_MC)
      {
        std::vector<TH1F*> truth;
        for (TString t : tru)
        { 
          TString condition = t + "_" + L + "_" + BN;
          std::map<TString, TH1F*>::iterator it;
          for (it = Data_MC.begin(); it != Data_MC.end(); it++)
          {
            TString name = it -> second -> GetName();
            if (name == condition)
            {
              truth.push_back(it -> second);
            }
          }
        }
        Closure_Dic[L + "_" + BN].push_back(truth); 
      }

      for (TString d : Constants::Data_Names)
      {
        TString condition = d + "_" + L + "_" + BN;
        std::map<TString, TH1F*>::iterator it;
        for (it = Data_MC.begin(); it != Data_MC.end(); it++)
        {
          TString name = it -> second -> GetName();
          if (name == condition)
          {
            Data_Dic[L + "_" + BN].push_back(it -> second);
          }
        }
      }
    }

    for (TString d : Constants::Data_Names)
    {
      TString condition = d + "_" + L;
      std::map<TString, TH1F*>::iterator it;
      for (it = Data_MC.begin(); it != Data_MC.end(); it++)
      {
        TString name = it -> second -> GetName();
        if (name == condition)
        {
          Data_Dic[L].push_back(it -> second); 
        }
      } 
    }

    for (std::vector<TString> d : Constants::All_MC)
    {
      std::vector<TH1F*> truth; 
      for (TString x : d)
      {
        TString condition = x + "_" + L;
        std::map<TString, TH1F*>::iterator it;
        for (it = Data_MC.begin(); it != Data_MC.end(); it++)
        {
          TString name = it -> second -> GetName();
          if (name == condition)
          {
            truth.push_back(it -> second); 
          }
        } 
      }
      Closure_Dic[L].push_back(truth); 
    }
  }


  TString Title_Params = "GlobalSettings"; Title_Params += ("-iter-"); Title_Params += (iter); Title_Params += ("-cor_loop-"); Title_Params += (cor_loop); 
   
  std::vector<TString> Name = {"m_s", "m_e", "s_s", "s_e"}; 
  for ( TString n : Name )
  {
    std::vector<float> t = Params[n]; 
    Title_Params += (n + ":");
    for (float x : t)
    {
      double h = (int)(10000*x);
      h = h/10000;
      Title_Params += (h); Title_Params += (",");
    }
    Title_Params += (".");
  }
  Title_Params += (".root");

  std::map<TString, std::vector<TH1F*>>::iterator it; 
  int z = 1; 
  for (it = Data_Dic.begin(); it != Data_Dic.end(); it++)
  {
    TString Test = it -> first; 
    std::cout << "Current Test Being Conducted: " << Test << " :: " << z++ << "/" << Data_Dic.size() << std::endl; 
    std::vector<TH1F*> ntrk_Data = it -> second; 
    std::vector<std::vector<TH1F*>> Truth_Sets = Closure_Dic[Test];

    // Here we check if the data histograms actually have enough events 
    bool run = true;  
    for (TH1F* H : ntrk_Data)
    {
      int e = H -> GetEntries();
      if (e == 0){run = false;}
    } 
    if (run == true)
    { 
      TFile Output(Title_Params, "UPDATE"); 
      std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> res =  DF.MainAlgorithm(ntrk_Data, Params, offset, iter, cor_loop, Truth_Sets); 
      
      Output.mkdir(Test);
      for (int p(0); p < res.size(); p++)
      {
        Output.cd(Test);   
        TString n = "Iteration_"; n+=(p); 
        gDirectory -> mkdir(n); 
        gDirectory -> cd(n); 
        TH1F* FLost = res[p].first;
        FLost -> Write(); 
        std::vector<TH1F*> Others = res[p].second; 
        for (TH1F* H_R : Others)
        {
          H_R -> Write(); 
        }
        gDirectory -> cd(); 
      }
      Output.Close(); 
    }
   
    else
    {
      std::cout << "skipped " << Test << std::endl;
    }
  } 
}

void Presentation::AlgorithmPlots(TString dir, int iter)
{
  auto Make =[](TCanvas* can, std::vector<TH1F*> trk_D, std::vector<std::vector<TH1F*>> trk_P, std::vector<std::vector<TH1F*>> trk_T, TString Name)
  {
    Plotting P; 
    can -> SetWindowSize(600, 1200); 
    gStyle -> SetOptStat(0); 
    P.PlotHists(trk_T, trk_P, trk_D, can); 
    can -> Draw();
    can -> Print(Name);
  };

  auto RatioPlot =[] (TH1F* H1, TH1F* H2, TString name)
  {
    
    
    auto rp1 = new TRatioPlot(H1); 

    TH1F* R = (TH1F*)H1 -> Clone(name); 
    R -> Add(H1, -1);  
    return R;
  };



  Plotting P; 
  
  std::vector<TString> Detector_Layer = {"IBL", "Blayer", "layer1", "layer2", "All"};
  std::vector<TString> Batch_Names = {"_200", "_200-600", "_600-1200", "_1200+", ""};
  std::vector<TString> BaseName = {"trk1", "trk2", "trk3", "trk4", "trk5"}; 
  std::vector<TString> BaseName2 = {"F1", "F2", "F3", "F4", "F5", "Pure"};  
  std::vector<TString> BaseName3 = {"ntru_1_", "ntru_2_", "ntru_3_", "ntru_4_", "ntru_5_", ""}; 
  std::vector<TString> BaseName4 = {"dEdx_ntrk_1_", "dEdx_ntrk_2_", "dEdx_ntrk_3_", "dEdx_ntrk_4_", "dEdx_ntrk_5_"}; 

  std::set<TString> Names; 
  Names.insert("FLost_T.at."); 
  Names.insert("FLost_P.at."); 
  for (TString T : BaseName)
  {
    for (TString B : BaseName2)
    {
      TString n = T + "_" + B; 
      Names.insert(n);  
    }
  }
  
  std::set<TString> DirNames; 
  for (TString Layer : Detector_Layer)
  {
    for (TString Ba : Batch_Names)
    {
      TString n = Layer + Ba; 
      DirNames.insert(n); 
    }
  }

  for (TString X : DirNames)
  {
    for (TString N : BaseName4)
    {
      for (TString G : BaseName3)
      {
        Names.insert(N+G+X); 
      } 
    }
  } 

  TFile* f = new TFile(dir);
  std::vector<TString> Dirs; 
  for (TString D : DirNames)
  {
    if (f -> cd(D))
    {
      for (int i(0); i < iter; i++)
      {
        TString it = "Iteration_"; it += (i); 
        f -> cd(D + "/" + it);
        Dirs.push_back(D + "/" + it);
      }
    }  
  }

  std::map<TString, std::vector<std::vector<TH1F*>>> Container;
  for (TString n : Dirs)
  {
    int chop = n.Last(*"_");
    TString r = n(chop+1, n.Length());
    std::vector<TH1F*> Hist_Iteration;  
    for (TString x : Names)
    {
      f -> cd(n);
      TString n_h = x+r;
      TH1F* H = (TH1F*)gDirectory -> Get(n_h); 
      if (H != 0){Hist_Iteration.push_back(H);} 
    }
    chop = n.Last(*"/"); 
    r = n(0, chop); 
    Container[r].push_back(Hist_Iteration);
  }

  for (TString N : DirNames)
  {
    TCanvas* can = new TCanvas(); 
    TString FileName = N + ".pdf";
    can -> Print(FileName + "["); 

    TCanvas* can_1 = new TCanvas(); 
    TString FileName_1 = N + "_1.pdf";  
    can_1 -> Print(FileName_1 + "["); 

    TCanvas* can_2 = new TCanvas(); 
    TString FileName_2 = N + "_2.pdf";  
    can_2 -> Print(FileName_2 + "["); 

    TCanvas* can_3 = new TCanvas(); 
    TString FileName_3 = N + "_3.pdf";  
    can_3 -> Print(FileName_3 + "["); 

    TCanvas* can_4 = new TCanvas(); 
    TString FileName_4 = N + "_4.pdf";  
    can_4 -> Print(FileName_4 + "["); 

    std::vector<std::vector<TH1F*>> Iterations = Container[N]; 
    for (int i(0); i < Iterations.size(); i++)
    {
      std::vector<TH1F*> Hists = Iterations[i]; 
      std::vector<TH1F*> trk1_PDF; 
      std::vector<TH1F*> trk2_PDF; 
      std::vector<TH1F*> trk3_PDF; 
      std::vector<TH1F*> trk4_PDF; 
      std::vector<TH1F*> ntrk;  
      std::vector<TH1F*> FLost; 
      std::vector<TH1F*> ntrk_data;  
      std::vector<TH1F*> trk1_tru; 
      std::vector<TH1F*> trk2_tru; 
      std::vector<TH1F*> trk3_tru;          
      std::vector<TH1F*> trk4_tru; 
      
      for (TH1F* H : Hists)
      {
        TString title = H -> GetName(); 
        if (title.Contains("trk1_F")){ trk1_PDF.push_back(H); }
        if (title.Contains("trk2_F")){ trk2_PDF.push_back(H); }
        if (title.Contains("trk3_F")){ trk3_PDF.push_back(H); }
        if (title.Contains("trk4_F")){ trk4_PDF.push_back(H); }
        if (title.Contains("trk1_Pure")){ ntrk.push_back(H); }
        if (title.Contains("trk2_Pure")){ ntrk.push_back(H); }
        if (title.Contains("trk3_Pure")){ ntrk.push_back(H); }
        if (title.Contains("trk4_Pure")){ ntrk.push_back(H); }
        if (title.Contains("FLost")){ FLost.push_back(H); }
        if (title.Contains("dEdx_ntrk_1_"+N)){ntrk_data.push_back(H);}
        if (title.Contains("dEdx_ntrk_2_"+N)){ntrk_data.push_back(H);}
        if (title.Contains("dEdx_ntrk_3_"+N)){ntrk_data.push_back(H);}
        if (title.Contains("dEdx_ntrk_4_"+N)){ntrk_data.push_back(H);} 
        if (title.Contains("dEdx_ntrk_1_ntru")){trk1_tru.push_back(H);}
        if (title.Contains("dEdx_ntrk_2_ntru")){trk2_tru.push_back(H);}
        if (title.Contains("dEdx_ntrk_3_ntru")){trk3_tru.push_back(H);}        
        if (title.Contains("dEdx_ntrk_4_ntru")){trk4_tru.push_back(H);}  
      }

      can -> Clear(); 
      can -> Divide(2,2); 
      Make(can, ntrk_data, {trk1_PDF, trk2_PDF, trk3_PDF, trk4_PDF}, {trk1_tru, trk2_tru, trk3_tru, trk4_tru}, FileName); 


      TString Name1 = N + "trk1"; Name1 += (i); 
      TString Name2 = N + "trk2"; Name2 += (i); 
      TString Name3 = N + "trk3"; Name3 += (i); 
      TString Name4 = N + "trk4"; Name4 += (i); 

      can_1 -> Clear(); 
      can_1 -> Divide(1,2); 
      can_1 -> cd(1); 
      Make(can_1, {ntrk_data[0]}, {{trk1_PDF[0]}}, {{trk1_tru[0]}}, FileName_1); 
      can_1 -> cd(2); 
      TH1F* R_1 = RatioPlot(trk1_PDF[0], trk1_tru[0], Name1); 
      R_1 -> Draw(); 
      can_1 -> Update(); 

      can_2 -> Clear(); 
      can_2 -> Divide(1,2); 
      can_2 -> cd(1); 
      Make(can_2, {ntrk_data[1]}, {{trk2_PDF[1]}}, {{trk2_tru[1]}}, FileName_2); 
      can_2 -> cd(2); 
      TH1F* R_2 = RatioPlot(trk2_PDF[0], trk2_tru[0], Name2); 
      R_2 -> Draw(); 
      can_2 -> Update(); 

      can_3 -> Clear(); 
      can_3 -> Divide(1,2); 
      can_3 -> cd(1); 
      Make(can_3, {ntrk_data[2]}, {{trk3_PDF[2]}}, {{trk3_tru[2]}}, FileName_3); 
      can_3 -> cd(2); 
      TH1F* R_3 = RatioPlot(trk3_PDF[0], trk3_tru[0], Name3); 
      R_3 -> Draw(); 
      can_3 -> Update(); 

      can_4 -> Clear(); 
      can_4 -> Divide(1,2); 
      can_4 -> cd(1); 
      Make(can_4, {ntrk_data[3]}, {{trk4_PDF[3]}}, {{trk4_tru[3]}}, FileName_4); 
      can_4 -> cd(2); 
      TH1F* R_4 = RatioPlot(trk4_PDF[0], trk4_tru[0], Name4); 
      R_4 -> Draw(); 
      can_4 -> Update(); 

     
      delete R_1; 
      delete R_2; 
      delete R_3; 
      delete R_4; 
       
    }
    can -> Print(FileName + ")"); 
    can_1 -> Print(FileName_1 + ")"); 
    can_2 -> Print(FileName_2 + ")"); 
    can_3 -> Print(FileName_3 + ")"); 
    can_4 -> Print(FileName_4 + ")"); 
    delete can;
    delete can_1;
    delete can_2;
    delete can_3;
    delete can_4;


  }
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























