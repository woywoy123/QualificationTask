#include<PostAnalysis/Functions.h>
#include<PostAnalysis/UnitClosures.h>

void UnitClosures::TestFit(std::vector<TH1F*> PDF, std::vector<TH1F*> Data, float min, float max, std::vector<std::vector<float>> Closure)
{
  Fit_Functions f;
  Benchmark B;
   
  // Get each of the trks  
  TH1F* trk1 = Data[0];
  TH1F* trk2 = Data[1];
  TH1F* trk3 = Data[2];
  TH1F* trk4 = Data[3];

  // Get the closure constants
  std::vector<float> t1_C = Closure[0];
  std::vector<float> t2_C = Closure[1];
  std::vector<float> t3_C = Closure[2];
  std::vector<float> t4_C = Closure[3];

  // Fit them to data
  std::vector<RooRealVar*> vars1 = f.FitPDFtoData(PDF, trk1, min, max);
  std::vector<RooRealVar*> vars2 = f.FitPDFtoData(PDF, trk2, min, max);
  std::vector<RooRealVar*> vars3 = f.FitPDFtoData(PDF, trk3, min, max);
  std::vector<RooRealVar*> vars4 = f.FitPDFtoData(PDF, trk4, min, max);

  // Get the fractions 
  std::vector<float> v1 = f.Fractionalizer(vars1, trk1);
  std::vector<float> v2 = f.Fractionalizer(vars2, trk2);
  std::vector<float> v3 = f.Fractionalizer(vars3, trk3);
  std::vector<float> v4 = f.Fractionalizer(vars4, trk4);
  std::vector<std::vector<float>> v = {v1, v2, v3, v4};

  // Comparing the difference between the prediction and closure
  float s1 = B.PythagoreanDistance(v1, t1_C);
  float s2 = B.PythagoreanDistance(v2, t2_C);
  float s3 = B.PythagoreanDistance(v3, t3_C);
  float s4 = B.PythagoreanDistance(v4, t4_C);
  std::vector<float> s = {s1, s2, s3, s4};

  for (int i(0); i < Data.size(); i++)
  {
    std::cout << "###################### TRK: " << i+1 << " ################### " << std::endl;
    std::cout << "Prediction: " << v[i][0] << " "  << v[i][1] << " " << v[i][2] << " " << v[i][3] << std::endl;
    std::cout << "Truth: " << Closure[i][0] << " "  << Closure[i][1] << " " << Closure[i][2] << " " << Closure[i][3] << std::endl;
    std::cout << "Error: " << s[i] << std::endl;
  }
}

void UnitClosures::TestTailAndDeconv(TH1F* trk1, TH1F* trk2, int iter, float min, float max)
{
  Fit_Functions f;
  Functions F;

  std::vector<float> Data_V = f.TH1FDataVector(trk2, 0.1);
  std::vector<float> deconv(Data_V.size(), 0.5);

  TCanvas* can = new TCanvas();
  can -> SetLogy();
  can -> Update();

  TH1F* Target = new TH1F("Target1", "Target1", deconv.size(), min, max); 
  TH1F* Closure = new TH1F("Closure1", "Closure1", deconv.size(), min, max);
  F.ExpandTH1F(trk1, Closure);
  Closure -> SetLineStyle(kDashed);
  Closure -> Draw("SAMEHIST");  
  
  Target -> SetLineColor(kRed);
  Target -> GetYaxis() -> SetRangeUser(1e-2, 1e6);
  
  for (int i(0); i < iter; i++)
  {
    deconv = f.LRDeconvolution(Data_V, deconv, deconv, 0.75);
    deconv = f.TailReplaceClosure(trk1, deconv);
    F.VectorToTH1F(deconv, Target);
    Target -> Draw("SAMEHIST");
    can -> Update();
  }
}

void UnitClosures::TestDeconvolution(TH1F* h1, TH1F* PSF, int iter)
{
  Fit_Functions f;
  Functions F;
 
  TH1F* Conv = (TH1F*)h1 -> Clone("Convolved");
  Conv -> Reset();
  Conv -> SetLineStyle(kDotted); 

  TH1F* Deconv = (TH1F*)h1 -> Clone("Deconv");
  Deconv -> Reset();
  Deconv -> SetLineColor(kRed);
  Deconv -> GetYaxis() -> SetRangeUser(1e-2, h1 -> Integral());

  h1 -> SetLineColor(kBlack);
  h1 -> SetLineStyle(kDotted);

  PSF -> SetLineColor(kGreen);

  f.ConvolveHists(h1, PSF, Conv);
  f.ArtifactRemove(Conv, "b");
  f.Normalizer(Conv);
  Conv -> Scale(PSF -> Integral());
  
  TCanvas* can = new TCanvas();
  can -> SetLogy();

  std::vector<float> H1 = f.TH1FDataVector(Conv, 0.1);
  std::vector<float> PSF_V = f.TH1FDataVector(PSF, 0.1);
  std::vector<float> deconv(H1.size(), 0.5);

  Conv -> Draw("SAMEHIST");
   
  for (int i(0); i < iter; i++)
  {
    deconv = f.LRDeconvolution(H1, PSF_V, deconv, 0.75);
    F.VectorToTH1F(deconv, Deconv);
    f.Normalizer(Deconv);
    Deconv -> Scale(PSF -> Integral());
    Deconv -> Draw("SAMEHIST");
    PSF -> Draw("SAMEHIST");
    h1 -> Draw("SAMEHIST");
    can -> Update();
  }
  
}

void UnitClosures::TestSubtraction(TH1F* Data, int trk, std::vector<TH1F*> PDFs, float min, float max, std::vector<float> Closure)
{
  Fit_Functions f;

  // ======= Closure plot 
  TCanvas* can = new TCanvas();
  gStyle -> SetOptStat(0);
  can -> Divide(2,1);
  can -> cd(1) -> SetLogy();
  can -> cd(1);
  
  TH1F* Data2 = (TH1F*)Data -> Clone("Data"); 
  TH1F* Clo = (TH1F*)Data -> Clone("Closure Figure");
  Clo -> SetLineColor(kBlack);
  Clo -> Draw("SAMEHIST*");
   
  TLegend* leg = new TLegend(0.9, 0.9, 0.75, 0.75);
  leg -> AddEntry(Clo, "Data");
   
  for (int i(0); i < PDFs.size(); i++)
  {
    f.Normalizer(PDFs[i]);
    PDFs[i] -> SetLineColor(Constants::Colors[i]);
    PDFs[i] -> Scale(Closure[i]*(Data -> Integral())); 
    PDFs[i] -> Draw("SAMEHIST"); 
    leg -> AddEntry(PDFs[i], PDFs[i] -> GetTitle());
  }
  leg -> Draw("SAME");

  // ======== Prediction Plot   
  std::vector<RooRealVar*> vars = f.FitPDFtoData(PDFs, Data, min, max);
  f.Subtraction(PDFs, Data2, trk, vars); 
 
  can -> cd(2) -> SetLogy();
  can -> cd(2);
  Data2 -> Draw("SAMEHIST"); 
  PDFs[trk -1] -> Draw("SAMEHIST*");
}

void Presentation::TestMinimalAlgorithm(std::vector<TH1F*> Data, float min, float max, float offset, std::vector<TH1F*> Pure, std::vector<std::vector<float>> Closure)
{
  Algorithms A;
  Functions F;
  Fit_Functions f; 
  Benchmark B; 

  // Define the histograms individually
  TH1F* trk1_Clone = (TH1F*)Data[0] -> Clone("1-Track");
  TH1F* trk2_Clone = (TH1F*)Data[1] -> Clone("2-Track");

  // Make PDF output histograms
  std::vector<TString> PDFNames = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> PDFs = F.MakeTH1F(PDFNames, trk1_Clone -> GetNbinsX(), min, max, "_PDF");

  for (int i(0); i < 10; i++)
  {
    A.MinimalAlgorithm(trk1_Clone, trk2_Clone, PDFs, min, max, offset, 25);
  }
  
  std::vector<std::vector<float>> Prediction;
  for (int i(0); i < PDFs.size(); i++)
  { 
    std::vector<RooRealVar*> vars = f.FitPDFtoData(PDFs, Data[i], min, max);
    std::vector<float> v = f.Fractionalizer(vars, Data[i]); 
    Prediction.push_back(v);
  }
  
  // ====== Plotting ====== //
  // === Truth Canvas 
  TCanvas* Truth = B.ClosurePlot("Truth", Data, Pure, Closure);
   
  // === Algorithm Canvas
  TCanvas* Algorithm = B.ClosurePlot("Algorithm", Data, PDFs, Prediction);

  Truth -> Draw();
  Algorithm -> Draw();
}

void Presentation::TestGaussianAlgorithm(std::vector<TH1F*> Data, float min, float max, float offset, std::vector<TH1F*> Pure, std::vector<std::vector<float>> Closure)
{
  Algorithms A;
  Functions F;
  Fit_Functions f; 
  Benchmark B; 

  // Define the histograms individually
  TH1F* trk1_Clone = (TH1F*)Data[0] -> Clone("1-Track");
  TH1F* trk2_Clone = (TH1F*)Data[1] -> Clone("2-Track");

  // Make PDF output histograms
  std::vector<TString> PDFNames = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> PDFs = F.MakeTH1F(PDFNames, trk1_Clone -> GetNbinsX(), min, max, "_PDF");
  
  // Gaussian algorithm
  float mean_s = -1;
  float mean_e = 1; 
  float stdev_s= 0.001; 
  float stdev_e= 1;

  std::vector<float> Parameters;
  for (int x(0); x < 1; x++)
  { 
    Parameters = A.GaussianAlgorithm(trk1_Clone, trk2_Clone, PDFs, min, max, offset, mean_s, mean_e, stdev_s, stdev_e, 25);
  }
   
  std::vector<std::vector<float>> Prediction;
  for (int i(0); i < PDFs.size(); i++)
  { 
    std::vector<RooRealVar*> vars = f.FitPDFtoData(PDFs, Data[i], min, max);
    std::vector<float> v = f.Fractionalizer(vars, Data[i]); 
    Prediction.push_back(v);
  }
  
  // ====== Plotting ====== //
  // === Truth Canvas 
  TCanvas* Truth = B.ClosurePlot("Truth", Data, Pure, Closure);
   
  // === Algorithm Canvas
  TCanvas* Algorithm = B.ClosurePlot("Algorithm", Data, PDFs, Prediction);

  Truth -> Draw();
  Algorithm -> Draw();

  for (float i : Parameters)
  {
    std::cout << "### Parameter value: " << i << std::endl;
  }
}


void Presentation::Threshold(TString DataDir)
{
  TFile* File = new TFile(DataDir);
  if (!File -> IsOpen()){ std::cout << "Failed to open file" << std::endl; };
 
  // I am using the 1-trk in hist for the 1 track distribution because in the paper this is how 
  // the 1 track template is sampled. The 2 track templates are generated within the jet code  
  int bins = 300;
  float min = -0.5;
  float max = 14.5;

  std::vector<TString> energies = {"/200_up_GeV", "/200_400_GeV", 
                                   "/400_600_GeV", "/600_800_GeV", 
                                   "/800_1000_GeV", "/1000_1200_GeV", 
                                   "/1200_1400_GeV", "/1400_1600_GeV", 
                                   "/1600_1800_GeV", "/1800_2000_GeV", 
                                   "/2000_2200_GeV", "/2200_2400_GeV", 
                                   "/2400_2600_GeV", "/2600_2800_GeV", 
                                   "/2800_3000_GeV", "/higher_GeV"}; 

  std::vector<std::vector<TString>> Batch = {{"/200_up_GeV"}, {"/200_400_GeV","/400_600_GeV"}, {"/600_800_GeV","/800_1000_GeV", "/1000_1200_GeV"}, {"/1200_1400_GeV", "/1400_1600_GeV","/1600_1800_GeV"}, {"/1800_2000_GeV","/2000_2200_GeV", "/2200_2400_GeV","/2400_2600_GeV", "/2600_2800_GeV","/2800_3000_GeV", "/higher_GeV"}};

  std::vector<TString> Titles = { "<200", "200-600", "600-1200", "12000-1800", "1800+"};
  std::vector<TString> Histograms = {"dEdx_out_ntrk1_calib", "dEdx_in_ntrk2_calib"};
  std::vector<TString> Layers = Constants::Detector;
   
  Functions F;
  Fit_Functions f;
  std::vector<TH1F*> trk2_Hists = F.MakeTH1F(Titles, bins, min, max, "_2trk");
  std::vector<TH1F*> trk1_Hists = F.MakeTH1F(Titles, bins, min, max, "_1trk");
  std::vector<TH1F*> trk2_trk1_C = F.MakeTH1F(Titles, bins, min, max, "_Convolved");

  TCanvas* can = new TCanvas();
  can -> Divide(2,1);
  can -> cd(1) -> SetLogy();
  can -> cd(2) -> SetLogy();
  gStyle -> SetOptStat(0);
 
  TLegend* leg1 = new TLegend(0.9, 0.9, 0.6, 0.75);
  TLegend* leg2 = new TLegend(0.9, 0.9, 0.6, 0.75);
  leg1 -> SetTextSize(.01);
  leg2 -> SetTextSize(.01);
    
  for (TString layer : Layers)
  {
    for (int i(0); i < Batch.size(); i++)
    {
      std::vector<TString> B = Batch[i];

      for (int x(0); x < B.size(); x++)
      {
           TString Dir = B[x];
           File -> cd(layer + Dir);
           trk2_Hists[i] -> Add((TH1F*)gDirectory -> Get(Histograms[1]));
           trk1_Hists[i] -> Add((TH1F*)gDirectory -> Get(Histograms[0]));        
           File -> cd(); 
      }  
      f.ConvolveHists(trk1_Hists[i], trk1_Hists[i], trk2_trk1_C[i]);
      f.Normalizer(trk1_Hists[i]);
      f.Normalizer(trk2_Hists[i]);  
      f.Normalizer(trk2_trk1_C[i]);  
      f.ArtifactRemove(trk2_trk1_C[i], "b");
      trk2_Hists[i] -> GetXaxis() -> SetTitle("dEdx (MeV g^{-1} cm^{2})");
      trk1_Hists[i] -> GetXaxis() -> SetTitle("dEdx (MeV g^{-1} cm^{2})");
      trk2_trk1_C[i] -> GetXaxis() -> SetTitle("dEdx (MeV g^{-1} cm^{2})");
      trk1_Hists[i] -> SetLineColor(Constants::Colors[i]);
      trk2_Hists[i] -> SetLineColor(Constants::Colors[i]);    
      trk2_trk1_C[i] -> SetLineColor(Constants::Colors[i]);   
      trk2_trk1_C[i] -> SetLineStyle(kDashed);
      leg1 -> AddEntry(trk2_trk1_C[i], trk2_trk1_C[i] -> GetTitle());
      leg1 -> AddEntry(trk2_Hists[i], trk2_Hists[i] -> GetTitle());
      leg2 -> AddEntry(trk1_Hists[i], trk1_Hists[i] -> GetTitle());
      
      can -> cd(1);
      trk1_Hists[i] -> Draw("SAMEHIST");
      leg2 -> Draw("SAME");

      can -> cd(2);
      trk2_Hists[i] -> Draw("SAMEHIST");
      trk2_trk1_C[i] -> Draw("SAMEHIST");
      leg2 -> Draw("SAME"); 
    } 
  }
  can -> Update();
  can -> Print("Threshold.png");
}


