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
  TString Energy_trk1 = "/200_up_GeV";
  TString Energy_trk2 = "/600_800_GeV";
  std::vector<TString> Histograms = {"dEdx_out_ntrk1_calib", "dEdx_in_ntrk2_calib"};
  std::vector<TString> Layers = Constants::Detector;
   
  TH1F* trk1_D = new TH1F("trk1_D", "trk1_D", bins, min, max);
  TH1F* trk2_D = new TH1F("trk2_D", "trk2_D", bins, min, max);
  TH1F* trk2_C = new TH1F("trk2_C", "trk2_C", bins, min, max);
 
  for (TString layer : Layers)
  {
    File -> cd(layer + Energy_trk1);
    trk1_D -> Add((TH1F*)gDirectory -> Get(Histograms[0]));
    File -> cd(layer + Energy_trk2);
    trk2_D -> Add((TH1F*)gDirectory -> Get(Histograms[1])); 
    File -> cd(); 
  }

  Fit_Functions f;
  f.ConvolveHists(trk1_D, trk1_D, trk2_C);
  f.ArtifactRemove(trk2_C, "b");

  f.Normalizer(trk2_C);
  trk2_C -> Scale(trk2_D -> Integral()); 

  trk1_D -> SetLineColor(kRed);
  trk2_D -> SetLineColor(kGreen);
  trk2_C -> SetLineColor(kBlack);
  trk1_D -> GetXaxis() -> SetTitle("dEdx (MeV g^{-1} cm^{2})");
  trk1_D -> SetTitle("Convolution of 1-Track data for 2-Track Production");

  TLegend* leg = new TLegend(0.9, 0.9, 0.75, 0.75);
  leg -> AddEntry(trk1_D, "1-Track Data" + Energy_trk1);
  leg -> AddEntry(trk2_D, "2-Track Data" + Energy_trk2);
  leg -> AddEntry(trk2_C, "Convolved 2-Track");

  TCanvas* can = new TCanvas();
  can -> SetLogy();
  gStyle -> SetOptStat(0);
  trk1_D -> Draw("SAMEHIST");
  trk2_D -> Draw("SAMEHIST");
  trk2_C -> Draw("SAMEHIST");
  leg -> Draw("SAME");
  can -> Update();
  can -> Print("Threshold.png");
}


