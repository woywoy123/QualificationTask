#include<PostAnalysis/Functions.h>
#include<PostAnalysis/UnitClosures.h>
#include<TF1.h>
#include<iostream>

using namespace RooFit;

void UnitClosures::TestFit(std::vector<std::vector<TH1F*>> PDF, std::vector<TH1F*> Data, float min, float max, std::vector<std::vector<float>> Closure)
{

  for (int i(0); i < PDF.size(); i++)
  {
    std::vector<TH1F*> PDFs = PDF[i]; 
    TH1F* trk = Data[i];
    std::vector<float> Clo_trk = Closure[i];
    std::cout << "###################### TRK: " << i+1 << " ################### " << std::endl; 
    ClosureBaseFit(PDFs, trk, Clo_trk, min, max); 
  }
}

void UnitClosures::TestFit(std::vector<TH1F*> PDF, std::vector<TH1F*> Data, float min, float max, std::vector<std::vector<float>> Closure)
{
   
  for (int i(0); i < Data.size(); i++)
  {
    TH1F* trk = Data[i];
    std::vector<float> Clo_trk = Closure[i];
    std::cout << "###################### TRK: " << i+1 << " ################### " << std::endl; 
    ClosureBaseFit(PDF, trk, Clo_trk, min, max); 
  }
}

void UnitClosures::ClosureBaseFit(std::vector<TH1F*> PDFs, TH1F* d_trk, std::vector<float> Clos, float min, float max)
{
  Fit_Functions f;
  Benchmark B;

  int fd = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);  
  dup2(nullfd, 1); 
  close(nullfd);

  std::vector<RooRealVar*> var = f.FitPDFtoData(PDFs, d_trk, min, max);
  dup2(fd, 1);
  close(fd);

  std::vector<float> v = f.Fractionalizer(var, d_trk);
  float s = B.PythagoreanDistance(v, Clos);
  std::cout << "Prediction: " << v[0] << " "  << v[1] << " " << v[2] << " " << v[3] << std::endl;
  std::cout << "Truth: " << Clos[0] << " "  << Clos[1] << " " << Clos[2] << " " << Clos[3] << std::endl;
  std::cout << "Error: " << s << std::endl;
  
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

  TH1F* Target = new TH1F("Target", "Target", deconv.size(), min, max); 
  TH1F* Closure = new TH1F("Closure", "Closure", deconv.size(), min, max);
  F.ExpandTH1F(trk1, Closure);
  Closure -> SetLineStyle(kDashed);
  Closure -> Draw("SAMEHIST");  
  
  Target -> SetLineColor(kRed);
  Target -> GetYaxis() -> SetRangeUser(1e-2, 1e6);
  
  for (int i(0); i < iter; i++)
  {
    deconv = f.LRDeconvolution(Data_V, deconv, deconv, 1);
    deconv = f.TailReplaceClosure(trk1, deconv);
    F.VectorToTH1F(deconv, Target);
     
    f.Normalizer(Target);
    Target -> Scale(trk1 -> Integral());
    Target -> Draw("SAMEHIST");
    can -> Update();
  }
}

void UnitClosures::TestDeconvolution(TH1F* h1, TH1F* PSF, int iter)
{
  Fit_Functions f;
  Functions F;
 
  TH1F* Conv = (TH1F*)h1 -> Clone("Convolved"); 
  TH1F* Deconv = (TH1F*)h1 -> Clone("Deconv");
  Conv -> Reset();
  Deconv -> Reset();

  f.ConvolveHists(h1, PSF, Conv);
  f.ArtifactRemove(Conv, "b");
  f.Normalizer(Conv);
  f.Normalizer(PSF);
  Conv -> Scale(h1 -> GetBinContent(h1 -> GetMaximumBin()));
 
  std::vector<float> H1 = f.TH1FDataVector(Conv, 0.1);
  std::vector<float> PSF_V = f.TH1FDataVector(PSF, 0.1);
  std::vector<float> deconv(H1.size(), 0.5);

  Conv -> SetLineStyle(kDotted); 
  Deconv -> SetLineColor(kRed);
  h1 -> SetLineColor(kBlack);
 
  TCanvas* can = new TCanvas();
  can -> SetLogy();
      
  for (int i(0); i < iter; i++)
  {
    deconv = f.LRDeconvolution(H1, PSF_V, deconv, 1);
    F.VectorToTH1F(deconv, Deconv);
    f.Normalizer(Deconv);
    Deconv -> Scale(h1 -> Integral());
  
    Deconv -> GetYaxis() -> SetRangeUser(1, h1 -> GetBinContent(h1 -> GetMaximumBin()));
    Deconv -> Draw("SAMEHIST");
    Conv -> Draw("SAMEHIST");
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
  Functions F;
  Benchmark B; 

  // Make PDF output histograms
  std::vector<TString> PDFNames = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> PDFs = F.MakeTH1F(PDFNames, Data[0] -> GetNbinsX(), min, max, "_DataPDF"); 
  std::vector<std::vector<float>> Prediction = MinimalAlgorithmBase(Data, PDFs, min, max, offset);
  
  // ====== Plotting ====== //
  // === Truth Canvas 
  TCanvas* Truth = B.ClosurePlot("Truth", Data, Pure, Closure);
   
  // === Algorithm Canvas
  TCanvas* Algorithm = B.ClosurePlot("Algorithm", Data, PDFs, Prediction);

  Truth -> Draw();
  Algorithm -> Draw();
}

void Presentation::TestMinimalAlgorithm(std::vector<TH1F*> Data, float min, float max, float offset, std::vector<std::vector<TH1F*>> Pure)
{
  Benchmark B; 
  Functions F;

  // Make PDF output histograms
  std::vector<TString> PDFNames = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> PDFs = F.MakeTH1F(PDFNames, Data[0] -> GetNbinsX(), min, max, "_DataPDF"); 
  std::vector<std::vector<float>> Prediction = MinimalAlgorithmBase(Data, PDFs, min, max, offset);
  
  // ====== Plotting ====== //
  // === Truth Canvas  
  TCanvas* Truth = B.ClosurePlot("Truth", Data, Pure);
   
  // === Algorithm Canvas
  TCanvas* Algorithm = B.ClosurePlot("Algorithm", Data, PDFs, Prediction);

  Truth -> Draw();
  Algorithm -> Draw();
}

std::vector<std::vector<float>> Presentation::MinimalAlgorithmBase(std::vector<TH1F*> Data, std::vector<TH1F*> PDFs, float min, float max, float offset)
{
  Algorithms A; 
  Fit_Functions f;
  
  for (int i(0); i < PDFs.size(); i++)
  {
    PDFs[i] -> SetLineStyle(kDashed);
    PDFs[i] -> SetLineStyle(Constants::Colors[i]);
  }
  TH1F* trk2 = (TH1F*)Data[1] -> Clone("2 Track");
  TH1F* trk1 = (TH1F*)Data[0] -> Clone("1 Track");
  for (int i(0); i < 2; i++){A.MinimalAlgorithm(trk1, trk2, PDFs, min, max, offset, 100);}

  std::vector<std::vector<float>> Prediction;
  for (int i(0); i < PDFs.size(); i++)
  { 
    std::vector<RooRealVar*> vars = f.FitPDFtoData(PDFs, Data[i], min, max);
    std::vector<float> v = f.Fractionalizer(vars, Data[i]); 
    Prediction.push_back(v);
  }

  return Prediction;
}

void Presentation::TestGaussianAlgorithm(std::vector<TH1F*> Data, float min, float max, float offset)
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
  float mean_s = -2;
  float mean_e = 2; 
  float stdev_s= 0.001; 
  float stdev_e= 1;

  std::vector<float> Parameters;
  for (int x(0); x < 1; x++)
  { 
    Parameters = A.GaussianAlgorithm(trk1_Clone, trk2_Clone, PDFs, min, max, offset, mean_s, mean_e, stdev_s, stdev_e, 150);
  }
   
  std::vector<std::vector<float>> Prediction;
  for (int i(0); i < PDFs.size(); i++)
  { 
    std::vector<RooRealVar*> vars = f.FitPDFtoData(PDFs, Data[i], min, max);
    std::vector<float> v = f.Fractionalizer(vars, Data[i]); 
    Prediction.push_back(v);
  }
  
  // === Algorithm Canvas
  TCanvas* Algorithm = B.ClosurePlot("Algorithm", Data, PDFs, Prediction);
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

  std::vector<std::vector<TString>> Batch = {{"/200_up_GeV"}, {"/200_400_GeV","/400_600_GeV", "/600_800_GeV","/800_1000_GeV", "/1000_1200_GeV", "/1200_1400_GeV", "/1400_1600_GeV","/1600_1800_GeV"}, {"/1800_2000_GeV","/2000_2200_GeV", "/2200_2400_GeV","/2400_2600_GeV", "/2600_2800_GeV","/2800_3000_GeV", "/higher_GeV"}};

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

void DataGeneration::MonteCarlo(std::vector<TH1F*> Hists, TString dir, TString Extension)
{
  Functions F;
  TFile* File = new TFile(dir);
  if (!File -> IsOpen()) 
  { 
    std::cout << "Error Reading the entered file" << std::endl; 
  }

  for (TString Layer : Constants::Detector)
  {
    F.FillTH1F_From_File(Hists, File, Layer, Extension);
  } 
}

void DataGeneration::IdealLandau(std::vector<TH1F*> Hists, std::vector<float> COMP, std::vector<float> LandauParams)
{
  auto FillHist = [](std::vector<TH1F*> Hists, std::vector<float> Comps, std::vector<float> Params)
  {
    TF1 Lan("Lan", "landau", 0, 20);
    for (int i(0); i < Params.size(); i++)
    {
      Lan.SetParameter(i, Params.at(i));
    }
  
    for ( int i(0); i < 500000; i++)
    {
      double r1 = Lan.GetRandom();
      double r2 = Lan.GetRandom();
      double r3 = Lan.GetRandom();
      double r4 = Lan.GetRandom();
      
      if ( Comps[0] > 0){ Hists[0] -> Fill(r1, Comps[0]); }
      if ( Comps[1] > 0){ Hists[1] -> Fill(r1+r2, Comps[1]); }
      if ( Comps[2] > 0){ Hists[2] -> Fill(r1 + r2 + r3, Comps[2]); }
      if ( Comps[3] > 0){ Hists[3] -> Fill(r1 + r2 + r3 + r4, Comps[3]); }
    }
  };
  FillHist(Hists, COMP, LandauParams);
}

std::vector<float>  DataGeneration::MergeToys(std::vector<TH1F*> Hists, TH1F* trk)
{
  for (TH1F* H : Hists)
  {
    trk -> Add(H);
  }

  float lumi = trk -> Integral();
  std::vector<float> Out; 
  for (TH1F* H : Hists)
  {
    float ratio = H -> Integral() / lumi;
    Out.push_back(ratio);
  }
  return Out;
}

void Debug::GaussianFlow(std::vector<TH1F*> Data, float min, float max, float offset)
{ 
  // Importing functions 
  Fit_Functions f; 
  Functions F;
  Benchmark B;

  // == Data:
  TH1F* trk1_D = (TH1F*)Data[0] -> Clone("trk1_D");
  TH1F* trk2_D = (TH1F*)Data[1] -> Clone("trk2_D"); 

  // === Bins
  int bins = trk1_D -> GetNbinsX();

  // ===== Define the names  
  std::vector<TString> trk_P_n = Constants::Pure_Names;
  std::vector<TString> trk_1_n = Constants::trk_1;
  std::vector<TString> trk_2_n = Constants::trk_2;
  std::vector<TString> trk_3_n = Constants::trk_3;
  std::vector<TString> trk_4_n = Constants::trk_4;

  // ===== Create the histograms
  std::vector<TH1F*> O_PDFs = F.MakeTH1F(trk_P_n, bins, min, max, "_P");
  std::vector<TH1F*> trk_1_V = F.MakeTH1F(trk_1_n, bins, min, max);
  std::vector<TH1F*> trk_2_V = F.MakeTH1F(trk_2_n, bins, min, max);
  std::vector<TH1F*> trk_3_V = F.MakeTH1F(trk_3_n, bins, min, max);
  std::vector<TH1F*> trk_4_V = F.MakeTH1F(trk_4_n, bins, min, max);
  std::vector<std::vector<TH1F*>> trk_P = {trk_1_V, trk_2_V, trk_3_V, trk_4_V};
   
  // ############ Constants and histograms ############  
  // Temporary constants  
  int iter = 25;
  
  // == Forward declarations 
  std::vector<RooRealVar*> var;
  std::vector<float> Fit_Var; 
 
  // == Data Vector
  std::vector<float> Data_Vector = f.TH1FDataVector(trk2_D, offset); 

  // == Prior: Flat
  std::vector<float> deconv(bins + offset*bins, offset);
  
  // ################################################## 
  for (int y(0); y < iter; y++)
  {
    // Deconvolution  
    deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 1); 
     
    // Tail Replace with a 1trk dataset       
    deconv = f.TailReplaceClosure(trk1_D, deconv);
  } 
 
  // === Start building the n-track histograms via convolution 
  // 1-track
  F.VectorToTH1F(deconv, O_PDFs[0]);
  f.ArtifactRemoveExperimental(O_PDFs[0]); 
  f.Normalizer(O_PDFs[0]);

  // 2-track
  f.ConvolveHists(O_PDFs[0], O_PDFs[0], O_PDFs[1], 0);
  f.ArtifactRemove(O_PDFs[1], "b");
  f.Normalizer(O_PDFs[1]);

  // 3-track 
  f.ConvolveHists(O_PDFs[1], O_PDFs[0], O_PDFs[2], 0);
  //f.CorrectArtifactShift(O_PDFs[2]);
  f.ArtifactRemove(O_PDFs[2], "b");
  f.Normalizer(O_PDFs[2]);
  
  // 4-track
  f.ConvolveHists(O_PDFs[1], O_PDFs[1], O_PDFs[3], 0);
  //f.CorrectArtifactShift(O_PDFs[3]);
  f.ArtifactRemove(O_PDFs[3], "b"); 
  f.Normalizer(O_PDFs[3]);      


  // Constants in the code 
  const float Padding = max / 2;
  const float ss = (max-min) / bins;
  const int Pad = std::round(std::abs(Padding) / ss);
  const float GaussianMean = 0; /// need to add later
  const float STDev = 0.2; // Need to add later
  const float Damp = 1;
  const int iteration = 25;

  // ===== Temp
  float mean_s = -1;
  float mean_e = 3;
  float stdev_s = 0.001;
  float stdev_e = 1;

  // Define the static Gaussian PSF
  TH1F* PSF_HL = new TH1F("PSF_HL", "PSF_HL", bins + 2*Pad, min - Padding, max + Padding);   
  f.GaussianGenerator(GaussianMean, STDev, 500000, PSF_HL);
  f.Normalizer(PSF_HL); 
  std::vector<float> PSF_V = F.TH1FToVector(PSF_HL);
  delete PSF_HL;

  std::vector<TH1F*> PDF_HLs;
  for (TH1F* trk : O_PDFs)
  {
    // Configure the histograms 
    TString name = trk -> GetName(); name+=("_PDF");
    TH1F* PDF = new TH1F(name, name, bins + 2*Pad, min - Padding, max + Padding);
    
    // Fill the histograms with an offset 
    std::vector<float> temp = f.TH1FDataVector(trk, offset);
    F.VectorToTH1F(temp, PDF, Pad);

    // Normalize the histograms 
    f.Normalizer(PDF);
    
    // Convert to vector 
    std::vector<float> PDF_V = F.TH1FToVector(PDF);

    // Start the deconvolution 
    std::vector<float> deconv(PSF_V.size(), 0.5);
    for (int i(0); i < iteration; i++)
    {
      deconv = f.LRDeconvolution(PDF_V, PSF_V, deconv, Damp);  
    }
    
    // Write deconv to a histogram 
    PDF -> Reset(); 
    F.VectorToTH1F(deconv, PDF, Pad);
    f.Normalizer(PDF);
    PDF_HLs.push_back(PDF);
  }

  // Define the range of the dEdx
  RooRealVar* x = new RooRealVar("x", "x", min, max); 

  // Define the Gaussian Parameter: Mean
  std::vector<TString> Means_String = { "m1", "m2", "m3", "m4" };
  std::vector<float> Means_Begin(Means_String.size(), mean_s);
  std::vector<float> Means_End(Means_String.size(), mean_e);
  std::vector<RooRealVar*> Means = f.GenerateVariables(Means_String, Means_Begin, Means_End);

  // Define the Gaussian Parameter: Standard Deviation
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4" };
  std::vector<float> Stdev_Begin(Stdev_String.size(), stdev_s);
  std::vector<float> Stdev_End(Stdev_String.size(), stdev_e);
  std::vector<RooRealVar*> Stdev = f.GenerateVariables(Stdev_String, Stdev_Begin, Stdev_End);

  // Define the Gaussian Variables
  std::vector<TString> Gaus_String = { "g1", "g2", "g3", "g4"};
  std::vector<RooGaussian*> G_Vars = f.GaussianVariables(Gaus_String, Means, Stdev, x);
 
  // Import the PDFs as a RooDataHist
  std::vector<RooHistPdf*> PDF_Vars = f.ConvertTH1FtoPDF(PDF_HLs, x);
   
  // Define the ntrack coefficients:
  float Lumi = trk2_D -> Integral();
  std::vector<TString> C_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4" };
  std::vector<float> C_Begin = { 0, 0, 0, 0 };
  std::vector<float> C_End = { Lumi, Lumi, Lumi, Lumi };
  std::vector<RooRealVar*> C_Vars = f.GenerateVariables(C_String, C_Begin, C_End);

  // Convolve the PDFs with the Gaussians
  std::vector<TString> Conv_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4" };
  std::vector<RooFFTConvPdf*> Conv_Vars = f.ConvolveVariables(Conv_String, PDF_Vars, G_Vars, x);
  
  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgList(*Conv_Vars[0], *Conv_Vars[1], *Conv_Vars[2], *Conv_Vars[3]), RooArgList(*C_Vars[0], *C_Vars[1], *C_Vars[2], *C_Vars[3]));   
  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk = new RooDataHist("trk", "trk", *x, trk2_D); 
  model.fitTo(*trk, Constrain(*Means[0]), Constrain(*Means[1]), Constrain(*Means[2]), Constrain(*Means[3]), Constrain(*Stdev[0]), Constrain(*Stdev[1]), Constrain(*Stdev[2]), Constrain(*Stdev[3]));


  std::vector<float> Scales(C_Vars.size()); 
  
  std::vector<TString> GxT_names = {"g_trk1", "g_trk2", "g_trk3", "g_trk4"}; 
  std::vector<TH1F*> GxT_P = F.MakeTH1F(GxT_names, bins + 2*Pad, min - Padding, max + Padding, "_PDF");
  std::vector<TH1F*> GxT_G = F.MakeTH1F(GxT_names, bins + 2*Pad, min - Padding, max + Padding, "_Gau"); 
  std::vector<TH1F*> GxT_O = F.MakeTH1F(GxT_names, bins + 2*Pad, min - Padding, max + Padding, "_Out");

  TCanvas* z = new TCanvas();
  z -> SetLogy(); 
  for (int i(0); i < C_Vars.size(); i++)
  {
    Scales[i] = C_Vars[i] -> getVal(); 
    f.GaussianGenerator(Means[i] -> getVal(), Stdev[i] -> getVal(), 500000, GxT_G[i]);
    f.Normalizer(GxT_G[i]); 
 
  
    F.ExpandTH1F(O_PDFs[i], GxT_P[i]);

    // ====
    f.ConvolveHists(GxT_G[i], GxT_P[i], GxT_O[i]);
    f.Normalizer(GxT_O[i]); 
 
    f.CorrectArtifactShift(GxT_O[i]);
    O_PDFs[i] -> Reset(); 
    F.CutTH1F(GxT_O[i], O_PDFs[i], Pad);

    GxT_O[i] -> SetLineColor(Constants::Colors[i]);
    GxT_O[i] -> SetLineStyle(kDashed);

    GxT_G[i] -> SetLineColor(Constants::Colors[i]);
    GxT_G[i] -> SetLineStyle(kDotted);

    GxT_P[i] -> SetLineColor(Constants::Colors[i]);
  
    O_PDFs[i] -> Scale(C_Vars[i] -> getVa 
    O_PDFs[i] -> Draw("SAMEHIST"); 
    z -> Update();
  }
  trk2_D -> Draw("SAMEHIST");






 // // Plotting 
 // RooPlot* xframe = x -> frame(RooFit::Title("Hello"));
 // trk -> plotOn(xframe, Name("Data"), DataError(RooAbsData::None), XErrorSize(0));
 // model.plotOn(xframe, Name("1trk"), Components(*Conv_Vars[0]), LineStyle(kDotted), LineColor(kBlue));
 // model.plotOn(xframe, Name("2trk"), Components(*Conv_Vars[1]), LineStyle(kDotted), LineColor(kCyan)); 
 // model.plotOn(xframe, Name("3trk"), Components(*Conv_Vars[2]), LineStyle(kDotted), LineColor(kOrange));
 // model.plotOn(xframe, Name("4trk"), Components(*Conv_Vars[3]), LineStyle(kDotted), LineColor(kGreen));

 // TCanvas* z = new TCanvas();
 // gPad -> SetLogy();
 // xframe -> SetMinimum(1);
 // xframe -> Draw();

 // // TLegend*
 // TLegend* le = new TLegend(0.9, 0.9, 0.75, 0.75);
 // le -> AddEntry("Data", "Data");
 // le -> AddEntry("1trk", "1trk");
 // le -> AddEntry("2trk", "2trk");
 // le -> AddEntry("3trk", "3trk");
 // le -> AddEntry("4trk", "4trk");
 // le -> Draw();

 // z -> Update();







//  TCanvas* can = new TCanvas();
//  can -> SetLogy();
//  trk2_D -> SetLineColor(kBlack);
//  trk2_D -> Draw("SAMEHIST");
//  for (int i(0); i < O_PDFs.size(); i++)
//  {
//    O_PDFs[i] -> SetLineStyle(kDashed);
//    O_PDFs[i] -> SetLineColor(Constants::Colors[i]);
//    O_PDFs[i] -> Scale(Scales[i] * trk2_D -> Integral());
//    O_PDFs[i] -> Draw("SAMEHIST");
//  }


}




















