#include<PostAnalysis/Verification.h>
#include<PostAnalysis/Functions.h>
#include<TF1.h>

using namespace RooFit;

void Verification::MainAlgorithm(std::vector<TH1F*> Data, TH1F* Target, std::vector<TH1F*> Closure)
{
  // Things that need an input 
  int trkn = 2; 
  float offset = 0.1;
  int bins = Target -> GetNbinsX();
  float min = 0;
  float max = 20;
  
  Fit_Functions f;
  Functions F;
  Benchmark B;

  TH1F* trk1 = Data.at(0);
  TH1F* trk2 = Data.at(1);
  TH1F* trk3 = Data.at(2);
  TH1F* trk4 = Data.at(3);

  TH1F* trk1_Clo = Closure.at(0);
  TH1F* trk2_Clo = Closure.at(1);
  TH1F* trk3_Clo = Closure.at(2);
  TH1F* trk4_Clo = Closure.at(3);
  trk2_Clo -> SetLineColor(kBlack);

  // Data we are using as the target 
  TH1F* Meas = (TH1F*)Target -> Clone("Target");

  // Now we copy the TH1F into a data vector and take the image of the tail 
  std::vector<float> Data_Vector(bins + bins*offset, 0);
  for (int i(0); i < bins; i++)
  {
    Data_Vector[i] = trk2 -> GetBinContent(i+1);
    if (i < bins*offset) { Data_Vector[i+bins] = trk2 -> GetBinContent(bins - i - 1); }
  }

  // Some histograms for debugging and a TCanvas 
  TCanvas* can = new TCanvas("can", "can", 800, 800);
  can -> Divide(2,1);
  can -> cd(1) -> SetLogy();
  can -> cd(2) -> SetLogy();

  TH1F* trk1_C = new TH1F("trk1_C", "trk1_C", bins, min, max);
  TH1F* trk2_C = new TH1F("trk2_C", "trk2_C", bins, min, max); 
  TH1F* trk3_C = new TH1F("trk3_C", "trk3_C", bins, min, max); 
  TH1F* trk4_C = new TH1F("trk4_C", "trk4_C", bins, min, max); 
  trk1_C -> SetLineStyle(kDashed);
  trk2_C -> SetLineStyle(kDashed);
  trk3_C -> SetLineStyle(kDashed);
  trk4_C -> SetLineStyle(kDashed);
  
  // Add some styles 
  trk1_C -> SetLineColor(kRed);
  trk2_C -> SetLineColor(kBlue);
  trk3_C -> SetLineColor(kOrange);
  trk4_C -> SetLineColor(kGreen); 
 
  std::vector<TH1F*> PDFs = {trk1_C, trk2_C, trk3_C, trk4_C};
  std::vector<RooRealVar*> var;
  std::vector<float> Fit_Var;
 
  std::vector<float> prediction(bins, 0);
  std::vector<float> closure(bins, 0); 

  // Assume flat prior for deconv
  std::vector<float> deconv(bins + offset*bins, 0.5);

  for (int i(0); i < 4; i++)
  {
    // Deconvolution process 
    for (int y(0); y < 100; y++)
    {
      deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 0.75); 
       
      // Tail Replace with a 1trk dataset       
      deconv = f.TailReplaceClosure(trk1, deconv); 

      can -> cd(1); 
      F.VectorToTH1F(deconv, trk1_C);  
      f.Normalizer(trk1_C); 
      f.ArtifactRemove(trk1_C);  
          
      trk1_C -> Scale(trk1_Clo -> Integral());
      trk1_Clo -> Draw("SAMEHIST");
      trk1_C -> Draw("SAMEHIST");
      can -> Update();
    }  
    can -> Update();
   
    // === Start building the n-track histograms via convolution 
    // 1-track
    F.VectorToTH1F(deconv, trk1_C);
    f.ArtifactRemove(trk1_C);
    f.Normalizer(trk1_C);

    // 2-track
    f.ConvolveHists(trk1_C, trk1_C, trk2_C, 0);
    f.ArtifactRemove(trk2_C);
    f.Normalizer(trk2_C);

    // 3-track 
    f.ConvolveHists(trk2_C, trk1_C, trk3_C, 0);
    f.ArtifactRemove(trk3_C);
    f.Normalizer(trk3_C);
   
    // 4-track
    f.ConvolveHists(trk2_C, trk2_C, trk4_C, 0);
    f.ArtifactRemove(trk4_C); 
    f.Normalizer(trk4_C);      

    // Perform the fit 
    var = f.FitPDFtoData(PDFs, Meas, 0, 20);
    Fit_Var = f.Fractionalizer(var, Meas); 
    
    std::cout << "Fraction of 1trk: " << Fit_Var[0] << std::endl;    
    std::cout << "Fraction of 2trk: " << Fit_Var[1] << std::endl;
    std::cout << "Fraction of 3trk: " << Fit_Var[2] << std::endl;
    std::cout << "Fraction of 4trk: " << Fit_Var[3] << std::endl;

    // Subtract the estimated cross contamination from the Target copy
    f.Subtraction(PDFs, Meas, trkn, var);
  
    float lumi = Target -> Integral();
    trk1_C -> Scale(lumi*Fit_Var.at(0));
    trk2_C -> Scale(lumi*Fit_Var.at(1));
    trk3_C -> Scale(lumi*Fit_Var.at(2));
    trk4_C -> Scale(lumi*Fit_Var.at(3));
   
    can -> cd(2);
    //Meas -> Draw("SAMEHIST");
    trk2_Clo -> Draw("SAMEHIST");
    trk1_C -> Draw("SAMEHIST");
    trk2_C -> Draw("SAMEHIST");
    trk3_C -> Draw("SAMEHIST");
    trk4_C -> Draw("SAMEHIST");
    can -> Update();

    for (int x(0); x < bins; x++)
    {
      Data_Vector[x] = Meas -> GetBinContent(x+1);
      if (x < bins*offset) { Data_Vector[x+bins] = Meas -> GetBinContent(bins - x - 1); }
    }

    for (int x(0); x < bins; x++)
    {
      prediction[x] = trk2_C -> GetBinContent(x+1);
      closure[x] = trk2_Clo -> GetBinContent(x+1);
    }
    std::cout << "Error normalized in the shape of the prediction: " << B.WeightedEuclidean(prediction, closure) << std::endl;
    std::cout << "Error total in the shape of the prediction: " << B.PythagoreanDistance(prediction, closure) << std::endl; 
  } 
}

void Verification::MainGaussianUnfolding(std::vector<TH1F*> Data, TH1F* Target, std::vector<TH1F*> Closure)
{
  
  // Parameters needed to define 
  float offset = 0.1;
  float min = 0;
  float max = 20;
  float mean_s = -5; 
  float mean_e = 5;
  float stdev_s = 0.001;
  float stdev_e = 1;
  int trkn = 2;
  int bins = Target -> GetNbinsX();

  Fit_Functions f;
  Functions F; 
  Benchmark B;

  // Define the datasets explicitly
  TH1F* trk1 = (TH1F*)Data[0] -> Clone("trk1_Copy");
  TH1F* trk2 = (TH1F*)Data[1] -> Clone("trk2_Copy");
  TH1F* trk3 = (TH1F*)Data[2] -> Clone("trk3_Copy");
  TH1F* trk4 = (TH1F*)Data[3] -> Clone("trk4_Copy");
  std::vector<TH1F*> Data_Copy = {trk1, trk2, trk3, trk4};
  
  // Create the trk convolved hists
  TH1F* trk1_C = new TH1F("trk1_C", "trk1_C", bins, min, max); 
  TH1F* trk2_C = new TH1F("trk2_C", "trk2_C", bins, min, max); 
  TH1F* trk3_C = new TH1F("trk3_C", "trk3_C", bins, min, max); 
  TH1F* trk4_C = new TH1F("trk4_C", "trk4_C", bins, min, max); 
   
  // ========= Debug Stuff ================== //
  TCanvas* can = new TCanvas();
  can -> SetLogy();
  trk1_C -> SetLineColor(kRed);
  trk2_C -> SetLineColor(kOrange);
  trk3_C -> SetLineColor(kGreen);
  trk4_C -> SetLineColor(kBlue);

  // Vector declaration: Data_Vector
  std::vector<float> Data_Vector = f.TH1FDataVector(trk2, offset);

  // Vector declarations: Priors 
  std::vector<float> deconv(bins + offset*bins, 0.5);
 
  for (int y(0); y < 3; y++)
  { 
    for (int z(0); z < 100; z++)
    {
      deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 0.75);
      deconv = f.TailReplaceClosure(trk1, deconv);

      F.VectorToTH1F(deconv, trk1_C);  
      f.ArtifactRemove(trk1_C); 
      f.Normalizer(trk1_C); 
    }

    // 1-Track Histogram  
    F.VectorToTH1F(deconv, trk1_C); 
    f.ArtifactRemove(trk1_C);
    f.Normalizer(trk1_C);
    
    // 2-Track Histogram
    f.ConvolveHists(trk1_C, trk1_C, trk2_C);
    f.ArtifactRemove(trk2_C);
    f.Normalizer(trk2_C);

    // 3-Track Histogram
    f.ConvolveHists(trk1_C, trk2_C, trk3_C);
    f.ArtifactRemove(trk3_C);
    f.Normalizer(trk3_C);
    
    // 4-Track Histogram
    f.ConvolveHists(trk2_C, trk2_C, trk4_C);
    f.ArtifactRemove(trk4_C);
    f.Normalizer(trk4_C); 

    // Vectors of the templates
    std::vector<TH1F*> PDFs = {trk1_C, trk2_C, trk3_C, trk4_C};

    std::vector<RooRealVar*> vars = f.GaussianConvolutionFit(PDFs, trk2, min, max, offset, stdev_s, stdev_e, mean_s, mean_e);
    std::vector<float> O = f.Fractionalizer(vars, trk2);

    std::cout << "Fraction of trk 1: " << O[0] << std::endl;
    std::cout << "Fraction of trk 2: " << O[1] << std::endl; 
    std::cout << "Fraction of trk 3: " << O[2] << std::endl;
    std::cout << "Fraction of trk 4: " << O[3] << std::endl;


    float Lumi = trk2 -> Integral();
    trk1_C -> Scale(Lumi*O[0]);    
    trk2_C -> Scale(Lumi*O[1]);
    trk3_C -> Scale(Lumi*O[2]); 
    trk4_C -> Scale(Lumi*O[3]);
    trk1_C -> Draw("SAMEHIST");
    trk2_C -> Draw("SAMEHIST");
    trk3_C -> Draw("SAMEHIST");
    trk4_C -> Draw("SAMEHIST");
    trk2 -> Draw("SAMEHIST*");
    can -> Update();
     

    // Subtract the estimated cross contamination from the Target copy
    f.Subtraction(PDFs, trk2, trkn, vars);    
    Data_Vector = f.TH1FDataVector(trk2, offset);   
    
  }
}

void Verification::DeconvolutionGaussianDebug(TH1F* trk1, TH1F* trk2)
{

  // Things that need an input 
  float offset = 0.1;
  int bins = trk1 -> GetNbinsX();
  float min = 0;
  float max = 20;
  float Padding = 10;
  
  Fit_Functions f;
  Functions F;

  // Conversion factors
  float ss = (max-min) / bins;
  int Pad = std::round(std::abs(Padding) / ss);

  // Some histograms for debugging and a TCanvas 
  TCanvas* can = new TCanvas("can", "can", 800, 800);
  can -> SetLogy();
 
  // Histograms  
  TH1F* PSF_HL = new TH1F("PSF_HL", "PSF_HL", bins + 2*Pad, min - Padding, max + Padding);  
  TH1F* Data_HL = new TH1F("Data_HL", "Data_HL", bins + 2*Pad, min - Padding, max + Padding);
  TH1F* Deconv_HL = new TH1F("Deconv_HL", "Deconv_HL", bins + 2*Pad, min - Padding, max + Padding);
  TH1F* Closure_HL = new TH1F("Closure_HL", "Closure_HL", bins + 2*Pad, min - Padding, max + Padding);
 
  // Colors 
  PSF_HL -> SetLineColor(kRed);
  Data_HL -> SetLineColor(kViolet);
  Deconv_HL -> SetLineColor(kBlack);
  Closure_HL -> SetLineColor(kCyan);
  // Note To Self:
  // RED = Gaus/PSF
  // BLACK = Deconv
  // VIOLET = Data
  // CYAN = Closure

  // Fill the PSF with a Gaussian
  f.GaussianGenerator(-0.2, 0.2, 500000, PSF_HL);

  // Fill the Data hist
  std::vector<float> temp = f.TH1FDataVector(trk1, offset);
  F.VectorToTH1F(temp, Data_HL, Pad);

  // Normalize the histograms 
  f.Normalizer(PSF_HL);
  f.Normalizer(Data_HL);

  // Populate the vectors 
  std::vector<float> PSF_V = F.TH1FToVector(PSF_HL);
  std::vector<float> Data_V = F.TH1FToVector(Data_HL);
  std::vector<float> deconv(PSF_V.size(), 0.5);

  PSF_HL -> Draw("SAMEHIST");
  Data_HL -> Draw("SAMEHIST");
  can -> Update();
  
  // Deconvolve the data with the PSF
  for (int i(0); i < 25; i++)
  {
    deconv = f.LRDeconvolution(Data_V, PSF_V, deconv, 1.);

  }
  F.VectorToTH1F(deconv, Deconv_HL, Pad);
  f.Normalizer(Deconv_HL);

  Deconv_HL -> Draw("SAMEHIST");
  can -> Update();   

  // Now we try to fit the data histogram to the deconvolution in RooFit 
  // ============================= RooFit =========================== //

  // Define the range of the dEdx
  RooRealVar* x = new RooRealVar("x", "x", min - Padding, max + Padding); 

  // Define the Gaussian Parameter: Mean
  std::vector<TString> Means_String = { "m1" };
  std::vector<float> Means_Begin = {-2};
  std::vector<float> Means_End = {2};
  std::vector<RooRealVar*> Means = f.GenerateVariables(Means_String, Means_Begin, Means_End);

  // Define the Gaussian Parameter: Standard Deviation
  std::vector<TString> Stdev_String = { "s1" };
  std::vector<float> Stdev_Begin = {0.001};
  std::vector<float> Stdev_End = {1};
  std::vector<RooRealVar*> Stdev = f.GenerateVariables(Stdev_String, Stdev_Begin, Stdev_End);

  // Define the Gaussian Variables
  std::vector<TString> Gaus_String = { "g1" };
  std::vector<RooGaussian*> G_Vars = f.GaussianVariables(Gaus_String, Means, Stdev, x);

  // Define the ntrack coefficients:
  float Lumi = Data_HL -> Integral();
  std::vector<TString> C_String = { "n_trk1" };
  std::vector<float> C_Begin = { 0 };
  std::vector<float> C_End = { Lumi };
  std::vector<RooRealVar*> C_Vars = f.GenerateVariables(C_String, C_Begin, C_End);

  // Import the PDFs as a RooDataHist
  std::vector<RooHistPdf*> PDF_Vars = f.ConvertTH1FtoPDF({Deconv_HL}, x);
 
  // Convolve the PDFs with the Gaussians
  std::vector<TString> Conv_String = { "P1xG1" };
  std::vector<RooFFTConvPdf*> Conv_Vars = f.ConvolveVariables(Conv_String, PDF_Vars, G_Vars, x);

  // Define the model we are using for the fit:
  RooAddPdf model("model", "model", RooArgList(*Conv_Vars[0]), RooArgList(*C_Vars[0]));   

  // Import the trk 2 data as a RooDataHist
  RooDataHist* trk2_D = new RooDataHist("trk2_D", "trk2_D", *x, Data_HL); 
  model.fitTo(*trk2_D, Constrain(*Means[0]), Constrain(*Stdev[0]));
  
  RooPlot* xframe = x -> frame(RooFit::Title("loL"));
  trk2_D -> plotOn(xframe);
  model.plotOn(xframe);
  xframe -> SetMinimum(1e-5);
  xframe -> Draw();
    
}


