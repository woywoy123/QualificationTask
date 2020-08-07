#include<PostAnalysis/Verification.h>
#include<PostAnalysis/Functions.h>
#include<TF1.h>

using namespace RooFit;

// Test case the Fitting
std::vector<float> TestFit(std::vector<TH1F*> PDF, TH1F* Data)
{
  Fit_Functions f;
  float Lumi = Data -> Integral();

  std::vector<RooRealVar*> Var = f.FitPDFtoData(PDF, Data, 0, 20); 
 
  std::vector<float> Results; 
  for ( int i(0); i < Var.size(); i++)
  {
    float e = Var[i] -> getVal(); 
    float frac = e/Lumi;
    Results.push_back(frac); 
    delete Var[i];
  }

  return Results;
 
}


void Verification::UnitTesting()
{
  Functions F;
  Fit_Functions f;

  std::vector<TString> Hist_Names = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> Hists = F.MakeTH1F(Hist_Names, 500, 0, 20);
  //Debug(Hists, {1, 0.9, 0.1}); 
 
  // Histograms - With coloring  
  TH1F* trk1 = Hists.at(0);
  TH1F* trk2 = Hists.at(1);
  TH1F* trk3 = Hists.at(2);
  TH1F* trk4 = Hists.at(3);
  trk1 -> SetLineColor(kRed);
  trk2 -> SetLineColor(kBlue);
  trk3 -> SetLineColor(kOrange);
  trk4 -> SetLineColor(kGreen); 
  
  // ==================== Testing Units: Uncomment test units ==================== //
  // Fit testing section
  //std::vector<TH1F*> PDF = {trk1, trk2, trk3, trk4};
  //TH1F* Data1 = (TH1F*)trk1 -> Clone("Data1");
  //TH1F* Data2 = (TH1F*)trk2 -> Clone("Data2");
  //TH1F* Data3 = (TH1F*)trk3 -> Clone("Data3");
  //TH1F* Data4 = (TH1F*)trk4 -> Clone("Data4");
  //std::vector<float> Fit_test1 = TestFit(PDF, Data1);
  //std::vector<float> Fit_test2 = TestFit(PDF, Data2);
  //std::vector<float> Fit_test3 = TestFit(PDF, Data3);
  //std::vector<float> Fit_test4 = TestFit(PDF, Data4);  

  // Test Tail replace
  // Create the deconv fake vector  
  //float offset = 0.1;
  //int nbins = trk2 -> GetNbinsX(); 
  //std::vector<float> deconv(nbins + nbins*offset, 0);
  //
  //for (int i(0); i < deconv.size(); i++)
  //{
  //  if (i < nbins){deconv[i] = trk2 -> GetBinContent(i+1);}
  //  else { deconv[i] = trk2 -> GetBinContent(2*nbins - i -1);}
  //}
  // 
  //deconv = f.TailReplace(trk1, deconv); 
   
  // Deconvolution Test  
  // Create the Data fake vector  
  //float offset = 0.1;
  //int nbins = trk2 -> GetNbinsX(); 
  //std::vector<float> Data(nbins + nbins*offset, 0);
  //std::vector<float> deconv(nbins + nbins*offset, 0.5);
  // 
  //for (int i(0); i < Data.size(); i++)
  //{
  //  if (i < nbins){Data[i] = trk4 -> GetBinContent(i+1);}
  //  else { Data[i] = trk4 -> GetBinContent(2*nbins - i -1);}
  //}

  //TH1F* Deconv = (TH1F*)trk2 -> Clone("Deconv");
  //TCanvas* can = new TCanvas;
  //can -> SetLogy(); 
  //Deconv -> GetYaxis() -> SetRangeUser(1e-2, 1e6);
  //Deconv -> SetLineColor(kGreen);
  //
  //for (int y(0); y < 25; y++)
  //{
  //  for (int x(0); x < 25; x++)
  //  {
  //    Deconv -> Reset();
  //    deconv = f.LRDeconvolution(Data, deconv, deconv, 0.75); 
  //    F.VectorToTH1F(deconv, Deconv);
  //    Deconv -> Draw("SAMEHIST");
  //    trk2 -> Draw("SAMEHIST");
  //    can -> Update();
  //  }
  //  deconv = f.TailReplace(trk1, deconv);
  //}
}

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
      deconv = f.TailReplace(trk1, deconv); 

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
  float mean = 0; 
  float stdev = 0.2;
  float GaussianOffSet = -2;
  int trkn = 2;

  Fit_Functions f;
  Functions F; 
  Benchmark B;

  // Define the datasets explicitly
  TH1F* trk1 = (TH1F*)Data[0] -> Clone("trk1_Copy");
  TH1F* trk2 = (TH1F*)Data[1] -> Clone("trk2_Copy");
  TH1F* trk3 = (TH1F*)Data[2] -> Clone("trk3_Copy");
  TH1F* trk4 = (TH1F*)Data[3] -> Clone("trk4_Copy");
  std::vector<TH1F*> Data_Copy = {trk1, trk2, trk3, trk4};

  // Conversion factors
  float bins = trk2 -> GetNbinsX(); 
  float ss = (max-min) / bins;
  int GausOff = std::round(std::abs(GaussianOffSet) / ss);

  // Create the trk convolved hists
  TH1F* trk1_C = new TH1F("trk1_C", "trk1_C", bins + offset*bins + GausOff, GaussianOffSet, max); 
  TH1F* trk2_C = new TH1F("trk2_C", "trk2_C", bins + offset*bins + GausOff, GaussianOffSet, max); 
  TH1F* trk3_C = new TH1F("trk3_C", "trk3_C", bins + offset*bins + GausOff, GaussianOffSet, max); 
  TH1F* trk4_C = new TH1F("trk4_C", "trk4_C", bins + offset*bins + GausOff, GaussianOffSet, max);  

  // Static Gaussian Histogram used as a PSF
  TH1F* GausStatic = new TH1F("Gaus_S", "Gaus_S", bins + offset*bins + GausOff, GaussianOffSet, max);

  // ========= Debug Stuff ================== //
  TCanvas* can = new TCanvas();
  can -> SetLogy();
  trk1_C -> SetLineColor(kRed);
  trk2_C -> SetLineColor(kOrange);
  trk3_C -> SetLineColor(kGreen);
  trk4_C -> SetLineColor(kBlue);
  trk1_C -> SetLineStyle(kDashed);
  trk2_C -> SetLineStyle(kDashed);
  trk3_C -> SetLineStyle(kDashed);
  trk4_C -> SetLineStyle(kDashed);
  GausStatic -> SetLineStyle(kDashed);
  GausStatic -> SetLineColor(kBlack);
  TH1F* GausXTrk1 = new TH1F("Gaus_S", "Gaus_S", bins + offset*bins + GausOff, GaussianOffSet, max);

  // Vector declarations: Priors 
  std::vector<float> deconv(bins + offset*bins, 0.5);
  std::vector<float> deconv1;
  std::vector<float> deconv2;
  std::vector<float> deconv3;
  std::vector<float> deconv4;
  
  // Vector declarations: Data Vectors
  std::vector<float> Data_Vector = f.TH1FDataVector(trk2, 0.1);
  std::vector<float> DV1;
  std::vector<float> DV2;
  std::vector<float> DV3;
  std::vector<float> DV4;

  // Constants
  const int len = trk1_C -> GetNbinsX();

  // Generate a static Gaussian for deconvolution
  f.GaussianGenerator(0, 0.2, 500000, GausStatic);
  f.Normalizer(GausStatic);
  std::vector<float> Gaus = F.TH1FToVector(GausStatic);


  for (int y(0); y < 5; y++)
  { 
    for (int i(0); i < 50; i++)
    {
      deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 0.75);
      deconv = f.TailReplace(trk1, deconv);
    }

    // 1-Track Histogram  
    F.VectorToTH1F(deconv, trk1_C, GausOff); 
    f.ArtifactRemove(trk1_C, "b");
    f.Normalizer(trk1_C);
    
    // 2-Track Histogram
    f.ConvolveHists(trk1_C, trk1_C, trk2_C, GausOff);
    f.ArtifactRemove(trk2_C);
    f.Normalizer(trk2_C);

    // 3-Track Histogram
    f.ConvolveHists(trk1_C, trk2_C, trk3_C, GausOff);
    f.ArtifactRemove(trk3_C);
    f.Normalizer(trk3_C);
    
    // 4-Track Histogram
    f.ConvolveHists(trk2_C, trk2_C, trk4_C, GausOff);
    f.ArtifactRemove(trk4_C);
    f.Normalizer(trk4_C); 
 
    // === Use the Gaussian as the PSF function and do the deconvolution 
    // Define the Data Vectors
    DV1 = f.TH1FDataVector(trk1_C, 0);
    DV2 = f.TH1FDataVector(trk2_C, 0);
    DV3 = f.TH1FDataVector(trk3_C, 0);
    DV4 = f.TH1FDataVector(trk4_C, 0);

    // Define the priors for the deconv
    deconv1 = std::vector<float>(len, 0.5);
    deconv2 = std::vector<float>(len, 0.5);
    deconv3 = std::vector<float>(len, 0.5);
    deconv4 = std::vector<float>(len, 0.5);

    for (int i(0); i < 50; i++)
    {
      deconv1 = f.LRDeconvolution(DV1, Gaus, deconv1, 0.75);
      deconv2 = f.LRDeconvolution(DV2, Gaus, deconv2, 0.75);
      deconv3 = f.LRDeconvolution(DV3, Gaus, deconv3, 0.75);
      deconv4 = f.LRDeconvolution(DV4, Gaus, deconv4, 0.75);
      
      F.VectorToTH1F(deconv1, trk1_C, GausOff);
      F.VectorToTH1F(deconv2, trk2_C, GausOff);
      F.VectorToTH1F(deconv3, trk3_C, GausOff);
      F.VectorToTH1F(deconv4, trk4_C, GausOff);
     
      f.Normalizer(trk1_C);   
      f.Normalizer(trk2_C);
      f.Normalizer(trk3_C);
      f.Normalizer(trk4_C);
    }

    std::vector<TH1F*> PDFs = {trk1_C, trk2_C, trk3_C, trk4_C};
    std::vector<RooRealVar*> vars = f.GaussianConvolutionFit(PDFs, Target, 0, 20, 0.01, 1, -1, 1);
    std::vector<float> O = f.Fractionalizer(vars, trk2);

    std::cout << "Fraction of trk 1: " << O[0] << std::endl;
    std::cout << "Fraction of trk 2: " << O[1] << std::endl; 
    std::cout << "Fraction of trk 3: " << O[2] << std::endl;
    std::cout << "Fraction of trk 4: " << O[3] << std::endl;

    trk2 -> Draw("SAMEHIST");
    trk1_C -> Scale(Target -> Integral()*O[0]);
    trk2_C -> Scale(Target -> Integral()*O[1]);
    trk3_C -> Scale(Target -> Integral()*O[2]);
    trk4_C -> Scale(Target -> Integral()*O[3]);
     
    trk1_C -> Draw("SAMEHIST");
    trk4_C -> Draw("SAMEHIST");
    trk2_C -> Draw("SAMEHIST");
    trk3_C -> Draw("SAMEHIST");
    can -> Update();

    // Subtract the estimated cross contamination from the Target copy
    f.Subtraction(PDFs, trk2, trkn, vars);   
 
    Data_Vector = f.TH1FDataVector(Target, 0.1);   
    
  }
}































