#include<PostAnalysis/PresentationFigures.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Convolution.h>
#include<PostAnalysis/DistributionGenerator.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>

void FigureCompiler(TFile* F)
{
  std::map<TString, std::vector<TH1F*>> Map = ReadEntries(F); 
  std::map<TString, std::vector<TH1F*>>::iterator M; 
  TCanvas* can = new TCanvas(); 
  TString filename = "Output.pdf"; 
  can -> Print("Output.pdf["); 

  for (M = Map.begin(); M != Map.end(); M++)
  {
    TString CompilerName = M ->first; 
    std::vector<TH1F*> Hist_V = M ->second;
    
    can -> Print(filename);  
    if (CompilerName == "TestLandau"){ PlotLandau(can, Hist_V, filename); }
    if (CompilerName == "TestGaussian"){ PlotGaussian(can, Hist_V, filename); } 
    if (CompilerName == "TestGaussianXGaussian"){ PlotGaussianXGaussian(can, Hist_V, filename); } 
    if (CompilerName == "TestLandauXLandau"){ PlotLandauXLandau(can, Hist_V, filename); }
    if (CompilerName == "TestLandauXGaussian"){ PlotLandauXGaussian(can, Hist_V, filename); }
    if (CompilerName == "TestDeconvGausXGaus"){ PlotDeconvGausXGaus(can, Hist_V, filename); }
    if (CompilerName == "TestDeconvLandauXLandau"){ PlotDeconvLandauXLandau(can, Hist_V, filename); }
    if (CompilerName == "TestDeconvLandauXGaussian"){ PlotDeconvLandauXGaussian(can, Hist_V, filename); }
    if (CompilerName == "TestGaussianDeconvolutionFit"){ PlotGaussianDeconvolutionFit(can, Hist_V, filename); }
    if (CompilerName == "TestLandauXGausFit"){ PlotLandauXGausFit(can, Hist_V, filename); }
    if (CompilerName == "TestNLandauXNGausFit"){ PlotNLandauXNGausFit(can, Hist_V, filename); }
    if (CompilerName == "TestDeconvolutionFit"){ PlotDeconvolutionFit(can, Hist_V, filename); }
    if (CompilerName == "TestComparisonBinCenteringLandauXLandau"){ PlotComparisonBinCenteringLandauXLandau(can, Hist_V, filename);}
    if (CompilerName == "TestOscillationLucyRichardson"){ PlotOscillationLucyRichardson(can, Hist_V, filename);}
    if (CompilerName == "TestAlgorithm"){ PlotAlgorithm(can, Hist_V, filename);}
    if (CompilerName == "TestReadFile"){ PlotTestReadFile(can, Hist_V, filename);}
    if (CompilerName == "TestReadFileTrackEnergy") {PlotReadFileTrackEnergy(can, Hist_V, filename); }
    if (CompilerName == "TestMonteCarloMatchConvolution") { PlotMonteCarloMatchConvolution(can, Hist_V, filename);}
    if (CompilerName == "TestMonteCarloFit"){PlotMonteCarloFit(can, Hist_V, filename); }
    can -> Clear(); 
    can -> SetLogy(0);
  }
  can -> Print("Output.pdf]");
}

void PlotLandau(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  can -> SetLogy(); 
  Hist_V[0] -> GetYaxis() -> SetRangeUser(1, 2*GetMaxValue(Hist_V[0])); 
  GeneratePlot(Hist_V[0], "Example Landau Distribution being Generated", can, kBlue, kSolid, "HIST", 0.5); 
  can -> Print(filename); 
}

void PlotGaussian(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  GeneratePlot(Hist_V[0], "Example Gaussian Distribution being Generated", can, kRed, kSolid, "HIST", 0.5); 
  can -> Print(filename); 
}

void PlotGaussianXGaussian(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  // Plot the two analytical gaussians we are going to use as examples 
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(0, 0.2); 
  GeneratePlot(Empty, "Two analytical Gaussians Shifted", can, kBlack, kSolid, "HIST", 1); 
  GeneratePlot(Hist_V[0], "", can, kOrange, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[1], "", can, kBlue, kDashed, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[0], Hist_V[1]}, can);  
  can -> Update(); 
  can -> Print(filename); 

  // Plot the target Gaussian with the solution gaussian 
  Hist_V[2] -> SetLineColor(kRed);
  Hist_V[2] -> GetYaxis() -> SetRangeUser(0, 0.12); 
  GenerateRatioPlot(Hist_V[2], Hist_V[1], can, "Convolved Gaussian overlayed with a Shifted Target Gaussian", ""); 
  can -> Update();
  can -> Print(filename); 
} 

void PlotLandauXLandau(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{

  can -> SetLogy(); 
  // First plot the different Landaus that were numerically generated 
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-6, 0.1); 
  GeneratePlot(Empty, "Numerically Generated Landau Distributions", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "", can, kGreen, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[3], "", can, kBlue, kSolid, "SAMEHIST", 0.5);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2], Hist_V[3]}, can);  
  can -> Update(); 
  can -> Print(filename);
  
  // Now draw Ratio plots of the individual Landaus that were convolved 
  Hist_V[1] -> GetYaxis() -> SetRangeUser(1e-8, 0.02); 
  Hist_V[4] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[1], Hist_V[4], can, "A Convolved Landau overlayed with a Numerically generated Landau",  "LOG"); 
  can -> Update();  
  can -> Print(filename);
  can -> Clear();
  
  Hist_V[2] -> GetYaxis() -> SetRangeUser(1e-8, 0.014); 
  Hist_V[5] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[2], Hist_V[5], can, "A Convolved Landau overlayed with a Numerically generated Landau",  "LOG"); 
  can -> Update();  
  can -> Print(filename);
  can -> Clear(); 
   
  Hist_V[3] -> GetYaxis() -> SetRangeUser(1e-8, 0.012); 
  Hist_V[6] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[3], Hist_V[6], can, "A Convolved Landau overlayed with a Numerically generated Landau",  "LOG"); 
  can -> Update();  
  can -> Print(filename);
}

void PlotLandauXGaussian(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  can -> SetLogy(); 
  Normalize(Hist_V); 
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-7, 0.5); 
  GeneratePlot(Empty, "Landau Distribution Convolved with a Gaussian", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[1], "", can, kRed, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "", can, kBlack, kSolid, "SAMEHIST", 0.5); 
  GenerateLegend(Hist_V, can); 
  can -> Update(); 
  can -> Print(filename); 
}

void PlotDeconvGausXGaus(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-7, 1.2); 
  
  // Simple overlayer plot 
  GeneratePlot(Empty, "Illustration of Recovering an Original Gaussian using Lucy Richardson Deconvolution", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "", can, kBlack, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[1], "", can, kRed, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GenerateLegend(Hist_V, can); 
  can -> Update(); 
  can -> Print(filename);
  
  // Plot a Ratio plot of the Gaussians 
  Hist_V[0] -> GetYaxis() -> SetRangeUser(0, 1.2);
  GenerateRatioPlot(Hist_V[0], Hist_V[2], can, "Recovery of Original Point Spread Function using Lucy Richardson Deconvolution", ""); 
  can -> Update(); 
  can -> Print(filename);
}

void PlotDeconvLandauXLandau(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  can -> SetLogy(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-2, 1e8); 
 
  // Plots the Landau first 
  GeneratePlot(Empty, "Numerically Generated Landau Distributions", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "", can, kGreen, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[3], "", can, kBlue, kSolid, "SAMEHIST", 0.5);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2], Hist_V[3]}, can);  
  can -> Update(); 
  can -> Print(filename);
  can -> Clear(); 
  
  // Now plot the Results
  GeneratePlot(Empty, "Results from Deconvolving the Landau4 distribution n-times with Landau1", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[4], "", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[5], "", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[6], "", can, kGreen, kSolid, "SAMEHIST", 0.5); 
  GenerateLegend({Hist_V[4],  Hist_V[5], Hist_V[6]}, can);  
  can -> Update(); 
  can -> Print(filename);
  can -> Clear(); 
 
  // Now we plot the comparison between them 
  // Landau 1
  Hist_V[4] -> GetYaxis() -> SetRangeUser(1e-2, 8e5); 
  Hist_V[4] -> SetLineStyle(kDashed); 
  Hist_V[4] -> SetLineColor(kBlue); 
  GenerateRatioPlot(Hist_V[4], Hist_V[0], can, "Comparison Between Numerically Generated Landau1 and Recovered through Deconvolution", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Landau 2
  Hist_V[5] -> GetYaxis() -> SetRangeUser(1e-2, 6e5); 
  Hist_V[5] -> SetLineStyle(kDashed); 
  Hist_V[5] -> SetLineColor(kBlue); 
  GenerateRatioPlot(Hist_V[5], Hist_V[1], can, "Comparison Between Numerically Generated Landau2 and Recovered through Deconvolution", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
 
  // Landau 3 
  Hist_V[6] -> GetYaxis() -> SetRangeUser(1e-2, 4e5); 
  Hist_V[6] -> SetLineStyle(kDashed); 
  Hist_V[6] -> SetLineColor(kBlue); 
  GenerateRatioPlot(Hist_V[6], Hist_V[2], can, "Comparison Between Numerically Generated Landau3 and Recovered through Deconvolution", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Plot the convergence of the deconvolution 
  delete Empty; 
  Empty = (TH1F*)Hist_V[7] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1e1, 1e7);  
  
  can -> SetLogy(); 
  GeneratePlot(Empty, "Convergence of the different Landau Deconvolutions", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[7], "Landau 3", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[8], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[9], "Landau 1", can, kGreen, kSolid, "SAMEHIST", 0.5); 
  GenerateLegend({Hist_V[7],  Hist_V[8], Hist_V[9]}, can);  
  can -> Update(); 
  can -> Print(filename);
  can -> Clear(); 
}

void PlotDeconvLandauXGaussian(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  can -> SetLogy(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-7, 1); 
 
  // Plotting the Landau and the used PSF to perform the convolutions
  GeneratePlot(Empty, "Numerically Generated Landau Distributions with the PSF being used for the Convolution", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kGreen, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[3], "Landau 4", can, kBlue, kSolid, "SAMEHIST", 0.5);
  GeneratePlot(Hist_V[4], "Gaussian 1", can, kRed, kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[5], "Gaussian 2", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[6], "Gaussian 3", can, kGreen, kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[7], "Gaussian 4", can, kBlue, kDashed, "SAMEHIST", 0.5);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2], Hist_V[3],  Hist_V[4], Hist_V[5], Hist_V[6], Hist_V[7]}, can);  
  can -> Update(); 
  can -> Print(filename);
  can -> Clear(); 
 
  // Plot the convolution plots that we will try to deconvolve 
  GeneratePlot(Empty, "Landau Distributions Convolved with Gaussians", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[8], "Landau 1 X Gaussian 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[9], "Landau 2 X Gaussian 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[10], "Landau 3 X Gaussian 3", can, kGreen, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[11], "Landau 4 X Gaussian 4", can, kBlue, kSolid, "SAMEHIST", 0.5);
  GenerateLegend({Hist_V[8],  Hist_V[9], Hist_V[10], Hist_V[11]}, can);  
  can -> Update(); 
  can -> Print(filename);
  can -> Clear(); 
 
  // Now plot the deconvolutions and see how well they match the originals 
  // Landau 1
  Hist_V[12] -> GetYaxis() -> SetRangeUser(1e-9, 0.5); 
  Hist_V[12] -> SetLineStyle(kDashed); 
  Hist_V[12] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[12], Hist_V[0], can, "Recovered Landau 1 compared to the Original Landau Distributions", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Landau 2
  Hist_V[13] -> GetYaxis() -> SetRangeUser(1e-9, 0.3);  
  Hist_V[13] -> SetLineStyle(kDashed); 
  Hist_V[13] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[13], Hist_V[1], can, "Recovered Landau 2 compared to the Original Landau Distributions", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Landau 3
  Hist_V[14] -> GetYaxis() -> SetRangeUser(1e-9, 0.3); 
  Hist_V[14] -> SetLineStyle(kDashed); 
  Hist_V[14] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[14], Hist_V[2], can, "Recovered Landau 3 compared to the Original Landau Distributions", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Landau 4
  Hist_V[15] -> GetYaxis() -> SetRangeUser(1e-9, 0.3); 
  Hist_V[15] -> SetLineStyle(kDashed); 
  Hist_V[15] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[15], Hist_V[3], can, "Recovered Landau 4 compared to the Original Landau Distributions", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Now check the recovery of the PSF
  // Gaussian 1:
  Hist_V[16] -> GetYaxis() -> SetRangeUser(1e-9, 0.25); 
  Hist_V[16] -> SetLineStyle(kDashed); 
  Hist_V[16] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[16], Hist_V[4], can, "Recovered PSF compared to Original PSF", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
 
  // Gaussian 2:
  Hist_V[17] -> GetYaxis() -> SetRangeUser(1e-9, 0.25); 
  Hist_V[17] -> SetLineStyle(kDashed); 
  Hist_V[17] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[17], Hist_V[5], can, "Recovered PSF compared to Original PSF", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
 
  // Gaussian 3:
  Hist_V[18] -> GetYaxis() -> SetRangeUser(1e-9, 0.25); 
  Hist_V[18] -> SetLineStyle(kDashed); 
  Hist_V[18] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[18], Hist_V[6], can, "Recovered PSF compared to Original PSF", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
 
  // Gaussian 4:
  Hist_V[19] -> GetYaxis() -> SetRangeUser(1e-9, 0.25); 
  Hist_V[19] -> SetLineStyle(kDashed); 
  Hist_V[19] -> SetLineColor(kBlack); 
  GenerateRatioPlot(Hist_V[19], Hist_V[7], can, "Recovered PSF compared to Original PSF", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
}

void PlotGaussianDeconvolutionFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  // Plot the objective functions; The Source and the Target 
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  //Empty -> Reset(); 
  can -> SetLogy(0);
  Empty -> GetYaxis() -> SetRangeUser(0, 1.5); 
  GeneratePlot(Empty, "Illustration of the Source Gaussian which will be fitted to the Target Gaussian through RooFit", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GenerateLegend({Hist_V[1],  Hist_V[0]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
 
  Hist_V[2] -> SetLineStyle(kDashed); 
  Hist_V[2] -> GetYaxis() -> SetRangeUser(0, 1.5); 
  GenerateRatioPlot(Hist_V[2], Hist_V[1], can, "Outcome of Fitting with Convolution in RooFit", ""); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
}

void PlotLandauXGausFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  
  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(); 
  Hist_V[4] -> Scale(1e4); 
  GeneratePlot(Empty, "PDFs and PSF being used for creating the Fake Data", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[4], "Gaussian 1", can, kBlue, kSolid, "SAMEHIST", 0.4); 
  GeneratePlot(Hist_V[5], "Fake Data; L1 X G1", can, kBlack, kDashed, "SAMEHIST", 0.8); 
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2],  Hist_V[3],  Hist_V[4], Hist_V[5]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Plot the output of the fit. It is expected to only return 1 Landau.
  can -> SetLogy(); 
  GeneratePlot(Empty, "RooFit's Prediction to the Fake Data", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[5], "", can, kGreen, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[6], "G x Landau 1", can, kRed, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[7], "G x Landau 2", can, kOrange, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[8], "G x Landau 3", can, kViolet, kDashed, "SAMEHIST", 1);
  GeneratePlot(Hist_V[9], "G x Landau 4", can, kBlue, kDashed, "SAMEHIST", 1);
  GenerateLegend({Hist_V[6],  Hist_V[7], Hist_V[8],  Hist_V[9],  Hist_V[5]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Plot a ratio plot to see how well RooFit has predicted the Fake data
  Hist_V[5] -> SetLineStyle(kSolid); 
  Hist_V[5] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[5], Hist_V[6], can, "Shape and Luminosity Predictions Compared to Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 
}

void PlotNLandauXNGausFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  
  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(); 
  GeneratePlot(Empty, "Original Landau Distributions from which Fake Data and PDFs are derived", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2],  Hist_V[3]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(0); 
  Empty -> GetYaxis() -> SetRangeUser(0, 1.2); 
  GeneratePlot(Empty, "Original Gaussian Distributions used to Convolve Ideal Landaus", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[4], "Gaussian 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[5], "Gaussian 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[6], "Gaussian 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[7], "Gaussian 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[4],  Hist_V[5], Hist_V[6],  Hist_V[7]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Plot the Fake data with its underlying distribution that we need to predict.
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  can -> SetLogy();
  GeneratePlot(Empty, "Fake Data that is Generated by Convolving Landaus with a PSF and Summing those distributions", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[8], "G1XL1", can, kRed, kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[9], "G2XL2", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[10], "G3XL3", can, kViolet, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[11], "G4XL4", can, kGreen, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[12], "Fake Data", can, kBlack, kSolid, "SAMEHIST", 0.5); 
  GenerateLegend({Hist_V[12], Hist_V[8], Hist_V[10], Hist_V[11], Hist_V[9]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Plot the Prediction of RooFit compared to the actual truth 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  can -> SetLogy();
  GeneratePlot(Empty, "Evaluation of RooFit's Convolution Fit compared to Truth Convolution", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[8], "G1XL1", can, kRed, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[9], "G2XL2", can, kOrange, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[10], "G3XL3", can, kViolet, kDashed, "SAMEHIST", 1);
  GeneratePlot(Hist_V[11], "G4XL4", can, kGreen, kDashed, "SAMEHIST", 1);

  GeneratePlot(Hist_V[13], "G1XL1_RooFit", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[14], "G2XL2_RooFit", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[15], "G3XL3_RooFit", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[16], "G4XL4_RooFit", can, kGreen, kSolid, "SAMEHIST", 0.8);

  GenerateLegend({Hist_V[8], Hist_V[9], Hist_V[10], Hist_V[11], Hist_V[13], Hist_V[14], Hist_V[15], Hist_V[16]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Draw the Ratio Plots 
  Hist_V[8] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[8], Hist_V[13], can, "Shape and Luminosity Predictions From RooFit compared to Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[9] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[9], Hist_V[14], can, "Shape and Luminosity Predictions Compared to compared to Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[10] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[10], Hist_V[15], can, "Shape and Luminosity Predictions Compared tocompared to  Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[11] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[11], Hist_V[16], can, "Shape and Luminosity Predictions Compared tocompared to  Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Write out the statistics or the hists 
  std::cout << "############ Statistics of the Fit #################" << std::endl;
  std::cout << "# Centered: " << std::endl;   
  Statistics(Hist_V[8], Hist_V[13], 1, 16); 
  Statistics(Hist_V[9], Hist_V[14], 1, 16); 
  Statistics(Hist_V[10], Hist_V[15], 1, 16);
  Statistics(Hist_V[11], Hist_V[16], 1, 16);
  std::cout << std::endl; 

}

void PlotDeconvolutionFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  
  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(); 
  GeneratePlot(Empty, "Original Landau Distributions from which Fake Data and PDFs are derived", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2],  Hist_V[3]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(0); 
  Empty -> GetYaxis() -> SetRangeUser(0, 1.2); 
  GeneratePlot(Empty, "Original Gaussian Distributions used to DECONVOLVE Ideal Landaus", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[4], "Gaussian 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[5], "Gaussian 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[6], "Gaussian 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[7], "Gaussian 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[4],  Hist_V[5], Hist_V[6],  Hist_V[7]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Plot the Fake data with its underlying distribution that we need to predict.
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  can -> SetLogy();
  GeneratePlot(Empty, "Fake Data Being Generated from Summing Individual Landaus", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[12], "Fake Data", can, kBlack, kSolid, "SAMEHIST", 0.5); 
  GenerateLegend({Hist_V[12], Hist_V[0], Hist_V[1], Hist_V[2], Hist_V[3]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Compare the Deconvolved PDFs with the Original to verify them 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  can -> SetLogy();
  GeneratePlot(Empty, "Original PDFs overlayed with the Deconvolved PDFs (Inputs to RooFit)", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[8], "L1_G1_De", can, kRed, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[9], "L2_G2_De", can, kOrange, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[10], "L3_G3_De", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[11], "L4_G4_De", can, kGreen, kSolid, "SAMEHIST", 0.8);
  
  GenerateLegend({Hist_V[0], Hist_V[1], Hist_V[2], Hist_V[3], Hist_V[8], Hist_V[9], Hist_V[10], Hist_V[11]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Plot the Prediction of RooFit compared to the actual truth 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e8); 
  can -> SetLogy();
  GeneratePlot(Empty, "Evaluation of RooFit's Convolution Fit compared to Truth Convolution", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kDashed, "SAMEHIST", 1);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kDashed, "SAMEHIST", 1);

  GeneratePlot(Hist_V[13], "L1_RooFit", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[14], "L2_RooFit", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[15], "L3_RooFit", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[16], "L4_RooFit", can, kGreen, kSolid, "SAMEHIST", 0.8);

  GenerateLegend({Hist_V[0], Hist_V[1], Hist_V[2], Hist_V[3], Hist_V[13], Hist_V[14], Hist_V[15], Hist_V[16]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Draw the Ratio Plots 
  Hist_V[0] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[0], Hist_V[13], can, "Shape and Luminosity Predictions From RooFit compared to Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[1] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[1], Hist_V[14], can, "Shape and Luminosity Predictions From RooFit compared to Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[2] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[2], Hist_V[15], can, "Shape and Luminosity Predictions From RooFit compared to  Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[3] -> GetYaxis() -> SetRangeUser(1, 1e8); 
  GenerateRatioPlot(Hist_V[3], Hist_V[16], can, "Shape and Luminosity Predictions From RooFit compared to  Truth", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Write out the statistics or the hists 
  std::cout << "############ Statistics of the Fit #################" << std::endl;
  std::cout << "# Centered: " << std::endl;   
  Statistics(Hist_V[0], Hist_V[13], 1, 16); 
  Statistics(Hist_V[1], Hist_V[14], 1, 16); 
  Statistics(Hist_V[2], Hist_V[15], 1, 16);
  Statistics(Hist_V[3], Hist_V[16], 1, 16);
  std::cout << std::endl; 

}

void PlotComparisonBinCenteringLandauXLandau(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  
  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(); 
  GeneratePlot(Empty, "Landau Distribution Generated from not Centering the Bins", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2],  Hist_V[3]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Plot the PDFs and the PSF used for this fit 
  GeneratePlot(Empty, "Landau Distribution Generated from Centering the Bins", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[8], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[9], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[10], "Landau 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[11], "Landau 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[8],  Hist_V[9], Hist_V[10],  Hist_V[11]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear();  

  // Landau Wrong
  can -> SetLogy();
  Hist_V[1] -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  GenerateRatioPlot(Hist_V[1], Hist_V[5], can, "Shape Difference Between Numerically Generated Landau and Convolved Landau (Non Centering)", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[2] -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  GenerateRatioPlot(Hist_V[2], Hist_V[6], can, "Shape Difference Between Numerically Generated Landau and Convolved Landau (Non Centering)", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[3] -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  GenerateRatioPlot(Hist_V[3], Hist_V[7], can, "Shape Difference Between Numerically Generated Landau and Convolved Landau (Non Centering)", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Landau Correct
  Hist_V[9] -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  GenerateRatioPlot(Hist_V[9], Hist_V[13], can, "Shape Difference Between Numerically Generated Landau and Convolved Landau (Centering)", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[10] -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  GenerateRatioPlot(Hist_V[10], Hist_V[14], can, "Shape Difference Between Numerically Generated Landau and Convolved Landau (Centering)", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Hist_V[11] -> GetYaxis() -> SetRangeUser(1e-8, 2); 
  GenerateRatioPlot(Hist_V[11], Hist_V[15], can, "Shape Difference Between Numerically Generated Landau and Convolved Landau (Centering)", "LOG"); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  // Write out the statistics or the hists 
  std::cout << "############ Landau Bin Centering vs Landau Non Bin Centering #################" << std::endl;
  std::cout << "# Not Centered: " << std::endl;   
  Statistics(Hist_V[1], Hist_V[5], 1, 16); 
  Statistics(Hist_V[2], Hist_V[6], 1, 16); 
  Statistics(Hist_V[3], Hist_V[7], 1, 16); 
  std::cout << std::endl;
  std::cout << "# Centered: " << std::endl;   
  Statistics(Hist_V[9], Hist_V[13], 1, 16); 
  Statistics(Hist_V[10], Hist_V[14], 1, 16); 
  Statistics(Hist_V[11], Hist_V[15], 1, 16);
  std::cout << std::endl; 

}

void PlotOscillationLucyRichardson(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  can -> SetLogy(); 
  GeneratePlot(Hist_V[0], "", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();

  GeneratePlot(Hist_V[1], "", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();
  
  GeneratePlot(Hist_V[2], "", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();

  can -> SetLogy(0); 
  TH1F* Empty = CloneTH1F(Hist_V[3], {"Empty"})[0]; 
  GeneratePlot(Empty, "Gaussians with different Standard Deviations", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[3], "", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[4], "", can, kBlue, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[5], "", can, kOrange, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[3],  Hist_V[4], Hist_V[5]}, can); 
  can -> Print(filename); 
  can -> Clear();

  GeneratePlot(Hist_V[6], "", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();

  GeneratePlot(Hist_V[7], "", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();

  can -> SetLogy(); 
  Empty -> GetYaxis() -> SetRangeUser(1, 1e6);
  GeneratePlot(Empty, "Deconvolving a Landau with different Gaussian Widths", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[8], "", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[9], "", can, kBlue, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[10], "", can, kOrange, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[8],  Hist_V[9], Hist_V[10]}, can); 
  can -> Print(filename); 
  can -> Clear();

  GeneratePlot(Empty, "Deconvolving a Landau with the same Gaussian and different number of Iterations", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[11], "", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[12], "", can, kBlue, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[13], "", can, kOrange, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[11],  Hist_V[12], Hist_V[13]}, can); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[14] -> GetYaxis() -> SetRangeUser(1, 1e6);
  GeneratePlot(Hist_V[14], "Deconvolving Landau with a Gaussian - 100 bins", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[15] -> GetYaxis() -> SetRangeUser(1, 1e6);
  GeneratePlot(Hist_V[15], "Deconvolving Landau with a Gaussian - 200 bins", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[16] -> GetYaxis() -> SetRangeUser(1, 1e6);
  GeneratePlot(Hist_V[16], "Deconvolving Landau with a Gaussian - 500 bins", can, kRed, kSolid, "HIST", 1); 
  can -> Print(filename); 
  can -> Clear();
}

void PlotAlgorithm(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{  
  for (int i(0); i < 8; i++){ Normalize(Hist_V[i]); }

  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty Space");
  Empty -> Reset(); 
  Empty -> GetYaxis() -> SetRangeUser(1e-6, 2); 
  
  // Plot the PDFs and the PSF used for this fit 
  can -> SetLogy(); 
  GeneratePlot(Empty, "Original Landau Distributions with the Gaussians", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Landau 1", can, kRed, kSolid, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[1], "Landau 2", can, kOrange, kSolid, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[2], "Landau 3", can, kViolet, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[3], "Landau 4", can, kGreen, kSolid, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[4], "Gaussian 1", can, kRed,    kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[5], "Gaussian 2", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[6], "Gaussian 3", can, kViolet, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[7], "Gaussian 4", can, kGreen,  kDashed, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[0],  Hist_V[1], Hist_V[2],  Hist_V[3], Hist_V[4],  Hist_V[5], Hist_V[6],  Hist_V[7] }, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  Empty -> GetYaxis() -> SetRangeUser(1, 1e7); 
  GeneratePlot(Empty, "Fake Data with Smearing applied", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[12], "", can, kBlack,  kSolid, "SAMEHIST", 0.5);
  GeneratePlot(Hist_V[8], "", can, kRed, kDashed, "SAMEHIST", 0.5); 
  GeneratePlot(Hist_V[9], "", can, kOrange, kDashed, "SAMEHIST", 0.8); 
  GeneratePlot(Hist_V[10], "", can, kViolet, kDashed, "SAMEHIST", 0.8);
  GeneratePlot(Hist_V[11], "", can, kGreen, kDashed, "SAMEHIST", 0.8);
  GenerateLegend({Hist_V[12],  Hist_V[8], Hist_V[9],  Hist_V[10], Hist_V[11]}, can); 
  can -> Update(); 
  can -> Print(filename); 
  can -> Clear(); 

  for (int i(13); i < Hist_V.size()-1; i++)
  {
    std::cout << "############ Algorithm Test #################" << std::endl;
    Statistics(Hist_V[i], Hist_V[8], 1, 16); 
    Hist_V[i] -> GetYaxis() -> SetRangeUser(1, 1e7); 
    GenerateRatioPlot(Hist_V[i], Hist_V[8], can, "Ratio Plot of the Truth and Prediction", "LOG"); 
    can -> Print(filename); 
    can -> Clear();   
    i++;

    Statistics(Hist_V[i], Hist_V[9], 1, 16); 
    Hist_V[i] -> GetYaxis() -> SetRangeUser(1, 1e7); 
    GenerateRatioPlot(Hist_V[i], Hist_V[9], can, "Ratio Plot of the Truth and Prediction", "LOG"); 
    can -> Print(filename); 
    can -> Clear(); 
    i++; 

    Statistics(Hist_V[i], Hist_V[10], 1, 16); 
    Hist_V[i] -> GetYaxis() -> SetRangeUser(1, 1e7); 
    GenerateRatioPlot(Hist_V[i], Hist_V[10], can, "Ratio Plot of the Truth and Prediction", "LOG"); 
    can -> Print(filename); 
    can -> Clear(); 
    i++; 

    Statistics(Hist_V[i], Hist_V[11], 1, 16); 
    Hist_V[i] -> GetYaxis() -> SetRangeUser(1, 1e7); 
    GenerateRatioPlot(Hist_V[i], Hist_V[11], can, "Ratio Plot of the Truth and Prediction", "LOG"); 
    can -> Print(filename); 
    can -> Clear();
    std::cout << std::endl;
  }
  can -> SetLogy();
  GeneratePlot(Hist_V[Hist_V.size()-1], "", can, kWhite, kSolid, "HIST*", 1);
  can -> Print(filename); 


}

void PlotTestReadFile(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarlo(Dir); 
  std::vector<TString> Names = {"trk1_IBL", "trk2_IBL", "trk3_IBL", "trk4_IBL", "trk1_Blayer", "trk2_Blayer", "trk3_Blayer", "trk4_Blayer", "trk1_layer1", "trk2_layer1", "trk3_layer1", "trk4_layer1", "trk1_layer2", "trk2_layer2", "trk3_layer2", "trk4_layer2", "trk1_All", "trk2_All", "trk3_All", "trk4_All", "dEdx_IBL", "dEdx_Blayer", "dEdx_layer1", "dEdx_layer2", "dEdx_ntrk_All"}; 
 
  
  can -> SetLogy(); 
  std::vector<Color_t> Col = {kRed, kGreen, kBlue, kOrange}; 
  for (TString N : Names)
  {
    std::vector<TH1F*> Hist_V = MC[N];
    TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("EMPTY"); 
    Empty -> Reset(); 
    Empty -> GetYaxis() -> SetRangeUser(1, 1e7);  
    GeneratePlot(Empty, N, can, kWhite, kSolid, "HIST", 0); 
    for (int i(0); i < Hist_V.size(); i++)
    {
      Hist_V[i] -> Scale(1e10); 
      GeneratePlot(Hist_V[i], "", can, Col[i], kSolid, "SAMEHIST", 1); 
    }
    can -> Print(filename); 
    can -> Clear(); 
  }
}

void PlotReadFileTrackEnergy(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  TString Dir = "Merged.root"; 
  std::map<TString, std::vector<TH1F*>> MC = MonteCarloLayerEnergy(Dir);  
  std::vector<TString> JetEnergy = {"200_up_GeV", "200_400_GeV", "400_600_GeV", "600_800_GeV", "800_1000_GeV", 
                                    "1000_1200_GeV", "1200_1400_GeV", "1400_1600_GeV", "1600_1800_GeV", "1800_2000_GeV", 
                                    "2000_2200_GeV", "2200_2400_GeV", "2400_2600_GeV", "2600_2800_GeV", "2800_3000_GeV", 
                                    "higher_GeV"}; 
  std::vector<TString> Track_Names = {"Track_1_", "Track_2_", "Track_3_", "Track_4_"}; 
  std::vector<Color_t> Col = {kRed, kGreen, kBlue, kOrange}; 

  can -> SetLogy(); 
  for (TString J : JetEnergy)
  {
    for (TString T : Track_Names)
    {
      TString name = T + J; 
      std::vector<TH1F*> Hist_V = MC[name]; 
      
      int count(0); 
      for (TH1F* H : Hist_V)
      {
        float L = H -> GetEntries(); 
        if ( L == 0 ) {count++;}
      }
      
      if (count == Hist_V.size()) {continue;}

      TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("EMPTY"); 
      Empty -> Reset(); 
      Empty -> GetYaxis() -> SetRangeUser(1e-6, 1e2);  
      GeneratePlot(Empty, name, can, kWhite, kSolid, "HIST", 0); 
      for (int i(0); i < Hist_V.size(); i++)
      {
        GeneratePlot(Hist_V[i], "", can, Col[i], kSolid, "SAMEHIST", 1); 
      }
      can -> Print(filename); 
      can -> Clear(); 
    }
  }
  
  std::map<TString, std::vector<TH1F*>> MC_Other = MonteCarlo(Dir); 

  std::vector<TH1F*> MC_Other_trk1 = MC_Other["trk1_All"]; 
  std::vector<TH1F*> MC_Other_trk2 = MC_Other["trk2_All"]; 
  std::vector<TH1F*> MC_Other_trk3 = MC_Other["trk3_All"]; 
  std::vector<TH1F*> MC_Other_trk4 = MC_Other["trk4_All"]; 
  Normalize(MC_Other_trk1); 
  Normalize(MC_Other_trk2); 
  Normalize(MC_Other_trk3); 
  Normalize(MC_Other_trk4); 

  std::vector<TH1F*> MC_trk1 = MC["Track_1_All"]; 
  std::vector<TH1F*> MC_trk2 = MC["Track_2_All"]; 
  std::vector<TH1F*> MC_trk3 = MC["Track_3_All"]; 
  std::vector<TH1F*> MC_trk4 = MC["Track_4_All"]; 
  Normalize(MC_trk1); 
  Normalize(MC_trk2); 
  Normalize(MC_trk3); 
  Normalize(MC_trk4); 

  for (int i(0); i < MC_trk1.size(); i++)
  {
    TH1F* H_Other_1 = MC_Other_trk1[i]; 
    TH1F* H_1 = MC_trk1[i]; 
    GenerateRatioPlot(H_Other_1, H_1, can, "Track 1 compare", "LOG"); 
    can -> Print(filename); 
    can -> Clear(); 

    TH1F* H_Other_2 = MC_Other_trk2[i]; 
    TH1F* H_2 = MC_trk2[i]; 
    GenerateRatioPlot(H_Other_2, H_2, can, "Track 2 compare", "LOG"); 
    can -> Print(filename); 
    can -> Clear(); 

    TH1F* H_Other_3 = MC_Other_trk3[i]; 
    TH1F* H_3 = MC_trk3[i]; 
    GenerateRatioPlot(H_Other_3, H_3, can, "Track 3 compare", "LOG"); 
    can -> Print(filename); 
    can -> Clear(); 

    TH1F* H_Other_4 = MC_Other_trk4[i]; 
    TH1F* H_4 = MC_trk4[i]; 
    GenerateRatioPlot(H_Other_4, H_4, can, "Track 4 compare", "LOG"); 
    can -> Print(filename); 
    can -> Clear(); 
  }
  






}

void PlotMonteCarloMatchConvolution(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  int x(0); 
  for (TH1F* H : Hist_V)
  {
    std::cout << H -> GetTitle() << " --> " << x << std::endl; 
    x++; 
  }
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty"); 
  Empty -> Reset(); 
 
  // Original Monte Carlo pure templates 
  can -> SetLogy();  
  Empty -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Empty -> GetXaxis() -> SetRangeUser(-2, 10); 
  GeneratePlot(Empty, "Monte Carlo Distributions of n-tracks, n-truth", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[0], "Track 1, Truth 1", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[1], "Track 2, Truth 2", can, kOrange, kSolid, "SAMEHIST", 1);  
  GeneratePlot(Hist_V[2], "Track 3, Truth 3", can, kViolet, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[3], "Track 4, Truth 4", can, kGreen, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[0], Hist_V[1], Hist_V[2], Hist_V[3]}, can);  
  can -> Print(filename); 
  can -> Clear(); 

  //==== Compare the convolved distributions with the pure Monte Carlo distributions 
  // Clone the originals 
  Normalize(Hist_V); 

  // Make a Ratio Plot
  Hist_V[0] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[0] -> GetXaxis() -> SetRangeUser(0, 10); 
  Hist_V[4] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[0], Hist_V[4], can, "Ratio Plot of 1 Track - 1 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[1] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[1] -> GetXaxis() -> SetRangeUser(0, 10); 
  Hist_V[5] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[1], Hist_V[5], can, "Ratio Plot of 2 Track - 2 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[2] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[2] -> GetXaxis() -> SetRangeUser(0, 10); 
  Hist_V[6] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[2], Hist_V[6], can, "Ratio Plot of 3 Track - 3 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[3] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[3] -> GetXaxis() -> SetRangeUser(0, 10); 
  Hist_V[7] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[3], Hist_V[7], can, "Ratio Plot of 4 Track - 4 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  //==== Compare the Deconvolved and Fitted Convolutions with the pure Monte Carlo distributions 
  // Plot the Gaussian being used for the deconvolution
  Empty -> GetYaxis() -> SetRangeUser(1e-9, 10);
  GeneratePlot(Empty, "The Gaussian + Distributions being used for Deconvolution", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[8], "", can, kBlack, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[9], "", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[10], "", can, kOrange, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[11], "", can, kViolet, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[12], "", can, kGreen, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[8], Hist_V[9], Hist_V[10], Hist_V[11], Hist_V[12]}, can); 
  can -> Print(filename); 
  can -> Clear();  
  
  // Make a Ratio Plot
  Hist_V[0] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[0] -> GetXaxis() -> SetRangeUser(0, 10);
  Hist_V[0] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[0], Hist_V[13], can, "Ratio Plot of 1 Track - 1 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[1] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[1] -> GetXaxis() -> SetRangeUser(0, 10);
  Hist_V[1] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[1], Hist_V[14], can, "Ratio Plot of 2 Track - 2 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[2] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[2] -> GetXaxis() -> SetRangeUser(0, 10);
  Hist_V[2] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[2], Hist_V[15], can, "Ratio Plot of 3 Track - 3 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  Hist_V[3] -> GetYaxis() -> SetRangeUser(1e-9, 10);
  Hist_V[3] -> GetXaxis() -> SetRangeUser(0, 10);
  Hist_V[3] -> SetLineStyle(kDashed); 
  GenerateRatioPlot(Hist_V[3], Hist_V[16], can, "Ratio Plot of 4 Track - 4 Truth Compared to Monte Carlo", "LOG"); 
  can -> Print(filename); 
  can -> Clear();

  // ==== Now we compare the distributions using a statistical analysis 
  std::cout << "#################################### Statistics ########################################" << std::endl;
  std::cout << "==== Just using convolution:" << std::endl;
  Statistics(Hist_V[0], Hist_V[4], 1, 10); 
  Statistics(Hist_V[1], Hist_V[5], 1, 10); 
  Statistics(Hist_V[2], Hist_V[6], 1, 10); 
  Statistics(Hist_V[3], Hist_V[7], 1, 10); 
  std::cout << std::endl;

  std::cout << "==== Using Deconvolution + Gaussian:" << std::endl;
  Statistics(Hist_V[0], Hist_V[13], 1, 10); 
  Statistics(Hist_V[1], Hist_V[14], 1, 10); 
  Statistics(Hist_V[2], Hist_V[15], 1, 10); 
  Statistics(Hist_V[3], Hist_V[16], 1, 10); 
  std::cout << std::endl;
}

void PlotMonteCarloFit(TCanvas* can, std::vector<TH1F*> Hist_V, TString filename)
{
  auto GenClosure = [](std::vector<int> Index, std::vector<TH1F*> Hist_V, TH1F* Empty, TString Name, TCanvas* can, TString filename)
  {
    // Original Monte Carlo pure templates 
    can -> SetLogy();  
    GeneratePlot(Empty, Name, can, kWhite, kSolid, "HIST", 0); 
    GeneratePlot(Hist_V[Index[0]], "Truth 1", can, kRed, kSolid, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[Index[1]], "Truth 2", can, kOrange, kSolid, "SAMEHIST", 1);  
    GeneratePlot(Hist_V[Index[2]], "Truth 3", can, kViolet, kSolid, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[Index[3]], "Truth 4", can, kGreen, kSolid, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[Index[4]], "Summed Truth (Data)", can, kBlack, kSolid, "SAMEHIST", 1); 
    GenerateLegend({Hist_V[Index[0]], Hist_V[Index[1]], Hist_V[Index[2]], Hist_V[Index[3]], Hist_V[Index[4]]}, can);  
    can -> Print(filename); 
    can -> Clear(); 
  };

  auto GenClosureRatio = [](std::vector<int> I1, std::vector<int> I2, std::vector<TH1F*> Hist_V, TString Name, TCanvas* can, TString filename)
  {
    for (int i(0); i < I2.size(); i++)
    {
      TString N = Name; N +=(i+1);
      Hist_V[I1[i]] -> GetYaxis() -> SetRangeUser(1e-1, 100000);
      Hist_V[I1[i]] -> GetXaxis() -> SetRangeUser(0, 10); 
      Hist_V[I2[i]] ->  SetLineStyle(kDashed); 
      GenerateRatioPlot(Hist_V[I1[i]], Hist_V[I2[i]], can, "Ratio Plot of " + N + " Compared to Monte Carlo", "LOG"); 
      can -> Print(filename); 
      can -> Clear();
    }
  };

  auto Closure = [](std::vector<int> I1, std::vector<int> I2, std::vector<TH1F*> Hist_V, TH1F* Empty, TString Name, TCanvas* can, TString filename)
  {
    // Original Monte Carlo pure templates 
    can -> SetLogy();  
    GeneratePlot(Empty, Name, can, kWhite, kSolid, "HIST", 0); 
    GeneratePlot(Hist_V[I1[0]], "Truth 1", can, kRed, kSolid, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[I1[1]], "Truth 2", can, kOrange, kSolid, "SAMEHIST", 1);  
    GeneratePlot(Hist_V[I1[2]], "Truth 3", can, kViolet, kSolid, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[I1[3]], "Truth 4", can, kGreen, kSolid, "SAMEHIST", 1); 
    
    GeneratePlot(Hist_V[I2[0]], "Fit 1", can, kRed, kDashed, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[I2[1]], "Fit 2", can, kOrange, kDashed, "SAMEHIST", 1);  
    GeneratePlot(Hist_V[I2[2]], "Fit 3", can, kViolet, kDashed, "SAMEHIST", 1); 
    GeneratePlot(Hist_V[I2[3]], "Fit 4", can, kGreen, kDashed, "SAMEHIST", 1); 
    
    GenerateLegend({Hist_V[I1[0]], Hist_V[I1[1]], Hist_V[I1[2]], Hist_V[I1[3]], Hist_V[I2[0]], Hist_V[I2[1]], Hist_V[I2[2]], Hist_V[I2[3]]}, can);  
    can -> Print(filename); 
    can -> Clear(); 
  };

  int x(0); 
  for (TH1F* H : Hist_V)
  {
    std::cout << H -> GetTitle() << " --> " << x << std::endl; 
    x++; 
  }
  TH1F* Empty = (TH1F*)Hist_V[0] -> Clone("Empty"); 
  Empty -> Reset();   

  // Original Monte Carlo pure templates 
  can -> SetLogy();  
  Empty -> GetYaxis() -> SetRangeUser(1e-1, 1e8);
  Empty -> GetXaxis() -> SetRangeUser(-2, 10); 
  GenClosure({0, 1, 2, 3, 16}, Hist_V, Empty, "Track-1 Truth Distribution", can, filename); 
  GenClosure({4, 5, 6, 7, 17}, Hist_V, Empty, "Track-2 Truth Distribution", can, filename); 
  GenClosure({8, 9, 10, 11, 18}, Hist_V, Empty, "Track-3 Truth Distribution", can, filename); 
  GenClosure({12, 13, 14, 15, 19}, Hist_V, Empty, "Track-4 Truth Distribution", can, filename); 

  // Ratio Plots 
  GenClosureRatio({0, 1, 2, 3},{32, 33, 34, 35}, Hist_V, "Track-1, Truth-", can, filename);  
  GenClosureRatio({4, 5, 6, 7},{36, 37, 38, 39}, Hist_V, "Track-2, Truth-", can, filename);  
  GenClosureRatio({8, 9, 10, 11},{40, 41, 42, 43}, Hist_V, "Track-3, Truth-", can, filename);  
  GenClosureRatio({12, 13, 14, 15},{44, 45, 46, 47}, Hist_V, "Track-4, Truth-", can, filename);  

  // Plot the Gaussian being used for the deconvolution
  Empty -> GetYaxis() -> SetRangeUser(1e-1, 1e8);
  GeneratePlot(Empty, "The Gaussian + Distributions being used for Deconvolution", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[24], "Gaussian", can, kBlack, kDashed, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[20], "Conv 1", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[21], "Conv 2", can, kOrange, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[22], "Conv 3", can, kViolet, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[23], "Conv 4", can, kGreen, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[24], Hist_V[20], Hist_V[21], Hist_V[22], Hist_V[23]}, can); 
  can -> Print(filename); 
  can -> Clear();  
 
  // Plot the Deconvolved Hists
  Empty -> GetYaxis() -> SetRangeUser(1e-1, 1e10);
  GeneratePlot(Empty, "The Deconvolved Histograms", can, kWhite, kSolid, "HIST", 0); 
  GeneratePlot(Hist_V[28], "Deconv 1", can, kRed, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[29], "Deconv 2", can, kOrange, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[30], "Deconv 3", can, kViolet, kSolid, "SAMEHIST", 1); 
  GeneratePlot(Hist_V[31], "Deconv 4", can, kGreen, kSolid, "SAMEHIST", 1); 
  GenerateLegend({Hist_V[28], Hist_V[29], Hist_V[30], Hist_V[31]}, can); 
  can -> Print(filename); 
  can -> Clear();  

  Closure({0, 1, 2, 3},{32, 33, 34, 35}, Hist_V, Empty, "Track 1 Monte Carlo Truth vs Reconstructed", can, filename); 
  Closure({4, 5, 6, 7},{36, 37, 38, 39}, Hist_V, Empty, "Track 2 Monte Carlo Truth vs Reconstructed", can, filename); 
  Closure({8, 9, 10, 11},{40, 41, 42, 43}, Hist_V, Empty, "Track 3 Monte Carlo Truth vs Reconstructed", can, filename); 
  Closure({12, 13, 14, 15},{44, 45, 46, 47}, Hist_V, Empty, "Track 4 Monte Carlo Truth vs Reconstructed", can, filename); 

  // ==== Now we compare the distributions using a statistical analysis 
  std::cout << "#################################### Statistics ########################################" << std::endl;
  std::cout << "=========Track-1: Statistics" << std::endl; 
  Statistics(Hist_V[0], Hist_V[32], 1, 10); 
  Statistics(Hist_V[1], Hist_V[33], 1, 10); 
  Statistics(Hist_V[2], Hist_V[34], 1, 10); 
  Statistics(Hist_V[3], Hist_V[35], 1, 10); 
  std::cout << std::endl;

  std::cout << "=========Track-2: Statistics" << std::endl; 
  Statistics(Hist_V[4], Hist_V[36], 1, 10); 
  Statistics(Hist_V[5], Hist_V[37], 1, 10); 
  Statistics(Hist_V[6], Hist_V[38], 1, 10); 
  Statistics(Hist_V[7], Hist_V[39], 1, 10); 
  std::cout << std::endl;

  std::cout << "=========Track-3: Statistics" << std::endl; 
  Statistics(Hist_V[8], Hist_V[40], 1, 10); 
  Statistics(Hist_V[9], Hist_V[41], 1, 10); 
  Statistics(Hist_V[10], Hist_V[42], 1, 10); 
  Statistics(Hist_V[11], Hist_V[43], 1, 10); 
  std::cout << std::endl;

  std::cout << "=========Track-4: Statistics" << std::endl; 
  Statistics(Hist_V[12], Hist_V[44], 1, 10); 
  Statistics(Hist_V[13], Hist_V[45], 1, 10); 
  Statistics(Hist_V[14], Hist_V[46], 1, 10); 
  Statistics(Hist_V[15], Hist_V[47], 1, 10); 
  std::cout << std::endl;

}

