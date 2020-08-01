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

void Verification::MainAlgorithm(std::vector<TH1F*> Data, TH1F* Target)
{
  // Things that need an input 
  int trkn = 3; 
  float offset = 0.1;
  int bins = Target -> GetNbinsX();
  
  Fit_Functions f;
  Functions F;

  TH1F* trk1 = Data.at(0);
  TH1F* trk2 = Data.at(1);
  TH1F* trk3 = Data.at(2);
  TH1F* trk4 = Data.at(3);

  // Data we are using as the target 
  TH1F* Meas = (TH1F*)Target -> Clone("Target");

  // ==================== Preparation ============================ //
  // We perform a preliminary fit  
  std::vector<RooRealVar*> var = f.FitPDFtoData(Data, Meas, 0, 20);
  std::vector<float> Fit_Var = f.Fractionalizer(var, Meas);  

  // Monitoring purpose. Delete after 
  std::cout << Fit_Var.at(0) << std::endl;// <<<<<< Delete me 
  std::cout << Fit_Var.at(1) << std::endl;// <<<<<< Delete me   
  std::cout << Fit_Var.at(2) << std::endl;// <<<<<< Delete me 
  std::cout << Fit_Var.at(3) << std::endl;// <<<<<< Delete me 

  // Subtract the estimated cross contamination from the Target copy
  f.Subtraction(Data, Meas, trkn, var);

  // Now we copy the TH1F into a data vector and take the image of the tail 
  std::vector<float> Data_Vector(bins + bins*offset, 0);
  for (int i(0); i < bins; i++)
  {
    Data_Vector[i] = Meas -> GetBinContent(i+1);
    if (i < offset) { Data_Vector[i] = Meas -> GetBinContent(bins - i - 1); }
  }

  // ==================== Preparation: End ========================= //
 
  // Some histograms for debugging and a TCanvas 
  TCanvas* can = new TCanvas("can", "can", 1200, 1200);
  can -> Divide(2,2);
  can -> cd(1) -> SetLogy();
  can -> cd(2) -> SetLogy();
  can -> cd(3) -> SetLogy(); 
  can -> cd(4) -> SetLogy();
  
  TH1F* trk1_C = (TH1F*)trk1 -> Clone("trk1_C"); 
  trk1_C -> GetYaxis() -> SetRangeUser(1, 1e12);
  TH1F* trk2_C = (TH1F*)trk2 -> Clone("trk2_C"); 
  TH1F* trk3_C = (TH1F*)trk3 -> Clone("trk3_C"); 
  TH1F* trk4_C = (TH1F*)trk4 -> Clone("trk4_C"); 
  trk1_C -> SetLineStyle(kDashed);
  trk2_C -> SetLineStyle(kDashed);
  trk3_C -> SetLineStyle(kDashed);
  trk4_C -> SetLineStyle(kDashed);
  
  // Add some styles 
  trk1_C -> SetLineColor(kRed);
  trk2_C -> SetLineColor(kBlue);
  trk3_C -> SetLineColor(kOrange);
  trk4_C -> SetLineColor(kGreen); 
 
  // ==================== Main Part ======================= //
  // Assume flat prior for deconv
  std::vector<float> deconv = std::vector<float>(bins + offset*bins, 0.5);
  std::vector<TH1F*> PDFs = {trk1_C, trk2_C, trk3_C, trk4_C};
  for (int i(0); i < 25; i++)
  {
    //trk1_C -> Reset();
    //trk2_C -> Reset();
    //trk3_C -> Reset();
    //trk4_C -> Reset();

    // Deconvolution process
    for (int y(0); y < 25; y++)
    {
      deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 0.75);
 
      // Tail Replace with a 1trk dataset
      deconv = f.TailReplace(trk1, deconv);
    }
  
    
    // === Start building the n-track histograms via convolution 
    // 1-track
    F.VectorToTH1F(deconv, trk1_C);

    // 2-track
    f.ConvolveHists(trk1_C, trk1_C, trk2_C, 0);

    // 3-track 
    f.ConvolveHists(trk1_C, trk2_C, trk3_C, 0);

    // 4-track
    f.ConvolveHists(trk2_C, trk2_C, trk4_C, 0);
      
    // Now normalize these hists 
    f.Normalizer(trk1_C);
    f.Normalizer(trk2_C);
    f.Normalizer(trk3_C);
    f.Normalizer(trk4_C);
 
    can -> cd(1); 
    trk1_C -> Scale(trk1 -> Integral());
    //trk1 -> Draw("SAMEHIST*");
    trk1_C -> Draw("SAMEHIST");
    
//    can -> cd(2); 
//    trk2_C -> Scale(trk2 -> Integral()); 
//    trk2_C -> Draw("SAMEHIST");
//    trk2 -> Draw("SAMEHIST*");
//     
//    can -> cd(3); 
//    trk3_C -> Scale(trk3 -> Integral());    
//    trk3_C -> Draw("SAMEHIST");
//    trk3 -> Draw("SAMEHIST*");
// 
//    can -> cd(4); 
//    trk4_C -> Scale(trk4 -> Integral()); 
//    trk4_C -> Draw("SAMEHIST");
//    trk4 -> Draw("SAMEHIST*");
 
    can -> Update();




 
    // Perform the fit 
    // We perform a preliminary fit 
  //  Meas = (TH1F*)Target -> Clone("Target");
  //  var = f.FitPDFtoData(PDFs, Meas, 0, 20);
  //  Fit_Var = f.Fractionalizer(var, Meas);  

  //  // Monitoring purpose. Delete after 
  //  std::cout << Fit_Var.at(0) << std::endl;// <<<<<< Delete me 
  //  std::cout << Fit_Var.at(1) << std::endl;// <<<<<< Delete me   
  //  std::cout << Fit_Var.at(2) << std::endl;// <<<<<< Delete me 
  //  std::cout << Fit_Var.at(3) << std::endl;// <<<<<< Delete me 

  //  // Subtract the estimated cross contamination from the Target copy
  //  f.Subtraction(PDFs, Meas, trkn, var);
  }
























 
}


























