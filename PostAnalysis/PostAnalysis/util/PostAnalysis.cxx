// Including commonly reused functions
#include<PostAnalysis/Functions.h>
#include<PostAnalysis/Verification.h>
#include<TF1.h>

// Including standard C++ libraries 
#include<iostream>

using namespace RooFit;

void PostAnalysis()
{
  auto FillHist = [](TH1F* Hist, std::vector<float> Comps, std::vector<float> Params)
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

      Hist -> Fill(r1, Comps[0]);
      Hist -> Fill(r1+r2, Comps[1]);
      Hist -> Fill(r1 + r2 + r3, Comps[2]);
      Hist -> Fill(r1 + r2 + r3 + r4, Comps[3]);
    }



  };

  auto LandauGenerator = [](std::vector<TH1F*> Hists, std::vector<float> Params)
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
      Hists.at(0) -> Fill(r1);  
      Hists.at(1) -> Fill(r1 + r2); 
      Hists.at(2) -> Fill(r1 + r2 + r3); 
      Hists.at(3) -> Fill(r1 + r2 + r3 + r4); 
    }
  };

  Verification V;
  
  //V.UnitTesting();
 
  Functions F;

  // ====================== Generate the data =================== //  
  // Histograms: Pure with no cross contamination
  std::vector<TString> Hist_Names = {"trk1_P", "trk2_P", "trk3_P", "trk4_P"};
  std::vector<TH1F*> Hists = F.MakeTH1F(Hist_Names, 500, 0, 20);
  LandauGenerator(Hists, {1, 0.9, 0.1}); 
 
  // Create some Datasets which have some contamination
  // === COMPN is the composition of cross contamination
  std::vector<float> COMP1 = {1.0,  0.0, 0.0, 0.0};  
  std::vector<float> COMP2 = {0.05, 0.8,  0.1,  0.05};  
  std::vector<float> COMP3 = {0.01, 0.2,  0.59,  0.2};  
  std::vector<float> COMP4 = {0.02,   0.2, 0.2,   0.58}; 
  std::vector<TString> Data_Names = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> Data_ntrk = F.MakeTH1F(Data_Names, 500, 0, 20);
  
  // === Create the data sets
  // State the names explicitly 
  TH1F* trk1 = Data_ntrk.at(0);
  TH1F* trk2 = Data_ntrk.at(1);
  TH1F* trk3 = Data_ntrk.at(2);
  TH1F* trk4 = Data_ntrk.at(3);
 
  // Fill the histograms with the composition fractions  
  FillHist(trk1, COMP1, {1, 0.9, 0.1});
  FillHist(trk2, COMP2, {1, 0.9, 0.1});
  FillHist(trk3, COMP3, {1, 0.9, 0.1});
  FillHist(trk4, COMP4, {1, 0.9, 0.1});

  // Add some styles 
  trk1 -> SetLineColor(kRed);
  trk2 -> SetLineColor(kBlue);
  trk3 -> SetLineColor(kOrange);
  trk4 -> SetLineColor(kGreen); 
 
  // ======================= End of Data Generation ================ //  

  TCanvas* can = new TCanvas("Cant", "Cant", 800, 800);
  can -> SetLogy();
  trk1 ->Draw("SAMEHIST"); 
  trk2 -> Draw("SAMEHIST");
  trk3 -> Draw("SAMEHIST");
  trk4 -> Draw("SAMEHIST");
  trk2 -> Draw("SAMEHIST*");
  can -> Update();

 
  V.MainAlgorithm(Data_ntrk, trk2, Hists); 
  
  
  
  /*

  // === 2. State the Data we will be using: 
  TH1F* Data_Copy = (TH1F*)trk2Data -> Clone("FAKE");
  std::vector<float> Settings = COMP2;
  int NumTrak = 2; // <--------- Need to Automate this 
  float Offset = 0.1; // <------------------------------------------- Deconvolution Param
  
  // === 3.: Copy the track data
  TH1F* trk1_C = (TH1F*)trk1Data -> Clone("trk1_C");
  TH1F* trk2_C = (TH1F*)trk2Data -> Clone("trk2_C");
  TH1F* trk3_C = (TH1F*)trk3Data -> Clone("trk3_C");
  TH1F* trk4_C = (TH1F*)trk4Data -> Clone("trk4_C");
  
  // === 3.1: Normalize the data 
  f.Normalizer(trk1_C);
  f.Normalizer(trk2_C);
  f.Normalizer(trk3_C);
  f.Normalizer(trk4_C);
  
  // === 3.2: Store in vectors
  std::vector<TH1F*> PDFs= {trk1_C, trk2_C, trk3_C, trk4_C};
  std::vector<TH1F*> DataHists = {trk1Data, trk2Data, trk3Data, trk4Data}; 

  // ===================== END OF DATA ========================== // 




  // =========== Initial Estimate of contamination ============== //
  // === 0. Perform a simple fit with the hists that are given. 
  std::vector<RooRealVar*> vg = f.FitPDFtoData(PDFs, Data_Copy, 0, 20, Constants::Variable_Names, Constants::Begin, Constants::End);
   
  // === 1. Subtract the tracks from the Data Copy 
  f.Subtraction(PDFs, Data_Copy, NumTrak, vg); 

  // === 3. Copy the contents of the TH1F into a vector and create the mirror
  int nBins = Data_Copy -> GetNbinsX();
  int OS = Offset*nBins;
  std::vector<float> Data_Vector(nBins + OS, 0);
  
  for (int i(0); i < nBins; i++)
  {
    Data_Vector[i] = Data_Copy -> GetBinContent(i+1);
    if ( i < OS) { Data_Vector[nBins + i] = Data_Copy -> GetBinContent(nBins -i -1); } 
  }

  // === 4. Delete pointers and Reset any useful Histograms and set their colors 
 
  trk1_C -> SetMarkerColor(kBlue);
  trk1_C -> SetLineColor(kBlue);
  trk1_C -> SetLineWidth(3);

  trk2_C -> SetMarkerColor(kRed);
  trk2_C -> SetLineColor(kRed);
  trk2_C -> SetLineWidth(3);

  trk3_C -> SetMarkerColor(kGreen);
  trk3_C -> SetLineColor(kBlue);
  trk3_C -> SetLineWidth(3);

  trk4_C -> SetMarkerColor(kYellow);
  trk4_C -> SetLineColor(kYellow);
  trk4_C -> SetLineWidth(3);

  // =========== END OF: Initial Estimate of contamination ============== //


  // =================== Main Algorithm ================================= // 
  // === 0. Prepare the deconv vector 
  // === 0.1: Assume a flat prior of deconv 
  std::vector<float> deconv = std::vector<float>(nBins + OS, 0.5);

  // === 0.2: Define some probes for checking the deconvolution
  TCanvas* can = new TCanvas();
  can -> Divide(2,2);
    
  // === 1. Start the deconvolution process
  for (int i(0); i < 1; i++)
  {
 
    trk1_C -> Reset("ICES");
    trk2_C -> Reset("ICES");
    trk3_C -> Reset("ICES");
    trk4_C -> Reset("ICES");
    
    for (int x(0); x < 25; x++)
    {
      deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 0.75);
    }
    // === 1.1: Replace the tail of the deconv with the 1trk data 
    //deconv = f.TailReplace(trk1Data, deconv, 0, 20);
    
    // === 2. Start building the 2-trk, 3-trk, 4-trk histograms 
    // === 2.1: 1-Track histo -> Name: trk1_C 
    F.VectorToTH1F(deconv, trk1_C);
    f.Normalizer(trk1_C);
   
    // === 2.2: 2-Track histo -> Name: trk2_C
    f.ConvolveHists(trk1_C, trk1_C, trk2_C, 0);
    //f.Normalizer(trk2_C);

    // === 2.3: 3-Track histo -> Name: trk3_C
    f.ConvolveHists(trk2_C, trk1_C, trk3_C, 0);
    //f.Normalizer(trk3_C);

    // === 2.4: 4-Track histo -> Name: trk4_C
    f.ConvolveHists(trk2_C, trk2_C, trk4_C, 0);
    //f.Normalizer(trk4_C);
    std::vector<TH1F*> PDFs = {trk1_C, trk2_C, trk3_C, trk4_C};  

    TH1F* trk2 = (TH1F*)trk2Data -> Clone("TRK2");
    std::vector<RooRealVar*> var = f.FitPDFtoData(PDFs, trk2, 0, 20, 
                                                  Constants::Variable_Names, 
                                                  Constants::Begin, Constants::End);  
     
    // === 3.: Subtract the tracks from the trk2
    std::vector<float> vars = f.Subtraction(PDFs, trk2, NumTrak, var);
   
    // === 4. Convert the new data estimate to deconv  
    int n = trk2 -> GetNbinsX(); 
    for (int i(0); i < n; i++)
    {
      float e = trk2 -> GetBinContent(i+1);
      deconv[i] = e; 
    }
     
    P.View(can, PDFs); 
    P.View(can, DataHists);
        
    // Perform the fits on all the track data
    std::vector<RooRealVar*> trk1_v = f.FitPDFtoData(PDFs, trk1Data, 0, 20, Constants::Variable_Names, Constants::Begin, Constants::End);
    std::vector<RooRealVar*> trk2_v = f.FitPDFtoData(PDFs, trk2Data, 0, 20, Constants::Variable_Names, Constants::Begin, Constants::End);
    std::vector<RooRealVar*> trk3_v = f.FitPDFtoData(PDFs, trk3Data, 0, 20, Constants::Variable_Names, Constants::Begin, Constants::End);
    std::vector<RooRealVar*> trk4_v = f.FitPDFtoData(PDFs, trk4Data, 0, 20, Constants::Variable_Names, Constants::Begin, Constants::End);
    
    std::vector<float> trk1_f = f.Fractionalizer(trk1_v, trk1Data);
    std::vector<float> trk2_f = f.Fractionalizer(trk2_v, trk2Data);
    std::vector<float> trk3_f = f.Fractionalizer(trk3_v, trk3Data);
    std::vector<float> trk4_f = f.Fractionalizer(trk4_v, trk4Data);
    std::vector<std::vector<float>> trk_f = {trk1_f, trk2_f,  trk3_f, trk4_f};
    
    can -> Draw(); 
    can -> Update();   

    for (int x(0); x < trk1_v.size(); x++)
    {
      delete trk1_v[x];
      delete trk2_v[x];     
      delete trk3_v[x];
      delete trk4_v[x];  
      std::cout << "######################### Fit for TRK-" << x+1 << " #########################" << std::endl; 
      std::cout << " Trk1: " << trk_f[x][0] << " Trk2: " << trk_f[x][1] << " Trk3: "<< trk_f[x][2] << " Trk4: "<< trk_f[x][3] << std::endl; 
      std::cout << "      " << std::endl; 
    }
  } 

  */
}

void StandaloneApplications(int argc, char**argv)
{
  PostAnalysis();
}

int main(int argc, char** argv)
{
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplications(app.Argc(), app.Argv());
  app.Run();

  return 0;
}
