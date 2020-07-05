#include<PostAnalysis/Verification.h>
#include<PostAnalysis/Functions.h>

void Verification::RecoverScaling(RooAddPdf model, 
                                  std::vector<TH1F*> Histograms, 
                                  std::vector<RooHistPdf*> PDF, 
                                  RooRealVar* range, 
                                  std::vector<RooRealVar*> Variables,
                                  float S1, float S2)
{
  // Init the function that hosts common tools
  Functions F;
  Fit_Functions Fit;
  Plot_Functions plot;

  // 1. For Verification of scaling and correct implementation of RooFit functions
  // === 1.1: Create Scaled histograms and sum them to produce toy data sample 
  TH1F* Data_TH = (TH1F*)Histograms.at(0) -> Clone("hnew"); 
  Data_TH -> Scale(S1);
  TH1F* Data_TH_2 = (TH1F*)Histograms.at(1) -> Clone("hnew");
  Data_TH_2 -> Scale(S2);
  Data_TH -> Add(Data_TH_2);

  // === 1.2: Convert to RooDataHist and fit 
  RooDataHist* Data = Fit.ConvertTH1toDataHist(Data_TH, range);
  model.fitTo(*Data);
  
  // === 1.3: Plot the output
  TCanvas* Plot = plot.GeneratePlot("Addition of two scales histograms", range, Data, model, PDF, Constants::Pure_Names);
  Plot -> Draw(); 

  // === 1.4: Calculate the Scalings
  float trk1_e = Histograms.at(0) -> GetEntries(); 
  float trk2_e = Histograms.at(1) -> GetEntries();
  float Pre1_e = Variables.at(0) -> getVal(); 
  float Pre2_e = Variables.at(1) -> getVal();

  std::cout << "Scale of the 1-Track: " << Pre1_e / trk1_e << std::endl;
  std::cout << "Scale of the 2-Track: " << Pre2_e / trk2_e << std::endl;
}

void Verification::Subtraction(std::vector<TH1F*> Histograms, 
                               RooRealVar* range, 
                               RooAddPdf model, 
                               std::vector<RooHistPdf*> PDF, 
                               std::vector<RooRealVar*> Variables)
{

  // Init the function that hosts common tools
  Fit_Functions Fit;
  Plot_Functions plot;

  // 2. Stack the pure-ntrk templates and subtract at each stage
  // === 2.1: Create the stacked histograms for each n-track
  TH1F* Unscaled_Stack = new TH1F("Unscaled", "Unscaled", 50, 0, 20);
  Unscaled_Stack -> Add(Histograms.at(0));
  Unscaled_Stack -> Add(Histograms.at(1)); 
  Unscaled_Stack -> Add(Histograms.at(2));
  Unscaled_Stack -> Add(Histograms.at(3)); 
  RooDataHist* Data_Stack = Fit.ConvertTH1toDataHist(Unscaled_Stack, range); 
  
  // Get the entries in the pure histograms  
  float trk1_mc = Histograms.at(0) -> GetEntries();
  float trk2_mc = Histograms.at(1) -> GetEntries();
  float trk3_mc = Histograms.at(2) -> GetEntries();
  float trk4_mc = Histograms.at(3) -> GetEntries();
  
  // === 2.2: Fit the model to the data. 
  model.fitTo(*Data_Stack);
  TCanvas* Plot_1 = plot.GeneratePlot("Histogram subtraction", range, Data_Stack, model, PDF, Constants::Pure_Names);

  // Get the predicted values
  float trk1_fit_1 = Variables.at(0) -> getVal();
  float trk2_fit_1 = Variables.at(1) -> getVal();
  float trk3_fit_1 = Variables.at(2) -> getVal();
  float trk4_fit_1 = Variables.at(3) -> getVal();

  // === 2.3: Subtract the 4-trk histogram from data and fit  
  Unscaled_Stack -> Add(Histograms.at(3), -1);
  RooDataHist* Data_123_H = Fit.ConvertTH1toDataHist(Unscaled_Stack, range);
  model.fitTo(*Data_123_H);
  TCanvas* Plot_2 = plot.GeneratePlot("Histogram subtraction", range, Data_123_H, model, PDF, Constants::Pure_Names);

  // Get the prediction values
  float trk1_fit_2 = Variables.at(0) -> getVal();
  float trk2_fit_2 = Variables.at(1) -> getVal();
  float trk3_fit_2 = Variables.at(2) -> getVal();
  float trk4_fit_2 = Variables.at(3) -> getVal();
 
  // === 2.5: Subtract the 3-trk histogram from data and fit
  Unscaled_Stack -> Add(Histograms.at(2), -1);
  RooDataHist* Data_12_H = Fit.ConvertTH1toDataHist(Unscaled_Stack, range);
  model.fitTo(*Data_12_H);
  TCanvas* Plot_3 = plot.GeneratePlot("Histogram subtraction", range, Data_12_H, model, PDF, Constants::Pure_Names);

  // Get the prediction values
  float trk1_fit_3 = Variables.at(0) -> getVal();
  float trk2_fit_3 = Variables.at(1) -> getVal();
  float trk3_fit_3 = Variables.at(2) -> getVal();
  float trk4_fit_3 = Variables.at(3) -> getVal();
 

  std::cout << "========================== No subtraction results ===============================" << std::endl; 
  std::cout << "Fraction of 1-Track in Data from original: " << trk1_fit_1 / trk1_mc << std::endl;
  std::cout << "Fraction of 2-Track in Data from original: " << trk2_fit_1 / trk2_mc << std::endl;
  std::cout << "Fraction of 3-Track in Data from original: " << trk3_fit_1 / trk3_mc << std::endl;
  std::cout << "Fraction of 4-Track in Data from original: " << trk4_fit_1 / trk4_mc << std::endl;
 
  std::cout << "========================== 4 track subtraction  ===============================" << std::endl; 
  std::cout << "Fraction of 1-Track in Data from original: " << trk1_fit_2 / trk1_mc << std::endl;
  std::cout << "Fraction of 2-Track in Data from original: " << trk2_fit_2 / trk2_mc << std::endl;
  std::cout << "Fraction of 3-Track in Data from original: " << trk3_fit_2 / trk3_mc << std::endl;
  std::cout << "Fraction of 4-Track in Data from original: " << trk4_fit_2 / trk4_mc << std::endl;  

  std::cout << "========================== 3 track subtraction  ===============================" << std::endl; 
  std::cout << "Fraction of 1-Track in Data from original: " << trk1_fit_3 / trk1_mc << std::endl;
  std::cout << "Fraction of 2-Track in Data from original: " << trk2_fit_3 / trk2_mc << std::endl;
  std::cout << "Fraction of 3-Track in Data from original: " << trk3_fit_3 / trk3_mc << std::endl;
  std::cout << "Fraction of 4-Track in Data from original: " << trk4_fit_3 / trk4_mc << std::endl;
}

void Verification::Reconstruction(std::vector<TH1F*> trk1, 
                                  std::vector<TH1F*> trk2,
                                  std::vector<TH1F*> trk3,
                                  std::vector<TH1F*> trk4,
                                  std::vector<RooRealVar*> Variables, 
                                  std::vector<RooHistPdf*> PDFs,
                                  RooAddPdf model, 
                                  RooRealVar* range)
{
  // Lambda functions 
  auto Truth_Numbers = [](std::vector<TH1F*> trk)
  {
    std::vector<float> MC;
    for (TH1F* hist : trk)
    {
      MC.push_back(hist -> GetEntries());
    }
    return MC;
  };

  auto Combine_Hist = [](std::vector<TH1F*> trk, TH1F* Input)
  {
    for (TH1F* hist : trk)
    {
      Input -> Add(hist);
    }
  };

  auto Scores = [](std::vector<RooRealVar*> Var, std::vector<float> MC)
  {
    for (unsigned i = 0; i < Var.size(); i++)
    {
      float mc = MC.at(i);
      float fit = Var.at(i) -> getVal();
      std::cout << "Events predicted by fit / mc for tracks " << i+1 << ": " << fit / mc << std::endl;
    }
  };

  // Count the number of truth each of the track classifications 
  std::vector<float> trk1_MC = Truth_Numbers(trk1);
  std::vector<float> trk2_MC = Truth_Numbers(trk2);
  std::vector<float> trk3_MC = Truth_Numbers(trk3);
  std::vector<float> trk4_MC = Truth_Numbers(trk4);

  // Create the histograms we are trying to test for data
  TH1F* trk1_Hist = new TH1F("1-track", "1-track", 50, 0, 20);
  TH1F* trk2_Hist = new TH1F("2-track", "2-track", 50, 0, 20);
  TH1F* trk3_Hist = new TH1F("3-track", "3-track", 50, 0, 20);
  TH1F* trk4_Hist = new TH1F("4-track", "4-track", 50, 0, 20); 
  Combine_Hist(trk1, trk1_Hist);
  Combine_Hist(trk2, trk2_Hist);
  Combine_Hist(trk3, trk3_Hist);
  Combine_Hist(trk4, trk4_Hist);

  Fit_Functions F; 
  RooDataHist* trk1_H = F.ConvertTH1toDataHist(trk1_Hist, range);
  RooDataHist* trk2_H = F.ConvertTH1toDataHist(trk2_Hist, range);
  RooDataHist* trk3_H = F.ConvertTH1toDataHist(trk3_Hist, range); 
  RooDataHist* trk4_H = F.ConvertTH1toDataHist(trk4_Hist, range);

  model.fitTo(*trk1_H);
  
  std::cout << "=================== 1 Track Data ===================" << std::endl;
  Scores(Variables, trk1_MC); 

  Plot_Functions P; 
  TCanvas* f = P.GeneratePlot("1-Track", range, trk1_H, model, PDFs, Constants::Pure_Names); 

}
                                  
