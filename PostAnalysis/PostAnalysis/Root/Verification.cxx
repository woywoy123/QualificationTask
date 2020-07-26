#include<PostAnalysis/Verification.h>
#include<PostAnalysis/Functions.h>
#include<TF1.h>

using namespace RooFit;

void Verification::RecoverScaling(RooAddPdf model, 
                                  std::vector<TH1F*> Histograms, 
                                  std::vector<RooHistPdf*> PDF, 
                                  RooRealVar* range, 
                                  std::vector<RooRealVar*> Variables,
                                  float S1, float S2)
{
  // Init the function that hosts common tools
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

void Verification::FLostLayer(std::map<TString, std::vector<TH1F*>> Layer_Track, float lower, float upper)
{
  // Lambda functions
  auto ntrk_ntru = [](std::vector<TH1F*> Histo)
  {
    Functions F; 
    std::vector<std::pair<int, int>> Output;
    for (TH1F* i : Histo)
    {
      TString name = i -> GetTitle();
      std::vector<TString> list = F.SplitString("_", name);
      int ntrk = list.at(2).Atoi();
      int ntru = list.at(4).Atoi(); 
      Output.push_back(std::pair<int, int>(ntrk, ntru)); 
    }
    return Output;
  };

  auto FLost = [](std::vector<TH1F*> Histo, std::vector<std::pair<int, int>> tktu, float lower, float upper)
  {
    // loss sum for each truth
    float loss2 = 0;
    float loss3 = 0;
    float loss4 = 0; 

    // Total tracks
    float tot2 = 0;
    float tot3 = 0;
    float tot4 = 0;

    for (unsigned i = 0; i < Histo.size(); i++)
    {
      std::pair<int, int> it = tktu.at(i);
      int trk = it.first;
      int tru = it.second;
     
      TH1F* h = Histo.at(i);
      float W = h -> GetBinWidth(1);
      float entries = h -> Integral(lower/W, upper/W);
      float loss;

      // case where ntrk < ntru i.e. lost truth
      loss = (tru - trk)*entries;
     
      if ( tru == 2 && trk <= tru )
      {
        loss2 = loss2 + loss;
        tot2 = entries*(tru-trk+1) + tot2;
      }
      if ( tru == 3 && trk <= tru )
      {
        loss3 = loss3 + loss;
        tot3 = entries*(tru-trk+1) + tot3;
      }
      if ( tru == 4 && trk <= tru )
      {
        loss4 = loss4 + loss;
        tot4 = entries*(tru-trk+1) + tot4;
      }
    }
    
    std::vector<float> flost = {loss2/tot2, loss3/tot3, loss4/tot4};
    return flost;
  };

  std::cout << "" << std::endl;
  for (std::pair<TString, std::vector<TH1F*>> i : Layer_Track)
  {
    std::vector<std::pair<int, int>> tktu = ntrk_ntru(i.second);
    std::vector<float> FLOST = FLost(i.second, tktu, lower, upper);
    
    std::cout << "    Calculated FLost for Layer: " << i.first  << std::endl; 
    std::cout << "    2-Track: " << FLOST.at(0) << "| 3-Track: " << FLOST.at(1) << "| 4-Track: " << FLOST.at(2) << std::endl; 
    std::cout << "" << std::endl;
  }
   
}

void Verification::NormalizedSubtraction(float lower, float upper, std::vector<TH1F*> Hist1, std::vector<TH1F*> Hist2, RooAddPdf model, std::vector<RooRealVar*> Variables, RooRealVar* Range)
{

  auto GenerateData = [](std::vector<TH1F*> Hist1, TH1F* Data_H1, std::vector<TH1F*> Hist2, TH1F* Data_H2)
  {
    // Create the data histograms
    for (unsigned int i = 0; i < Hist1.size(); i++)
    {
      Data_H1 -> Add(Hist1.at(i));
      Data_H2 -> Add(Hist2.at(i));
    }
  };
 
  auto CompleteSubtraction = [](TH1F* Hist1, TH1F* Hist2, float lower, float upper)
  {
    float W = Hist1 -> GetBinWidth(1);
    // Assume Hist1 is the one we want to subtract from: e.g. 3-Track
    // Hist 2 is the one we are subtracting: e.g. 4-Track
    float E_H1 = Hist1 -> Integral(lower/W, upper/W);
    float E_H2 = Hist2 -> Integral(lower/W, upper/W);
    
    // Scale the histogram we are going to subtract
    Hist2-> Scale(E_H1/E_H2); 
    Hist1 -> Add(Hist2, -1);
  };
  
  auto AlternativeSubtraction = [](TH1F* Hist1, TH1F* Hist2, float lower, float upper)
  {
    float W = Hist1 -> GetBinWidth(1);

    // Assume Hist1 is the one we want to subtract from: e.g. 3-Track
    // Hist 2 is the one we are subtracting: e.g. 4-Track
    float E_H1 = Hist1 -> Integral(lower/W, upper/W);
    float E_H2 = Hist2 -> Integral(lower/W, upper/W);

    float Total = E_H1 + E_H2;
    std::cout << E_H1/Total << ":::" << E_H2/Total << std::endl;

    Hist2 -> Scale((E_H1/E_H2)*(E_H2/Total));
    Hist1 -> Add(Hist2, -1);

    return (E_H1/Total);
  };
  
   
  // Clone the TH1F structure and generate the toy data
  TH1F* Data_H1 = (TH1F*)Hist1.at(0) -> Clone("Data1");
  Data_H1 -> Reset();
  Data_H1 -> SetTitle("Remaining Distribution");

  TH1F* Data_H2 = (TH1F*)Data_H1 -> Clone("Data2");
  GenerateData(Hist1, Data_H1, Hist2, Data_H2);

  TH1F* Original1 = (TH1F*)Data_H1 -> Clone("Data1");
  TH1F* Original2 = (TH1F*)Data_H2 -> Clone("Data2");

  // This causes most of the data in the domain to be removed. Try another method 
  //CompleteSubtraction(Data_H1, Data_H2, lower, upper);

  // Different approach
  float Scale = AlternativeSubtraction(Data_H1, Data_H2, lower, upper);

  // Perform the model fit. 
  Fit_Functions F;
  RooDataHist* Data = F.ConvertTH1toDataHist(Data_H1, Range);
  model.fitTo(*Data, SumW2Error(kTRUE));
 
  // Predicted events  
  float c1 = Variables.at(0) -> getVal();
  float c2 = Variables.at(1) -> getVal();
  float c3 = Variables.at(2) -> getVal();
  float c4 = Variables.at(3) -> getVal();

  // Closure test: Get all the mc entries for each class. 
  float W = Hist1.at(0) -> GetBinWidth(1); 
  float t1_H1 = Hist1.at(0) -> GetEntries(); // Truth was 1 trk
  float t2_H1 = Hist1.at(1) -> GetEntries(); // Truth was 2 trk
  float t3_H1 = Hist1.at(2) -> GetEntries(); // Truth was 3 trk
  float t4_H1 = Hist1.at(3) -> GetEntries(); // Truth was 4 trk
  float t1_H2 = Hist2.at(0) -> GetEntries(); // Truth was 1 trk
  float t2_H2 = Hist2.at(1) -> GetEntries(); // Truth was 2 trk
  float t3_H2 = Hist2.at(2) -> GetEntries(); // Truth was 3 trk
  float t4_H2 = Hist2.at(3) -> GetEntries(); // Truth was 4 trk

  // Truth scaled with the numbers being subtracted
  float t1 = t1_H1 - Scale*t1_H2;
  float t2 = t2_H1 - Scale*t2_H2;
  float t3 = t3_H1 - Scale*t3_H2;
  float t4 = t4_H1 - Scale*t4_H2;

  std::cout << t1 << "___" << t2 << "___" << t3 << "___" << t4 <<  std::endl;
  std::cout << c1/t1 << " ::::: " << c2/t2 << " ::::: " << c3/t3 << " ::::: " << c4/t4 << std::endl;

  // Ignore the title of the histogram
  std::vector<TString> Legend = {"Final Dist", "Subtract from", "Subtract"};
  std::vector<TH1F*> H = {Data_H1, Original1, Original2};

  auto f = new TCanvas();
  gPad -> SetLogy();
  gStyle -> SetOptStat(0);
  
  Data_H1 -> SetLineColor(kRed);
  Original1 -> SetLineColor(kBlue);
  Original2 -> SetLineColor(kOrange);

  TLegend* Legend_L = new TLegend(0.9, 0.9, 0.75, 0.75);
  for (unsigned int i = 0; i < Legend.size(); i++)
  {
    Legend_L -> AddEntry(H.at(i), Legend.at(i));
    H.at(i) -> Draw("HISTSAME");
  }
  Legend_L -> Draw();

  f -> Draw();
}

void Verification::Debug(std::vector<TH1F*> Hist, std::vector<std::vector<float>> Params)
{
  for (int x(0); x < Hist.size(); x++)
  {   
    TH1F* H = Hist.at(x); 
    std::vector<float> Par = Params.at(x);

    TF1 Lan("Lan", "landau", 0, 20);
    for (int i(0); i < Par.size(); i++)
    {
      Lan.SetParameter(i, Par.at(i));
    }
  
    for ( int i(0); i < 1000000; i++)
    {
      double r1 = Lan.GetRandom(); 
      H -> Fill(r1);  
    } 
  } 
}

