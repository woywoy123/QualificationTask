// Including commonly reused functions
#include<PostAnalysis/Functions.h>
#include<PostAnalysis/Verification.h>

// Including standard C++ libraries 
#include<iostream>

using namespace RooFit;

void PostAnalysis()
{
  auto FillHist = [](TH1F* Hist, std::vector<TH1F*> Hists, std::vector<float> Comps)
  {
    for ( int i(0); i < Comps.size(); i++)
    {
      Hist -> Add(Hists.at(i), Comps.at(i));
    }
  };


  // Init the function that hosts common tools
  Functions F;
  Fit_Functions f;
  Verification V;
  Benchmark B; 
  Plot_Functions P;

  // ====================== Generate the data =================== // 
  // === 0. Parameters for the N-Tracks and fill the generated Hists 
  std::vector<float> Param1 = {1, 1, 0.5};
  std::vector<float> Param2 = {1, 2, 1};
  std::vector<float> Param3 = {1, 5, 10};
  std::vector<float> Param4 = {1, 6, 20};
  std::vector<std::vector<float>> Params = {Param1, Param2, Param3, Param4}; 
  std::vector<TString> nTrk = {"1trk", "2trk", "3trk", "4trk"}; 
  std::vector<TString> trk = {"trk1", "trk2", "trk3", "trk4"};
   
  // Generate the toys 
  std::vector<TH1F*> nTrkHist = F.MakeTH1F(nTrk, 500, 0, 20);
  std::vector<TH1F*> ntrk = F.MakeTH1F(trk, 500, 0, 20);
  V.Debug(nTrkHist, Params);
   
  // === 1.: Generate the fake n-Track datasets with closure  
  // 1-Track Data 
  std::vector<float> COMP1 = {0.8, 0.1, 0.075, 0.025}; 
  FillHist(ntrk.at(0), nTrkHist, COMP1);
  TH1F* trk1Data = ntrk.at(0); 
 
  // 2-Track Data 
  std::vector<float> COMP2 = {0.1, 0.7, 0.1, 0.1}; 
  FillHist(ntrk.at(1), nTrkHist, COMP2);
  TH1F* trk2Data = ntrk.at(1); 
   
  // 3-Track Data 
  std::vector<float> COMP3 = {0.01, 0.1, 0.79, 0.1}; 
  FillHist(ntrk.at(2), nTrkHist, COMP3);
  TH1F* trk3Data = ntrk.at(2); 
  
  // 4-Track Data 
  std::vector<float> COMP4 = {0., 0.02, 0.2, 0.78}; 
  FillHist(ntrk.at(3), nTrkHist, COMP4);
  TH1F* trk4Data = ntrk.at(3);  

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
  for (int i(0); i < 10; i++)
  {
 
    trk1_C -> Reset("ICES");
    trk2_C -> Reset("ICES");
    trk3_C -> Reset("ICES");
    trk4_C -> Reset("ICES");
    
    for (int x(0); x < 10; x++)
    {
      deconv = f.LRDeconvolution(Data_Vector, deconv, deconv, 0.75);
    }
    // === 1.1: Replace the tail of the deconv with the 1trk data 
    deconv = f.TailReplace(trk1Data, deconv, 0, 20);
    
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
