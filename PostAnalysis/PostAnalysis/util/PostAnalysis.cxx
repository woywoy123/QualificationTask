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
  std::vector<TH1F*> nTrkHist = F.MakeTH1F(nTrk, 50, 0, 20);
  std::vector<TH1F*> ntrk = F.MakeTH1F(trk, 50, 0, 20);
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
  // === 3. Combine all the data hists together for easy access
  std::vector<TH1F*> DataHists = {trk1Data, trk2Data, trk3Data, trk4Data}; 

  // ===================== END OF DATA ========================== // 




  // =========== Initial Estimate of contamination ============== //
  // === 0. Perform a simple fit with the hists that are given. 
  // === 0.1 Create the variables  
  RooRealVar* x = new RooRealVar("x", "x", 0, 20); 
  std::vector<RooRealVar*> var = f.GenerateVariables(Constants::Variable_Names, Constants::Begin, Constants::End);
  
  // === 0.2 Create the PDFs with the measured n-tracks
  std::vector<RooHistPdf*> pdf = f.ConvertTH1FtoPDF(DataHists, x);  
  RooAddPdf fit("fit", "fit", f.VectorToArgList(pdf), f.VectorToArgList(var));
  
  // === 0.3 Import the data and fit the model   
  RooDataHist* Data_Fit = f.ConvertTH1toDataHist(Data_Copy, x);
  fit.fitTo(*Data_Fit);
  
  // === 1. Check the varibles and verify that they match the set scales
  // === 1.1 Collect from variables  
  float p1 = var.at(0) -> getVal();
  float p2 = var.at(1) -> getVal();
  float p3 = var.at(2) -> getVal();
  float p4 = var.at(3) -> getVal(); 
  
  // === 2. Use the pn values to make a rough estimate of the number of entries are due to cross contamination 
  // === 2.1: Scale the number of entries according to Trk numbers and scale with Data Copy lumi 
  float Data_Lumi = Data_Copy -> Integral();
  float p1s = p1*(1/Data_Lumi);
  float p2s = 2*p2*(1/Data_Lumi);
  float p3s = 3*p3*(1/Data_Lumi);
  float p4s = 4*p4*(1/Data_Lumi);
  std::vector<float> Scales = {p1s, p2s, p3s, p4s};

  // === 2.2: Copy the track data
  TH1F* trk1_C = (TH1F*)trk1Data -> Clone("trk1_C");
  TH1F* trk2_C = (TH1F*)trk2Data -> Clone("trk2_C");
  TH1F* trk3_C = (TH1F*)trk3Data -> Clone("trk3_C");
  TH1F* trk4_C = (TH1F*)trk4Data -> Clone("trk4_C");
  std::vector<TH1F*> Hist_Copy = {trk1_C, trk2_C, trk3_C, trk4_C};

  // === 2.3: Subtract the scaled tracks from the Data Copy 
  for (int i(0); i < Scales.size(); i++)
  {
    if ( i == NumTrak -1 ) { continue; }
    else
    {
      Hist_Copy.at(i) -> Scale(Scales[i]);
      Data_Copy -> Add(Hist_Copy.at(i), -1);
    }
  }

  // === 3. Copy the contents of the TH1F into a vector and create the mirror
  int nBins = Data_Copy -> GetNbinsX();
  int OS = Offset*nBins;
  std::vector<float> Data_Vector(nBins + OS, 0);
  
  for (int i(0); i < nBins; i++)
  {
    Data_Vector[i] = Data_Copy -> GetBinContent(i+1);
    if ( i < OS) { Data_Vector[nBins + i] = Data_Copy -> GetBinContent(nBins -i -1); } 
  }

  // === 4. Reset any useful Histograms and set their colors 
  trk1_C -> Reset();
  trk2_C -> Reset();
  trk3_C -> Reset();
  trk4_C -> Reset();
  
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




  // ===================== Start of Deconvolution =============== //
   






  //// ======================================== Experimental =================================== // 
  //int nBins = nTrk_H.at(0) -> GetNbinsX(); 
  //int Offset = nBins*0.1;

  ////nTrk_H.at(1) -> Draw("SAMEHIST");

  //std::vector<float> data(nBins + Offset, 0);
  //std::vector<float> deconv(nBins + Offset, 1);

  //// Here he gets the entries from hist and stores into a vector.
  //for (size_t i(0); i < nBins; i++)
  //{
  //  data[i] = nTrk_H.at(1) -> GetBinContent(i+1);
  //  ProbeV.push_back(nTrk_H.at(1) -> GetBinContent(i+1)); 
  //}

  //// Does this weird tail inversion. Not sure why. 
  //for (int i(0); i < Offset; i++)
  //{
  //  data[ nBins + i ] = data[ nBins -1 -i ];
  //}


  //// *** Main Algorithm Part ================== //
  //TGraph* gr = new TGraph(); 
  //for (int i(0); i < 100; i++)
  //{
  //  deconv = f.LRDeconvolution(data, deconv, deconv, 0.75); 
  //  F.VectorToTH1F(deconv, Probe); 
  //  f.ConvolveHists(Probe, Probe, trk2, Offset);
  //  //Probe -> Draw("SAMEHIST*"); 

  //  std::vector<float> Prog;
  //  for (int i(0); i < trk2 -> GetNbinsX(); i++)
  //  {
  //    Prog.push_back(trk2 -> GetBinContent(i+1));
  //  }
  //   
  //  float dist = B.WeightedEuclidean(ProbeV, Prog);  
  //  gr -> SetPoint(i, i, dist);
  //  gr -> Draw("*ap");
  //  can -> Update();
  //}

 













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
