#include<PostAnalysis/ReportFigures.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/IO.h>
#include<TGraphSmooth.h>

void Figure_Proxy()
{
  Figure_N_TrackTemplates();
  Figure_Trk1_Subtraction_Justification(); 

}

void Figure_N_TrackTemplates()
{

  std::map<TString, std::map<TString, std::vector<TH1F*>>> M = ReadCTIDE("Merged_MC.root");  
  TH1F* Blayer_trk1_outside = M[example]["ntrk_1_M_O"][0]; 
  TH1F* Blayer_trk2_outside = M[example]["ntrk_2_M_O"][0]; 
  TH1F* Blayer_trk3_outside = M[example]["ntrk_3_M_O"][0]; 
  TH1F* Blayer_trk4_outside = M[example]["ntrk_4_M_O"][0]; 
  
  std::vector<TH1F*> Blayer_trk1_O_M = {Blayer_trk1_outside, Blayer_trk2_outside, Blayer_trk3_outside, Blayer_trk4_outside};
  SubtractData(Blayer_trk1_O_M, Blayer_trk1_outside, 0, false); 
  std::vector<TH1F*> ntrk = BuildNtrkMtru(4, Blayer_trk1_outside, "T")[0];  

  for (int i(0); i < ntrk.size(); i++)
  {
    TString na; na += (i+1); na += ("-Track Template"); 
    ntrk[i] -> SetTitle(na);
  }
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  GenerateNiceStacks(ntrk, "Example n-Track Template Distributions from Blayer at Jet PT 1200 - 1400 GeV", can, "dE/dx [MeV g^{-1} cm^{2}]", "", "nostack hist"); 
  can -> Print("ExampleTemplates.pdf"); 
}

void Figure_Trk1_Subtraction_Justification()
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> M = ReadCTIDE("Merged_MC.root");  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  TString name = "Outside_Subtraction_Justification.pdf"; 
  can -> Print(name + "["); 

  // per layer error maps 
  std::map<TString, float> IBL_Error_GS; 
  std::map<TString, float> IBL_Error_GNS; 
  std::map<TString, float> IBL_Error_Smo; 
  
  std::map<TString, float> Blayer_Error_GS; 
  std::map<TString, float> Blayer_Error_GNS; 
  std::map<TString, float> Blayer_Error_Smo; 

  std::map<TString, float> Layer1_Error_GS; 
  std::map<TString, float> Layer1_Error_GNS; 
  std::map<TString, float> Layer1_Error_Smo; 

  std::map<TString, float> Layer2_Error_GS; 
  std::map<TString, float> Layer2_Error_GNS; 
  std::map<TString, float> Layer2_Error_Smo; 

  std::map<TString, float> Combined_Error_GS; 
  std::map<TString, float> Combined_Error_GNS; 
  std::map<TString, float> Combined_Error_Smo; 

  for (MMVi h = M.begin(); h != M.end(); h++)
  {
    // Outside 1-track: We use this to initialize the fitting procedure 
    std::map<TString, std::vector<TH1F*>> Hists = M[h -> first]; 
    TH1F* ntrk1_O = Hists["ntrk_1_M_O"][0]; 
    TH1F* ntrk2_O = Hists["ntrk_2_M_O"][0]; 
    TH1F* ntrk3_O = Hists["ntrk_3_M_O"][0]; 
    TH1F* ntrk4_O = Hists["ntrk_4_M_O"][0]; 
    std::vector<TH1F*> Outside_v = {ntrk1_O, ntrk2_O, ntrk3_O, ntrk4_O}; 
    
    TString n_ = "No Subtraction ntrk-1"; 
    ntrk1_O -> SetTitle(n_); 

    TString n_temp = "Subtraction - ntrk 1"; 
    TH1F* ntrk1_O_S = (TH1F*)ntrk1_O -> Clone(n_temp); 
    ntrk1_O_S -> SetTitle(n_temp); 

    TH1F* ntrk1_O_Smo = (TH1F*)ntrk1_O -> Clone("No Subtract ntrk-1 SMOOTH"); 
    ntrk1_O_Smo -> SetTitle("No Subtract ntrk-1 SMOOTH");

    // Inside 1-track 1-Truth: Check how well the histograms match 
    TH1F* trk1_tru1 = Hists["ntrk_1_T_I"][0]; 
    if (trk1_tru1 -> Integral() == 0) { continue; }

    // Subtraction of 1-track outside 
    SubtractData(Outside_v, ntrk1_O_S, 0); 
    
    // Scale the outside measured histograms to the same integral as truth 
    TH1F* tmp = (TH1F*)ntrk1_O_S -> Clone("TMP"); 
    Normalize(tmp); 
    ntrk1_O_S -> Reset(); 
    ntrk1_O_S -> FillRandom(tmp, trk1_tru1 -> GetEntries()); 
    delete tmp; 

    tmp = (TH1F*)ntrk1_O -> Clone("TMP"); 
    Normalize(tmp); 
    ntrk1_O -> Reset(); 
    ntrk1_O -> FillRandom(tmp, trk1_tru1 -> GetEntries()); 
    delete tmp; 

    tmp = (TH1F*)ntrk1_O_Smo -> Clone("TMP"); 
    Normalize(tmp); 
    ntrk1_O_Smo -> Reset(); 
    ntrk1_O_Smo -> FillRandom(tmp, trk1_tru1 -> GetEntries()); 
    Smooth(ntrk1_O_Smo, 0.1); 
    delete tmp; 

    PlotHists(trk1_tru1, {ntrk1_O, ntrk1_O_S, ntrk1_O_Smo}, can); 
    can -> Print(name); 
    can -> Clear(); 

    // Normalize hists 
    Normalize(trk1_tru1); 
    Normalize(ntrk1_O_S); 
    Normalize(ntrk1_O); 
    Normalize(ntrk1_O_Smo); 

    // Compare subtracted vs non subtracted to truth 
    float err_s = ChiSquareError(trk1_tru1, ntrk1_O_S) * 100; 
    float err_ns = ChiSquareError(trk1_tru1, ntrk1_O) * 100;
    float err_ns_smo = ChiSquareError(trk1_tru1, ntrk1_O_Smo) * 100; 
   
    TString L_JPT = h -> first;
    if (L_JPT.Contains("IBL"))
    {
      TString l = L_JPT.ReplaceAll("IBL_", "");
      IBL_Error_GS[l] = err_s; 
      IBL_Error_GNS[l] = err_ns; 
      IBL_Error_Smo[l] = err_ns_smo; 
      continue; 
    }

    if (L_JPT.Contains("Blayer"))
    {
      TString l = L_JPT.ReplaceAll("Blayer_", ""); 
      Blayer_Error_GS[l] = err_s; 
      Blayer_Error_GNS[l] = err_ns; 
      Blayer_Error_Smo[l] = err_ns_smo; 
      continue; 
    }

    if (L_JPT.Contains("layer1"))
    {
      TString l = L_JPT.ReplaceAll("layer1_", ""); 
      Layer1_Error_GS[l] = err_s; 
      Layer1_Error_GNS[l] = err_ns; 
      Layer1_Error_Smo[l] = err_ns_smo; 
      continue; 
    }

    if (L_JPT.Contains("layer2"))
    {
      TString l = L_JPT.ReplaceAll("layer2_", ""); 
      Layer2_Error_GS[l] = err_s; 
      Layer2_Error_GNS[l] = err_ns; 
      Layer2_Error_Smo[l] = err_ns_smo; 
      continue; 
    }

    if (L_JPT.Contains("_"))
    {
      TString l = L_JPT; //.ReplaceAll("_", " "); 
      Combined_Error_GS[l] = err_s; 
      Combined_Error_GNS[l] = err_ns; 
      Combined_Error_Smo[l] = err_ns_smo; 
      continue; 
    }
  }

  float min = 1e-4;
  float max = 1000; 
  can -> Print(name + "]"); 
  can -> Clear();
  can -> SetLogy(false);
  can -> Print("ShapeToInsideTruth.pdf[");
  gStyle -> SetOptStat(false); 

  TGraph* IBL = GenerateGraph(IBL_Error_GNS, "IBL No Subtraction"); 
  TGraph* IBL_sub = GenerateGraph(IBL_Error_GS, "IBL Subtraction"); 
  TGraph* IBL_Smo = GenerateGraph(IBL_Error_Smo, "IBL No Subtraction Smooth"); 
  std::vector<TGraph*> ng = {IBL, IBL_sub, IBL_Smo}; 
  GeneratePerformanceGraphs(ng, "IBL Shape Performance", "Jet Energy (GeV)", "Error (%)", min, max, can);
  can -> Print("ShapeToInsideTruth.pdf");
  can -> Clear(); 
  
  can -> SetLogy(false);
  TGraph* Blayer = GenerateGraph(Blayer_Error_GNS, "BLayer No Subtraction"); 
  TGraph* Blayer_sub = GenerateGraph(Blayer_Error_GS, "BLayer Subtraction"); 
  TGraph* Blayer_Smo = GenerateGraph(Blayer_Error_Smo, "BLayer No Subtraction Smooth"); 
  ng = {Blayer, Blayer_sub, Blayer_Smo}; 
  GeneratePerformanceGraphs(ng, "BLayer Shape Performance", "Jet Energy (GeV)", "Error (%)", min, max, can);
  can -> Print("ShapeToInsideTruth.pdf");
  can -> Clear(); 

  can -> SetLogy(false);
  TGraph* Layer1 = GenerateGraph(Layer1_Error_GNS, "Layer1 No Subtraction"); 
  TGraph* Layer1_sub = GenerateGraph(Layer1_Error_GS, "Layer1 Subtraction"); 
  TGraph* Layer1_Smo = GenerateGraph(Layer1_Error_Smo, "Layer1 No Subtraction Smooth"); 
  ng = {Layer1, Layer1_sub, Layer1_Smo}; 
  GeneratePerformanceGraphs(ng, "Layer1 Shape Performance", "Jet Energy (GeV)", "Error (%)", min, max, can);
  can -> Print("ShapeToInsideTruth.pdf");
  can -> Clear(); 

  can -> SetLogy(false);
  TGraph* Layer2 = GenerateGraph(Layer2_Error_GNS, "Layer2 No Subtraction"); 
  TGraph* Layer2_sub = GenerateGraph(Layer2_Error_GS, "Layer2 Subtraction"); 
  TGraph* Layer2_Smo = GenerateGraph(Layer2_Error_Smo, "Layer2 No Subtraction Smooth"); 
  ng = {Layer2, Layer2_sub, Layer2_Smo}; 
  GeneratePerformanceGraphs(ng, "Layer2 Shape Performance", "Jet Energy (GeV)", "Error (%)", min, max, can);
  can -> Print("ShapeToInsideTruth.pdf");

  can -> SetLogy(false);
  TGraph* Comb = GenerateGraph(Combined_Error_GNS, "Combined Layers No Subtraction"); 
  TGraph* Comb_sub = GenerateGraph(Combined_Error_GS, "Combined Layers Subtraction"); 
  TGraph* Comb_Smo = GenerateGraph(Combined_Error_Smo, "Combined Layers No Subtraction Smooth"); 
  ng = {Comb, Comb_sub, Comb_Smo}; 
  GeneratePerformanceGraphs(ng, "Combined Shape Performance", "Jet Energy (GeV)", "Error (%)", min, max, can);
  can -> Print("ShapeToInsideTruth.pdf");

  can -> Print("ShapeToInsideTruth.pdf]"); 
}
