#include<PostAnalysis/Debug.h>
#include<PostAnalysis/DistributionGenerator.h>

void Proxy(TString sample, TString Files, TString Mode)
{
  
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(Files); 

  TFile* File = new TFile(Mode + ".root", "RECREATE"); 
  File -> mkdir(sample); 
  File -> cd(sample); 

  MVT M = F[sample]; 
  std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
  TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
  TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
  TH1F* trk2_start = M["ntrk_2_M_O"][0];
  TH1F* trk3_start = M["ntrk_3_M_O"][0];
  TH1F* trk4_start = M["ntrk_4_M_O"][0];
  std::vector<TH1F*> starter = {trk1_start, trk2_start, trk3_start, trk4_start};  
  if (Mode.Contains("_Subtract")) {SubtractData(starter, trk1_start, 0, false);}
  if (Mode.Contains("_Smooth")) { Smooth(trk1_start, 0.1); }

  std::vector<TH1F*> trk1_templates = BuildNtrkMtru(1, trk1_start, "_template", 2)[0]; 
    
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print("Package.pdf["); 
  ntrk_1_M -> SetTitle(sample); 
  PlotHists(ntrk_1_M, ntrk_1_T, can); 
  can -> Print("Package.pdf"); 
  
  TH1F* titlePage = (TH1F*)trk1_templates[0] -> Clone("example"); 
  titlePage -> Clear();  
  titlePage -> SetTitle("Example Templates for ntruth fits"); 
  PlotHists(titlePage, trk1_templates, can); 
  can -> Print("Package.pdf");
  can -> Clear(); 

  std::vector<TString> c_out; 
  for (int i(0); i < FitRanges_Names.size(); i++)
  {
    TString c = FitRanges_Names[i]; c += (" -> start: "); c += (Ranges[i][0]); c += (" end: "); c+= (Ranges[i][1]); 
    c_out.push_back(c); 
  }
  c_out.push_back(""); 

  std::vector<TString> c_normal; 
  std::pair<std::vector<MVF>, std::vector<TString>> Normal; 
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShift; 
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFT; 
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFTWidth; 
  
  Normal = StepFit(trk1_start, ntrk_1_T, delta, can, "NormalCase1", sample); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 
  Normal = StepFit(trk1_start, ntrk_1_T, delta, can, "NormalCase2", sample); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 
  Normal = StepFit(trk1_start, ntrk_1_T, delta, can, "NormalCase3", sample); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
  
  
  // ======================== Normalization + Shift x-a ====================== //
  NormalShift = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalCase1", sample); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 

  NormalShift = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalCase2", sample); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 

  NormalShift = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalCase3", sample); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 

  // ======================== Normalization + Shift FFT ====================== //
  NormalShiftFFT = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalFFTCase1", sample); 
  FitParameterError(NormalShiftFFT.first, &c_out, {"Normalization", "Mean"}, NormalShiftFFT.second);  
  CompareNormalization(Normal.first, NormalShiftFFT.first, NormalShiftFFT.second, &c_normal); 

  NormalShiftFFT = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalFFTCase2", sample); 
  FitParameterError(NormalShiftFFT.first, &c_out, {"Normalization", "Mean"}, NormalShiftFFT.second);  
  CompareNormalization(Normal.first, NormalShiftFFT.first, NormalShiftFFT.second, &c_normal); 

  NormalShiftFFT = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalFFTCase3", sample); 
  FitParameterError(NormalShiftFFT.first, &c_out, {"Normalization", "Mean"}, NormalShiftFFT.second);  
  CompareNormalization(Normal.first, NormalShiftFFT.first, NormalShiftFFT.second, &c_normal); 

  // ======================== Normalization + Shift FFT Width ====================== //
  NormalShiftFFTWidth = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalWidthFFTCase1", sample); 
  FitParameterError(NormalShiftFFTWidth.first, &c_out, {"Normalization", "Mean", "Stdev"}, NormalShiftFFTWidth.second);  
  CompareNormalization(Normal.first, NormalShiftFFTWidth.first, NormalShiftFFTWidth.second, &c_normal); 
 
  NormalShiftFFTWidth = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalWidthFFTCase2", sample); 
  FitParameterError(NormalShiftFFTWidth.first, &c_out, {"Normalization", "Mean", "Stdev"}, NormalShiftFFTWidth.second);  
  CompareNormalization(Normal.first, NormalShiftFFTWidth.first, NormalShiftFFTWidth.second, &c_normal); 
 
  NormalShiftFFTWidth = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalWidthFFTCase3", sample); 
  FitParameterError(NormalShiftFFTWidth.first, &c_out, {"Normalization", "Mean", "Stdev"}, NormalShiftFFTWidth.second);  
  CompareNormalization(Normal.first, NormalShiftFFTWidth.first, NormalShiftFFTWidth.second, &c_normal); 

  gDirectory -> mkdir("Results"); 
  gDirectory -> cd("Results"); 
  
  BulkWrite(Graphs_1Tru); 
  BulkWrite(Graphs_2Tru); 

  PlotGraphs(Graphs_1Tru, "Performance of Prediction vs Truth of 1-Track", can); 
  can -> Print("Package.pdf"); 

  PlotGraphs(Graphs_2Tru, "Performance of Prediction vs Truth of 2-Track", can); 
  can -> Print("Package.pdf"); 


  can -> Print("Package.pdf]");
  
  std::ofstream myfile; 
  myfile.open("NormalizationWithRanges.txt"); 
  for (TString S : c_out)
  {
    std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 

  std::ofstream myNormal; 
  myNormal.open("NormalizationComparison.txt"); 
  for (TString S : c_normal)
  {
    std::cout << S << std::endl;
    myNormal << S << "\n"; 
  }
  myNormal.close(); 

  File -> Close();
}

void SplitTest(TString Sample)
{

  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE(Sample); 

  TString name = "SplitAnalysis.pdf"; 
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print(name + "["); 
  for (MMVi x = F.begin(); x != F.end(); x++)
  {
    auto Plotter =[&] (std::vector<TH1F*> Fits, std::vector<TH1F*> Truth, TString current, TString Mode)
    {
      TString name = "SplitAnalysis.pdf"; 
      TCanvas* can = new TCanvas(); 
      can -> SetLogy(); 
      PlotHists(Fits, Truth, current + "_" + Mode, can); 
      can -> Print(name); 
      can -> Clear();
    };

    std::map<TString, std::vector<TH1F*>> M = F[x -> first]; 
    std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
    std::vector<TH1F*> ntrk_2_T = M["ntrk_2_T_I"]; 
    std::vector<TH1F*> ntrk_3_T = M["ntrk_3_T_I"]; 
    std::vector<TH1F*> ntrk_4_T = M["ntrk_4_T_I"]; 
    std::vector<std::vector<TH1F*>> TruthVector = { ntrk_1_T,  ntrk_2_T, ntrk_3_T, ntrk_4_T };

    std::vector<TH1F*> ntrk_1_T_IsSplit = M["ntrk_1_T_I_IsSplit"]; 
    std::vector<TH1F*> ntrk_2_T_IsSplit = M["ntrk_2_T_I_IsSplit"]; 
    std::vector<TH1F*> ntrk_3_T_IsSplit = M["ntrk_3_T_I_IsSplit"]; 
    std::vector<TH1F*> ntrk_4_T_IsSplit = M["ntrk_4_T_I_IsSplit"]; 
    std::vector<std::vector<TH1F*>> TruthVector_Split = { ntrk_1_T_IsSplit,  ntrk_2_T_IsSplit, ntrk_3_T_IsSplit, ntrk_4_T_IsSplit };

    std::vector<TH1F*> ntrk_1_T_NotSplit = M["ntrk_1_T_I_NotSplit"]; 
    std::vector<TH1F*> ntrk_2_T_NotSplit = M["ntrk_2_T_I_NotSplit"]; 
    std::vector<TH1F*> ntrk_3_T_NotSplit = M["ntrk_3_T_I_NotSplit"]; 
    std::vector<TH1F*> ntrk_4_T_NotSplit = M["ntrk_4_T_I_NotSplit"]; 
    std::vector<std::vector<TH1F*>> TruthVector_NoSplit = { ntrk_1_T_NotSplit,  ntrk_2_T_NotSplit, ntrk_3_T_NotSplit, ntrk_4_T_NotSplit };
    
    bool Skip = false; 
    for (int i(0); i < TruthVector.size(); i++)
    {
      for (int j(0); j < TruthVector[i].size(); j++)
      {
        if (TruthVector[i][j] -> GetEntries() < 1000){Skip = true;}
      }
      if (Skip){continue;}
      Plotter(TruthVector_Split[i], TruthVector_NoSplit[i], x -> first, "Split vs Not Split"); 
      Plotter(TruthVector[i], TruthVector_NoSplit[i], x -> first, "Nominal vs Not Split"); 
      Plotter(TruthVector[i], TruthVector_Split[i], x -> first, "Nominal vs Split"); 
    }
  }
  can -> Print(name + "]"); 

}

