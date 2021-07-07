#include<PostAnalysis/Debug.h>
#include<PostAnalysis/DistributionGenerator.h>

void Proxy(TString sample)
{

  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE("Merged_MC.root"); 
  MVT M = F[sample]; 
  
  std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
  TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
  TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
  TH1F* trk2_start = M["ntrk_2_M_O"][0];
  TH1F* trk3_start = M["ntrk_3_M_O"][0];
  TH1F* trk4_start = M["ntrk_4_M_O"][0];
  std::vector<TH1F*> starter = {trk1_start, trk2_start, trk3_start, trk4_start};  
  SubtractData(starter, trk1_start, 0, false); 

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
}

void SmoothingTest()
{
  TString sample = "Blayer_1400_1600_GeV"; 
  TString run = "Normalization"; 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE("Merged_MC.root"); 
  MVT M = F[sample]; 
  TH1F* trk1_start = M["ntrk_1_M_O"][0]; 
  TH1F* trk1 = (TH1F*)trk1_start -> Clone("Clone"); 

  TH1F* Small = Snipping(trk1, 0.4, 8); 
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print("Smooth.pdf["); 
  Small -> Draw("HIST");  
  can -> Print("Smooth.pdf");
  can -> Clear(); 
  trk1 -> Draw("HIST"); 
  can -> Print("Smooth.pdf");
  
  TH1F* Ori = (TH1F*)trk1 -> Clone("Ori"); 
  Ori -> SetLineColor(kBlack); 
  RevertSnipping(Small, Ori);  
  can -> Clear(); 
  trk1 -> Draw("HIST"); 
  Ori -> Draw("SAMEHIST");  
  can -> Print("Smooth.pdf");

  can -> Clear(); 
  can -> Print("Smooth.pdf]");
  
  for (int i(0); i < Ori -> GetNbinsX(); i++)
  {
    float e = Ori -> GetBinContent(i+1); 
    float f = trk1 -> GetBinContent(i+1); 
    std::cout << e-f << std::endl;
  }





  delete can; 

}

