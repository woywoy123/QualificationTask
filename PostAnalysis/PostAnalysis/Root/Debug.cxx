#include<PostAnalysis/Debug.h>

void Proxy()
{

  TString sample = "Blayer_1400_1600_GeV"; 
  TString run = "Normalization"; 
  std::map<TString, std::map<TString, std::vector<TH1F*>>> F = ReadCTIDE("Merged_MC.root"); 
  MVT M = F[sample]; 
  
  std::vector<TH1F*> ntrk_1_T = M["ntrk_1_T_I"]; 
  TH1F* ntrk_1_M = M["ntrk_1_M_I"][0]; 
  TH1F* trk1_start = M["ntrk_1_M_O"][0]; 

  std::vector<TH1F*> trk1_templates = BuildNtrkMtru(1, trk1_start, "_template", 2)[0]; 
  //SmoothHist(trk1_start, 0); 
  
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
  std::pair<std::vector<MVF>, std::vector<TString>> Normal = StepFit(trk1_start, ntrk_1_T, 10, can, "Normal"); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 
  // ======================== Normalization + Shift x-a ====================== //
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShift = StepFit(trk1_start, ntrk_1_T, 10, can, "ShiftNormal"); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 

  // ======================== Normalization + Shift FFT ====================== //
  //std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFT = StepFit(trk1_start, ntrk_1_T, 10, can, "ShiftNormalFFT"); 
  //FitParameterError(NormalShiftFFT.first, &c_out, {"Normalization", "Mean"}, NormalShiftFFT.second);  
  //CompareNormalization(Normal.first, NormalShiftFFT.first, NormalShiftFFT.second, &c_normal); 

  //// ======================== Normalization + Shift FFT Width ====================== //
  //std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFTWidth = StepFit(trk1_start, ntrk_1_T, 10, can, "ShiftNormalWidthFFT"); 
  //FitParameterError(NormalShiftFFTWidth.first, &c_out, {"Normalization", "Mean", "Stdev"}, NormalShiftFFTWidth.second);  
  //CompareNormalization(Normal.first, NormalShiftFFTWidth.first, NormalShiftFFTWidth.second, &c_normal); 

  can -> Print("Package.pdf]");
  
  std::ofstream myfile; 
  myfile.open("NormalizationWithRanges.txt"); 
  for (TString S : c_out)
  {
    //std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 

  std::ofstream myNormal; 
  myNormal.open("NormalizationComparison.txt"); 
  for (TString S : c_normal)
  {
    //std::cout << S << std::endl;
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
  TH1F* H = (TH1F*)trk1_start -> Clone("Clone"); 
  SmoothHist(trk1_start, 0); 

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  trk1_start -> GetYaxis() -> SetRangeUser(1, 1e5); 
  trk1_start -> Draw("HIST");
  H -> Draw("SAMEHIST"); 
  can -> Print("Smooth.pdf"); 
  delete can; 

}

