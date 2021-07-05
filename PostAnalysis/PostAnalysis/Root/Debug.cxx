#include<PostAnalysis/Debug.h>
#include<PostAnalysis/DistributionGenerator.h>

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
  //ReplaceTail(trk1_templates[0], "./Smoothed.root", sample, "H");
    
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
  float delta = 10;  

  std::vector<TString> c_normal; 
  std::pair<std::vector<MVF>, std::vector<TString>> Normal = StepFit(trk1_start, ntrk_1_T, delta, can, "NormalCase1", sample); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 
  // ======================== Normalization + Shift x-a ====================== //
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShift = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalCase1", sample); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 


  Normal = StepFit(trk1_start, ntrk_1_T, delta, can, "NormalCase2", sample); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 
  // ======================== Normalization + Shift x-a ====================== //
  NormalShift = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalCase2", sample); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 

  Normal = StepFit(trk1_start, ntrk_1_T, delta, can, "NormalCase3", sample); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 
  // ======================== Normalization + Shift x-a ====================== //
  NormalShift = StepFit(trk1_start, ntrk_1_T, delta, can, "ShiftNormalCase3", sample); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  
  CompareNormalization(Normal.first, NormalShift.first, NormalShift.second, &c_normal); 

  //// ======================== Normalization + Shift FFT ====================== //
  //std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFT = StepFit(trk1_start, ntrk_1_T, 100, can, "ShiftNormalFFT", sample); 
  //FitParameterError(NormalShiftFFT.first, &c_out, {"Normalization", "Mean"}, NormalShiftFFT.second);  
  //CompareNormalization(Normal.first, NormalShiftFFT.first, NormalShiftFFT.second, &c_normal); 

  //// ======================== Normalization + Shift FFT Width ====================== //
  //std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFTWidth = StepFit(trk1_start, ntrk_1_T, 100, can, "ShiftNormalWidthFFT", sample); 
  //FitParameterError(NormalShiftFFTWidth.first, &c_out, {"Normalization", "Mean", "Stdev"}, NormalShiftFFTWidth.second);  
  //CompareNormalization(Normal.first, NormalShiftFFTWidth.first, NormalShiftFFTWidth.second, &c_normal); 
 

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
   
  //Smooth1Trk(trk1_start, 200); 
  //WriteHistsToFile({trk1_start}, sample); 
  //FS -> Close(); 
  
  //SmoothHist(trk1, 10, 0.25); 

  //float min = trk1 -> GetXaxis() -> GetXmin(); 
  //float max = trk1 -> GetXaxis() -> GetXmax();

  //TGraph* gr = new TGraph(trk1);  
  //TGraphSmooth* g = new TGraphSmooth("Smoother"); 

  //TGraph* h = g -> SmoothKern(gr, "normal", 0.1); 

  //int n = h -> GetN(); 
  //for (int i(0); i < n; i++)
  //{
  //  double x, y; 
  //  h -> GetPoint(i, x, y); 
  //  trk1 -> SetBinContent(i+1, y);
  //}
  //Average(trk1);

  //
  ////ReplaceTail(trk1, trk1_start);
  // 
  ////ConvolutionFFT(trk1, {Lan}, Params_WidthFFT);
  //trk1 -> SetLineColor(kBlack);



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

