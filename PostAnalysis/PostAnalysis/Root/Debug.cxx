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
  SmoothHist(trk1_start, 0); 
  
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

  VT trk1_Normal = BuildNtrkMtru(1, trk1_start, "_Normal", 2)[0]; 
  MVF N_Res = Normalization(ntrk_1_M, trk1_Normal, Params_N); 
  
  ntrk_1_M -> SetTitle("Fit of Templates using only Normalization Fit on 1-Track Data"); 
  PlotHists(ntrk_1_M, ntrk_1_T, trk1_Normal, can); 
  can -> Print("Package.pdf"); 

  // Now we perform a stress test of when the fit will fail 
  // Store original parameters
  // ======================== Normalization only fit ====================== //
  auto StepFit =[&](std::vector<TH1F*> templ, std::vector<TH1F*> truth, float delta, float S_L, TCanvas* can)
  {
    int steps = 100 / delta; 
    float cu = 0; 

    std::vector<MVF> out; 
    std::vector<TString> setting; 
    for (int i(0); i < 2*steps; i++)
    {
      cu += delta; 
      TString title = "Track 1 Normalization Fit with Truth-2 being "; 
      title += (cu); title += "% of Original Normalization";
      TString name = "_ntru_2_"; name += (cu); name += "%"; 

      std::cout << "=========== >" + title << std::endl;
      // Scale ntrk1_tru2 by 1%
      Normalize(templ); 
      Normalize(truth[1]); 
      truth[1] -> Scale(S_L*(cu / 100)); 
      TH1F* ntrk_1 = SumHists(truth, title); 
      ntrk_1 -> SetTitle(title);
      ntrk_1 -> SetLineColor(kBlack);
      MVF N_Res = Normalization(ntrk_1, templ, Params_N); 
      PlotHists(ntrk_1, {truth[0], truth[1]}, templ, can); 
      can -> Print("Package.pdf"); 
      can -> Clear(); 
      
      out.push_back(N_Res); 
      setting.push_back(title);

      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
    }
    return std::pair<std::vector<MVF>, std::vector<TString>>(out, setting); 
  }; 

  float S_L = ntrk_1_T[1] -> Integral(); 
  Normalize(ntrk_1_T[1]); 
  
  std::pair<std::vector<MVF>, std::vector<TString>> Normal = StepFit(trk1_Normal, ntrk_1_T, 10, S_L, can); 
  std::vector<TString> c_out; 
  
  for (int i(0); i < FitRanges_Names.size(); i++)
  {
    TString c = FitRanges_Names[i]; c += (" -> start: "); c += (Ranges[i][0]); c += (" end: "); c+= (Ranges[i][1]); 
    c_out.push_back(c); 
  }
  
  c_out.push_back(""); 
  FitParameterError(Normal.first, &c_out, {"Normalization"}, Normal.second);  
 

  // ======================== Normalization + Shift x-a ====================== //
  auto StepFitShift =[&](TH1F* start, std::vector<TH1F*> truth, float delta, float S_L, TCanvas* can)
  {
    int steps = 200 / delta; 
    float cu = 0; 

    std::vector<MVF> out; 
    std::vector<TString> setting; 
    for (int i(0); i < steps; i++)
    {

      cu += delta; 
      TString title = "Track 1 Normalization + Shifting Fit with Truth-2 being "; 
      title += (cu); title += "% of Original Normalization";
      TString name = "_ntru_2_"; name += (cu); name += "%"; 
      VT templ = BuildNtrkMtru(1, start, "_Norm_Shift", 2)[0];
   
      std::cout << "=========== >" + title << std::endl;
      // Scale ntrk1_tru2 by 1%
      Normalize(templ); 
      Normalize(truth[1]); 
      truth[1] -> Scale(S_L*(cu / 100)); 
      TH1F* ntrk_1 = SumHists(truth, name); 
      ntrk_1 -> SetTitle(title);
      ntrk_1 -> SetLineColor(kBlack);
      MVF N_Res = NormalizationShift(ntrk_1, templ, Params_NS); 
      PlotHists(ntrk_1, {truth[0], truth[1]}, templ, can); 
      can -> Print("Package.pdf"); 
      can -> Clear(); 
      out.push_back(N_Res); 
      setting.push_back(title);

      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;

    }
    return std::pair<std::vector<MVF>, std::vector<TString>>(out, setting); 
  }; 

  Normalize(ntrk_1_T[1]); 
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShift = StepFitShift(trk1_start, ntrk_1_T, 10, S_L, can); 
  FitParameterError(NormalShift.first, &c_out, {"Normalization", "Shift"}, NormalShift.second);  


  // ======================== Normalization + Shift FFT ====================== //
  auto StepFitFFTShift =[&](TH1F* start, std::vector<TH1F*> truth, float delta, float S_L, TCanvas* can)
  {
    int steps = 100 / delta; 
    float cu = 0; 

    std::vector<MVF> out; 
    std::vector<TString> setting; 
    for (int i(0); i < 2*steps; i++)
    {

      cu += delta; 
      TString title = "Track 1 Normalization + Shifting using FFT in Fit with Truth-2 being "; 
      title += (cu); title += "% of Original Normalization";
      TString name = "_ntru_2_"; name += (cu); name += "%"; 
      VT templ = BuildNtrkMtru(1, start, "_Norm_Shift", 2)[0];
   
      std::cout << "=========== >" + title << std::endl;
      // Scale ntrk1_tru2 by 1%
      Normalize(templ); 
      Normalize(truth[1]); 
      truth[1] -> Scale(S_L*(cu / 100)); 
      TH1F* ntrk_1 = SumHists(truth, name); 
      ntrk_1 -> SetTitle(title);
      ntrk_1 -> SetLineColor(kBlack);
      MVF N_Res = ConvolutionFFT(ntrk_1, templ, Params_FFT); 
      PlotHists(ntrk_1, {truth[0], truth[1]}, templ, can); 
      can -> Print("Package.pdf"); 
      can -> Clear(); 
      
      out.push_back(N_Res); 
      setting.push_back(title);

      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;

    }
    return std::pair<std::vector<MVF>, std::vector<TString>>(out, setting); 
  }; 

  Normalize(ntrk_1_T[1]); 
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFT = StepFitFFTShift(trk1_start, ntrk_1_T, 10, S_L, can); 
  FitParameterError(NormalShiftFFT.first, &c_out, {"Normalization", "Mean"}, NormalShiftFFT.second);  

  // ======================== Normalization + Shift FFT Width ====================== //
  auto StepFitFFTShiftWidth =[&](TH1F* start, std::vector<TH1F*> truth, float delta, float S_L, TCanvas* can)
  {
    int steps = 100 / delta; 
    float cu = 0; 

    std::vector<MVF> out; 
    std::vector<TString> setting; 
    for (int i(0); i < 2*steps; i++)
    {

      cu += delta; 
      TString title = "Track 1 Normalization + Shifting + Width using FFT in Fit with Truth-2 being "; 
      title += (cu); title += "% of Original Normalization";
      TString name = "_ntru_2_"; name += (cu); name += "%"; 
      VT templ = BuildNtrkMtru(1, start, "_Norm_Shift_Width", 2)[0];
   
      std::cout << "=========== >" + title << std::endl;
      // Scale ntrk1_tru2 by 1%
      Normalize(templ); 
      Normalize(truth[1]); 
      truth[1] -> Scale(S_L*(cu / 100)); 
      TH1F* ntrk_1 = SumHists(truth, name); 
      ntrk_1 -> SetTitle(title);
      ntrk_1 -> SetLineColor(kBlack);
      MVF N_Res = ConvolutionFFT(ntrk_1, templ, Params_WidthFFT); 
      PlotHists(ntrk_1, {truth[0], truth[1]}, templ, can); 
      can -> Print("Package.pdf"); 
      can -> Clear(); 
      
      out.push_back(N_Res); 
      setting.push_back(title);

      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "" << std::endl;

    }
    return std::pair<std::vector<MVF>, std::vector<TString>>(out, setting); 
  }; 

  Normalize(ntrk_1_T[1]); 
  std::pair<std::vector<MVF>, std::vector<TString>> NormalShiftFFTWidth = StepFitFFTShiftWidth(trk1_start, ntrk_1_T, 10, S_L, can); 
  FitParameterError(NormalShiftFFTWidth.first, &c_out, {"Normalization", "Mean", "Stdev"}, NormalShiftFFTWidth.second);  


  can -> Print("Package.pdf]");
  
  std::ofstream myfile; 
  myfile.open("NormalizationWithRanges.txt"); 
  for (TString S : c_out)
  {
    //std::cout << S << std::endl;
    myfile << S << "\n"; 
  }
  myfile.close(); 

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

