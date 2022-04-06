#include "../PlottingCode/Debugging.h"

void CaseAnalysis(maps debug, _TS LE, float delta, float start, float end, bool Plotter, 
                  std::map<_TS, std::vector<float>>* Stepper1, 
                  std::map<_TS, std::vector<float>>* Stepper2,
                  std::map<_TS, std::vector<float>>* Stepper3
)
{
  std::vector<_TS> v = Split(LE, "-"); 
  _TS Layer = v[0]; 
  _TS Energy = v[1];
  LE.ReplaceAll("-", "_");
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy();
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(3); 
  can -> SetTopMargin(0.1); 
  for (int n(start); n < end;)
  {
    n = n + delta;
    TString step = ""; step += (n);

    // Truth is the same for all cases except the steps
    TH1F* ntrk1_ntru1_T = debug[LE][step]["NormalCase1"]["ntrk1-ntru1-T"]; 
    TH1F* ntrk1_ntru2_T = debug[LE][step]["NormalCase1"]["ntrk1-ntru2-T"];
    std::vector<TH1F*> Truth = {ntrk1_ntru1_T, ntrk1_ntru2_T}; 

    // ====== Case1: Use 1-Truth to construct templates
    // -> Normalization Fits
    TH1F* ntrk1_ntru1_D_N1 = debug[LE][step]["NormalCase1"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_N1 = debug[LE][step]["NormalCase1"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> N1 = {ntrk1_ntru1_D_N1, ntrk1_ntru2_D_N1};
    
    // -> Normalization Shifting Fits
    TH1F* ntrk1_ntru1_D_NS1 = debug[LE][step]["ShiftNormalCase1"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_NS1 = debug[LE][step]["ShiftNormalCase1"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> NS1 = {ntrk1_ntru1_D_NS1, ntrk1_ntru2_D_NS1};

    // -> Normalization Shifting FFT Fits
    TH1F* ntrk1_ntru1_D_NSFFT1 = debug[LE][step]["ShiftNormalFFTCase1"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_NSFFT1 = debug[LE][step]["ShiftNormalFFTCase1"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> NSFFT1 = {ntrk1_ntru1_D_NSFFT1, ntrk1_ntru2_D_NSFFT1}; 

    // ====== Case2: Use 2-Truth to perform the fit
    // -> Normalization Fits
    TH1F* ntrk1_ntru1_D_N2 = debug[LE][step]["NormalCase2"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_N2 = debug[LE][step]["NormalCase2"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> N2 = {ntrk1_ntru1_D_N2, ntrk1_ntru2_D_N2};

    // -> Normalization Shifting Fits
    TH1F* ntrk1_ntru1_D_NS2 = debug[LE][step]["ShiftNormalCase2"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_NS2 = debug[LE][step]["ShiftNormalCase2"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> NS2 = {ntrk1_ntru1_D_NS2, ntrk1_ntru2_D_NS2};

    // -> Normalization Shifting FFT Fits
    TH1F* ntrk1_ntru1_D_NSFFT2 = debug[LE][step]["ShiftNormalFFTCase2"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_NSFFT2 = debug[LE][step]["ShiftNormalFFTCase2"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> NSFFT2 = {ntrk1_ntru1_D_NSFFT2, ntrk1_ntru2_D_NSFFT2}; 

    // ====== Case3: Use normal templates to perform the fit
    // -> Normalization Fits
    TH1F* ntrk1_ntru1_D_N3 = debug[LE][step]["NormalCase3"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_N3 = debug[LE][step]["NormalCase3"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> N3 = {ntrk1_ntru1_D_N3, ntrk1_ntru2_D_N3};

    // -> Normalization Shifting Fits
    TH1F* ntrk1_ntru1_D_NS3 = debug[LE][step]["ShiftNormalCase3"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_NS3 = debug[LE][step]["ShiftNormalCase3"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> NS3 = {ntrk1_ntru1_D_NS3, ntrk1_ntru2_D_NS3};

    // -> Normalization Shifting FFT Fits
    TH1F* ntrk1_ntru1_D_NSFFT3 = debug[LE][step]["ShiftNormalFFTCase3"]["ntrk1-ntru1-D"]; 
    TH1F* ntrk1_ntru2_D_NSFFT3 = debug[LE][step]["ShiftNormalFFTCase3"]["ntrk1-ntru2-D"]; 
    std::vector<TH1F*> NSFFT3 = {ntrk1_ntru1_D_NSFFT3, ntrk1_ntru2_D_NSFFT3}; 

    (*Stepper1)[LE + "-Normal"].push_back((ntrk1_ntru2_D_N1 -> Integral())/ (ntrk1_ntru2_T -> Integral())); 
    (*Stepper1)[LE + "-ShiftNormalFFT"].push_back((ntrk1_ntru2_D_NSFFT1 -> Integral())/ (ntrk1_ntru2_T -> Integral())); 

    (*Stepper2)[LE + "-Normal"].push_back((ntrk1_ntru2_D_N2 -> Integral())/ (ntrk1_ntru2_T -> Integral())); 
    (*Stepper2)[LE + "-ShiftNormalFFT"].push_back((ntrk1_ntru2_D_NSFFT2 -> Integral())/ (ntrk1_ntru2_T -> Integral())); 

    (*Stepper3)[LE + "-Normal"].push_back((ntrk1_ntru2_D_N3 -> Integral())/ (ntrk1_ntru2_T -> Integral())); 
    (*Stepper3)[LE + "-ShiftNormalFFT"].push_back((ntrk1_ntru2_D_NSFFT3 -> Integral())/ (ntrk1_ntru2_T -> Integral())); 

    if (!Plotter){continue;}
    Printer(can, Layer + "/" + Energy + "/Normal/Case1/" + step + "_per.png", N1, Truth, "Case-1: Normalization at " + step + "%"); 
    Printer(can, Layer + "/" + Energy + "/ShiftNormal/Case1/" + step + "_per.png", NS1, Truth, "Case-1: Normalization and Shifting at " + step + "%");
    Printer(can, Layer + "/" + Energy + "/ShiftNormalFFT/Case1/" + step + "_per.png", NSFFT1, Truth, "Case-1: Normalization and Shifting (FFT) at " + step + "%");

    Printer(can, Layer + "/" + Energy + "/Normal/Case2/" + step + "_per.png", N2, Truth, "Case-2: Normalization at " + step + "%"); 
    Printer(can, Layer + "/" + Energy + "/ShiftNormal/Case2/" + step + "_per.png", NS2, Truth, "Case-2: Normalization and Shifting at " + step + "%");
    Printer(can, Layer + "/" + Energy + "/ShiftNormalFFT/Case2/" + step + "_per.png", NSFFT2, Truth, "Case-2: Normalization and Shifting (FFT) at " + step + "%");

    Printer(can, Layer + "/" + Energy + "/Normal/Case3/" + step + "_per.png", N3, Truth, "Case-3: Normalization at " + step + "%"); 
    Printer(can, Layer + "/" + Energy + "/ShiftNormal/Case3/" + step + "_per.png", NS3, Truth, "Case-3: Normalization and Shifting at " + step + "%");
    Printer(can, Layer + "/" + Energy + "/ShiftNormalFFT/Case3/" + step + "_per.png", NSFFT3, Truth, "Case-3: Normalization and Shifting (FFT) at " + step + "%");
  }
}



int main(int argc, char* argv[])
{
  TString Directory = argv[1];
  std::vector<_TS> Layer_Energy = Layer_Energy_H();
  maps debug = ReadDebugging(Directory);
  
  float steps_size = 20; 
  float start = 0; 
  float end = 200;
  bool Plot = true;
  int atstep1 = 40; 
  int atstep2 = 60; 
  int atstep3 = 80;

  std::map<_TS, std::vector<float>> Stepper1; 
  std::map<_TS, std::vector<float>> Stepper2; 
  std::map<_TS, std::vector<float>> Stepper3;

  for (_TS LH : Layer_Energy)
  {
    CaseAnalysis(debug, LH, steps_size, start, end, Plot, &Stepper1, &Stepper2, &Stepper3); 
  }
  

  for (_TS L : Layer_H)
  {
    std::vector<float> Case1_N_L;
    std::vector<float> Case2_N_L;
    std::vector<float> Case3_N_L;
    
    std::vector<float> Case1_NSFFT_L;
    std::vector<float> Case2_NSFFT_L;
    std::vector<float> Case3_NSFFT_L;

    for (_TS E : Energy_H)
    {
      TCanvas* can = new TCanvas();
      ConformCanvas(can);
      std::vector<float> Case1_N = Stepper1[L + "_" + E + "-Normal"]; 
      std::vector<float> Case2_N = Stepper2[L + "_" + E + "-Normal"]; 
      std::vector<float> Case3_N = Stepper3[L + "_" + E + "-Normal"]; 
      
      std::vector<float> Case1_NSFFT = Stepper1[L + "_" + E + "-ShiftNormalFFT"]; 
      std::vector<float> Case2_NSFFT = Stepper2[L + "_" + E + "-ShiftNormalFFT"]; 
      std::vector<float> Case3_NSFFT = Stepper3[L + "_" + E + "-ShiftNormalFFT"]; 
    
      TGraph* g1_N = MakeGraphSteps(Case1_N, "Normal-Case1", steps_size, start); 
      TGraph* g2_N = MakeGraphSteps(Case2_N, "Normal-Case2", steps_size, start);
      TGraph* g3_N = MakeGraphSteps(Case3_N, "Normal-Case3", steps_size, start);

      TGraph* g1_NSFFT = MakeGraphSteps(Case1_NSFFT, "ShiftNormalWidthFFT-Case1", steps_size, start); 
      TGraph* g2_NSFFT = MakeGraphSteps(Case2_NSFFT, "ShiftNormalWidthFFT-Case2", steps_size, start);
      TGraph* g3_NSFFT = MakeGraphSteps(Case3_NSFFT, "ShiftNormalWidthFFT-Case3", steps_size, start);
      
      std::vector<TGraph*> GR = {g1_N, g2_N, g3_N, g1_NSFFT, g2_NSFFT, g3_NSFFT}; 
      _TS Title = "Integral Ratio Between Fit Prediction and Truth at " +  L + " " + E.ReplaceAll("_GeV", " GeV").ReplaceAll("_", "-"); 
      _TS YTitle = "Prediction/Truth";
      _TS XTitle = "Scaling of Original 2-Truth Distribution (%)";

      CombineStepGraphs(GR, Title, can, "Ratio_2Truth/" + L + "_" + E + ".png", XTitle, YTitle, 0.1, 10);

      
      BulkDelete(GR);
      delete can;

      Case1_N_L.push_back(Case1_N[int(atstep1/steps_size)-1]); 
      Case2_N_L.push_back(Case2_N[int(atstep2/steps_size)-1]);
      Case3_N_L.push_back(Case3_N[int(atstep3/steps_size)-1]);

      Case1_NSFFT_L.push_back(Case1_NSFFT[int(atstep1/steps_size)-1]); 
      Case2_NSFFT_L.push_back(Case2_NSFFT[int(atstep2/steps_size)-1]);
      Case3_NSFFT_L.push_back(Case3_NSFFT[int(atstep3/steps_size)-1]); 
    }
 
    TCanvas* can = new TCanvas();
    ConformCanvas(can);

    _TS atstep1_s = ""; atstep1_s += (atstep1);  
    _TS atstep2_s = ""; atstep2_s += (atstep2);  
    _TS atstep3_s = ""; atstep3_s += (atstep3);  
    
    TGraph* g1_N_L = MakeGraphJets(Case1_N_L, "Normal-Case1-" + atstep1_s + "%"); 
    TGraph* g2_N_L = MakeGraphJets(Case2_N_L, "Normal-Case2-" + atstep2_s + "%");
    TGraph* g3_N_L = MakeGraphJets(Case3_N_L, "Normal-Case3-" + atstep3_s + "%");

    TGraph* g1_NSFFT_L = MakeGraphJets(Case1_NSFFT_L, "ShiftNormalWidthFFT-Case1-" + atstep1_s + "%"); 
    TGraph* g2_NSFFT_L = MakeGraphJets(Case2_NSFFT_L, "ShiftNormalWidthFFT-Case2-" + atstep2_s + "%");
    TGraph* g3_NSFFT_L = MakeGraphJets(Case3_NSFFT_L, "ShiftNormalWidthFFT-Case3-" + atstep3_s + "%");
    
    std::vector<TGraph*> GRL = {g1_N_L, g2_N_L, g3_N_L, g1_NSFFT_L, g2_NSFFT_L, g3_NSFFT_L}; 
    _TS Title = "Integral Ratio Between Fit Prediction and Truth in the " +  L; 
    _TS YTitle = "Prediction/Truth";
    _TS XTitle = "Jet Energy (GeV)";

    CombineStepGraphs(GRL, Title, can, "Ratio_2Truth_Layer/" + L + ".png", XTitle, YTitle, 0.1, 10);
    
    BulkDelete(GRL);
    delete can;
  }






  return 0; 
} 
