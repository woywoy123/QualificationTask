#ifndef DEBUG_H
#define DEBUG_H

#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/AlgorithmFunctions.h>

void Proxy(); 
void SmoothingTest();

const static std::vector<float> k1 = {0.4, 10.0}; 
const static std::vector<float> k2 = {0.4, 10.0}; 
const static std::vector<float> k3 = {0.4, 10.0}; 
const static std::vector<float> k4 = {0.4, 10.0}; 

const static std::vector<std::vector<float>> Ranges = {k1, k2, k3, k4}; 

const float m = 0.4; 
// Normalization parameters
const MVF Params_N = { 
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"Minimizer", {10000}}, 
  //{"Print", {1}},
}; 

// Normalization Shift parameters
const MVF Params_NS = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"dx_s", {-m, -m, -m, -m}},
  {"dx_e", {m, m, m, m}},
  //{"dx_G", {0, 0, 0, 0}}, 
  //{"dx_C", {1, 1, 1, 1}},
  //{"fft_cache", {10000}},  
  {"Minimizer", {10000}}, 
  //{"Print", {1}},
};

// Normalization Shift FFT parameters
const MVF Params_FFT = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"m_s", {-m, -m, -m, -m}}, 
  {"m_e", {m, m, m, m}}, 
  {"m_G", {0, 0, 0, 0}}, 
  {"s_C", {1, 1, 1, 1}},
  {"fft_cache", {50000}},  
  {"Minimizer", {10000}},  
  //{"Print", {1}},
}; 

// Normalization Shift Width FFT parameters
const MVF Params_WidthFFT = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"m_s", {-m, -m, -m, -m}}, 
  {"m_e", {m, m, m, m}}, 
  {"m_G", {0, 0, 0, 0}}, 
  {"s_s", {0.001, 0.001, 0.001, 0.001}}, 
  {"s_e", {0.1, 0.1, 0.1, 0.1}}, 
  {"fft_cache", {50000}},  
  {"Minimizer", {10000}},  
  //{"Print", {1}},
};

static void FitParameterError(std::vector<MVF> Params, std::vector<TString>* out, std::vector<TString> keys, std::vector<TString> Setting)
{
  
  for (int j(0); j < Params.size(); j++)
  {
    MVF x = Params[j];
    (*out).push_back(Setting[j]); 
    
    std::map<TString, TString> Temp; 
    for (int l(0); l < keys.size(); l++)
    {
      TString k = keys[l];
      std::vector<float> res = x[k]; 
      std::vector<float> res_er = x[k + "_Error"];
      
      for (int t(0); t < res.size(); t++)
      {
        TString trk = "trk_1_ntru_"; trk += (t+1); 
        TString res_s = PrecisionString(res[t], 4, true); 
        float res_er_perc = std::abs(res_er[t] / res[t]) * 100;
        TString res_er_s = PrecisionString(res_er_perc, 4, false); 
        
        TString info = " | " + k + "-> " + res_s + " +- " + PrecisionString(res_er[t], 4, true) + " (" + res_er_s + "%)"; 
        Temp[trk] += info; 
      }
    }
    
    typedef std::map<TString, TString>::iterator TTi; 
    for (TTi m = Temp.begin(); m != Temp.end(); m++)
    {
      TString name = m -> first + " ::" + m -> second;  
      (*out).push_back(name);
    
    }
    TString status = "-----> Fit Status: "; status += (x["fit_status"][0]); 
    (*out).push_back(status); 
    (*out).push_back("");  
    (*out).push_back("===========> HERE COMES THE ERROR MATRIX"); 
      
    for (MVFi t = x.begin(); t != x.end(); t++)
    {
      if ((t -> first).Contains("Covar_M_i"))
      {
        TString u = " | "; 
        for (float g : t -> second){ u += (PrecisionString(g, 4, true)); u += (" | "); }
        (*out).push_back(u); 
      }
    }
    (*out).push_back("");  
    (*out).push_back("");  
    (*out).push_back("");  
  }

}

static std::pair<std::vector<MVF>, std::vector<TString>> StepFit(TH1F* start, std::vector<TH1F*> truth, float delta, TCanvas* can, TString alg)
{
  int steps = 100 / delta; 
  float cu = 0; 
  
  std::vector<MVF> out; 
  std::vector<TString> setting; 
  float S_L = truth[1] -> Integral(); 
  
  Normalize(truth[1]); 
  for (int i(0); i < steps; i++)
  {
    cu += delta; 
    TString title = "Fit Algo: " +  alg + " Truth-2 Normalization being "; 
    title += (cu); title += "% of Original Normalization";
    TString name = "_ntru_2_"; name += (cu); name += "%"; 
    VT templ = BuildNtrkMtru(1, start, "_" + alg, 2)[0];
    std::cout << "=========== >" + title << std::endl;

    truth[1] -> Scale(S_L*(cu / 100)); 
    TH1F* ntrk_1 = SumHists(truth, title); 
    ntrk_1 -> SetTitle(title);
    ntrk_1 -> SetLineColor(kBlack);
   
    MVF N_Res; 
    if (alg == "Normal"){ N_Res = Normalization(ntrk_1, templ, Params_N); }
    if (alg == "ShiftNormal"){ N_Res = NormalizationShift(ntrk_1, templ, Params_NS); }
    if (alg == "ShiftNormalFFT"){ N_Res = ConvolutionFFT(ntrk_1, templ, Params_FFT); }
    if (alg == "ShiftNormalWidthFFT"){ N_Res = ConvolutionFFT(ntrk_1, templ, Params_WidthFFT); }
   
    //Normalization(ntrk_1, templ, Params_N);

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

    Normalize(truth[1]); 
    
    
  }
  
  truth[1] -> Scale(S_L); 
  return std::pair<std::vector<MVF>, std::vector<TString>>(out, setting); 
}

static void CompareNormalization(std::vector<MVF> Normal, std::vector<MVF> Algo, std::vector<TString> Setting, std::vector<TString>* c_out)
{
  

  for (int x(0); x < Normal.size(); x++)
  {
    MVF N = Normal[x]; 
    MVF A = Algo[x];

    std::vector<float> N_V = N["Normalization"]; 
    std::vector<float> N_E = N["Normalization_Error"]; 
    std::vector<float> A_V = A["Normalization"]; 
    std::vector<float> A_E = A["Normalization_Error"]; 

    c_out -> push_back(Setting[x]); 
    for (int i(0); i < N_V.size(); i++)
    {
      float n = N_V[i]; 
      float n_e = N_E[i]; 

      float a = A_V[i]; 
      float a_e = A_E[i]; 
    
      float div = a - n; 
      float r = a/n; 

      float div_er = div/n; 

      TString t = PrecisionString(n, 4, true) + " +- " + 
                  PrecisionString(n_e, 4, true) + " | " + 
                  PrecisionString(a, 4, true) + " +- " + 
                  PrecisionString(a_e, 4, true) + " | " + 
                  PrecisionString(div, 4, true) + " +- " +
                  PrecisionString(div_er*100, 4, true) + " | " +
                  PrecisionString(r, 4, false) + " |"; 
      
      TString g = "=> tru-"; g += (i+1); g += (" "); g += t; 
      c_out -> push_back(g); 

    }

  }

}









#endif
