#ifndef DEBUG_H
#define DEBUG_H

#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/AlgorithmFunctions.h>

void Proxy(); 
void SmoothingTest();

const static std::vector<float> k1 = {0.1, 10.0}; 
const static std::vector<float> k2 = {0.2, 10.0}; 
const static std::vector<float> k3 = {0.2, 10.0}; 
const static std::vector<float> k4 = {0.2, 10.0}; 

const static std::vector<std::vector<float>> Ranges = {k1, k2, k3, k4}; 

const float m = 0.4; 
// Normalization parameters
const MVF Params_N = { 
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"Minimizer", {10000}}, 
  {"Print", {-1}},
}; 

// Normalization Shift parameters
const MVF Params_NS = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"dx_s", {-m, -m, -m, -m}},
  {"dx_e", {m, m, m, m}},
  //{"dx_G", {0, 0, 0, 0}}, 
  {"fft_cache", {10000}},  
  {"Minimizer", {100000}}, 
  {"Print", {-1}},
}; 

// Normalization Shift FFT parameters
const MVF Params_FFT = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"m_s", {-m, -m, -m, -m}}, 
  {"m_e", {m, m, m, m}}, 
  {"m_G", {0, 0, 0, 0}}, 
  {"s_C", {1, 1, 1, 1}},
  {"fft_cache", {10000}},  
  {"Minimizer", {10000}},  
  {"Print", {-1}},
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
  {"fft_cache", {10000}},  
  {"Minimizer", {10000}},  
  {"Print", {-1}},
};

// Simultaneous Fitting method 
const MVF Params_Sim = {
  {"Range_ntrk_1", Ranges[0]},  
  {"Range_ntrk_2", Ranges[1]},    
  {"m_e", {m, m, m, m}}, 
  {"m_s", {-m, -m, -m, -m}}, 
  {"s_s", {0.001, 0.001, 0.001, 0.001}}, 
  {"s_e", {0.075, 0.05, 0.075, 0.075}}, 
  {"fft_cache", {10000}},  
  {"Minimizer", {50000}},  
  {"Print", {1}}  
}; 

// Experimental Fitting method 
const MVF Params_Exp = { 
  {"Range_ntrk_1", Ranges[0]}, 
  {"Range_ntrk_2", Ranges[1]}, 
  {"m_e", {m, m, m, m}}, 
  {"m_s", {-m, -m, -m, -m}}, 
  {"s_s", {0.001, 0.001, 0.001, 0.001}}, 
  {"s_e", {0.05, 0.05, 0.05, 0.05}}, 
  {"fft_cache", {10000}},  
  {"Minimizer", {10000}},  
  {"Print", {-1}},  
  {"G_Mean", {0, 0, 0, 0}}, 
  {"G_Stdev", {0.01, 0.01, 0.01, 0.01}}, 
  {"LR", {20}}, 
  {"Print", {1}}  
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




#endif
