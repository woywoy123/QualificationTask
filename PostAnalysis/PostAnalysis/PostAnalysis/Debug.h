#ifndef DEBUG_H
#define DEBUG_H

#include<PostAnalysis/RooFitBaseFunctions.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/IO.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<TGraphSmooth.h>
#include<TMultiGraph.h>

void Proxy(); 
void SmoothingTest();

const static std::vector<float> k1 = {0.1, 8.0}; 
const static std::vector<float> k2 = {0.1, 8.0}; 
const static std::vector<float> k3 = {0.1, 8.0}; 
const static std::vector<float> k4 = {0.1, 8.0}; 
const static std::vector<std::vector<float>> Ranges = {k1, k2, k3, k4}; 
static std::vector<TH1F*> Graphs_1Tru; 
static std::vector<TH1F*> Graphs_2Tru; 

const float m = 0.5; 
// Normalization parameters
const MVF Params_N = { 
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"Minimizer", {100000}}, 
  {"seek", {1}},
  //{"GSL", {1}},
  //{"Print", {3}},
}; 

// Normalization Shift parameters
const MVF Params_NS = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"dx_s", {-m, -m, -m, -m}},
  {"dx_e", {m, m, m, m}},
  {"dx_G", {0, 0, 0, 0}}, 
  {"dx_C", {0, 0, 0, 0}},
  //{"fft_cache", {10000}},  
  {"Minimizer", {10000}}, 
  {"seek", {1}},
  //{"GSL", {1}},
  //{"Strategy", {2}},
  //{"Print", {3}},
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
  {"s_e", {0.5, 0.5, 0.5, 0.5}}, 
  {"fft_cache", {10000}},  
  {"Minimizer", {10000}},
  //{"seek", {1}},
  //{"GSL", {1}},

  //{"Print", {1}},
};

const MVF Params_WidthDeFFT = {
  //{"Range_ntrk_1", Ranges[0]},  
  //{"Range_ntrk_2", Ranges[1]},    
  {"m_s", {-m, -m, -m, -m}}, 
  {"m_e", {m, m, m, m}}, 
  {"m_G", {0, 0, 0, 0}}, 
  {"s_s", {0.001, 0.001, 0.001, 0.001}}, 
  {"s_e", {0.5, 0.5, 0.5, 0.5}}, 
  {"fft_cache", {10000}},  
  {"Minimizer", {10000}},
  //{"seek", {1}},
  //{"GSL", {1}},
  {"G_Mean", {0., 0., 0., 0.}}, 
  {"G_Stdev", {0.5, 0.5, 0.5, 0.5}}, 
  {"LR", {100}}, 

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

static std::pair<std::vector<MVF>, std::vector<TString>> StepFit(TH1F* start_P, std::vector<TH1F*> truth_P, float delta, TCanvas* can, TString alg, TString sample)
{
 
  int min = truth_P[0] -> GetXaxis() -> GetXmin(); 
  int max = truth_P[0] -> GetXaxis() -> FindBin(5.)+1; 
  float max_n = truth_P[0] -> GetBinCenter(max); 
  float min_n = truth_P[0] -> GetXaxis() -> GetXmin();

  std::vector<TH1F*> truth = MakeTH1F({"ntrk1-tru1", "ntrk1-tru2", "ntrk1-tru3", "ntrk1-tru4"}, max-min, min_n, max_n);
  TH1F* start = MakeTH1F({"start"}, max-min, min_n, max_n)[0]; 

  std::vector<TH1F*> T = {truth[0], truth[1], truth[2], truth[3], start};  
  std::vector<TH1F*> Fillers = {truth_P[0], truth_P[1], truth_P[2], truth_P[3], start_P}; 
  
  for (int i(0); i < Fillers.size(); i++)
  {
    int k = 0; 
    for (int j(min); j < max; j++)
    {
      T[i] -> SetBinContent(k+1, Fillers[i] -> GetBinContent(j+1)); 
      k++; 
    }
  }
  
  int steps = 200 / delta; 
  float cu = 0; 
  TH1F* G = (TH1F*)start -> Clone("x");  
  std::vector<MVF> out; 
  std::vector<TString> setting; 
  float S_L = truth[1] -> Integral(); 
  Normalize(truth[1]); 
  VT temp; 

  // Create the output graph
  TH1F* P1 = new TH1F(alg + "tru1", alg, steps, 0, 200); 
  Graphs_1Tru.push_back(P1); 

  TH1F* P2 = new TH1F(alg + "tru2", alg, steps, 0, 200); 
  Graphs_2Tru.push_back(P2);

  // Use the 1-Track, 1-Truth as starting point
  if ( alg.Contains("Case1") )
  {
    temp = BuildNtrkMtru(1, truth[0], "_" + alg, 2)[0];
    temp[0] -> Reset(); 
    temp[0] -> Add(truth[0]); 
    P1 -> SetTitle("C1: 1-Trk, 1-Tru Starting Point " + alg); 
    P2 -> SetTitle("C1: 1-Trk, 1-Tru Starting Point " + alg); 
  }
  // Case 2: Use the 2-Truth as the template  
  else if (alg.Contains("Case2"))
  {
    temp = BuildNtrkMtru(1, start, "_" + alg, 2)[0]; 
    temp[1] -> Reset(); 
    temp[1] -> Add(truth[1]); 
    P1 -> SetTitle("C2: 1-Trk 0.1 < R, 2-Tru Replaced by Truth " + alg); 
    P2 -> SetTitle("C2: 1-Trk 0.1 < R, 2-Tru Replaced by Truth " + alg); 
  }
  else if (alg.Contains("Case3"))
  {
    temp = BuildNtrkMtru(1, start, "_" + alg, 2)[0]; 
    P1 -> SetTitle("C3: Convolution " + alg);
    P2 -> SetTitle("C3: Convolution " + alg);
  }
  alg = alg(0, alg.Length() -5); 
  
   
  for (int i(0); i < steps; i++)
  {
    cu += delta; 
    TString title = "Fit Algo: " +  alg + " Truth-2 Normalization being "; 
    title += (cu); title += "% of Original Normalization";
    TString name = "_ntru_2_"; name += (cu); name += "%"; 
    
    VT templ; 
    for (int k(0); k < temp.size(); k++)
    {  
      TH1F* H = temp[k]; 
      TString n = H -> GetTitle(); n += (k+1); 
      templ.push_back((TH1F*)H -> Clone(n)); 
    }

    std::cout << "=========== >" + title << std::endl;

    truth[1] -> Scale(S_L*(cu / 100)); 
    TH1F* ntrk_1 = SumHists({truth[0], truth[1]}, title); 
    ntrk_1 -> SetTitle(title);
    ntrk_1 -> SetLineColor(kBlack);
    MVF N_Res; 
    if (alg == "Normal"){ N_Res = Normalization(ntrk_1, templ, Params_N, "Normal"); }
    if (alg == "ShiftNormal"){ N_Res = NormalizationShift(ntrk_1, templ, Params_NS, "NormalShift"); }
    if (alg == "ShiftNormalFFT"){ N_Res = ConvolutionFFT(ntrk_1, templ, Params_FFT); }
    if (alg == "ShiftNormalWidthFFT"){ N_Res = ConvolutionFFT(ntrk_1, templ, Params_WidthFFT); }
    if (alg == "ShiftNormalWidthDeFFT"){ N_Res = DeConvolutionFFT(ntrk_1, templ, Params_WidthDeFFT); } 

    PlotHists(ntrk_1, {truth[0], truth[1]}, templ, can); 
    can -> Print("Package.pdf"); 
    can -> Clear(); 
    
    N_Res["Truth_Int"].push_back(truth[0] -> Integral()); 
    N_Res["Truth_Int"].push_back(truth[1] -> Integral()); 
    N_Res["Truth_Int"].push_back(truth[2] -> Integral()); 
    N_Res["Truth_Int"].push_back(truth[3] -> Integral()); 
  
    out.push_back(N_Res); 
    setting.push_back(title);

    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    Normalize(truth[1]); 

    // Collect Data Point for the performance plots 
    float r1 = (N_Res["Normalization"][0] / N_Res["Truth_Int"][0])*100; 
    float r2 = (N_Res["Normalization"][1] / N_Res["Truth_Int"][1])*100; 
    P1 -> SetBinContent(i+1, r1); 
    P2 -> SetBinContent(i+1, r2); 
    
    delete ntrk_1; 
    BulkDelete(templ); 
  }
  
  truth[1] -> Scale(S_L); 
  BulkDelete(truth); 
  delete start; 
  BulkDelete(temp); 
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

    std::vector<float> T_V = A["Truth_Int"]; 

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

      TString t = " Norm: " + PrecisionString(n, 4, true) + " +- " + 
                  PrecisionString(n_e, 4, true) + " | " + 
                  " Algo: " + PrecisionString(a, 4, true) + " +- " + 
                  PrecisionString(a_e, 4, true) + " | " + 
                  " A-N: " + PrecisionString(div, 4, true) + " +- " +
                  PrecisionString(div_er*100, 4, true) + " | " +
                  " A/N: " + PrecisionString(r, 4, false) + " |"; 
      
      TString g = "=> tru-"; g += (i+1); g += (" | "); g += t; 
      c_out -> push_back(g); 


      n = T_V[i]; 
      a = A_V[i]; 
      r = a/n;

      t = " Truth: " + PrecisionString(n, 4, true) + " | " + 
          " Algo: " + PrecisionString(a, 4, true) + " | " + 
          " A/T: " + PrecisionString(r, 4, false) + " |" +
          " N/T: " + PrecisionString(N_V[i] / n, 4, false) + " |"; 

      g = "=> tru-(Truth)"; g += (i+1); g += (" | "); g += t; 
      c_out -> push_back(g); 
    }
    c_out -> push_back(""); 

  }

}









#endif
