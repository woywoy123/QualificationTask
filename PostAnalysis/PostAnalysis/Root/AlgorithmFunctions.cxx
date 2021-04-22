#include<PostAnalysis/AlgorithmFunctions.h>

// basic n-track, n-truth fits 
std::vector<TH1F*> Normalization_Fit(std::vector<TH1F*> Data, TH1F* trk1_Start, std::map<TString, std::vector<float>> Params, TString JE)
{
  std::vector<TString> ntrk_ntru_names; 
  for (int i(0); i < Data.size(); i++)
  {
    TString base = "dEdx_ntrk_"; base += (i+1); base += ("_ntru_"); base += (i+1); 
    ntrk_ntru_names.push_back(base); 
  }
  TString ext = "_" + JE + "_Normal"; 
  std::vector<TH1F*> ntrk_ntru_H = ConvolveNTimes(trk1_Start, ntrk_ntru_names.size(), ntrk_ntru_names, ext); 

  for (int i(0); i < ntrk_ntru_H.size(); i++)
  {
    TString r_name = "Range_ntrk_"; r_name += (i+1); 
    Params["Range"] = Params[r_name]; 

    TH1F* ntrk = ntrk_ntru_H[i];
    TH1F* ntrk_ntru_T = Data[i]; 

    TString base = "Fit_"; base += (i+1); base += (ext); 
    Normalization(ntrk_ntru_T, {ntrk}, Params, base); 
  }
  return ntrk_ntru_H; 
}

// Continue here!!! Build all the algo wrappers






















