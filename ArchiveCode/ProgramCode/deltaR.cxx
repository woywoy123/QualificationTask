#include "BaseFunctions.h"
#include "IO.h"
#include "Plotting.h"

int main(int argc, char* argv[])
{
  TString dir = argv[1]; 
  
  MMTH3F F = ReadDeltaR(dir);

  
  TCanvas* can = new TCanvas();
  can -> SetLogy();
  gStyle -> SetOptStat(0);
  gStyle -> SetImageScaling(3);


  TH1D* dEdx_dR_tru1;
  TH1D* dEdx_dR_tru2;
  TH1D* dEdx_dR_me; 

  TH1D* dEdx_je_tru1;
  TH1D* dEdx_je_tru2;
  TH1D* dEdx_je_me; 

  std::map<TString, std::map<TString, float>> trk1_tru1_peak_con_je;
  std::map<TString, std::map<TString, float>> trk1_tru2_peak_con_je;
  std::map<TString, std::map<TString, float>> trk1_gr2_peak_con_je;

  std::map<TString, std::map<TString, float>> trk1_tru1_peak_con_dR;
  std::map<TString, std::map<TString, float>> trk1_tru2_peak_con_dR;
  std::map<TString, std::map<TString, float>> trk1_gr2_peak_con_dR;
  
  float dx_min = 1.8; 
  float dx_max = 3.0; 

  can -> Print("Hists/Plots.pdf["); 

  for (MMTH3Fi i = F.begin(); i != F.end(); i++)
  {
    TString layer = i -> first; 
    
    TH3F* ntrk1_ntru1 = F[layer]["dEdx_dR_JetEnergy_ntrk1_ntru1"]; 
    TH3F* ntrk1_ntru2 = F[layer]["dEdx_dR_JetEnergy_ntrk1_ntru2"]; 
    TH3F* ntrk1_ntru2_gr = F[layer]["dEdx_dR_JetEnergy_ntrk1_gr_ntru2"]; 
    TH3F* ntrk1 = F[layer]["dEdx_dR_JetEnergy_ntrk1"]; 
    
    int delR = ntrk1_ntru1 -> GetYaxis() -> GetNbins(); 
    int JE = ntrk1_ntru1 -> GetXaxis() -> GetNbins();
    int dEdx = ntrk1_ntru1 -> GetZaxis() -> GetNbins(); 
    
    float dr_step = 0.01; 
    float dr_start = 0.00; 
      
    for (float x = dr_start; x < 0.4; x+= dr_step)
    {
      int dr_i = 0; 
      for (int dr(0); dr < delR; dr++)
      {
        float tick_dr = ntrk1_ntru1 -> GetYaxis() -> GetBinCenter(dr); 
        if (x >= tick_dr){dr_i = dr;}
      }
      dr_i += 1;
      
      TH1D* dEdx_dR_tru1 = ntrk1_ntru1 -> ProjectionZ("1-tru", 0, JE, dr_i, delR); 
      TH1D* dEdx_dR_tru2 = ntrk1_ntru2 -> ProjectionZ("2-tru", 0, JE, dr_i, delR);
      TH1D* dEdx_dR_gr_tru2 = ntrk1_ntru2_gr -> ProjectionZ("> 2-tru", 0, JE, dr_i, delR); 
      TH1D* dEdx_dR_me = ntrk1 -> ProjectionZ("1-trk", 0, JE, dr_i, delR); 

      int bin_min = dEdx_dR_tru1 -> GetXaxis() -> FindBin(dx_min); 
      int bin_max = dEdx_dR_tru1 -> GetXaxis() -> FindBin(dx_max); 
      
      float tru1 = dEdx_dR_tru1 -> Integral(bin_min, bin_max); 
      float tru2 = dEdx_dR_tru2 -> Integral(bin_min, bin_max); 
      float gr2 = dEdx_dR_gr_tru2 -> Integral(bin_min, bin_max);
      float me = dEdx_dR_me -> Integral(bin_min, bin_max); 

      float frac1 = tru1/me; 
      float frac2 = tru2/me;
      float frac_gr_2 = gr2/me;

      TString p = ""; p += (x);
      trk1_tru1_peak_con_dR[layer][p] = frac1; 
      trk1_tru2_peak_con_dR[layer][p] = frac2;
      trk1_gr2_peak_con_dR[layer][p] = frac_gr_2; 

      TString f = PrecisionString(x, 2, false);
      TString Title = layer + "_dR_gr_"; Title += (f);;
      dEdx_dR_tru1 -> SetTitle("1-Truth"); 
      dEdx_dR_tru2 -> SetTitle("2-Truth"); 
      dEdx_dR_gr_tru2 -> SetTitle("> 2-Truth");
      dEdx_dR_me -> SetTitle("Observed"); 
      TString Add = " dR > "; Add += (f);
      GenerateNiceStacks(dEdx_dR_me, {dEdx_dR_tru1, dEdx_dR_tru2, dEdx_dR_gr_tru2}, Add + ": Two or more Truth Contributions in 1-Track Monte Carlo Sample in "+layer , can, "dE/dx [MeV g^{-1} cm^]", "#Clusters", "HIST");
      can  -> Print("./Hists/" + Title + ".png");
      can -> Print("Hists/Plots.pdf");     

      delete dEdx_dR_tru1; 
      delete dEdx_dR_tru2;
      delete dEdx_dR_me; 
      delete dEdx_dR_gr_tru2;
    }

    float je_step = 100; 
    float je_start = 0; 
      
    for (float x = je_start; x <= 2500; x+= je_step)
    {
      int je_i = 0; 
      for (int je(0); je < JE; je++)
      {
        float tick_je = ntrk1_ntru1 -> GetXaxis() -> GetBinCenter(je); 
        if (x >= tick_je){je_i = je;}
        else{break;}
      }
      je_i += 1;

      TH1D* dEdx_je_tru1 = ntrk1_ntru1 -> ProjectionZ("1-tru", je_i, JE, 0, delR); 
      TH1D* dEdx_je_tru2 = ntrk1_ntru2 -> ProjectionZ("2-tru", je_i, JE, 0, delR); 
      TH1D* dEdx_je_gr_tru2 = ntrk1_ntru2_gr -> ProjectionZ("> 2-tru", je_i, JE, 0, delR);
      TH1D* dEdx_je_me = ntrk1 -> ProjectionZ("1-trk", je_i, JE, 0, delR); 

      int bin_min = dEdx_je_tru1 -> GetXaxis() -> FindBin(dx_min); 
      int bin_max = dEdx_je_tru1 -> GetXaxis() -> FindBin(dx_max); 

      float tru1 = dEdx_je_tru1 -> Integral(bin_min, bin_max); 
      float tru2 = dEdx_je_tru2 -> Integral(bin_min, bin_max);
      float gr2 = dEdx_je_gr_tru2 -> Integral(bin_min, bin_max);
      float me = dEdx_je_me -> Integral(bin_min, bin_max); 

      float frac1 = tru1/me; 
      float frac2 = tru2/me;
      float frac_gr_2 = gr2/me;
      
      TString p = ""; p += (x);
      trk1_tru1_peak_con_je[layer][p] = frac1; 
      trk1_tru2_peak_con_je[layer][p] = frac2;
      trk1_gr2_peak_con_je[layer][p] = frac_gr_2; 

      TString f = PrecisionString(x, 2, false);
      TString Title = layer + "_je_gr_"; Title += (f);;
      dEdx_je_tru1 -> SetTitle("1-Truth"); 
      dEdx_je_tru2 -> SetTitle("2-Truth"); 
      dEdx_je_gr_tru2 -> SetTitle("> 2-Truth");
      dEdx_je_me -> SetTitle("Observed"); 
      TString Add = " Jet Energy > "; Add += (f); Add += ("GeV");
      GenerateNiceStacks(dEdx_je_me, {dEdx_je_tru1, dEdx_je_tru2, dEdx_je_gr_tru2}, Add + ": Two or more Truth Contributions in 1-Track Monte Carlo Sample in "+layer , can, "dE/dx [MeV g^{-1} cm^]", "#Clusters", "HIST");
      can  -> Print("./Hists/" + Title + ".png");
      can -> Print("Hists/Plots.pdf");     

      delete dEdx_je_tru1; 
      delete dEdx_je_tru2;
      delete dEdx_je_me; 
      delete dEdx_je_gr_tru2;
    }

    TGraph* gr_je_tru1 = MapGraph(trk1_tru1_peak_con_je[layer], "1-Truth", "Jet Energy (GeV)", "fraction");
    TGraph* gr_je_tru2 = MapGraph(trk1_tru2_peak_con_je[layer], "2-Truth", "Jet Energy (GeV)", "fraction");
    TGraph* gr_je_gr2 = MapGraph(trk1_gr2_peak_con_je[layer], "> 2-Truth", "Jet Energy (GeV)", "fraction");
    
    TString add = ""; add += (PrecisionString(dx_min, 2, false)); add += ("->"); add += (dx_max); 
    CombinedGraphs({gr_je_tru1, gr_je_tru2, gr_je_gr2}, layer + ": Fraction of Truth Contribution in 1-Track Monte Carlo dE/dx range: "+add, 0.001, 2, 100, 2500, "fraction", "Jet-Energy (GeV)", can, true); 
    can -> Print("./Hists/" + layer + "_jet_energy.png"); 
    can -> Print("Hists/Plots.pdf");     
    can -> Clear(); 

    TGraph* gr_dR_tru1 = MapGraph(trk1_tru1_peak_con_dR[layer], "1-Truth", "dR", "fraction");
    TGraph* gr_dR_tru2 = MapGraph(trk1_tru2_peak_con_dR[layer], "2-Truth", "dR", "fraction");
    TGraph* gr_dR_gr2 = MapGraph(trk1_gr2_peak_con_dR[layer], "> 2-Truth", "dR", "fraction");

    CombinedGraphs({gr_dR_tru1, gr_dR_tru2, gr_dR_gr2}, layer + ": Fraction of Truth Contribution in 1-Track Monte Carlo dE/dx range: "+add, 0.0001, 10, 0., 0.39, "fraction", "dR", can, true); 
    can -> Print("./Hists/" + layer + "_dR.png"); 
    can -> Print("Hists/Plots.pdf");     
    can -> Clear(); 
  }
  can -> Print("Hists/Plots.pdf]");
}
