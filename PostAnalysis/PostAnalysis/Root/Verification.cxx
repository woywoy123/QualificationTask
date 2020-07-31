#include<PostAnalysis/Verification.h>
#include<PostAnalysis/Functions.h>
#include<TF1.h>

using namespace RooFit;

std::vector<float> LR( const std::vector<float>& data, const std::vector<float>& current, const std::vector<float>& psf, float dampen)
{
  size_t offset = (data.size() - current.size())/2;
  std::vector<float> next(current.size(), 0);

  for ( size_t i(0); i < current.size(); i++)
  {
    float sum_j = 0;
    size_t upperLimitJ = psf.size();
    for ( size_t j(i); j < upperLimitJ; j++)
    {
      float c_j = 0;
      size_t upperLimitK = j;
      for ( size_t k(0); k <= upperLimitK; k++ )
      {
        c_j += psf[j-k]*current[k];
      }
      if (c_j != 0)
      {
        sum_j += data[j] / c_j * psf[j-i];
      }
    }
    next[i] = current[i] * (1+dampen*(sum_j - 1));
    if (next[i] < 0. || std::isnan(next[i]) || std::isinf(next[i]))
    {
      next[i] = 0;
    }
  }
  return next;
};

void Verification::Debug(std::vector<TH1F*> Hist, std::vector<float> Params)
{
  TF1 Lan("Lan", "landau", 0, 20);
  for (int i(0); i < Params.size(); i++)
  {
    Lan.SetParameter(i, Params.at(i));
  }

  for ( int i(0); i < 500000; i++)
  {
    double r1 = Lan.GetRandom();
    double r2 = Lan.GetRandom();
    double r3 = Lan.GetRandom();
    double r4 = Lan.GetRandom();
    Hist.at(0) -> Fill(r1);  
    Hist.at(1) -> Fill(r1 + r2); 
    Hist.at(2) -> Fill(r1 + r2 + r3); 
    Hist.at(3) -> Fill(r1 + r2 + r3 + r4); 
  } 
}

// Test case the Fitting
std::vector<float> TestFit(std::vector<TH1F*> PDF, TH1F* Data)
{
  Fit_Functions f;
  float Lumi = Data -> Integral();

  std::vector<RooRealVar*> Var = f.FitPDFtoData(PDF, Data, 0, 20); 
 
  std::vector<float> Results; 
  for ( int i(0); i < Var.size(); i++)
  {
    float e = Var[i] -> getVal(); 
    float frac = e/Lumi;
    Results.push_back(frac); 
    delete Var[i];
  }

  return Results;
 
}



std::vector<float> TailReplace(TH1F* hist, std::vector<float> deconv)
{

  // Get basic information 
  float bins = hist -> GetNbinsX(); // Number of bins in hist
  float n_D = deconv.size(); // Number of bins for deconv 
  int delta = n_D - bins; // Difference in length 

  // Make a copy/create new histograms 
  TH1F* Deconv = (TH1F*)hist -> Clone("Deconv");
  Deconv -> Reset();
  Deconv -> SetLineColor(kOrange); // <------ delete after
  TH1F* Hist = (TH1F*)hist -> Clone("Hist");

  // Fill the Deconv
  for (int i(0); i < bins; i++)
  {
    Deconv -> SetBinContent(i+1, deconv[i]);
  }

  // Find the peaks of the two histograms 
  int n_max_D = Deconv -> GetMaximumBin();
  int n_max_H = Hist -> GetMaximumBin();
  float max_D = Deconv -> GetBinContent(n_max_D);
  float max_H = Hist -> GetBinContent(n_max_H);
  int delta_peak = n_max_D - n_max_H;

  for (int i(n_max_D+1); i < bins; i++)
  {
    Deconv -> SetBinContent(i+1, (Hist -> GetBinContent(i+1 - delta_peak))*(max_D/max_H));
  }

  TCanvas* can = new TCanvas("can", "can", 1600, 800);
  can -> SetLogy();
  Hist -> Draw("SAMEHIST");
  Deconv -> Draw("SAMEHIST");

  
  std::vector<float> De;
  for (int i(0); i < n_D; i++)
  {
    if ( i < bins ) { De.push_back(Deconv -> GetBinContent(i+1)); }
    else { De.push_back(deconv[i]*(max_D/max_H)); }
  }
  
  return De;

}


void Verification::UnitTesting()
{
  Functions F;
  Fit_Functions f;

  std::vector<TString> Hist_Names = {"trk1", "trk2", "trk3", "trk4"};
  std::vector<TH1F*> Hists = F.MakeTH1F(Hist_Names, 500, 0, 20);
  Debug(Hists, {1, 0.9, 0.1}); 
 
  // Histograms - With coloring  
  TH1F* trk1 = Hists.at(0);
  TH1F* trk2 = Hists.at(1);
  TH1F* trk3 = Hists.at(2);
  TH1F* trk4 = Hists.at(3);
  trk1 -> SetLineColor(kRed);
  trk2 -> SetLineColor(kBlue);
  trk3 -> SetLineColor(kOrange);
  trk4 -> SetLineColor(kGreen); 

  //// Fit testing section
  //std::vector<TH1F*> PDF = {trk1, trk2, trk3, trk4};
  //TH1F* Data1 = (TH1F*)trk1 -> Clone("Data1");
  //TH1F* Data2 = (TH1F*)trk2 -> Clone("Data2");
  //TH1F* Data3 = (TH1F*)trk3 -> Clone("Data3");
  //TH1F* Data4 = (TH1F*)trk4 -> Clone("Data4");
  //std::vector<float> Fit_test1 = TestFit(PDF, Data1);
  //std::vector<float> Fit_test2 = TestFit(PDF, Data2);
  //std::vector<float> Fit_test3 = TestFit(PDF, Data3);
  //std::vector<float> Fit_test4 = TestFit(PDF, Data4);  

  // Test Tail replace
  // Create the deconv fake vector  
  float offset = 0.1;
  int nbins = trk2 -> GetNbinsX(); 
  std::vector<float> deconv(nbins + nbins*offset, 0);
  
  for (int i(0); i < deconv.size(); i++)
  {
    if (i < nbins){deconv[i] = trk2 -> GetBinContent(i+1);}
    else { deconv[i] = trk2 -> GetBinContent(2*nbins - i -1);}
  }
  
 
  std::vector<float> deconv = TailReplace(trk1, deconv); 




 
}
