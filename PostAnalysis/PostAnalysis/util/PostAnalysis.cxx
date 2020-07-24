// Including commonly reused functions
#include<PostAnalysis/Functions.h>
#include<PostAnalysis/Verification.h>

// Including standard C++ libraries 
#include<iostream>

using namespace RooFit;

void PostAnalysis()
{
  // Init the function that hosts common tools
  Functions F;
  Fit_Functions f;
  Verification V;
  Benchmark B; 

  // Probe for checking hist steps in code  
  TH1F* Probe = new TH1F("Probe", "probe", 500, 0, 20);
  TH1F* trk2 = new TH1F("trk2", "trk2", 500, 0, 20); 
  TCanvas* can = new TCanvas();
  can -> SetLogy();

  // Create the histograms for toy. 
  std::vector<TString> nTrk = {"1trk", "2trk", "3trk", "4trk"};
  std::vector<TString> nTrkS = {"1trkS", "2trkS", "3trkS", "4trkS"};
  std::vector<TString> nTrkG = {"1trkG", "2trkG", "3trkG", "4trkG"};
  std::vector<TH1F*> nTrk_H = F.MakeTH1F(nTrk, 500, 0, 20);
  std::vector<TH1F*> nTrk_S = F.MakeTH1F(nTrkS, 500, 0, 20);
  std::vector<TH1F*> nTrk_G = F.MakeTH1F(nTrkG, 500, 0, 20);

  // Fill with Landau
  std::vector<float> Params1 = {1, 0.9, 0.1};
  V.Debug(nTrk_H, Params1); 
  std::vector<float> Params2 = {1, 0.91, 0.075};
  V.Debug(nTrk_G, Params2); 

  // ======================================== Experimental =================================== // 
  int nBins = nTrk_H.at(0) -> GetNbinsX(); 
  float Offset = nBins*0.1;

  //nTrk_H.at(1) -> Draw("SAMEHIST");

  std::vector<float> data(nBins + Offset, 0);
  std::vector<float> deconv(nBins + Offset, 1);

  // Here he gets the entries from hist and stores into a vector.
  for (size_t i(0); i < nBins; i++)
  {
    data[i] = nTrk_H.at(1) -> GetBinContent(i+1);
  }

  // Does this weird tail inversion. Not sure why. 
  for (int i(0); i < Offset; i++)
  {
    data[ nBins + i ] = data[ nBins -1 -i ];
  }

  // *** Main Algorithm Part ================== //
  for (int i(0); i < 100; i++)
  {
    deconv = f.LRDeconvolution(data, deconv, deconv, 0.75); 
    F.VectorToTH1F(deconv, Probe); 
    f.ConvolveHists(Probe, Probe, trk2, 0, 20);
    trk2 -> Draw("SAMEHIST*");  
    can -> Update();
  }




 
//  for (std::vector<float> prg : Progress)
//  { 
//    float dist = B.WeightedEuclidean(data, prg);  
//    std::cout << dist << std::endl; 
//  }















}

void StandaloneApplications(int argc, char**argv)
{
  PostAnalysis();
}

int main(int argc, char** argv)
{
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplications(app.Argc(), app.Argv());
  app.Run();

  return 0;
}
