#include<TH1F.h>
#include<TCanvas.h>


using namespace RooFit;

void FitTail(TH1F* histA, std::vector<float> deconv, float min, float max)
{
  
  int Abins = histA -> GetNbinsX();
  int Bbins = deconv.size(); 

  // MeV to bins;
  float ss = (min - max)/Abins;
  float Pad = (Bbins - Abins)/2; 

  TH1F* Deconv = new TH1F("Deconv", "Deconv", Bbins, min, max + Pad*2);
  TH1F* TRK1 = (TH1F*)Deconv -> Clone("TRK1");

  // Fill the histograms 
  for (int i(0); i < Bbins; i++)
  {
    Deconv -> SetBinContent(i+1, deconv[i]);
    if ( i < Abins)
    {
      TRK1 -> SetBinContent(i+1, histA -> GetBinContent(i+1));
    }
    else
    {
      TRK1 -> SetBinContent(i+1, histA -> GetBinContent(Abins - i - 1));
    }
  }


  // Get the sum of the hists for normalizing them 
  float sum_A(0);
  float sum_B(0); 
  for (int i(0); i < Abins; i++)
  {
    sum_A += TRK1 -> GetBinContent(i+1);
    sum_B += Deconv -> GetBinContent(i+1);
  }

  // Normalize them 
  TRK1 -> Scale(1/sum_A);
  Deconv -> Scale(1/sum_B);
  
  // Fitting range
  float a = 3;
  float b = 15;

  // Fitting being performed here 
  RooRealVar x("x", "x", min - Pad * ss, max + Pad * ss);
  RooRealVar s("s", "s", 0., -6, 6);
  RooFormulaVar delta("delta", "x-s", RooArgSet(x, s));

  RooDataHist D("D", "D", x, Deconv);
  RooHistPdf model("model", "model", delta, x, D, 1);

  RooDataHist trk1("trk1", "trk1", x, TRK1);
  RooFitResult* Result = model.fitTo(trk1, RooFit::Range(a, b));

  RooPlot* xframe = x.frame(RooFit::Title("Shift Test"));
  D.plotOn(xframe, RooFit::LineColor(kAzure));
  trk1.plotOn(xframe, RooFit::LineColor(kRed));
  model.plotOn(xframe, RooFit::LineColor(kBlue));

  TCanvas* can = new TCanvas();
  can -> SetLogy();
  xframe -> SetMinimum(1);
  xframe -> Draw();

  float shift = s.getVal();
  float S = std::round(shift/ss);

  std::cout << "Shift: " << shift << std::endl;
  std::cout << "bins shift: " << S << std::endl;
   

}



void CustomShift(int bins, TH1F* hist)
{
  std::vector<float> temp(hist -> GetNbinsX(), 0);
  for (int i(0); i < hist -> GetNbinsX(); i++)
  {
    temp[i] = hist -> GetBinContent(i+1);
  }

  hist -> Reset();
  for (int i(0); i < temp.size(); i++)
  {
    hist -> SetBinContent(i+1 + bins, temp[i]);
  }
}



void Debug()
{
  float nBins = 500;
  float Padding = 25;
  float min = 0; 
  float max = 20;
  float step = (20 - 0)/nBins;
  
  // Histograms to hold information
  TH1F* histA  = new TH1F("histA", ";x;y;", nBins,0,20); 
  TH1F* histB  = new TH1F("histB", ";x;y;", nBins + 2*Padding, 0 - Padding*step, 20 + Padding*step);

  // True PDF of the energy deposition
  TF1 Landau1("Landau1","landau",0,20);
  Landau1.SetParameter(0,1);
  Landau1.SetParameter(1,0.9);
  Landau1.SetParameter(2,0.1);

  TF1 Landau2("Landau2","landau",0,20);
  Landau2.SetParameter(0,0.6);
  Landau2.SetParameter(1,1.2);
  Landau2.SetParameter(2,0.075);

  // Fill some distributions
  for( int i(0); i < 500000; ++i){
    double r1 = Landau1.GetRandom();
    double r2 = Landau2.GetRandom();
    histA->Fill(r1);
    histB->Fill(r1);
  }

  // HistB is the Deconv vector 
  // HistA is the trk1 distribution - i.e. target
  for (int i(0); i < Padding*2; i++)
  {
     histB -> SetBinContent(Padding + nBins + i + 1, histB -> GetBinContent(Padding + nBins - i -1)); 
  }

  // make histB the shifted histogram 
  CustomShift(25, histB);
  histB -> Scale(2);

  std::vector<float> deconv(nBins + Padding*2, 0);
  for (int i(0); i < histB -> GetNbinsX(); i++)
  {
    deconv[i] = histB -> GetBinContent(i+1);
  }
 
  TCanvas* can = new TCanvas();
  can -> SetLogy();
  histA -> GetYaxis() -> SetRangeUser(1, 1e6);
  histB -> Draw("SAMEHIST");
  histA -> SetLineColor(kRed);
  histA -> Draw("SAMEHSIT");
  can -> Update();

  FitTail(histA, deconv, min, max);

}
