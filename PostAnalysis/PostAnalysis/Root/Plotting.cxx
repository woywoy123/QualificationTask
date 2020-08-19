#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/Constants.h>

TCanvas* Plotting::SimplePlot(TH1F* Hist)
{
  TCanvas* can = new TCanvas();
  gStyle -> SetOptStat(0);
  can -> SetLogy();
  Hist -> Draw("SAMEHIST");
  
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  len -> AddEntry(Hist, Hist -> GetTitle()); 
  len -> Draw("SAME");
  can -> Update(); 

  return can; 
}

void Plotting::Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len)
{
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = Hists[i];
    H -> SetLineColor(Constants::Colors[i]);
    H -> SetLineStyle(kDashed);
    H -> Draw("SAMEHIST");
    len -> AddEntry(H, H -> GetTitle());
    len -> Draw("SAME");
    can -> Update(); 
  }
}

TCanvas* Plotting::PlotHists(std::vector<TH1F*> Hists)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  gStyle -> SetOptStat(0); 
  Populate(Hists, can, len);
  return can;
}

TCanvas* Plotting::PlotHists(std::vector<TH1F*> Hists, TH1F* Data)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  gStyle -> SetOptStat(0); 
  Data -> SetLineColor(kBlack);
  Data -> Draw("SAMEHIST");
  Populate(Hists, can, len);
  len -> AddEntry(Data, Data -> GetTitle()); 
  return can; 
}

TCanvas* Plotting::PlotHists(TH1F* Hists, TH1F* Data)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  gStyle -> SetOptStat(0); 
  Data -> SetLineColor(kBlack);
  Data -> Draw("SAMEHIST");
  Populate({Hists}, can, len);
  len -> AddEntry(Data, Data -> GetTitle()); 
  return can; 
}

TCanvas* Plotting::PlotHists(std::vector<std::vector<TH1F*>> Hists, std::vector<TH1F*> Data)
{
  TCanvas* can = new TCanvas();  
  can -> SetLogy(); 
  gStyle -> SetOptStat(0);
  can -> Divide(Data.size()); 
  for (int i(0); i < Hists.size(); i++)
  {
    TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
    can -> cd(i+1) -> SetLogy(); 
    can -> cd(i+1);
    Data[i] -> SetLineColor(kBlack);
    Data[i] -> Draw("SAMEHIST"); 
    Populate(Hists[i], can, len);
    len -> AddEntry(Data[i], Data[i] -> GetTitle()); 
  }
  return can; 
}

TCanvas* Plotting::PlotHists(std::vector<std::vector<TH1F*>> Hists)
{
  TCanvas* can = new TCanvas();  
  can -> SetLogy(); 
  gStyle -> SetOptStat(0);
  can -> Divide(Hists.size()); 
  for (int i(0); i < Hists.size(); i++)
  {
    TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
    can -> cd(i+1) -> SetLogy(); 
    can -> cd(i+1);
    Populate(Hists[i], can, len);
  }
  return can; 
}

TCanvas* Plotting::PlotHists(RooAddPdf model, RooRealVar* Domain, std::vector<RooHistPdf*> PDFs, RooDataHist* Data)
{
  RooPlot* xframe = Domain -> frame(RooFit::Title("Figure"));
  Data -> plotOn(xframe, RooFit::Name("Data"));
  for (int i(0); i < PDFs.size(); i++)
  {
    TString name = "trk-"; name+=(i+1);
    model.plotOn(xframe, RooFit::Name(name), RooFit::Components(*PDFs[i]), RooFit::LineStyle(kDotted), RooFit::LineColor(Constants::Colors[i]));  
  }
  TCanvas* can = new TCanvas();
  gPad -> SetLogy();
  xframe -> SetMinimum(1); 
  xframe -> Draw();
  can -> Update();
  return can;
}

// ============================================= Generators ======================================= //
// === Public 
void DistributionGenerators::Landau(std::vector<TH1F*> Hists, std::vector<float> COMP, std::vector<float> Parameters, int Number, float min, float max)
{
  // Define the generator and initialize the parameters 
  TF1 Lan("Lan", "landau", min, max);
  for (int i(0); i < Parameters.size(); i++){Lan.SetParameter(i, Parameters[i]);}
 
  // Fill the hist vector
  for (int i(0); i < Number; i++)
  {
    for (int y(0); y < Hists.size(); y++)
    {
      float r(0);
      for (int x(0); x < y+1; x++){r = r + Lan.GetRandom();}
      if (COMP[y] > 0){ Hists[y] -> Fill(r, COMP[y]); } 
    }   
  }
}

void DistributionGenerators::Gaussian(float mean, float stdev, int Number, TH1F* Hist)
{
  gRandom = new TRandom(); 
  for (int i(0); i < Number; i++){ Hist -> Fill(gRandom -> Gaus(mean, stdev)); }
  delete gRandom; 
}

std::vector<TH1F*> DistributionGenerators::FillTH1F(std::vector<TString> Names, TString dir)
{
  TFile* File = new TFile(dir); 
  std::vector<TH1F*> Hists; 
  for (TString name : Names)
  {
    TH1F* Hist; 
    for (int i(0); i < Constants::Detector.size(); i++)
    {
      if (i == 0){ Hist = TH1FFromFile(name, Constants::Detector[i], File); }
      else { Hist -> Add(TH1FFromFile(name, Constants::Detector[i], File)); }
    }
    Hists.push_back(Hist); 
  }

  return Hists; 
}


// === Private 
TH1F* DistributionGenerators::TH1FFromFile(TString Name, TString Layer, TFile* file)
{
  file -> cd(Layer);
  TH1F* Hist = (TH1F*)gDirectory -> Get(Name);
  file -> cd();
  return Hist; 
}




