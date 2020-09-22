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

void Plotting::Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style)
{
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = Hists[i];
    H -> SetLineColor(Constants::Colors[i]);
    H -> SetLineStyle(Style);
    H -> SetLineWidth(1);
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

void Plotting::PlotHists(std::vector<TH1F*> Hists, TCanvas* can)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  can -> SetLogy(); 
  gStyle -> SetOptStat(0); 
  Populate(Hists, can, len);
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

void Plotting::PlotHists(std::vector<TH1F*> Hists, TH1F* Data, TCanvas* can)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  can -> SetLogy(); 
  gStyle -> SetOptStat(0); 
  Data -> SetLineColor(kBlack);
  Data -> SetMinimum(1e-8);
  Data -> Draw("SAMEHIST");
  Populate(Hists, can, len);
  len -> AddEntry(Data, Data -> GetTitle()); 
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
    can -> Update();
  }
  return can;
}

void Plotting::PlotHists(std::vector<std::vector<TH1F*>> Hists, std::vector<std::vector<TH1F*>> Closure, std::vector<TH1F*> Data, TCanvas* can)
{
  can -> SetWindowSize(2400, 1200); 
  for (int i(0); i < Hists.size(); i++)
  {
    TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
    can -> cd(i+1) -> SetLogy(); 
    can -> cd(i+1);
    Data[i] -> SetLineColor(kBlack);
    Data[i] -> Draw("SAMEHIST"); 
    Populate(Hists[i], can, len, kSolid);
    Populate(Closure[i], can, len, kDashed); 
    len -> AddEntry(Data[i], Data[i] -> GetTitle()); 
    can -> Update();
  }
}

void Plotting::PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Closure, TH1F* Data, TCanvas* can)
{
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  can -> cd(1) -> SetLogy(); 
  Data -> SetLineColor(kBlack);
  Data -> Draw("SAMEHIST"); 
  Populate(Hists, can, len, kSolid);
  Populate(Closure, can, len, kDotted); 
  len -> AddEntry(Data, Data -> GetTitle()); 
  can -> Update();
}

void Plotting::PlotHists(std::vector<std::vector<TH1F*>> Hists, std::vector<TH1F*> Data, TCanvas* can)
{
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
}

TCanvas* Plotting::PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Data)
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
    Populate({Hists[i]}, can, len);
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

void Plotting::PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Subtract, std::vector<TH1F*> Closure, TCanvas* can)
{
  auto Fill =[](TH1F* H, Color_t col, TCanvas* can, TLegend* len)
  {
    H -> SetLineColor(col); 
    H -> SetLineWidth(1); 
    H -> SetAxisRange(1, H -> Integral(), "Y");
    H -> Draw("SAMEHIST");
  
    len -> AddEntry(H, H -> GetTitle());
    len -> Draw("SAME");
    can -> Update(); 
  }; 

  gStyle -> SetOptStat(0);
  for (int i(0); i < Hists.size(); i++)
  {
    TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
    can -> cd(i+1) -> SetLogy(); 
    can -> cd(i+1);
    Fill(Hists[i], kRed, can, len);
    Fill(Subtract[i], kOrange, can, len);
    Fill(Closure[i], kViolet, can, len);
  }
}

void Plotting::PlotHists(std::vector<TH1F*> Hists, std::vector<TH1F*> Subtract, TCanvas* can)
{
  auto Fill =[](TH1F* H, Color_t col, TCanvas* can, TLegend* len)
  {
    H -> SetLineColor(col); 
    H -> SetLineWidth(1); 
    H -> SetAxisRange(1e-9, 1, "Y");
    H -> Draw("SAMEHIST");
  
    len -> AddEntry(H, H -> GetTitle());
    len -> Draw("SAME");
    can -> Update(); 
  }; 

  gStyle -> SetOptStat(0);
  for (int i(0); i < Hists.size(); i++)
  {
    TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
    can -> cd(i+1) -> SetLogy(); 
    can -> cd(i+1);
    Fill(Hists[i], kRed, can, len);
    Fill(Subtract[i], kBlack, can, len);
  }
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

TCanvas* Plotting::PlotHists(RooAddPdf model, RooRealVar* Domain, std::vector<RooFFTConvPdf*> PDFs, RooDataHist* Data)
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

TCanvas* Plotting::PlotHists(RooHistPdf model, RooRealVar Domain, RooDataHist Data)
{
  RooPlot* xframe = Domain.frame(RooFit::Title("Figure"));
  Data.plotOn(xframe, RooFit::Name("Data"));
  model.plotOn(xframe);  
  TCanvas* can = new TCanvas();
  gPad -> SetLogy();
  xframe -> SetMinimum(1); 
  xframe -> Draw();
  can -> Update();
  return can;
}

void Plotting::DifferencePlot(TH1F* H1, TH1F* H2, TCanvas* can)
{

  TH1F* Ratio = (TH1F*)H1 -> Clone("Difference");
  Ratio -> Clear();  
  for (int i(0); i < H1 -> GetNbinsX(); i++)
  {
    float e1 = H1 -> GetBinContent(i+1); 
    float e2 = H2 -> GetBinContent(i+1); 
    float diff = e1 - e2; 
    Ratio -> SetBinContent(i+1, diff); 
  }

  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1.0);
  P1 -> Draw(); 
  P1 -> cd(); 
  P1 -> SetLogy(); 
  H1 -> Draw("SAMEHIST"); 
  H1 -> SetStats(0);  
  H2 -> SetLineColor(kRed);
  H2 -> Draw("SAMEHIST"); 
  
  can -> cd();

  TPad *P2 = new TPad("P2", "P2", 0.0, 0.05, 1, 0.3); 
  P2 -> Draw();
  P2 -> cd(); 
  Ratio -> SetStats(0); 
  Ratio -> Draw("SAMEHIST"); 
}

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

TH1F* DistributionGenerators::FillTH1F(TString name, std::vector<TString> SubDir, TString dir, std::vector<TString> Detector)
{
  TFile* File = new TFile(dir); 
  TH1F* Hist; 
  int i=0;
  for (TString Layer : Detector)
  {
    for (TString sub : SubDir)
    {
      TString directory = Layer + sub;
      if (i == 0){ Hist = TH1FFromFile(name, directory, File); }
      else { Hist -> Add(TH1FFromFile(name, directory, File)); } 
      i++;
    }
  }

  return Hist; 
}


// === Private 
TH1F* DistributionGenerators::TH1FFromFile(TString Name, TString Layer, TFile* file)
{
  file -> cd(Layer);
  TH1F* Hist = (TH1F*)gDirectory -> Get(Name);
  file -> cd();
  return Hist; 
}





