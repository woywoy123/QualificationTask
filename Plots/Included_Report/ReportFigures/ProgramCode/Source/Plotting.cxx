#include "IO.h"
#include "Plotting.h"

float GetMaxValue(TH1F* H)
{
  int bin_m = H -> GetMaximumBin(); 
  float v = H -> GetBinContent(bin_m+1); 
  return v; 
}

void Populate(std::vector<TH1F*> Hists, TCanvas* can, TLegend* len, ELineStyle Style)
{
  std::vector<Color_t> Colors = {kRed, kGreen, kBlue, kCyan, kViolet, kOrange, kCoffee, kAurora}; 
  for (int i(0); i < Hists.size(); i++)
  {
    TH1F* H = Hists[i];
    
    TString name = H -> GetTitle(); 
    if (name.Contains("005_ntru_1")){ H -> SetTitle("Truth-1 (Monte Carlo)"); }
    if (name.Contains("005_ntru_2")){ H -> SetTitle("Truth-2 (Monte Carlo)"); }
    if (name.Contains("005_ntru_3")){ H -> SetTitle("Truth-3 (Monte Carlo)"); }
    if (name.Contains("005_ntru_4")){ H -> SetTitle("Truth-4 (Monte Carlo)"); }

    if (!name.Contains("rless") && name.Contains("ntru_1")){ H -> SetTitle("Truth-1 (Fit)"); }
    if (!name.Contains("rless") && name.Contains("ntru_2")){ H -> SetTitle("Truth-2 (Fit)"); }
    if (!name.Contains("rless") && name.Contains("ntru_3")){ H -> SetTitle("Truth-3 (Fit)"); }
    if (!name.Contains("rless") && name.Contains("ntru_4")){ H -> SetTitle("Truth-4 (Fit)"); }


    H -> SetLineColor(Colors[i]);
    H -> SetFillColorAlpha(Colors[i], 0.2); 
    H -> SetLineStyle(Style);
    H -> SetLineWidth(4);
    H -> Draw("SAMEHIST");
    len -> AddEntry(H, H -> GetTitle());
    len -> Draw("SAME");
    can -> Update(); 
  }
}

void PlotHists(TH1F* Hist, TCanvas* can)
{
  can -> SetLogy(); 
  Hist -> Draw("HIST"); 
  can -> Update(); 
}

void PlotHists(std::vector<TH1F*> Hists, std::vector<TString> Legend_Titles, TCanvas* can)
{
  gStyle -> SetOptStat(0); 
  for (int i(0); i < Hists.size(); i++){Hists[i] -> SetTitle(Legend_Titles[i]); }
  
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  can -> cd(1); 
  Populate(Hists, can, len, kSolid);
}

void PlotHists(std::vector<TH1F*> Hists, TCanvas* can)
{
  float m; 
  for (TH1F* H : Hists)
  {
    int b = H -> GetMaximumBin(); 
    float m_n = H -> GetBinContent(b+1); 
    if (m_n > m){ m = m_n; }
  }
  
  can -> Clear();
  Hists[0] -> Draw("HIST");
  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75);
  Hists[0] -> GetYaxis() -> SetRangeUser(1, m*2); 
  Populate(Hists, can, len, kSolid); 
  can -> Update();
}

void PlotHists(TH1F* Data, std::vector<TH1F*> Hists, TCanvas* can)
{
  int bin = Data -> GetMaximumBin(); 
  float m = Data -> GetBinContent(bin+1); 
  
  can -> Clear();
  gStyle -> SetOptStat(0); 
  Data -> GetYaxis() -> SetRangeUser(1, m);  
  //Data -> GetXaxis() -> SetRangeUser(0, 10);
  Data -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(Hists, can, len, kSolid); 
}

void PlotHists(TH1F* Data, std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TCanvas* can)
{
  int bin = Data -> GetMaximumBin(); 
  float m = Data -> GetBinContent(bin+1); 
  
  can -> Clear();
  gStyle -> SetOptStat(0); 
  Data -> GetYaxis() -> SetRangeUser(1, m*2);
  Data -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  Populate(truth, can, len, kSolid); 
  Populate(prediction, can, len, kDashed); 
}

void PlotHists(std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TCanvas* can)
{
  float sum = 0; 
  for (int i(0); i < prediction.size(); i++)
  {
    int bin = prediction[i] -> GetMaximumBin(); 
    float m = prediction[i] -> GetBinContent(bin+1); 
    sum += m;  
  }

  can -> Clear();
  gStyle -> SetOptStat(0); 
  truth[0] -> GetYaxis() -> SetRangeUser(1, sum*2);
  truth[0] -> GetXaxis() -> SetRangeUser(0, 10);
  truth[0] -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  len -> SetBorderSize(0); 
  Populate(truth, can, len, kDashed); 
  Populate(prediction, can, len, kSolid); 

}

void PlotHists(std::vector<TH1F*> prediction, std::vector<TH1F*> truth, TString title, TCanvas* can)
{
  float sum = 0; 
  for (int i(0); i < prediction.size(); i++)
  {
    int bin = prediction[i] -> GetMaximumBin(); 
    float m = prediction[i] -> GetBinContent(bin+1); 
    sum += m;  
  }

  can -> Clear();
  gStyle -> SetOptStat(0); 
  gStyle -> SetOptTitle(1); 
  TH1F* empty = (TH1F*)truth[0] -> Clone(title); 
  empty -> SetTitle(title); 
  empty -> Reset();

  empty -> GetYaxis() -> SetRangeUser(1, sum*2);
  empty -> GetYaxis() -> SetTitle("Clusters");
  empty -> GetYaxis() -> CenterTitle();

  empty -> GetXaxis() -> SetRangeUser(0, 10);
  empty -> GetXaxis() -> CenterTitle();

  empty -> Draw("HIST"); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
  len -> SetBorderSize(0); 

  Populate(truth, can, len, kDashed); 
  Populate(prediction, can, len, kSolid); 
}

void PlotHists(TH1F* Data, std::vector<TH1F*> Prediction, std::vector<TH1F*> Truth, TString Title, float FLost_P, float FLost_T, TCanvas* can)
{
  can -> SetLogy(); 
  Data -> SetTitle(Title);
 	Data -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^2]");
  if (Truth.size() == 0){PlotHists(Data, Prediction, can);}
  else{PlotHists(Data, Prediction, Truth, can);}
  
  TPaveText* pt = new TPaveText(.7, .5, 0.9, .7, "NDCNB"); 
  pt -> AddText("#int TrackN:");
  pt -> SetFillStyle(4000);
  for (int i(0); i < Prediction.size(); i++)
  {
    float x = Prediction[i] -> Integral(); 
    std::stringstream stream; 
    stream << std::fixed << std::setprecision(2) << x; 
    TString Line = "Track-"; Line += (i+1); Line += (": "); Line += (stream.str()); 
    pt -> AddText(Line); 
  }
  std::stringstream stream; 
  stream << std::fixed << std::setprecision(4) << FLost_P; 
  TString l = "FLost: "; l += (stream.str()); 
  pt -> AddText(l); 
  pt -> Draw();

  if(FLost_T != 0)
  {
    TPaveText* t = new TPaveText(.5, .5, 0.7, .7, "NDCNB"); 
    t -> AddText("#int TrackN:");
    t -> SetFillStyle(4000);
    for (int i(0); i < Truth.size(); i++)
    {
      float x = Truth[i] -> Integral(); 
      std::stringstream stream; 
      stream << std::fixed << std::setprecision(2) << x; 
      TString Line = "Track-"; Line += (i+1); Line += (": "); Line += (stream.str()); 
      t -> AddText(Line); 
    }
    std::stringstream stream; 
    stream << std::fixed << std::setprecision(4) << FLost_T; 
    TString l = "FLost T: "; l += (stream.str()); 
    t -> AddText(l); 
    t -> Draw();
  }
  can -> Update(); 

}

void RatioPlot(TH1F* H1, TH1F* H2, TCanvas* can)
{
  can -> Clear(); 
  
  TString name = "Ratio Plot: "; name += (H1 -> GetTitle()); name += ("/"); name += H2 -> GetTitle(); 
  TH1F* Ratio = (TH1F*)H1 -> Clone(name);
  Ratio -> SetTitle(name);
  Ratio -> Divide(H2); 
  
  for (int i(0); i < Ratio -> GetNbinsX(); i++)
  {
    if (std::isnan(Ratio -> GetBinContent(i+1)) || std::isnan(H1 -> GetBinContent(i+1)))
    {
      Ratio -> SetBinContent(i+1, 1e-10);
      H1 -> SetBinContent(i+1, 1e-10); 
    }
  }
  
  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.6, 0.75); 
	H1 -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^2]");

  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1.0);
  P1 -> Draw(); 
  P1 -> cd(); 
  H1 -> SetLineColor(kBlack); 
  H1 -> GetYaxis() -> SetRangeUser(1, 2*GetMaxValue(H1));
  H1 -> Draw("HIST"); 
  len -> AddEntry(H1, H1 -> GetTitle() ); 
  len -> Draw("SAME");
  H1 -> SetStats(0);  
  H2 -> SetLineColor(kRed);
  H2 -> Draw("SAMEHIST"); 
  len -> AddEntry(H2, H2 -> GetTitle() );
  len -> Draw("SAME");  
  P1 -> SetLogy();  
  can -> cd();

  TPad *P2 = new TPad("P2", "P2", 0.0, 0.05, 1, 0.3);  
  P2 -> Draw();
  P2 -> cd(); 
	Ratio -> SetStats(0); 
  Ratio -> SetLineColor(kWhite);
  Ratio -> GetYaxis() -> SetRangeUser(0.5, 1.5);
  Ratio -> SetMarkerSize(1); 
  Ratio -> Draw("epl"); 
}

// ====================== Proper histogram plotting for the presentation ======================== //
void GeneratePlot(TH1F* H, TString Title, TCanvas* can, Color_t color, ELineStyle style, TString DrawOption, float Intensity)
{
  gStyle -> SetOptStat(0); 
  if (Title != "")
  {
    H -> SetTitle(Title); 
  }
  H -> SetLineColor(color); 
  H -> SetLineColorAlpha(color, Intensity); 
  H -> SetLineStyle(style); 
  H -> SetLineWidth(1); 
  H -> Draw(DrawOption);
  can -> Update();  
}

TLegend* GenerateLegend(std::vector<TH1F*> Hist_V, TCanvas* can, float x_size, float x_position, float box_size, float size)
{
  TLegend* len = new TLegend(x_size, box_size, x_position, size); 
  len -> SetBorderSize(0); 
  for (TH1F* H : Hist_V)
  {
    len -> AddEntry(H); 
  }
  len -> Draw("SAME");
  can -> Update(); 
  return len; 
}

TLegend* GenerateLegend(std::vector<TGraph*> Hist_V, TCanvas* can, float x_size, float x_position, float box_size, float size)
{
  TLegend* len = new TLegend(x_size, box_size, x_position, size); 
  len -> SetBorderSize(0); 
  for (TGraph* H : Hist_V)
  {
    len -> AddEntry(H);
    std::cout << H -> GetTitle() << std::endl;
  }
  len -> Draw("SAME");
  can -> Update(); 
  return len; 
}

void GenerateRatioPlot(TH1F* H1, TH1F* H2, TCanvas* can, TString Title, TString Xaxis)
{
  TString name = "Ratio Plot: "; name += (H1 -> GetTitle()); name += (" / "); name += H2 -> GetTitle(); 
  TH1F* Ratio = (TH1F*)H1 -> Clone(name);
  Ratio -> SetTitle(name);
  Ratio -> Divide(H2); 

  gStyle -> SetOptStat(0); 
  TLegend* len = new TLegend(0.9, 0.9, 0.5, 0.8); 

  TH1F* Empty = (TH1F*)H1 -> Clone("empty"); 
  Empty -> Reset();
  Empty -> SetTitle(Title);
  
  len -> AddEntry(H1); 
  len -> AddEntry(H2);  
  
  TPad *P1 = new TPad("P1", "P1", 0, 0.3, 1, 1.0);

  if (Xaxis == "LOG")
  {
    P1 -> SetLogy(); 
    Empty -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
  }
  else { Empty -> GetXaxis() -> SetTitle(Xaxis);}

  P1 -> Draw(); 
  P1 -> cd();
  Empty -> Draw("HIST");   
   
  H1 -> Draw("SAMEHIST"); 
  H2 -> Draw("SAMEHIST"); 
  len -> Draw("SAME");  

  can -> cd();

  TPad *P2 = new TPad("P2", "P2", 0.0, 0.05, 1, 0.3);  
  P2 -> Draw();
  P2 -> cd(); 
  Ratio -> SetLineColor(kWhite);
  Ratio -> GetYaxis() -> SetRangeUser(0.25, 1.75);
  Ratio -> SetMarkerSize(1); 
  Ratio -> Draw("epl"); 
}

void PlotGraphs(std::vector<TH1F*> Hists, TString Title, TCanvas* can)
{

  can -> Clear(); 
  can -> SetLogy(false); 
  gStyle -> SetOptStat(0);
  TMultiGraph* mg = new TMultiGraph(); 
  TLegend* len = new TLegend(0.9, 0.9, 0.5, 0.8); 
  mg -> SetTitle(Title); 
  TH1F* H = (TH1F*)Hists[0] -> Clone(""); 
  
  for (int i(0); i < H -> GetNbinsX(); i++){ H -> SetBinContent(i+1, 110); }
  H -> SetLineColor(kWhite); 
  H -> GetYaxis() -> SetRangeUser(60, 200); 
  TGraph* gr = new TGraph(H); 
  mg -> Add(gr); 
  for (int i(0); i < Hists.size(); i++)
  {
    TGraph* G1 = new TGraph(Hists[i]); 
    G1 -> SetLineColor(i+1); 
    mg -> Add(G1); 
    len -> AddEntry(G1);
  }
  mg -> Draw("SAMEAC.");
  mg -> GetXaxis() -> SetTitle("Percentage of Original 2-Truth Normalization (%)"); 
  mg -> GetYaxis() -> SetTitle("Prediction / Truth (%)"); 

  len -> Draw("SAME"); 

  can -> Update(); 
}

void GenerateNiceStacks(std::vector<TH1F*> vec, TString Title, TCanvas* can, TString x_axis, TString y_axis, TString options)
{
  can -> Clear(); 

  THStack* he = new THStack("hs", Title); 
  std::vector<Color_t> col = {kRed, kBlack, kBlue, kGreen}; 
  for (int i(0); i < vec.size(); i++)
  {
    vec[i] -> SetFillColorAlpha(col[i], 0.2); 
    vec[i] -> SetLineColor(col[i]); 
    vec[i] -> GetYaxis() -> SetRangeUser(0.00001, 0.1); 
    vec[i] -> GetXaxis() -> SetTitle(x_axis); 
    vec[i] -> SetLineWidth(1); 
    he -> Add(vec[i]); 
  }

  TH1F* Empty = (TH1F*)vec[0] -> Clone(Title); 
  Empty -> SetTitle(Title); 
  Empty -> GetXaxis() -> SetTitle(x_axis); 
  Empty -> GetYaxis() -> SetTitle(y_axis); 
  Empty -> GetYaxis() -> SetRangeUser(0.0001, 1); 
  Empty -> SetLineColor(kWhite); 
  Empty -> SetFillColorAlpha(kWhite, 0); 
  Empty -> Draw("HIST"); 

  he -> Draw(options);
  he -> GetYaxis() -> SetRangeUser(0.0001, 1); 
  he -> GetXaxis() -> SetTitle(x_axis); 
  he -> Draw(options); 
  can -> Update(); 

  TLegend* le = GenerateLegend(vec, can, 0.89, 0.6, 0.89, 0.81); 
  le -> Draw("SAME"); 
}

TGraph* GenerateGraph(std::vector<float> Input, TString name)
{
  Double_t x[Input.size()]; 
  Double_t y[Input.size()]; 

  for (int i(0); i < Input.size(); i++){ y[i] = Input[i]; x[i] = i; }
  TGraph* Gr = new TGraph(Input.size(), x, y); 

  return Gr; 
}

TGraph* GenerateGraph(std::map<TString, float> Input, TString name, TString yname)
{
  TGraph* Gr = new TGraph(); 
  int i = 0;
  for (TString x : JetEnergy)
  {
    if (Input[x] == 0){continue;}
    Gr -> SetPoint(i, i+0.5, float(Input[x]));
    i++;
  }


  TH1F* H = new TH1F(name, name, Input.size()-2, 0, Input.size()-2); 
  Gr -> SetHistogram(H);
  Gr -> SetTitle(name);
  H -> SetTitle(name);
  H -> GetYaxis() -> SetTitle(yname); 
  H -> GetYaxis() -> CenterTitle(yname); 

  i = 0; 
  for (TString x : JetEnergy)
  {
    if (Input[x] == 0){continue;}
    Gr -> GetXaxis() -> SetBinLabel(i+1, x.ReplaceAll("_GeV", "").ReplaceAll("_", "-"));
    i++;
  }

  return Gr; 
}

std::vector<TGraph*> GenerateMultiGraph(MMF ntru, TString Title, TCanvas* can, float min, float max, TString yname)
{
  std::vector<Color_t> col = {kRed, kYellow, kBlue, kGreen}; 
  int i = 0; 
  std::vector<TGraph*> gr; 
  for (MMFi p = ntru.begin(); p != ntru.end(); p++)
  {
    TString tru = p -> first; 
    MF JE = p -> second;
    TGraph* g = GenerateGraph(JE, tru, yname); 
    g -> GetYaxis() -> SetRangeUser(min, max); 
    g -> GetXaxis() -> SetTitle("Jet Energy (GeV)");
    g -> GetXaxis() -> CenterTitle();

    g -> SetLineColor(col[i]);  
    g -> SetLineWidth(4.);
    gr.push_back(g); 
    
    if (i == 0){ g -> Draw("ALP"); }
    g -> Draw("SAME");
    i++; 
  }
  GenerateLegend(gr, can); 
  gPad -> SetLogy();
  can -> Update();

  return gr; 
}





