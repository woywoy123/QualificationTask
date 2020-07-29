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

std::vector<float> TailReplace(TH1F* Hist, std::vector<float> deconv, float min, float max)
{

  // Check that the vector and the hist have the same length 
  float n_H = Hist -> GetNbinsX();
  float n_V = deconv.size(); 
  float ss = (max - min)/n_H;
  int offset = 0;

  // Define the histograms used in this code 
  TH1F* Deconv = new TH1F("Deconv", "Deconv", n_H + 2*offset, min - ss*offset, max + ss*offset);  
  TH1F* hist = (TH1F*)Deconv -> Clone("hist");
  for (int i(0); i < n_H + 2*offset; i++)
  {
    if ( i < offset || n_H + offset <= i) 
    {
      Deconv -> SetBinContent(i+1, 0);
      hist -> SetBinContent(i+1, 0);   
    }
    else
    {
      Deconv -> SetBinContent(i+1, deconv[i-offset]);
      hist -> SetBinContent(i+1, Hist -> GetBinContent(i +1 - offset));
    }
  }
  
  // Find the min max positions of the hist    
  float max_b = Deconv -> GetMaximumBin();
  float min_D = Deconv -> GetBinContent(max_b);
  int min_b(0);
  for (int i(max_b); i < n_H; i++)
  {
    float e = Deconv -> GetBinContent(i+1);
    if (e < min_D)
    {
      min_D = e;
      min_b = i+1;
    } 
  }

  // Define the step size being used for the dEdx
  float start = ss*(max_b + offset);
  float end = ss*(min_b + offset); 


  // ========================= Perform the fitting =============================== //
  // Define the variables for RooFit 
  RooRealVar x("x", "x", min, max); // Used for fitting range
  RooRealVar s("s", "s", -20, 20); // Determine the shifts between histos
  RooFormulaVar delta("delta", "x-s", RooArgSet(x,s)); // Relation between x and s

  // Define the fitting range
  x.setRange("Signal", 2, 10);

  // Create the PDF and the model for fitting 
  RooDataHist D("Decon", "Decon", x, Deconv);
  RooHistPdf model("model", "model", delta, x, D, 1);

  // Define the data hist for the fit
  RooDataHist trk1("trk1", "trk1", x, hist);

  // Perform the fit 
  RooFitResult* Result = model.fitTo(trk1, RooFit::Save(), RooFit::Range("Signal"), SumW2Error(true));

  
  // ============== Debug/Experimental ========================== // 

  float shift = s.getVal();
  float S = shift/ss;

  RooPlot* xframe = x.frame(RooFit::Title("Shift Test"));
  trk1.plotOn(xframe, RooFit::LineColor(kOrange), RooFit::LineStyle(kSolid), RooFit::YErrorSize(0));
  D.plotOn(xframe, RooFit::LineColor(kRed), RooFit::LineStyle(kSolid), RooFit::YErrorSize(0));
  model.plotOn(xframe, RooFit::LineColor(kGreen), RooFit::LineStyle(kSolid));

  TCanvas* can = new TCanvas();
  can -> SetLogy();
  xframe -> SetMinimum(1);
  xframe -> Draw("HIST");

 
//  for (int i(0); i < n_H; i++)
//  {
//    Deconv -> SetBinContent(i+1, deconv[i]);
//  }
//
//  float l_H = Hist -> Integral();
//  float l_D = Deconv -> Integral();




  std::vector<float> dec;
  return dec;
}

std::vector<float> TestShift(std::vector<TH1F*> Hists, TH1F* reference)
{
  Fit_Functions f;

  std::vector<float> Results = {};
  for (TH1F* hist : Hists)
  {
    std::vector<float> deconv(hist -> GetNbinsX(), 0);
    for (int i(0); i < hist -> GetNbinsX(); i++)
    { 
      deconv[i] = hist -> GetBinContent(i+1);  
    }

    std::vector<float> something = TailReplace(reference, deconv, 0, 20);
   
  }  
  return Results; 
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

  // Fit testing section
  std::vector<TH1F*> PDF = {trk1, trk2, trk3, trk4};
  TH1F* Data1 = (TH1F*)trk1 -> Clone("Data1");
  TH1F* Data2 = (TH1F*)trk2 -> Clone("Data2");
  TH1F* Data3 = (TH1F*)trk3 -> Clone("Data3");
  TH1F* Data4 = (TH1F*)trk4 -> Clone("Data4");
  std::vector<float> Fit_test1 = TestFit(PDF, Data1);
  std::vector<float> Fit_test2 = TestFit(PDF, Data2);
  std::vector<float> Fit_test3 = TestFit(PDF, Data3);
  std::vector<float> Fit_test4 = TestFit(PDF, Data4);  

  // Test Tail replace
  std::vector<TString> Shift_N = {"S1"};
  std::vector<TH1F*> Shift_Hist = F.MakeTH1F(Shift_N, 500, 0, 20);
  std::vector<float> Shift_Numbers = {5};

  //TCanvas* can = new TCanvas();
  //can -> SetLogy();
  for (int i(0); i < Shift_Numbers.size(); i++)
  {
    for (int x(0); x < Shift_Hist[i] -> GetNbinsX(); x++)
    {
      Shift_Hist[i] -> SetBinContent(x + 1 + Shift_Numbers[i], trk1 -> GetBinContent(x+1)); 
    }
    //Shift_Hist[i] -> SetLineColor(Constants::Colors[i]);
    //Shift_Hist[i] -> Draw("SAMEHIST");
    //can -> Update();
  }
  
  std::vector<float> Shift_Results = TestShift(Shift_Hist, trk1); 
  


}
