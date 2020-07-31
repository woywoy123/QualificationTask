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


std::vector<float> TailReplace(TH1F* hist, std::vector<float> deconv, float min, float max)
{
  float n_H = hist -> GetNbinsX();
  float n_D = deconv.size();
  float Pad = (n_D - n_H) / 2.;
  float ss = (max - min)/n_D;

  TH1F* Deconv = new TH1F("Deconv", "Deconv", n_D, min, max);
  TH1F* Hist = (TH1F*)Deconv -> Clone("Hist"); 
  
  // Fill the histogram with deconv data
  for (int i(0); i < n_D; i++)
  {
    Deconv -> SetBinContent(i+1, deconv[i]);
    if (i < Pad || n_H + Pad < i)
    {
      Hist -> SetBinContent(i + 1, 0);
    }
    else
    {
      Hist -> SetBinContent(i + 1, hist -> GetBinContent(i + 1 - Pad));
    }
  }

  // Get the domain of the fit by checking the peak and the tail length
  float max_H = Hist -> GetMaximumBin();
  float max_D = Deconv -> GetMaximumBin();
  float min_H = max_H;
  float min_D = max_D;
  float val_H = Hist -> GetBinContent(max_H);
  float val_D = Deconv -> GetBinContent(max_D);

  for (int i(0); i < n_D; i++)
  {
    float H = Hist -> GetBinContent(i+1);
    float D = Deconv -> GetBinContent(i+1);

    if (i > max_H && H < val_H)
    {
      val_H = H;
      min_H = i+1;
    }

    if (i > max_D && D < val_D)
    {
      val_D = D;
      min_D = i+1;
    }

  }

  float peak_delta = max_H - max_D;

//  float start(0);
//  float end(0);
//  if (peak_delta > 0)  
//  {  
//    start = (max_D - Pad +1)*ss;
//    end = (min_D - Pad - peak_delta)*ss;
//  }
//  else
//  {
  double b = 2; //min + (max_H)*ss; <<< ------ This part has a weird bug. When I change this number, so does the variable e not sure why...
  double e = 18;//min + (min_H -10)*ss;
  

  //TCanvas* can = new TCanvas();
  //can -> SetLogy();
  //Hist -> Draw("SAMEHIST");
  //Deconv -> Draw("SAMEHIST");

  RooRealVar x("x", "x", min, max);

  RooRealVar s("s", "s", 0., -6, 6);
  RooFormulaVar delta("delta", "x-s", RooArgSet(x,s));
  
  RooDataHist D("Decon", "Decon", x, Deconv); 
  RooHistPdf model("Model", "Model", delta, x, D, 1);

  RooDataHist trk1("trk1", "trk1", x, Hist);

  RooFitResult* Result = model.fitTo(trk1, Range(1.3, 18));

  RooPlot* xframe = x.frame(RooFit::Title("Shift Test"));
  D.plotOn(xframe, RooFit::Name("Decon"), RooFit::LineColor(kAzure)); 
  trk1.plotOn(xframe, RooFit::Name("1Track"), RooFit::LineColor(kRed)); 
  model.plotOn(xframe, RooFit::Name("Model"), RooFit::LineColor(kBlue)); 
  
  TCanvas* can = new TCanvas();
  can -> SetLogy();
  xframe -> SetMinimum(1);
  xframe -> Draw();

  float shift = s.getVal();
  float S = std::round(shift/ss);
  std::cout << "###################### " << S << std::endl; 
  std::cout << "Start" << b << " ::: " << e << std::endl;
  std::cout << peak_delta << std::endl;

  std::vector<float> dec;
  return dec;

}


void ShiftDetect(TH1F* hist, std::vector<float> deconv, float min, float max)
{

  float n_H = hist -> GetNbinsX();
  float n_D = deconv.size();
  float Pad = (n_D - n_H) / 2.;
  float ss = (max - min)/n_D;

  TH1F* Deconv = new TH1F("Deconv", "Deconv", n_D, min, max);
  TH1F* Hist = (TH1F*)Deconv -> Clone("Hist"); 
  
  // Fill the histogram with deconv data
  for (int i(0); i < n_D; i++)
  {
    Deconv -> SetBinContent(i+1, deconv[i]);
    if (i < Pad || n_H + Pad < i)
    {
      Hist -> SetBinContent(i + 1, 0);
    }
    else
    {
      Hist -> SetBinContent(i + 1, hist -> GetBinContent(i + 1 - Pad));
    }
  }

  // Get the domain of the fit by checking the peak and the tail length
  float max_H = Hist -> GetMaximumBin();
  float max_D = Deconv -> GetMaximumBin();
  float min_H = max_H;
  float min_D = max_D;
  float val_H = Hist -> GetBinContent(max_H);
  float val_D = Deconv -> GetBinContent(max_D);

  for (int i(0); i < n_D; i++)
  {
    float H = Hist -> GetBinContent(i+1);
    float D = Deconv -> GetBinContent(i+1);

    if (i > max_H && H < val_H)
    {
      val_H = H;
      min_H = i+1;
    }

    if (i > max_D && D < val_D)
    {
      val_D = D;
      min_D = i+1;
    } 
  }

  TH1F* Hist_Clip = (TH1F*)Hist -> Clone("Clip_H");
  TH1F* Deco_Clip = (TH1F*)Deconv -> Clone("Clip_D");
  Deco_Clip -> Reset();
  Hist_Clip -> Reset();
  Deco_Clip -> SetLineColor(kYellow);  
   
  // Peak and minimum positions; max_H, min_H :: max_D, min_D
  //Clip the pre-tail of the histograms 
  for (int i(0); i < n_D; i++)
  {
    if ( i >= max_D )
    {
      Deco_Clip -> SetBinContent(i, Deconv -> GetBinContent(i));
    }
    
    if ( i >= max_H )
    {
      Hist_Clip -> SetBinContent(i, Hist -> GetBinContent(i));
    }
  }

 
  TH1F* temp = (TH1F*)Hist -> Clone("Temp");  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy();
  for (int i(0); i < n_D; i++)
  {
    temp -> Reset();
    for (int y(i); y < n_D - i; y++)
    {
      float e = Hist_Clip -> GetBinContent(y+1);
      temp -> SetBinContent(i+1+y, e);
    } 
    Deco_Clip -> Draw("SAMEHIST");
    temp -> Draw("SAMEHIST");
    can -> Update();       
    
  }
  


  std::cout << "fin" << std::endl;



}












































// Test case the shift test and tail replace
std::vector<float> TestShift(std::vector<TH1F*> Hists, TH1F* reference, float min, float max)
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

    //std::vector<float> something = TailReplace(reference, deconv, min, max);
    ShiftDetect(reference, deconv, min, max); 
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
  std::vector<TString> Shift_N = {"S1"};
  float Padding = 100;
  float step_size = (20.0-0.0)/500.;
  float max = 20 + step_size*Padding;
  float min =  - step_size*Padding;
  
  std::vector<TH1F*> Shift_Hist = F.MakeTH1F(Shift_N, 500 + 2*Padding, min, max);
  std::vector<float> Shift_Numbers = {50};

  // TCanvas* can = new TCanvas();
  // can -> SetLogy();
  for (int i(0); i < Shift_Numbers.size(); i++)
  {
    for (int x(0); x < Shift_Hist[i] -> GetNbinsX(); x++)
    {
      if ( x < Padding + Shift_Numbers[i] || Shift_Hist[i] -> GetNbinsX() - Padding + Shift_Numbers[i] <= x)
      { 
        Shift_Hist[i] -> SetBinContent(x + 1, 0);
      }
      else
      {
        Shift_Hist[i] -> SetBinContent(x + 1, trk1 -> GetBinContent(x+1 - Padding - Shift_Numbers[i])); 
      }
    }
    Shift_Hist[i] -> SetLineColor(Constants::Colors[i]);
    // Shift_Hist[i] -> Draw("SAMEHIST");
    // can -> Update();
  }
  
  std::vector<float> Shift_Results = TestShift(Shift_Hist, trk1, min, max); 
  
}
