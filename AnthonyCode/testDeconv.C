#include <vector>
#include <iterator>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TComplex.h>
#include <TVirtualFFT.h>
#include <TFractionFitter.h>
#include <TLegend.h>

#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"


std::vector<float>  convolveWithAwithB( std::vector<float>::iterator beginA,
                                        std::vector<float>::iterator endA,
                                        std::vector<float>::iterator beginB,
                                        std::vector<float>::iterator endB)
{
 
  int nBins = std::distance( beginA, endA);
  int nBinsB = std::distance( beginB, endB);
  std::vector<float>  output(nBins,0);
  if(nBinsB != nBins){
    std::cout << "ERROR bin size don't match  --  Update code if you want to do this" <<  std::endl;
    output.resize(0);
    return output;
  }
  
  TVirtualFFT* fftr2cA = TVirtualFFT::FFT(1, &nBins, "R2C K P");
  TVirtualFFT* fftr2cB = TVirtualFFT::FFT(1, &nBins, "R2C K P");
  for (Int_t i=0 ; i<nBins ; ++i, ++beginA, ++beginB) {
    fftr2cA->SetPoint(i,*beginA,0);
    fftr2cB->SetPoint(i,*beginB,0);
  }
  fftr2cA->Transform();
  fftr2cB->Transform();

  TVirtualFFT* fftc2r  = TVirtualFFT::FFT(1, &nBins, "C2R K P");
  // Loop over first half +1 of complex output results, multiply
  // and set as input of reverse transform
  for (Int_t i=0 ; i<nBins/2+1 ; i++) {
    Double_t re1,re2,im1,im2 ;
    fftr2cA->GetPointComplex(i,re1,im1) ;
    fftr2cB->GetPointComplex(i,re2,im2) ;
    Double_t re = re1*re2 - im1*im2 ;
    Double_t im = re1*im2 + re2*im1 ;
    TComplex t(re,im) ;
    fftc2r->SetPointComplex(i,t) ;
  }

  // Reverse Complex->Real FFT transform product
  fftc2r->Transform() ;
  for (Int_t i=0 ; i<nBins ; i++) {
    Double_t re1,im1;
    fftc2r->GetPointComplex(i,re1,im1);
    output[i] = re1;
  }

  delete fftr2cA;
  delete fftr2cB;
  delete fftc2r;

  return output;
}

std::vector<float>  convolveWithOnesSelf( std::vector<float>::iterator begin, std::vector<float>::iterator end)
{
  return convolveWithAwithB( begin, end, begin, end);
}


std::vector<float> RichardsonLucy1DIter( const std::vector<float>& data,  //Target
                                       const std::vector<float>& current, //Current underlying image
                                       const std::vector<float>& psf, // Point spread function
                                       float dampen){

  //I think this could be written entirely in Fouriour space (and would be more efficient) but for now do it in image space
  // I guess you could also write this in Matrix notation;
  
  // Current Limitation  -- PSF is always assumed shift things to the right
  // Need to update code to allow double sided movement
  size_t offset = (data.size() - current.size() ) /2;
  std::vector<float>  next(current.size(),0);
  //  For each bin
  for( size_t i(0); i < current.size();++i){
    float sum_j = 0;
    size_t upperLimitJ = psf.size();
    for( size_t j(i); j < upperLimitJ; ++j ){
      float c_j = 0;
      size_t upperLimitK = j;
      // Mi
      for( size_t k(0); k <= upperLimitK; ++k ){
        c_j += psf[j-k]*current[k]; //  Nominal point spead function * the image guess to estimate how much flows out
      }
      if(c_j!=0)
        sum_j += data[j] / c_j  * psf[j-i]; // Flipped point spread function to estimate how much flows in
    }
    // Slow any ossiclations by dampening the correction a little
    next[i] = current[i] * ( 1 + dampen * (sum_j - 1) );
    if(next[i] < 0. ||  std::isnan(next[i]) || std::isinf(next[i]) )
      next[i] = 0.;
  }
  return next;
}

std::vector<float> RichardsonLucy1DIter( const std::vector<float>& data,
                                      const std::vector<float>& current,
                                      float dampen)
{
  // V special case of Richardson-Lucy algorithm for when the point spread function and the underlying image are one in the same thing.
  return RichardsonLucy1DIter( data, current, current, dampen);
}




double fitFraction(TH1* data,
            TH1* template2Part,
            TH1* template3Part,
            TCanvas* can)
{
  // Create WS
  RooWorkspace ws("ws");
  // Define observable
  RooRealVar* x = (RooRealVar*) ws.factory("x[0.5,20]");
  // Import data
  RooDataHist dataHist_Data( "dataHist_Data",
                             "dataHist_Data",
                             RooArgList(*x),
                             data );
  // Import templates
  RooDataHist dataHist_2( "dataHist_2",
                          "dataHist_2",
                          RooArgList(*x),
                          template2Part );
  
  RooDataHist dataHist_3( "dataHist_3",
                           "dataHist_#",
                           RooArgList(*x),
                          template3Part );
  
  RooHistPdf pdfHist_2( "pdfHist_2",
                        "pdfHist_2",
                        RooArgList(*x),
                        dataHist_2);
  RooHistPdf pdfHist_3( "pdfHist_3",
                        "pdfHist_3",
                        RooArgList(*x),
                        dataHist_3);
  
  ws.import(pdfHist_2);
  ws.import(pdfHist_3);

  // Build fit model
  ws.factory("SUM::pdf_Total(coeff_2[0.8,0.1,1]*pdfHist_2,pdfHist_3)");
  auto pdfTotal = ws.pdf("pdf_Total");
  // Fit data
  auto fitOk = pdfTotal->fitTo( dataHist_Data,
                                RooFit::PrintLevel(0),
                                RooFit::SumW2Error(true),
                                RooFit::Save() );
  if( !fitOk || fitOk->status() > 1 )
  {
    return -1;
  }
  
  
  // Draw some stuff
  if(!can){
    delete fitOk;
    return ws.var("coeff_2")->getVal();
  }
  can->cd();
  auto frameFit = x->frame();
  dataHist_Data.plotOn(frameFit);
  pdfTotal->plotOn(frameFit, RooFit::LineColor(kBlue));
  pdfTotal->plotOn(frameFit, RooFit::Components("pdfHist_2") ,
                                                RooFit::LineColor(kViolet),
                                                RooFit::LineStyle(kDashed));
  pdfTotal->plotOn(frameFit, RooFit::Components("pdfHist_3") ,
                                                RooFit::LineColor(kRed),
                                                RooFit::LineStyle(kDashed));
  frameFit->SetMinimum(data->Integral()/data->GetNbinsX()*1e-4);
  
  can->SetLogy();
  frameFit->Draw();
  
  TLatex p;
  p.SetNDC();
  p.SetTextFont(42);
  p.SetTextSize(0.05);
  p.SetTextColor(kBlack);
  p.DrawLatex(0.55,0.85,Form("Fraction 2: %0.3f", ws.var("coeff_2")->getVal())) ;
  
  can->Print("RooFitCan.pdf");
  delete frameFit;

  delete fitOk;
  return ws.var("coeff_2")->getVal();
}

void testDeconvolution(TH1* trueTarget, TH1* trueConv1, TH1* trueConv2,
                       TH1* resultTarget, TH1* resultConv1, TH1* resultConv2)
{
  size_t nBins = trueConv1->GetNbinsX();
  //  Store histogram data in vector
  size_t upperMirrorSize = nBins*0.1;
  std::vector<float>  data(nBins+upperMirrorSize,0);
  std::vector<float>  deconv(nBins+upperMirrorSize,1);
  for (size_t i=0 ; i <nBins ; i++) {
    data[i] = trueConv1->GetBinContent(i+1);
  }
  for (size_t i=0 ; i < upperMirrorSize ; i++) {
    data[ nBins + i ] = data[ nBins - 1 - i ];
  }
   
   
  int iter(0);
  // Perform deconvolution  and store residual in graph to track progress
  TGraph* gr = new TGraph();
  while( iter < 100 ){
    deconv = RichardsonLucy1DIter( data, deconv, 0.75 );
    // Check to make deconvolution is converging but convolving the obtained PDF with its self
    auto testConv  = convolveWithOnesSelf( deconv.begin(), deconv.end()-upperMirrorSize );
   
    double sum2(0);
    double normTest(0);
    double normData(0);
    for( size_t i(0); i<nBins;++i ){
      normTest += testConv[i];
      normData += data[i];
    }
    normTest = 1/normTest;
    normData = 1/normData;
    for( size_t i(0); i<nBins;++i ){
      sum2 += pow( testConv[i] * normTest - data[i] * normData , 2 );
    }
    gr->SetPoint(iter,iter,sum2);
    ++iter;
  }
  
  auto testConv1 = convolveWithOnesSelf( deconv.begin(),
                                         deconv.end()-upperMirrorSize );
  
  auto testConv2 = convolveWithAwithB( testConv1.begin(),
                                       testConv1.end(),
                                       deconv.begin(),
                                       deconv.end()-upperMirrorSize );
 
  for (size_t i=0 ; i<nBins ; i++) {
    resultTarget->SetBinContent(i+1, deconv[i]);
    resultConv1->SetBinContent(i+1, testConv1[i]);
    resultConv2->SetBinContent(i+1, testConv2[i]);
  }

  
  TCanvas* grCan =  new TCanvas("grCan","grCan",1200,600);
  grCan->Divide(2,1);
  grCan->cd(1)->SetLogy();
  auto mainHist = resultTarget->DrawNormalized();
  mainHist->GetYaxis()->SetRangeUser(1e-6,1e-1);
  trueTarget->DrawNormalized("same");
  resultConv1->DrawNormalized("same");
  trueConv1->DrawNormalized("same");
  resultConv2->DrawNormalized("same");
  trueConv2->DrawNormalized("same");
  
  grCan->cd(2)->SetLogy();
  gr->Draw("ap");
  gr->GetXaxis()->SetTitle("Iteration");
}


std::vector<float> ShiftReplace(TH1D* trk1, std::vector<float> deconv)
{
  // The deconv vector is of 550 length and the trk1 TH1D has a length of 500.
  // Need to standardize the two lengths. 
  // Create two hists of nbins 550, iterate all entries in trk1 until exhausted and replace excess with deconv 

  float Max = 20; 
  float Min = 0;  
  float nBins = deconv.size();
  TH1D* Decon = new TH1D("Decon", "Decon", nBins, Min, Max);
  TH1D* TRK1 = new TH1D("TRK1", "TRK1", nBins, Min, Max);

  // Get length of trk1 histogram 
  int len = trk1 -> GetNbinsX();

  // Convert the std::vector to a TH1D for fitting purposes. 
  for (int i(0); i < nBins; i++)
  {
    float e = deconv.at(i);
    if(i <= len){TRK1 -> SetBinContent(i+1, trk1 -> GetBinContent(i+1));}
    else{TRK1 -> SetBinContent(i+1, e);} 
    Decon -> SetBinContent(i+1, e);
  }
 
  // Get the post peak bin
  int m_b_de = Decon -> GetMaximumBin()+1;
  int m_b_tr = TRK1 -> GetMaximumBin()+1;
  
  // Derive the conversion factor bin -> dEdx
  float ss = (Max - Min)/nBins;

  // ================================ Here we perform the fitting of the tail ============================= //
  // Define the variables being used for the fit
  RooRealVar x("x", "x", 0, 20);
  RooRealVar s("s", "s", 0., -4, 4);
  RooFormulaVar delta("delta", "x-s", RooArgSet(x,s));

  // Define the fitting range 
  x.setRange("FinalFit", (m_b_de)*ss, 19);
  
  // RooData Hists used for the fit 
  RooDataHist original("original", "original", x, Decon);
  RooDataHist Data("Data", "Data", x, TRK1);
  
  // Create the PDF for the modelling of the shift in the distributions 
  RooHistPdf model("model", "model", delta, x, original, 1);
  
  // Fit the model 
  RooFitResult* Fitress = model.fitTo(Data, RooFit::Save(), RooFit::Range("FinalFit")); 
 
  // ==== Plot related stuff ^(0_o)^ ====== // 
  RooPlot* xframe = x.frame(RooFit::Title("Shift Test"));
  original.plotOn(xframe, RooFit::Name("Decon"), RooFit::LineColor(kAzure)); 
  Data.plotOn(xframe, RooFit::Name("1Track"), RooFit::LineColor(kRed)); 
  model.plotOn(xframe, RooFit::Name("Model"), RooFit::LineColor(kBlue)); 

  float shift = s.getVal();
  int Shift = std::round((shift)/((Max-Min)/nBins)); 

  // ============== Debug/Experimental ========================== // 

  TString name = "Shifted: "; name += (Shift); name += ("-"); name +=(std::rand());
  TCanvas* Cans = new TCanvas(name,name, 800, 800);
  Cans -> Divide(2,1);
  Cans -> cd(1) -> SetLogy();
  xframe -> SetMinimum(1); 
  xframe -> Draw();
  Cans -> cd(2) -> SetLogy();
  Decon -> SetLineColor(kGreen);
  TRK1 -> SetLineColor(kBlue);
  Decon -> Draw("SAMEHIST");
  TRK1 -> Draw("SAMEHIST"); 
  Cans -> Print("out.pdf");

  // ================================= Tail Replacement =================================================== //
  int Repl = m_b_tr + 5;
  for (int i(Repl); i < nBins; i++)
  {
    float e = trk1 -> GetBinContent(i+1);
    Decon -> SetBinContent(i-Shift, e);
  } 
 
  // Extract the entries from Decon and place inside a vector 
  std::vector<float> deconv_raw; 
  for (int i(0); i < nBins; i++)
  {
    deconv_raw.push_back(Decon -> GetBinContent(i+1));
  }
  
  delete TRK1; 

  return deconv_raw;

}


void testDeconv()
{
  size_t nBins = 500;
  // Histograms to hold information
  TH1D* histA  = new TH1D("histA", ";x;y;", nBins,0,20);
  TH1D* histB  = new TH1D("histB", ";x;y;", nBins,0,20);
  TH1D* histC  = new TH1D("histC", ";x;y;", nBins,0,20);
  TH1D* histD  = new TH1D("histD", ";x;y;", nBins,0,20);
  
  TH1D* histAs = new TH1D("histAs",";x;y;", nBins,0,20);
  TH1D* histBs = new TH1D("histBs",";x;y;", nBins,0,20);
  TH1D* histCs = new TH1D("histCs",";x;y;", nBins,0,20);
  TH1D* histDs = new TH1D("histDs",";x;y;", nBins,0,20);
  
  TH1D* histAg  = new TH1D("histAg", ";x;y;", nBins,0,20);
  TH1D* histBg  = new TH1D("histBg", ";x;y;", nBins,0,20);
  TH1D* histCg  = new TH1D("histCg", ";x;y;", nBins,0,20);

  // True PDF of the energy deposition
  TF1 myfunc("mylandau","landau",0,20);
  myfunc.SetParameter(0,1);
  myfunc.SetParameter(1,0.9);
  myfunc.SetParameter(2,0.1);
  TF1 myDummyFunc("myGaus","landau",0,20);
  myDummyFunc.SetParameter(0,1);
  myDummyFunc.SetParameter(1,0.91);
  myDummyFunc.SetParameter(2,0.075);

  
  
  // Fill some distributions
  for( int i(0); i < 500000; ++i){
    double ran1 = myfunc.GetRandom();
    double ran2 = myfunc.GetRandom();
    double ran3 = myfunc.GetRandom();
    double ran4 = myfunc.GetRandom();
    histA->Fill(ran1);
    histB->Fill(ran1+ran2);
    histC->Fill(ran1+ran2+ran3);
    histD->Fill(ran1+ran2+ran3+ran4);

    ran1 = myDummyFunc.GetRandom();
    ran2 = myDummyFunc.GetRandom();
    ran3 = myDummyFunc.GetRandom();    
    histAg->Fill( ran1 );
    histBg->Fill( ran1 + ran2);
    histCg->Fill( ran1 + ran2 + ran3);

  }
  
  // Draw the true PDFs
  TCanvas* can =  new TCanvas("can","can",800,800);
  can->SetLogy();
 
  histA->SetMarkerSize(0);
  histB->SetMarkerSize(0);
  histC->SetMarkerSize(0);
  histD->SetMarkerSize(0);
  histA->SetFillColorAlpha(kBlue,0.5);
  histB->SetFillColorAlpha(kRed,0.5);
  histC->SetFillColorAlpha(kGreen,0.5);
  histD->SetFillColorAlpha(kOrange,0.5);
 
  
  // Run a test to make sure things working ok
  testDeconvolution( histA,  histB,  histC,
                     histAs, histBs, histCs);
  
 
  //  Create some dummy data
  double faction3 = 0.2;
  size_t upperMirrorSize = nBins*0.1;
  TH1D* dummyDataHist = (TH1D*)histB->Clone("twoWithThree");
  dummyDataHist->Add( histC, faction3);
  
  can->Print("RooFitCan.pdf[");
  //  Check what it should really look like by performing the fit on data first
  double rooFitFraction = fitFraction( dummyDataHist,
                                        histB,
                                        histC,
                                        can);

  //  Now fit the data with our first guess of the 2 part and 3 part pdfs
  rooFitFraction = fitFraction( dummyDataHist,
                                       histBg,
                                       histCg,
                                       can);
  
  //  Lets subtract away the estimated  3 part contribution and go from there
  double dataNorm = dummyDataHist->Integral();
  double nThreePart = histCg->Integral();
  
  // Store data in a vector
  std::vector<float> data(nBins+upperMirrorSize,0);
  for (size_t i=0 ; i <nBins ; i++) {
    data[i] = dummyDataHist->GetBinContent(i+1);
  }
  for (size_t i=0 ; i < upperMirrorSize ; i++) {
    data[ nBins + i ] = data[ nBins - i - 1];
  }
  
  auto dataCopy =  data;
  // Subtract away the 3 part from the 2 Track and mirror the result
  // the factot of 0.4 is a complete fudge but is basically there to make suer
  // that there isn't a giant hole in the distribution
  for (size_t i=0 ; i <nBins ; i++) {
    dataCopy[i] -= 0.4 * dataNorm/nThreePart
                       * (1-rooFitFraction)*histCg->GetBinContent(i+1);
  }
  for (size_t i=0 ; i < upperMirrorSize ; i++) {
    dataCopy[ nBins + i ] = dataCopy[ nBins - i - 1];
  }

  // Histograms to store the results.
  TH1D* testConv0Hist = (TH1D*)histAs->Clone("testConv0");
  TH1D* testConv1Hist = (TH1D*)histBs->Clone("testConv1");
  TH1D* testConv2Hist = (TH1D*)histCs->Clone("testConv2");
  testConv0Hist->SetMarkerColor(kBlue);
  testConv1Hist->SetMarkerColor(kRed);
  testConv2Hist->SetMarkerColor(kGreen);
  testConv0Hist->SetLineColor(kBlue);
  testConv1Hist->SetLineColor(kRed);
  testConv2Hist->SetLineColor(kGreen);
  testConv0Hist->SetLineWidth(3);
  testConv1Hist->SetLineWidth(3);
  testConv2Hist->SetLineWidth(3);
 
  TH1D* trk1 = (TH1D*)histAs -> Clone("trk1"); 
  //  Ok lets do some iterations of the unfolding
  TGraph  graphFrac2;
  int iterFit = 0;
  //  Lets do the unfolding --  starting with a flat prior
  std::vector<float> deconv = std::vector<float>(nBins+upperMirrorSize,0.5);
  while( iterFit < 25 ){
    int iterDeconv = 0;
    // Need to develop some early termination criteria
    std::vector<float> deconvCopy = deconv;
    for(size_t j(0); j < deconv.size(); ++j )
    {
      double nEle(0);
      deconv[j] = 0;
      size_t lowerBound = j-10 > 0 ? j-10 : 0;
      size_t upperBound = j+10 > deconv.size() ? deconv.size() : j+10;
      for(size_t k = lowerBound; k< upperBound; ++k ){
        deconv[j] += deconvCopy[k];
        nEle += 1.;
      }
      if(nEle>0)
        deconv[j] /= nEle;
      
    }
    while( iterDeconv < 40 ){
      deconv = RichardsonLucy1DIter( dataCopy, deconv, 0.75 );
      ++iterDeconv;
    }

    // Tail replacement of deconv with the one from histA (i.e. the pure 1 trk histogram) 
    deconv = ShiftReplace(trk1, deconv); // <------------------------------- Mine 


    auto testConv2 = convolveWithOnesSelf( deconv.begin(),
                                           deconv.end()-upperMirrorSize );
       
    auto testConv3 = convolveWithAwithB( testConv2.begin(),
                                         testConv2.end(),
                                         deconv.begin(),
                                         deconv.end()-upperMirrorSize );
    
    
    //  Fill normlized hists
    double normTest1(0);
    double normTest2(0);
    double normTest3(0);
    for( size_t i(0); i<nBins;++i ){
      normTest1 += deconv[i];
      normTest2 += testConv2[i];
      normTest3 += testConv3[i];
    }
    normTest1 = 1/normTest1;
    normTest2 = 1/normTest2;
    normTest3 = 1/normTest3;
    for( size_t i(0); i<nBins;++i ){
      testConv0Hist->SetBinContent(i+1,deconv[i]*normTest1);
      testConv0Hist->SetBinError(i+1,1e-12);
      testConv1Hist->SetBinContent(i+1,testConv2[i]*normTest2);
      testConv1Hist->SetBinError(i+1,1e-12);
      testConv2Hist->SetBinContent(i+1,testConv3[i]*normTest3);
      testConv2Hist->SetBinError(i+1,1e-12);
    }
        
    double rooFitFraction = fitFraction( dummyDataHist,
                                         testConv1Hist,
                                         testConv2Hist,
                                         can);
    
    //  Something needs to be done here force convergence here
    if (rooFitFraction > 0) {                       // check on fit status
      //  Create new 2 paricle data by  subtracting away 3 particle
      for (size_t i=0 ; i <nBins ; i++) {
        dataCopy[i] =  data[i] -  dataNorm * (1-rooFitFraction)* testConv2Hist->GetBinContent(i+1);
        if(dataCopy[i] < 0)
          dataCopy[i] = 0;
      }
      for (size_t i=0 ; i < upperMirrorSize ; i++) {
        dataCopy[ nBins + i ] = data[ nBins - i - 1 ];
      }
      graphFrac2.SetPoint( iterFit, iterFit+1, rooFitFraction);
    }
    ++iterFit;
  }
  can->SetLogy(false);
  graphFrac2.Draw("ap");
  graphFrac2.GetXaxis()->SetTitle("Iteration");
  can->Print("RooFitCan.pdf");

  can->SetLogy(true);
  histA->SetLineWidth(0);
  histB->SetLineWidth(0);
  histC->SetLineWidth(0);
  histA->DrawNormalized("h");
  histB->DrawNormalized("sameh");
  histC->DrawNormalized("sameh");
  testConv0Hist->DrawNormalized("samelc");
  testConv1Hist->DrawNormalized("samelc");
  testConv2Hist->DrawNormalized("samelc");

  
  can->Print("RooFitCan.pdf)");
}
