#include<PostAnalysis/DerivedFunctions.h>

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, 
                                                     TH1F* Data, 
                                                     RooRealVar* Domain, 
                                                     std::vector<float> Begin, 
                                                     std::vector<float> End, 
                                                     std::vector<TString> Names)
{
  BaseFunctions B;
  std::vector<RooRealVar*> Var = B.RooVariables(Names, Begin, End);
  std::vector<RooHistPdf*> Hist = B.RooPDF(Hists, Domain);
  RooDataHist* data = B.RooData(Data, Domain);

  RooAddPdf model("Model", "Model", B.RooList(Hist), B.RooList(Var));
  RooFitResult* stat = model.fitTo(*data, RooFit::Save(), RooFit::SumW2Error(true));
   
  for (int i(0); i < Var.size(); i++)
  {
    delete Hist[i];
  }
  delete data;
	delete stat; 
  
  return Var;
} 

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, TH1F* Data, float min, float max)
{
  RooRealVar* x = new RooRealVar("x", "x", min, max);
  std::vector<float> Begin(Hists.size(), 0.);
  std::vector<float> End(Hists.size(), 1);
  std::vector<RooRealVar*> Vars = FitToData(Hists, Data, x, Begin, End, Constants::FitNames);
	delete x;
  return Vars; 
}

void DerivedFunctions::SafeScale(std::vector<TH1F*> PDFs, TH1F* Data)
{
  
  int bins = Data -> GetNbinsX(); 
  for (int i(0); i < bins; i++)
  {
    // Get the sum of the PDFs for the bin we are looking at 
    float sum_pdfs = 0; 
    for (int x(0); x < PDFs.size(); x++)
    {
      sum_pdfs += PDFs[x] -> GetBinContent(i+1); 
    }
   
    // Get the number of entries within the data for that bin 
    float data_e = Data -> GetBinContent(i+1); 
    
    // Get the ratio between the data and the PDFs predictions 
    float ratio = data_e/sum_pdfs;  

    // Some fail safe conditions 
    if (sum_pdfs == 0) {ratio = 1e-9;}
    if (std::isnan(ratio)) {ratio = 1e-9;}

    // Update the bins of the PDFs
    for (int x(0); x < PDFs.size(); x++)
    {
      float e = PDFs[x] -> GetBinContent(i+1); 
      PDFs[x] -> SetBinContent(i+1, ratio*e); 
    }
  }
}

void DerivedFunctions::SafeScaleNew(std::vector<TH1F*> PDFs, TH1F* Data)
{
  int bins = Data -> GetNbinsX();
  for (int i(0); i < bins; i++)
  {
    float e = Data -> GetBinContent(i+1);

    float sum = 0;
    for (TH1F* H : PDFs)
    {
      float f = H -> GetBinContent(i+1);
      sum = sum +f;
    }
    if (sum == 0) { sum = 1; } 
    if ( sum > e )
    {
      float ratio = e/sum;
      
      for (TH1F* H : PDFs)
      {
        float f = H -> GetBinContent(i+1);
        H -> SetBinContent(i+1, f*ratio);  
      }
    }
		else if ( sum < e ) 
		{
      float ratio = e/sum;
      
      for (TH1F* H : PDFs)
      {
        float f = H -> GetBinContent(i+1);
        H -> SetBinContent(i+1, f*ratio);  
      }
		}
  }
}






int DerivedFunctions::NumericalShift(TH1F* H1, TH1F* H2)
{
  BaseFunctions B; 

  int bin1 = H1 -> GetNbinsX();
  int bin2 = H2 -> GetNbinsX();
  int delta = bin1 - bin2;
  TH1F* Hist;
  TH1F* Clone;
  int max;
  
  // Means bin1 is bigger  
  if (delta > 0)
  {
    max = bin1;
    Hist = (TH1F*)H1 -> Clone("Hist");
    Hist -> Reset(); 
    B.ShiftExpandTH1F(H2, Hist, 0);
    Clone = (TH1F*)H1 -> Clone("Clone");
  }
  else
  {
    max = bin2;
    Hist = (TH1F*)H2 -> Clone("Hist"); 
    Hist -> Reset();
    B.ShiftExpandTH1F(H1, Hist, 0);
    Clone = (TH1F*)H2 -> Clone("Clone");
  }

  B.Normalize({Hist, Clone});
  TH1F* Temp = (TH1F*)Clone -> Clone("Temp");
  Temp -> Reset(); 
  
  float T = 100; 
  int shifted;
  for (int i(-max); i < max; i++)
  {
    B.ShiftExpandTH1F(Hist, Temp, i);
   
    float e(0); 
    for (int x(0); x < Temp -> GetNbinsX(); x++)
    {
      e += std::abs(Temp -> GetBinContent(x+1) - Clone -> GetBinContent(x+1));
    }
    
    if (T > e)
    {
      T = e; 
      shifted = i+1;
    }
  } 
  if (delta < 0){return -shifted;}

  delete Temp; 
  delete Clone; 
  delete Hist;  
  return shifted;
}

void DerivedFunctions::ReplaceShiftTail(TH1F* Source, TH1F* Target, float offset)
{
  BaseFunctions B;

  int bin_s = Target -> GetNbinsX(); 
  float Lumi = Target -> Integral();
  float Lumi_S = Source -> Integral();
  B.Normalize({Source, Target});

  // Extend the source and target hist to avoid nan 
  std::vector<float> Source_V = B.TH1FDataVector(Source, offset);
  std::vector<float> Target_V = B.TH1FDataVector(Target, offset);  
  TH1F* Source_E = new TH1F("Source_E", "Source_E", Source_V.size(), 0, Source_V.size());
  TH1F* Target_E = new TH1F("Target_E", "Target_E", Target_V.size(), 0, Target_V.size()); 
  B.ToTH1F(Source_V, Source_E); 
  B.ToTH1F(Target_V, Target_E);

  // Define the tail to be replaced and move away from peak 
  int max_bin = Target -> GetMaximumBin() + bin_s*0.01;
  int max_bin_s = Source -> GetMaximumBin();
  TH1F* Temp = (TH1F*)Source -> Clone("Temp"); 
  
  for (int i(0); i < max_bin_s; i++)
  {
    Temp -> SetBinContent(i+1, 0);
  }
  int shift = NumericalShift(Temp, Target); 
  
  float t_S = Source_E -> GetBinContent(max_bin + 1 - shift -1);
  float t_T = Target_E -> GetBinContent(max_bin + 1);
  float sum_S = 0;
  float sum_T = 0;  
  for (int i(max_bin); i < bin_s; i++) 
  {
    float S_e = Source_E -> GetBinContent(i + 1 - shift -1);
    float T_e = Target_E -> GetBinContent(i + 1); 
    sum_S += S_e;
    sum_T += T_e;  
    Target -> SetBinContent(i+1, S_e); //*(t_T/t_S));
  }
 
  Target -> Scale(Lumi);
  Source -> Scale(Lumi_S);
  delete Temp; 
  delete Source_E;
  delete Target_E;
}

std::vector<TH1F*> DerivedFunctions::nTRKGenerator(TH1F* trk1, TH1F* trk2, float offset, int iter)
{
  BaseFunctions B;
  std::vector<float> trk2_V = B.TH1FDataVector(trk2, offset);
  std::vector<float> deconv(trk2_V.size(), 0.5);

  // Generate the histograms as output
  std::vector<TString> Names = {"TRK1", "TRK2", "TRK3", "TRK4", "TRK5"};  
  std::vector<TH1F*> PDFs = B.MakeTH1F(Names, trk2);
 
  // Deconvolve the 2trk data into 1 trk  
  for (int i(0); i < iter; i++)
  {
    deconv = B.LucyRichardson(trk2_V, deconv, deconv); 
  }
  B.ToTH1F(deconv, PDFs[0]);

  // === TRK1
  B.Normalize(PDFs); 
  TH1F* temp = (TH1F*)PDFs[0] -> Clone("temp"); 
  //ReplaceShiftTail(trk1, PDFs[0], offset); 

  // === TRK2
  B.ConvolveHists(PDFs[0], temp, PDFs[1]); 
  B.Normalize(PDFs[1]);

  // === TRK3
  B.ConvolveHists(PDFs[1], temp, PDFs[2]);  
  B.Normalize(PDFs[2]);

  // === TRK4 
  B.ConvolveHists(PDFs[2], temp, PDFs[3]);  
  B.Normalize(PDFs[3]);

  // === TRK5 
  B.ConvolveHists(PDFs[3], temp, PDFs[4]);  
  B.Normalize(PDFs[4]);

  delete temp;
  return PDFs;
}

TH1F* DerivedFunctions::GaussianConvolve(TH1F* Hist, float mean, float stdev, int Toys)
{
  DistributionGenerators DG; 
  BaseFunctions B;

	int bins = Hist -> GetNbinsX(); 
	
	TH1F* H1 = new TH1F("H1", "H1", 2*bins+1, -bins-1, bins); 

	B.ShiftExpandTH1F(Hist, H1, 0.5*bins); 

	TH1F* Gaus = new TH1F("Gaus", "Gaus", 2*bins+1, -bins-0.5, bins+0.5); 
	DG.Gaussian(mean, stdev, Toys, Gaus); 
	B.Normalize(Gaus); 

	TString name = Hist -> GetTitle(); name += ("_C"); 
	TH1F* H_Clone = (TH1F*)H1 -> Clone(name); 
	H_Clone -> Reset();

	TH1F* Out = (TH1F*)Hist -> Clone("Out");
	Out -> Reset();	
	
	B.ConvolveHists(H1, Gaus, H_Clone); 
	B.ShiftExpandTH1F(H_Clone, Out, 0.5*bins+1);	

	B.Normalize(H_Clone);	
	B.Normalize(Out); 

	delete H1;
	delete Gaus;	
	delete H_Clone;
	return Out; 
}

std::map<TString, float> DerivedFunctions::ConvolveFit(TH1F* GxTrk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> Params, float offset, int iter)
{
  auto Deconvolve =[] (TH1F* Hist, TH1F* PDF, TH1F* Gaussian, TH1F* Output, float offset, int iter)
  {
    BaseFunctions B; 
    int bins = PDF -> GetNbinsX();  
   
    // Import the PDFs into the Temp Hist 
    B.ShiftExpandTH1F(PDF, Hist, 0.5*bins); 
		B.Normalize(Hist); 
		B.Normalize(Gaussian);

    // Variables for the LR deconvolution part
    std::vector<float> PDF_V = B.TH1FDataVector(Hist, offset); 
    std::vector<float> PSF_V = B.TH1FDataVector(Gaussian, offset); 
    std::vector<float> deconv(PSF_V.size(), 0.5); 
    
    // Lucy Richardson deconvolution of a gaussian with the PDF
    for (int i(0); i < iter; i++){deconv = B.LucyRichardson(PDF_V, PSF_V, deconv);}
    B.ToTH1F(deconv, Output); 
		B.Normalize(Output);
  }; 

  BaseFunctions B;  
  DistributionGenerators DG; 
   
  // Get the number of bins in the hists. This will be our bin domain 
  int bins = GxTrk -> GetNbinsX(); 

  // Create the temp hists that will be used to do the deconvolution on 
  std::vector<TString> Names;  
  for (int i(0); i < PDFs.size(); i++)
  {
    TString n = PDFs[i] -> GetTitle(); n +=(i); 
    Names.push_back(n); 
  }
  std::vector<TH1F*> Hists = B.MakeTH1F(Names, 2*bins, -bins, bins); 

  // Define an output histogram for the deconvoluted histograms with length equivalent to the bins
  std::vector<TH1F*> PDFs_L = B.MakeTH1F(Names, bins, 0, bins, "_L"); 
  B.ShiftExpandTH1F(Hists, PDFs_L, bins); 

  // Generate the Gaussian Distribution
  TH1F* Gaussian = (TH1F*)Hists[0] -> Clone("Gaussian");
  DG.Gaussian(Params["Gaussian"][0], Params["Gaussian"][1], 500000, Gaussian);

  // Deconvolve the above Gaussian with the track PDF using a multithreaded approach 
  std::vector<std::thread> th; 
  for (int i(0); i < PDFs.size(); i++){th.push_back(std::thread(Deconvolve, Hists[i], PDFs[i], Gaussian, PDFs_L[i], offset, iter));}
  for (std::thread &t : th){ t.join(); } 
	
  // Convert the data TH1F to a common length as the new PDFs
  TH1F* Data_L = new TH1F("Data_L", "Data_L", bins, 0, bins); 
  B.ShiftExpandTH1F(GxTrk, Data_L); 

	// Normalize all the PDFs and set the error for the bins 
	B.Normalize(PDFs_L); 
	B.Normalize(Data_L); 

	for (TH1F* H : PDFs_L)
	{
		for (int i(0); i < bins; i++){H -> SetBinError(i+1, 1e-32);}
	}
	for (int i(0); i < bins; i++){Data_L -> SetBinError(i+1, 1e-32);}

  // Define all the names of the variables we will be needing 
  std::vector<TString> Means_String = { "m1", "m2", "m3", "m4", "m5"};
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4", "s5"};
  std::vector<TString> Gaus_String = { "g1", "g2", "g3", "g4", "g5"};
  std::vector<TString> N_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4", "n_trk5"};
  std::vector<TString> GxT_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4", "P5xG5"};

  // Free defined variables
  std::vector<float> N_Begin(PDFs.size(), 0);
  std::vector<float> N_End(PDFs.size(), 1); 

  // ====== Define the RooRealVars we will need for the fit  
  // Define the range of the fit 
  RooRealVar* x = new RooRealVar("x", "x", 0, bins); 
  std::vector<RooRealVar*> Means = B.RooVariables(Means_String, Params["m_s"], Params["m_e"]); 
  std::vector<RooRealVar*> Stdev = B.RooVariables(Stdev_String, Params["s_s"], Params["s_e"]); 
  std::vector<RooRealVar*> N_Vars = B.RooVariables(N_String, N_Begin, N_End); 
  std::vector<RooGaussian*> G_Vars = B.RooVariables(Gaus_String, Means, Stdev, x); 
  std::vector<RooHistPdf*> PDF_Vars = B.RooPDF(PDFs_L, x);
  std::vector<RooFFTConvPdf*> Conv_Vars = B.RooVariables(GxT_String, PDF_Vars, G_Vars, x); 
 
  // Import the data 
  RooDataHist* trk = new RooDataHist("trk", "trk", *x, Data_L);

  // Create the RooArgSet to import into the model 
  RooArgSet Conv; 
  RooArgSet N;  
  for (int i(0); i < Conv_Vars.size(); i++)
  {
    Conv.add(*Conv_Vars[i]); 
    N.add(*N_Vars[i]);  
  }
  
  RooAddPdf model("model", "model", Conv, N);
  RooFitResult* stat = model.fitTo(*trk, RooFit::Save()); 
   
  // Store the output in a map and delete pointers 
  std::map<TString, float> Out;
  if (!stat){Out["Status"] = -1;}
  else{Out["Status"] = stat -> status();}

  for (int i(0); i < PDFs_L.size(); i++)
  {
    Out[Means_String[i]] = Means[i] -> getVal();
    Out[Stdev_String[i]] = Stdev[i] -> getVal();
    Out[N_String[i]] = N_Vars[i] -> getVal();

    delete Conv_Vars[i]; 
    delete Hists[i];
    delete PDF_Vars[i];
    delete PDFs_L[i];
  }
  
  for (int i(0); i < Stdev.size(); i++)
  {
    delete Means[i]; 
    delete G_Vars[i]; 
    delete Stdev[i]; 
    delete N_Vars[i]; 
  }


  delete trk; 
  delete x;
  delete Data_L;
  delete stat;

  return Out;    

}

std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> DerivedFunctions::MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure)
{
  // Lambda Function: This does the convolution and takes care of the labels 
  auto MakeGaussianConvoluted = [](TString n, std::vector<TH1F*> PDFs, std::map<TString, float> trkn_Params, TH1F* trkn_L, float Gamma)
  {
    BaseFunctions B;
    DerivedFunctions DF;
    const std::vector<TString> Names = {"n_trk1", "n_trk2", "n_trk3", "n_trk4", "n_trk5"}; 
    const std::vector<TString> Stdev = {"s1", "s2", "s3", "s4", "s5"}; 
    const std::vector<TString> Mean = {"m1", "m2", "m3", "m4", "m5"};
    
    std::vector<TH1F*> GxTrkN(Names.size());
    for (int i(0); i < Names.size(); i++)
    {
      TString Name = n; Name += Names[i];
      TH1F* trknXg = (TH1F*)trkn_L -> Clone(Name);
      trknXg -> SetTitle(Name);
      trknXg -> Reset();
      TH1F* trkn_H = DF.GaussianConvolve(PDFs[i], 0, trkn_Params[Stdev[i]]);
      B.ShiftExpandTH1F(trkn_H, trknXg);
      delete trkn_H;

      float lumi = trkn_L -> Integral();
      
      // Check for nan
      if (std::isnan(trkn_Params[Names[i]])){Gamma = 0;}
      trknXg -> Scale(Gamma*trkn_Params[Names[i]]*lumi);
      GxTrkN[i] = trknXg;
    }
    return GxTrkN;
  };

  auto DynamicVariables = [] (std::map<TString, std::vector<float>> *Params, std::map<TString, float> FitPredictions)
  {  
    std::vector<TString> Means_String = {"m1", "m2", "m3", "m4", "m5"};
    std::vector<TString> Stdev_String = {"s1", "s2", "s3", "s4", "s5"};
   
    std::map<TString, std::vector<float>> &P = *Params;
    for ( int i(0); i < Means_String.size(); i++ )
    {
      TString M = Means_String[i]; 
      TString S = Stdev_String[i]; 
      
      float P_M = FitPredictions[M];  
      float P_S = FitPredictions[S]; 
 
      if (P["m_s"][i] == P_M){P["m_s"][i] = P_M++;}
      if (P["m_e"][i] == P_M){P["m_e"][i] = P_M++;}
      if (P["s_s"][i] == P_S){P["s_s"][i] = P_S++;}
      if (P["s_e"][i] == P_S){P["s_e"][i] = P_S*1.5;}
    } 
  };

  // Lambda Function: Condensed version of the above. Wanted to use this as a multithread function but RooFit doesnt seem to be multithreading safe 
  auto Parallel = [&MakeGaussianConvoluted, &DynamicVariables](TH1F* trk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> *Params, float offset, int iter, std::vector<TH1F*> GxTrk, TString name)
  {
    DerivedFunctions DF; 
    std::map<TString, float> trk_Param = DF.ConvolveFit(trk, PDFs, *Params, offset, iter);
    std::vector<TH1F*> Out = MakeGaussianConvoluted(name, PDFs, trk_Param, trk, 1);
    DynamicVariables(Params, trk_Param);
    for (int i(0); i < Out.size(); i++){GxTrk.push_back(Out[i]);}
    return GxTrk;
  };

  BaseFunctions B;
  Plotting P; 
  TH1::AddDirectory(false);
  
  // Create histograms for the ntrack data
  TH1F* trk1_L = new TH1F("trk1_L", "trk1_L", ntrk[0] -> GetNbinsX(), 0, ntrk[0] -> GetNbinsX());
  TH1F* trk2_L = new TH1F("trk2_L", "trk2_L", ntrk[1] -> GetNbinsX(), 0, ntrk[1] -> GetNbinsX());
  TH1F* trk3_L = new TH1F("trk3_L", "trk3_L", ntrk[2] -> GetNbinsX(), 0, ntrk[2] -> GetNbinsX());
  TH1F* trk4_L = new TH1F("trk4_L", "trk4_L", ntrk[3] -> GetNbinsX(), 0, ntrk[3] -> GetNbinsX());
  TH1F* trk5_L = new TH1F("trk5_L", "trk5_L", ntrk[4] -> GetNbinsX(), 0, ntrk[4] -> GetNbinsX());

  // Fill those histograms with the ntrack data
  B.ShiftExpandTH1F(ntrk[0], trk1_L); 
  B.ShiftExpandTH1F(ntrk[1], trk2_L);  
  B.ShiftExpandTH1F(ntrk[2], trk3_L); 
  B.ShiftExpandTH1F(ntrk[3], trk4_L); 
  B.ShiftExpandTH1F(ntrk[4], trk5_L); 

  // Make clones of these for restoring the state after each iteration 
  TH1F* trk1_L_C = (TH1F*)trk1_L -> Clone("trk1_L_C");
  TH1F* trk2_L_C = (TH1F*)trk2_L -> Clone("trk2_L_C");
  TH1F* trk3_L_C = (TH1F*)trk3_L -> Clone("trk3_L_C");
  TH1F* trk4_L_C = (TH1F*)trk4_L -> Clone("trk4_L_C");
  TH1F* trk5_L_C = (TH1F*)trk5_L -> Clone("trk5_L_C");
      
  // Forward declaration 
  std::map<TString, std::vector<float>> Params1;
  std::map<TString, std::vector<float>> Params2;
  std::map<TString, std::vector<float>> Params3;
  std::map<TString, std::vector<float>> Params4;
  std::map<TString, std::vector<float>> Params5;

  std::vector<TH1F*> GxTrk1_C;
  std::vector<TH1F*> GxTrk2_C;
  std::vector<TH1F*> GxTrk3_C;
  std::vector<TH1F*> GxTrk4_C;
  std::vector<TH1F*> GxTrk5_C;

  std::vector<TH1F*> PDFs1;
  std::vector<TH1F*> PDFs2;
  std::vector<TH1F*> PDFs3;
  std::vector<TH1F*> PDFs4;
  std::vector<TH1F*> PDFs5;


  // ============= Output Variables ============= //
  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> Output; 
  TH1F* FLOST_Prediction = new TH1F("FLost_Pred", "FLost_Pred", cor_loop, 0, cor_loop); 
  TH1F* FLOST_Truth = new TH1F("FLost_Truth", "FLost_Truth", cor_loop, 0, cor_loop); 

  TCanvas* can_HD = new TCanvas();
  TString Title_HD = "out.pdf";
  for (int x(0); x < cor_loop; x++)
  { 
    // Forward declaration 
    std::vector<TH1F*> GxTrk1;   
    std::vector<TH1F*> GxTrk2;
    std::vector<TH1F*> GxTrk3;
    std::vector<TH1F*> GxTrk4; 
    std::vector<TH1F*> GxTrk5; 
    std::vector<TH1F*> PDFs;  

    // Create separate Params containers so that we can add some dynamics 
    if ( x == 0)
    {
      Params1 = Params; 
      Params2 = Params; 
      Params3 = Params; 
      Params4 = Params; 
      Params5 = Params; 
    }

    if (x <= 2)
    {
      // Generate the minimal version of the PDFs without any Gaussian convolution 
      PDFs = nTRKGenerator(trk1_L, trk2_L, offset, 100);

      // Generate the Gaussian Parameters for the convolution from a fit
      GxTrk1 = Parallel(trk1_L, PDFs, &Params1, offset, 100, GxTrk1, "GxTrk1");  
      GxTrk2 = Parallel(trk2_L, PDFs, &Params2, offset, 100, GxTrk2, "GxTrk2"); 
      GxTrk3 = Parallel(trk3_L, PDFs, &Params3, offset, 100, GxTrk3, "GxTrk3");
      GxTrk4 = Parallel(trk4_L, PDFs, &Params4, offset, 100, GxTrk4, "GxTrk4"); 
      GxTrk5 = Parallel(trk5_L, PDFs, &Params5, offset, 100, GxTrk5, "GxTrk5"); 

      for (int v(0); v < GxTrk1_C.size(); v++)
      {
        delete GxTrk1_C[v];
        delete GxTrk2_C[v];
        delete GxTrk3_C[v];
        delete GxTrk4_C[v];
        delete GxTrk5_C[v];
      }

      for (int i(0); i < PDFs.size(); i++)
      {
        delete PDFs[i];
      }
    }
    else
    {
      PDFs = nTRKGenerator(trk1_L, trk2_L, offset, 100);

      // Generate the Gaussian Parameters for the convolution from a fit
      GxTrk1 = Parallel(trk1_L, PDFs, &Params1, offset, 100, GxTrk1, "GxTrk1");  
      GxTrk2 = Parallel(trk2_L, PDFs, &Params2, offset, 100, GxTrk2, "GxTrk2"); 
      GxTrk3 = Parallel(trk3_L, PDFs, &Params3, offset, 100, GxTrk3, "GxTrk3");
      GxTrk4 = Parallel(trk4_L, PDFs, &Params4, offset, 100, GxTrk4, "GxTrk4"); 
      GxTrk5 = Parallel(trk5_L, PDFs, &Params5, offset, 100, GxTrk5, "GxTrk5"); 
      
      for (int v(0); v < GxTrk1_C.size(); v++)
      {
        delete GxTrk1_C[v];
        delete GxTrk2_C[v];
        delete GxTrk3_C[v];
        delete GxTrk4_C[v];
        delete GxTrk5_C[v];
        //delete PDFs[v];
      }

    }

    GxTrk1_C = B.CopyTH1F(GxTrk1, "_C");
    GxTrk2_C = B.CopyTH1F(GxTrk2, "_C");
    GxTrk3_C = B.CopyTH1F(GxTrk3, "_C");
    GxTrk4_C = B.CopyTH1F(GxTrk4, "_C");
    GxTrk5_C = B.CopyTH1F(GxTrk5, "_C");

    // Reset the state of the data histograms 
    trk1_L -> Reset();
    trk1_L -> Add(trk1_L_C);     
    trk2_L -> Reset();
    trk2_L -> Add(trk2_L_C); 
    trk3_L -> Reset();
    trk3_L -> Add(trk3_L_C);
    trk4_L -> Reset();
    trk4_L -> Add(trk4_L_C);
    trk5_L -> Reset();
    trk5_L -> Add(trk5_L_C);

    // Scale the PDFs according to the number of events in a bin 
    SafeScale(GxTrk1, trk1_L); 
    SafeScale(GxTrk2, trk2_L); 
    SafeScale(GxTrk3, trk3_L); 
    SafeScale(GxTrk4, trk4_L); 
    SafeScale(GxTrk5, trk5_L); 

    // Testing if temperature helps....
    float heat = 1; //float(x)/float(cor_loop);
    std::cout << heat << std::endl;

    // Do the subtraction  
    trk1_L -> Add(GxTrk1[1], -heat); 
    trk1_L -> Add(GxTrk1[2], -heat); 
    trk1_L -> Add(GxTrk1[3], -heat); 
    trk1_L -> Add(GxTrk1[4], -heat); 

    trk2_L -> Add(GxTrk2[0], -heat);
    trk2_L -> Add(GxTrk2[2], -heat);
    trk2_L -> Add(GxTrk2[3], -heat);
    trk2_L -> Add(GxTrk2[4], -heat);

    trk3_L -> Add(GxTrk3[0], -heat);
    trk3_L -> Add(GxTrk3[1], -heat);
    trk3_L -> Add(GxTrk3[3], -heat);
    trk3_L -> Add(GxTrk3[4], -heat);

    trk4_L -> Add(GxTrk4[0], -heat);
    trk4_L -> Add(GxTrk4[1], -heat);
    trk4_L -> Add(GxTrk4[2], -heat);
    trk4_L -> Add(GxTrk4[4], -heat);



    std::cout << "################### " << x << std::endl;
   
    // ==== Trk1 
    std::vector<TString> Names1 = {"trk1_F1", "trk2_F1", "trk3_F1", "trk4_F1", "trk5_F1"}; 
    std::vector<TH1F*> Trk1_PDFs = B.MakeTH1F(Names1, ntrk[0]); 
    B.ShiftExpandTH1F(GxTrk1, Trk1_PDFs);
    
    // ==== Trk2
    std::vector<TString> Names2 = {"trk1_F2", "trk2_F2", "trk3_F2", "trk4_F2", "trk5_F2"}; 
    std::vector<TH1F*> Trk2_PDFs = B.MakeTH1F(Names2, ntrk[1]); 
    B.ShiftExpandTH1F(GxTrk2, Trk2_PDFs);

    // ==== Trk3 
    std::vector<TString> Names3 = {"trk1_F3", "trk2_F3", "trk3_F3", "trk4_F3", "trk5_F3"}; 
    std::vector<TH1F*> Trk3_PDFs = B.MakeTH1F(Names3, ntrk[2]); 
    B.ShiftExpandTH1F(GxTrk3, Trk3_PDFs);

    // ==== Trk4 
    std::vector<TString> Names4 = {"trk1_F4", "trk2_F4", "trk3_F4", "trk4_F4", "trk5_F4"}; 
    std::vector<TH1F*> Trk4_PDFs = B.MakeTH1F(Names4, ntrk[3]); 
    B.ShiftExpandTH1F(GxTrk4, Trk4_PDFs);

    // Create the variabls needed for the Plotting below  
    std::vector<TString> Names_Sub = {"trk1_Sub", "trk2_Sub", "trk3_Sub", "trk4_Sub"}; 
    std::vector<TH1F*> Tracks = B.MakeTH1F(Names_Sub, ntrk[2]); 
    B.ShiftExpandTH1F({trk1_L, trk2_L, trk3_L, trk4_L}, Tracks);
    std::vector<TH1F*> Truth = {Closure[0][0], Closure[1][1], Closure[2][2], Closure[3][3]};


    can_HD -> Clear();
    TH1F* Temp = (TH1F*)Tracks[1] -> Clone("Lol"); 
    can_HD -> SetWindowSize(1200, 1200); 
    can_HD -> Divide(2,2);
    P.PlotHists({Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs}, Closure, Tracks, can_HD);


    // ======== Section for the output ========== //
    // === Save the PDFs and the subtracted Hists

    // === Save the subtracted Hists
    std::vector<TString> Names_Pure = {"trk1_Pure", "trk2_Pure", "trk3_Pure", "trk4_Pure"};    
    std::vector<TH1F*> Track_Out = B.MakeTH1F(Names_Pure, ntrk[2]); 
    B.ShiftExpandTH1F({trk1_L, trk2_L, trk3_L, trk4_L}, Track_Out);

    std::vector<TH1F*> Pred_PDF;
    std::vector<std::vector<TH1F*>> Set = {Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs, Track_Out, ntrk, Closure[0], Closure[1], Closure[2], Closure[3]}; 


    for (std::vector<TH1F*> PDF_S : Set)
    {
      for (int g(0); g < PDF_S.size(); g++)
      {
        TString name_pdf = PDF_S[g] -> GetName(); name_pdf += (x); 
        Pred_PDF.push_back((TH1F*)PDF_S[g] -> Clone(name_pdf));
      }
    }

    
    // === Save the FLost progress
    float FLost_Pred = B.FLost(ntrk, {Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs}); 
    float FLost_MC = B.FLost({ntrk[0], ntrk[1], ntrk[2], ntrk[3]}, Closure);   
    FLOST_Prediction -> SetBinContent(x+1, FLost_Pred); 
    FLOST_Truth -> SetBinContent(x+1, FLost_MC);   
    
    TString name_out =  "FLost_P.at."; name_out += (x); 
    TH1F* FL_P_Copy = (TH1F*)FLOST_Prediction -> Clone(name_out);

    name_out =  "FLost_T.at."; name_out += (x); 
    TH1F* FL_T_Copy = (TH1F*)FLOST_Truth -> Clone(name_out);
    Pred_PDF.push_back(FL_T_Copy); 

    Output[x] = std::make_pair(FL_P_Copy, Pred_PDF);
 
    // Clean up memory  
    for (int i(0); i < GxTrk2.size(); i++)
    {
       
      delete GxTrk1[i]; 
      delete GxTrk2[i];      
      delete GxTrk3[i];
      delete GxTrk4[i];

      delete Trk1_PDFs[i];
      delete Trk2_PDFs[i];
      delete Trk3_PDFs[i];
      delete Trk4_PDFs[i];   
    }

   // delete Temp; 
  }
  return Output;
}


