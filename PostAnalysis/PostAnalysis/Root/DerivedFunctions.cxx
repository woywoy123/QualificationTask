#include<PostAnalysis/DerivedFunctions.h>

std::vector<RooRealVar*> DerivedFunctions::FitToData(std::vector<TH1F*> Hists, 
                                                     TH1F* Data, 
                                                     RooRealVar* Domain, 
                                                     std::vector<float> Begin, 
                                                     std::vector<float> End, 
                                                     std::vector<TString> Names)
{
  BaseFunctions B;
	Plotting P; 
	 
  std::vector<RooRealVar*> Var = B.RooVariables(Names, Begin, End);
  std::vector<RooHistPdf*> Hist = B.RooPDF(Hists, Domain);
  RooDataHist* data = B.RooData(Data, Domain);

  RooAddPdf model("Model", "Model", B.RooList(Hist), B.RooList(Var));
  RooFitResult* stat = model.fitTo(*data, RooFit::SumW2Error(true), RooFit::Save());
 	
  for (int i(0); i < Hist.size(); i++)
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
  std::vector<float> Begin(Hists.size(), 0);
  std::vector<float> End(Hists.size(), Data-> Integral());
  std::vector<RooRealVar*> Vars = FitToData(Hists, Data, x, Begin, End, Constants::FitNames);
	delete x;
  return Vars; 
}

void DerivedFunctions::FitToData(std::vector<TH1F*> PDFs, std::vector<TH1F*> Data, float min, float max)
{
  for (int i(0); i < Data.size(); i++)
  {
    RooRealVar* x = new RooRealVar("x", "x", min, max); 
    std::vector<RooRealVar*> v = FitToData({PDFs[i]}, Data[i], x, {0}, {Data[i] -> Integral()}, {PDFs[i] -> GetTitle()}); 
    PDFs[i] -> Scale(v[0] -> getVal());  
  }
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
    if ( sum > e )
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
  int max_bin = Target -> GetMaximumBin() + bin_s*0.05;
  int max_bin_s = Source -> GetMaximumBin();
  TH1F* Temp = (TH1F*)Source -> Clone("Temp"); 
  
  for (int i(0); i < max_bin_s; i++)
  {
    Temp -> SetBinContent(i+1, 0);
  }
  int shift = NumericalShift(Temp, Target); 
  
  float t_S = Source_E -> GetBinContent(max_bin + 1 - shift);
  float t_T = Target_E -> GetBinContent(max_bin + 1);
  float sum_S = 0;
  float sum_T = 0;  
  for (int i(max_bin); i < bin_s; i++) 
  {
    float S_e = Source_E -> GetBinContent(i - 2 - shift);
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
  std::vector<TString> Names = {"TRK1", "TRK2", "TRK3", "TRK4", "TRK5", "TRK6"};  
  std::vector<TH1F*> PDFs = B.MakeTH1F(Names, trk2);
 
  // Deconvolve the 2trk data into 1 trk  
  //for (int i(0); i < iter; i++)
  //{
  //  deconv = B.LucyRichardson(trk2_V, deconv, deconv, 1); 
  //}
  //B.ToTH1F(deconv, PDFs[0]);

  // === TRK1
  B.Normalize(PDFs); 
  B.ShiftExpandTH1F(trk1, PDFs[0]);  
  //ReplaceShiftTail(trk1, PDFs[0], offset); 
	//B.ResidualRemove(PDFs[0]); 

  // === TRK2
  B.ConvolveHists(PDFs[0], PDFs[0], PDFs[1]); 
  B.Normalize(PDFs[1]);
	B.ResidualRemove(PDFs[1]); 

  // === TRK3
  B.ConvolveHists(PDFs[1], PDFs[0], PDFs[2]);  
  B.Normalize(PDFs[2]);
	B.ResidualRemove(PDFs[2]); 

  // === TRK4 
  B.ConvolveHists(PDFs[2], PDFs[0], PDFs[3]);  
  B.Normalize(PDFs[3]);
	B.ResidualRemove(PDFs[3]); 

  // === TRK5 
  B.ConvolveHists(PDFs[3], PDFs[0], PDFs[4]);  
  B.Normalize(PDFs[4]);
	B.ResidualRemove(PDFs[4]); 

  // === TRK6 
  B.ConvolveHists(PDFs[4], PDFs[0], PDFs[5]);  
  B.Normalize(PDFs[5]);
	B.ResidualRemove(PDFs[5]); 

  // === TRK7 
  //B.ConvolveHists(PDFs[5], PDFs[0], PDFs[6]);  
  //B.Normalize(PDFs[6]);
	//B.ResidualRemove(PDFs[4]); 
	
  return PDFs;
}

std::vector<TH1F*> DerivedFunctions::nTRKFrom1Trk(TH1F* trk1)
{
  BaseFunctions B; 

  std::vector<TString> Names = {"TRK1", "TRK2", "TRK3", "TRK4", "TRK5", "TRK6"};  
  std::vector<TH1F*> PDFs = B.MakeTH1F(Names, trk1);
  PDFs[0] -> Add(trk1);  

  for (int i(0); i < Names.size()-1; i++)
  {
    B.ConvolveHists(trk1, PDFs[i], PDFs[i+1]); 
    B.Normalize(PDFs);
    B.ResidualRemove(PDFs[i+1]);
  }

  return PDFs; 
}

TH1F* DerivedFunctions::GaussianConvolve(TH1F* Hist, float mean, float stdev, int Toys)
{
  DistributionGenerators DG; 
  BaseFunctions B;

	int bins = Hist -> GetNbinsX(); 
	
	TH1F* H1 = new TH1F("H1", "H1", 2*bins, -bins-1, bins); 

	B.ShiftExpandTH1F(Hist, H1, 0.5*bins); 

	TH1F* Gaus = new TH1F("Gaus_T", "Gaus_T", 2*bins, -bins-1, bins); 
	DG.Gaussian(mean, stdev, Toys, Gaus); 
	B.Normalize(Gaus); 

	TString name = Hist -> GetTitle(); name += ("_C"); 
	TH1F* H_Clone = (TH1F*)H1 -> Clone(name); 
	H_Clone -> Reset();

	TH1F* Out = (TH1F*)Hist -> Clone("Out");
	Out -> Reset();	
	
	B.ConvolveHists(H1, Gaus, H_Clone); 
	B.ShiftExpandTH1F(H_Clone, Out, 0.5*bins);	

	B.Normalize(H_Clone);	
	B.Normalize(Out); 

	delete H1;
	delete Gaus;	
	delete H_Clone;
	return Out; 
}

std::vector<float> DerivedFunctions::Deconvolve(TH1F* Hist, TH1F* PDF, TH1F* Gaussian, TH1F* Output, float offset, int iter)
{
  BaseFunctions B; 
  int bins = PDF -> GetNbinsX();  
 
  // Import the PDFs into the Temp Hist 
  B.ShiftExpandTH1F(PDF, Hist, 0.5*bins); 
	B.Normalize(Hist); 

  // Variables for the LR deconvolution part
  std::vector<float> PDF_V = B.TH1FDataVector(Hist, offset); 
  std::vector<float> PSF_V = B.TH1FDataVector(Gaussian, offset); 
  std::vector<float> deconv(PSF_V.size(), 0.5); 
  std::vector<float> Converge_Vector; 
   
  // Lucy Richardson deconvolution of a gaussian with the PDF
  for (int i(0); i < iter; i++)
  {
    std::vector<float> deconv_old = deconv; 
    deconv = B.LucyRichardson(PDF_V, PSF_V, deconv, 1);
   // calculate the difference between deconv old and new 
   float diff = 0; 
   for (int x(0); x < deconv.size(); x++)
   {
     float e = deconv[x];
     float f = deconv_old[x];
     diff = diff + f-e;
   }
   Converge_Vector.push_back(diff); 
   if (std::abs(diff) < 1e-9){break;}     
  }
 
  B.ToTH1F(deconv, Output);  
  return Converge_Vector;
}


std::vector<TH1F*>  DerivedFunctions::ConvolveFit(TH1F* GxTrk, std::vector<TH1F*> PDFs, std::map<TString, std::vector<float>> Params, float offset, int iter)
{
  auto Deconvolve =[] (TH1F* Hist, TH1F* PDF, TH1F* Gaussian, TH1F* Output, float offset, int iter)
  {
    BaseFunctions B; 
    int bins = PDF -> GetNbinsX();  
   
    // Import the PDFs into the Temp Hist 
    B.ShiftExpandTH1F(PDF, Hist, 0.5*bins); 
		B.Normalize(Hist); 

    // Variables for the LR deconvolution part
    std::vector<float> PDF_V = B.TH1FDataVector(Hist, offset); 
    std::vector<float> PSF_V = B.TH1FDataVector(Gaussian, offset); 
    std::vector<float> deconv(PSF_V.size(), 0.5); 
    
    // Lucy Richardson deconvolution of a gaussian with the PDF
     
    for (int i(0); i < iter; i++)
    {
      std::vector<float> deconv_old = deconv; 
      deconv = B.LucyRichardson(PDF_V, PSF_V, deconv, 1);
     // calculate the difference between deconv old and new 
     float diff = 0; 
     for (int x(0); x < deconv.size(); x++)
     {
       float e = deconv[x];
       float f = deconv_old[x];
       diff = diff + f-e;
     }
     if (std::abs(diff) < 1e-8){break;}     
    }
   
    B.ToTH1F(deconv, Output); 
  }; 

  BaseFunctions B;  
  DistributionGenerators DG; 
  
	// Capture the integrals of the PDFs and the data inserted 
	float Data_Lumi = GxTrk -> Integral(); 
	 
  // Get the number of bins in the hists. This will be our bin domain 
  int bins = GxTrk -> GetNbinsX(); 

  // Create the temp hists that will be used to do the deconvolution on 
	std::vector<TString> Names;
  for (int i(0); i < PDFs.size(); i++)
  {
    TString n = PDFs[i] -> GetTitle(); n +=(i); 
 		Names.push_back(n);  
	}

  std::vector<TH1F*> Hists = B.MakeTH1F(Names, 2*bins, -bins, bins, "_H"); 

  // Define an output histogram for the deconvoluted histograms with length equivalent to the bins
  std::vector<TH1F*> PDFs_L = B.MakeTH1F(Names, bins, 0, bins, "_L"); 
  //B.ShiftExpandTH1F(Hists, PDFs_L, bins); 

  // Generate the Gaussian Distribution
  TH1F* Gaus = new TH1F("Gaus", "Gaus", 2*bins, -bins-1, bins);
  DG.Gaussian(Params["Gaussian"][0], Params["Gaussian"][1], 0, Gaus);
	B.Normalize(Gaus);
   
  // Deconvolve(Hists[0], PDFs[0], Gaus, PDFs_L[0], offset, iter);  
  std::vector<std::thread> th; 
  for (int i(0); i < PDFs.size(); i++){th.push_back(std::thread(Deconvolve, Hists[i], PDFs[i], Gaus, PDFs_L[i], offset, iter));}
  for (std::thread &t : th){ t.join(); } 
  	
  // Convert the data TH1F to a common length as the new PDFs
  TH1F* Data_L = new TH1F("Data_L", "Data_L", bins, 0, bins); 
  B.ShiftExpandTH1F(GxTrk, Data_L); 

  //B.Normalize(Data_L);
  B.Normalize(PDFs_L);

  // Define all the names of the variables we will be needing 
  std::vector<TString> Means_String = { "m1", "m2", "m3", "m4", "m5", "m6"};
  std::vector<TString> Stdev_String = { "s1", "s2", "s3", "s4", "s5", "s6"};
  std::vector<TString> Gaus_String = { "g1", "g2", "g3", "g4", "g5", "g6"};
  std::vector<TString> N_String = { "n_trk1", "n_trk2", "n_trk3", "n_trk4", "n_trk5", "n_trk6"};
  std::vector<TString> GxT_String = { "P1xG1", "P2xG2", "P3xG3", "P4xG4", "P5xG5", "P6xG6"};

  // Free defined variables
  std::vector<float> N_Begin(PDFs.size(), 0.);
  std::vector<float> N_End(PDFs.size(), Data_L -> Integral()); 

  // ====== Define the RooRealVars we will need for the fit  
  // Define the range of the fit 
  RooRealVar* x = new RooRealVar("x", "x", 0, bins); 
  std::vector<RooRealVar*> Means = B.RooVariables(Means_String, Params["m_s"], Params["m_e"]); 
  std::vector<RooRealVar*> Stdev = B.RooVariables(Stdev_String, Params["s_s"], Params["s_e"]); 
  std::vector<RooRealVar*> N_Vars = B.RooVariables(N_String, N_Begin, N_End); 
  std::vector<RooGaussian*> G_Vars = B.RooVariables(Gaus_String, Means, Stdev, x); 
  std::vector<RooHistPdf*> PDF_Vars = B.RooPDF(PDFs_L, x);
  std::vector<RooFFTConvPdf*> Conv_Vars = B.RooVariables(GxT_String, G_Vars, PDF_Vars, x); 
 
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
  model.fitTo(*trk, RooFit::SumW2Error(true)); 	

  //Plotting P; 
  //P.PlotHists(model, x, Conv_Vars, trk);  
 
   
  std::vector<TH1F*> Out;	 
  for (int i(0); i < PDFs_L.size(); i++)
  {
    float M = Means[i] -> getVal();
    float S = Stdev[i] -> getVal();
    float N = N_Vars[i] -> getVal(); 

		PDFs[i] -> Reset(); 
		TH1F* GxT = GaussianConvolve(PDFs_L[i], M*0.5, S);

    GxT -> Scale(N); 
    Out.push_back(GxT);

    delete Conv_Vars[i]; 
    delete Hists[i];
    delete PDF_Vars[i];
    delete Means[i]; 
    delete G_Vars[i]; 
    delete Stdev[i]; 
    delete N_Vars[i]; 
		delete PDFs_L[i];
	}
  for (int i(0); i < Out.size(); i++)
  {
    PDFs[i] -> Reset(); 
		B.ShiftExpandTH1F(Out[i], PDFs[i], 0); 
		delete Out[i]; 
  }

  delete trk; 
  delete x;
	delete Data_L;
  delete Gaus;

	return PDFs;    
}

std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> DerivedFunctions::MainAlgorithm(std::vector<TH1F*> ntrk, std::map<TString, std::vector<float>> Params, float offset, int iter, int cor_loop, std::vector<std::vector<TH1F*>> Closure)
{

	Plotting P; 
	BaseFunctions B;

	// Forward Declaration: Distributions after each convolution.
	std::vector<TH1F*> GxTrk1;   
  std::vector<TH1F*> GxTrk2;
  std::vector<TH1F*> GxTrk3;
  std::vector<TH1F*> GxTrk4; 
  std::vector<TH1F*> PDFs;  

  //// ============= Output Variables ============= //
  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> Output; 
  
	// FLost iteration progress	
	TH1F* FLOST_Prediction = new TH1F("FLost_Pred", "FLost_Pred", cor_loop, 0, cor_loop); 
  TH1F* FLOST_Truth = new TH1F("FLost_Truth", "FLost_Truth", cor_loop, 0, cor_loop); 

	// Create clones of the original data that will be used to restore the estimate from data
	TH1F* trk1 = (TH1F*)ntrk[0] -> Clone("trk1");
	TH1F* trk2 = (TH1F*)ntrk[1] -> Clone("trk2");
	TH1F* trk3 = (TH1F*)ntrk[2] -> Clone("trk3");
	TH1F* trk4 = (TH1F*)ntrk[3] -> Clone("trk4");

	// Some plotting variables for checking 	
  TCanvas* can_HD = new TCanvas();
  TString Title_HD = "out.pdf";

	for (int i(0); i < cor_loop; i++)
	{
		
		if (i < 1)
		{
      // Generate the minimal version of the PDFs without any Gaussian convolution 
      PDFs = nTRKGenerator(trk1, trk2, offset, iter);
			
			// Now we clone these PDFs into their own sets
			GxTrk1 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk1");
			GxTrk2 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk2");			
			GxTrk3 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk3");		
			GxTrk4 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk4");		

			// Now we run these through the convolution 
			GxTrk1 = ConvolveFit(trk1, GxTrk1, Params, offset, iter); 
			GxTrk2 = ConvolveFit(trk2, GxTrk2, Params, offset, iter); 
			GxTrk3 = ConvolveFit(trk3, GxTrk3, Params, offset, iter); 
			GxTrk4 = ConvolveFit(trk4, GxTrk4, Params, offset, iter); 
			
			for (int c(0); c < PDFs.size(); c++)
			{
				delete PDFs[c];
			}

		}
		else
		{
      // Generate the minimal version of the PDFs without any Gaussian convolution 
      PDFs = nTRKGenerator(trk1, trk2, offset, iter);
			
			// Now we clone these PDFs into their own sets
			//GxTrk1 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk1");
			//GxTrk2 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk2");			
			//GxTrk3 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk3");		
			//GxTrk4 = B.CopyTH1F({PDFs[0], PDFs[1], PDFs[2], PDFs[3], PDFs[4], PDFs[5]}, "_trk4");		

      // Reset the state of the data histograms 
      trk1 -> Reset();
      trk2 -> Reset();
      trk3 -> Reset();
      trk4 -> Reset();
     
  	 	// Restore state from original data 
  		trk1 -> Add(ntrk[0], 1);
  		trk2 -> Add(ntrk[1], 1);
  		trk3 -> Add(ntrk[2], 1);
  		trk4 -> Add(ntrk[3], 1);

			// Now we run these through the convolution 
			GxTrk1 = ConvolveFit(trk1, GxTrk1, Params, offset, iter); 
			GxTrk2 = ConvolveFit(trk2, GxTrk2, Params, offset, iter); 
			GxTrk3 = ConvolveFit(trk3, GxTrk3, Params, offset, iter); 
			GxTrk4 = ConvolveFit(trk4, GxTrk4, Params, offset, iter); 
	
			for (int c(0); c < PDFs.size(); c++)
			{
				delete PDFs[c];
			}	
		
		}

		//SafeScale(GxTrk1, trk1);
		//SafeScale(GxTrk2, trk2);
		//SafeScale(GxTrk3, trk3);
		//SafeScale(GxTrk4, trk4);


    // Reset the state of the data histograms 
    trk1 -> Reset();
    trk2 -> Reset();
    trk3 -> Reset();
    trk4 -> Reset();
   
	 	// Restore state from original data 
		trk1 -> Add(ntrk[0], 1);
		trk2 -> Add(ntrk[1], 1);
		trk3 -> Add(ntrk[2], 1);
		trk4 -> Add(ntrk[3], 1);

    // Testing if temperature helps....
    float heat = 1; //float(i)/float(cor_loop);
    
    // Do the subtraction  
    trk1 -> Add(GxTrk1[1], -heat); 
    trk1 -> Add(GxTrk1[2], -heat); 
    trk1 -> Add(GxTrk1[3], -heat); 
    trk1 -> Add(GxTrk1[4], -heat);
    trk1 -> Add(GxTrk1[5], -heat);

    trk2 -> Add(GxTrk2[0], -heat);
    trk2 -> Add(GxTrk2[2], -heat);
    trk2 -> Add(GxTrk2[3], -heat);
    trk2 -> Add(GxTrk2[4], -heat);
    trk2 -> Add(GxTrk2[5], -heat);

    trk3 -> Add(GxTrk3[0], -heat);
    trk3 -> Add(GxTrk3[1], -heat);
    trk3 -> Add(GxTrk3[3], -heat);
    trk3 -> Add(GxTrk3[4], -heat);
    trk3 -> Add(GxTrk3[5], -heat);

    trk4 -> Add(GxTrk4[0], -heat);
    trk4 -> Add(GxTrk4[1], -heat);
    trk4 -> Add(GxTrk4[2], -heat);
    trk4 -> Add(GxTrk4[4], -heat);
    trk4 -> Add(GxTrk4[5], -heat);

    std::cout << "################### " << i << std::endl;
   
    // ==== Trk1 
    std::vector<TString> Names1 = {"trk1_F1", "trk2_F1", "trk3_F1", "trk4_F1", "trk5_F1", "trk6_F1"}; 
    std::vector<TH1F*> Trk1_PDFs = B.MakeTH1F(Names1, ntrk[0]); 
    B.ShiftExpandTH1F(GxTrk1, Trk1_PDFs);
    
    // ==== Trk2
    std::vector<TString> Names2 = {"trk1_F2", "trk2_F2", "trk3_F2", "trk4_F2", "trk5_F2", "trk6_F2"}; 
    std::vector<TH1F*> Trk2_PDFs = B.MakeTH1F(Names2, ntrk[1]); 
    B.ShiftExpandTH1F(GxTrk2, Trk2_PDFs);

    // ==== Trk3 
    std::vector<TString> Names3 = {"trk1_F3", "trk2_F3", "trk3_F3", "trk4_F3", "trk5_F3", "trk6_F3"}; 
    std::vector<TH1F*> Trk3_PDFs = B.MakeTH1F(Names3, ntrk[2]); 
    B.ShiftExpandTH1F(GxTrk3, Trk3_PDFs);

    // ==== Trk4 
    std::vector<TString> Names4 = {"trk1_F4", "trk2_F4", "trk3_F4", "trk4_F4", "trk5_F4", "trk6_F4"}; 
    std::vector<TH1F*> Trk4_PDFs = B.MakeTH1F(Names4, ntrk[3]); 
    B.ShiftExpandTH1F(GxTrk4, Trk4_PDFs);

    // Create the variabls needed for the Plotting below  
    std::vector<TString> Names_Sub = {"trk1_Sub", "trk2_Sub", "trk3_Sub", "trk4_Sub"}; 
    std::vector<TH1F*> Tracks = B.MakeTH1F(Names_Sub, ntrk[2]); 
    B.ShiftExpandTH1F({trk1, trk2, trk3, trk4}, Tracks);
    std::vector<TH1F*> Truth = {Closure[0][0], Closure[1][1], Closure[2][2], Closure[3][3]};

    can_HD -> Clear();
    TH1F* Temp = (TH1F*)Tracks[1] -> Clone("Lol"); 
    can_HD -> SetWindowSize(1200, 1200); 
    can_HD -> Divide(2,2);
    P.PlotHists({Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs}, Closure, Tracks, can_HD);
		can_HD -> Print("Out.pdf");

    // ======== Section for the output ========== //
    // === Save the PDFs and the subtracted Hists

    // === Save the subtracted Hists
    std::vector<TString> Names_Pure = {"trk1_Pure", "trk2_Pure", "trk3_Pure", "trk4_Pure"};    
    std::vector<TH1F*> Track_Out = B.MakeTH1F(Names_Pure, ntrk[2]); 
    B.ShiftExpandTH1F({trk1, trk2, trk3, trk4}, Track_Out);

    std::vector<TH1F*> Pred_PDF;
    std::vector<std::vector<TH1F*>> Set = {Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs, Track_Out, ntrk, Closure[0], Closure[1], Closure[2], Closure[3]}; 


    for (std::vector<TH1F*> PDF_S : Set)
    {
      for (int g(0); g < PDF_S.size(); g++)
      {
        TString name_pdf = PDF_S[g] -> GetName(); name_pdf += (i); 
        Pred_PDF.push_back((TH1F*)PDF_S[g] -> Clone(name_pdf));
      }
    }

    
    // === Save the FLost progress
    float FLost_Pred = B.FLost(ntrk, {Trk1_PDFs, Trk2_PDFs, Trk3_PDFs, Trk4_PDFs}); 
    float FLost_MC = B.FLost({ntrk[0], ntrk[1], ntrk[2], ntrk[3]}, Closure);   
    FLOST_Prediction -> SetBinContent(i+1, FLost_Pred); 
    FLOST_Truth -> SetBinContent(i+1, FLost_MC);   
    
    TString name_out =  "FLost_P.at."; name_out += (i); 
    TH1F* FL_P_Copy = (TH1F*)FLOST_Prediction -> Clone(name_out);

    name_out =  "FLost_T.at."; name_out += (i); 
    TH1F* FL_T_Copy = (TH1F*)FLOST_Truth -> Clone(name_out);
    Pred_PDF.push_back(FL_T_Copy); 

    Output[i] = std::make_pair(FL_P_Copy, Pred_PDF);
 
    // Clean up memory  
    for (int v(0); v < Trk1_PDFs.size(); v++)
    {
      delete Trk1_PDFs[v];
      delete Trk2_PDFs[v];
      delete Trk3_PDFs[v];
      delete Trk4_PDFs[v];   
    }



    delete Temp; 

	}
  return Output;
}


