#include "Evaluate.h"

void ShapePerformance(TString direct, TString Mode)
{
  MMVTH1F F = ReadAlgorithmResults(direct);
  
  // Layer - Algo - Energy
  MMMF Trk1_ShapePerformance; 
  MMMF Trk2_ShapePerformance; 
  MMMF Trk3_ShapePerformance; 
  MMMF Trk4_ShapePerformance; 

  int start; 
  int end; 
  if (Mode == "Test")
  {
    start = 0; 
    end = 4; 
  }
  if (Mode == "Template")
  {
    start = 4; 
    end = 8; 
  }

  TString PDF = Mode + "Fits.pdf"; 

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  can -> Print(PDF + "["); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(4);
  can -> SetTopMargin(0.1); 
  
  for (MMVTH1Fi x = F.begin(); x != F.end(); x++)
  {
   
    TString L_JE = x -> first; 
    TString Layer = "";  
    if (L_JE.Contains("IBL")){Layer = "IBL";}
    if (L_JE.Contains("Blayer")){Layer = "Blayer";}
    if (L_JE.Contains("layer1")){Layer = "layer1";}
    if (L_JE.Contains("layer2")){Layer = "layer2";}
    L_JE = L_JE.ReplaceAll(Layer + "_", ""); 
    
    if (!L_JE.Contains("GeV") || Layer == ""){continue;}


    MVTH1F M = F[x -> first];
    
    //Truth Histograms 
    VTH1F ntrk1_T = M["ntrk_1_Truth"];
    VTH1F ntrk2_T = M["ntrk_2_Truth"];
    VTH1F ntrk3_T = M["ntrk_3_Truth"];
    VTH1F ntrk4_T = M["ntrk_4_Truth"];

    if (ntrk1_T.size() == 0){continue;}
    int a = 0; 

    for (int k(0); k < Algos.size(); k++)
    {
      VTH1F np = GetSubVector(M[Algos[k]+"_ntrk_1"], start, end); 
      if (np.size() == 0){ continue;}
      a = k;
      break; 
    }

    std::cout << "-> " << L_JE << " -> " << Layer << " -> " << Algos[a] << std::endl;
    // Normalization Histograms 
    VTH1F ntrk1_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
    float trk1_Normal_Test = WeightedComparisonToTruth(ntrk1_Normal_Test, ntrk1_T); 
    float trk2_Normal_Test = WeightedComparisonToTruth(ntrk2_Normal_Test, ntrk2_T); 
    float trk3_Normal_Test = WeightedComparisonToTruth(ntrk3_Normal_Test, ntrk3_T); 
    float trk4_Normal_Test = WeightedComparisonToTruth(ntrk4_Normal_Test, ntrk4_T); 
    Trk1_ShapePerformance[Layer]["Normalization"][L_JE] = trk1_Normal_Test;
    Trk2_ShapePerformance[Layer]["Normalization"][L_JE] = trk2_Normal_Test;
    Trk3_ShapePerformance[Layer]["Normalization"][L_JE] = trk3_Normal_Test;
    Trk4_ShapePerformance[Layer]["Normalization"][L_JE] = trk4_Normal_Test;
    
    // ShiftNormal Histograms 
    a++; 
    VTH1F ntrk1_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
    float trk1_ShiftNormal_Test = WeightedComparisonToTruth(ntrk1_ShiftNormal_Test, ntrk1_T); 
    float trk2_ShiftNormal_Test = WeightedComparisonToTruth(ntrk2_ShiftNormal_Test, ntrk2_T); 
    float trk3_ShiftNormal_Test = WeightedComparisonToTruth(ntrk3_ShiftNormal_Test, ntrk3_T); 
    float trk4_ShiftNormal_Test = WeightedComparisonToTruth(ntrk4_ShiftNormal_Test, ntrk4_T); 
    Trk1_ShapePerformance[Layer]["ShiftNormal"][L_JE] = trk1_ShiftNormal_Test;
    Trk2_ShapePerformance[Layer]["ShiftNormal"][L_JE] = trk2_ShiftNormal_Test;
    Trk3_ShapePerformance[Layer]["ShiftNormal"][L_JE] = trk3_ShiftNormal_Test;
    Trk4_ShapePerformance[Layer]["ShiftNormal"][L_JE] = trk4_ShiftNormal_Test;

    // ShiftNormalFFT Histograms
    a++;
    VTH1F ntrk1_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
    float trk1_ShiftNormalFFT_Test = WeightedComparisonToTruth(ntrk1_ShiftNormalFFT_Test, ntrk1_T); 
    float trk2_ShiftNormalFFT_Test = WeightedComparisonToTruth(ntrk2_ShiftNormalFFT_Test, ntrk2_T); 
    float trk3_ShiftNormalFFT_Test = WeightedComparisonToTruth(ntrk3_ShiftNormalFFT_Test, ntrk3_T); 
    float trk4_ShiftNormalFFT_Test = WeightedComparisonToTruth(ntrk4_ShiftNormalFFT_Test, ntrk4_T); 
    Trk1_ShapePerformance[Layer]["ShiftNormalFFT"][L_JE] = trk1_ShiftNormalFFT_Test;
    Trk2_ShapePerformance[Layer]["ShiftNormalFFT"][L_JE] = trk2_ShiftNormalFFT_Test;
    Trk3_ShapePerformance[Layer]["ShiftNormalFFT"][L_JE] = trk3_ShiftNormalFFT_Test;
    Trk4_ShapePerformance[Layer]["ShiftNormalFFT"][L_JE] = trk4_ShiftNormalFFT_Test;

    // ShiftNormalWidthFFT Histograms
    a++;
    VTH1F ntrk1_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
    float trk1_ShiftNormalWidthFFT_Test = WeightedComparisonToTruth(ntrk1_ShiftNormalWidthFFT_Test, ntrk1_T); 
    float trk2_ShiftNormalWidthFFT_Test = WeightedComparisonToTruth(ntrk2_ShiftNormalWidthFFT_Test, ntrk2_T); 
    float trk3_ShiftNormalWidthFFT_Test = WeightedComparisonToTruth(ntrk3_ShiftNormalWidthFFT_Test, ntrk3_T); 
    float trk4_ShiftNormalWidthFFT_Test = WeightedComparisonToTruth(ntrk4_ShiftNormalWidthFFT_Test, ntrk4_T); 
    Trk1_ShapePerformance[Layer]["ShiftNormalWidthFFT"][L_JE] = trk1_ShiftNormalWidthFFT_Test;
    Trk2_ShapePerformance[Layer]["ShiftNormalWidthFFT"][L_JE] = trk2_ShiftNormalWidthFFT_Test;
    Trk3_ShapePerformance[Layer]["ShiftNormalWidthFFT"][L_JE] = trk3_ShiftNormalWidthFFT_Test;
    Trk4_ShapePerformance[Layer]["ShiftNormalWidthFFT"][L_JE] = trk4_ShiftNormalWidthFFT_Test;

    // ShiftNormalWidthFFT Histograms
    a++;
    VTH1F ntrk1_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
    float trk1_Incremental_Test = WeightedComparisonToTruth(ntrk1_Incremental_Test, ntrk1_T); 
    float trk2_Incremental_Test = WeightedComparisonToTruth(ntrk2_Incremental_Test, ntrk2_T); 
    float trk3_Incremental_Test = WeightedComparisonToTruth(ntrk3_Incremental_Test, ntrk3_T); 
    float trk4_Incremental_Test = WeightedComparisonToTruth(ntrk4_Incremental_Test, ntrk4_T); 
    Trk1_ShapePerformance[Layer]["Incremental"][L_JE] = trk1_Incremental_Test;
    Trk2_ShapePerformance[Layer]["Incremental"][L_JE] = trk2_Incremental_Test;
    Trk3_ShapePerformance[Layer]["Incremental"][L_JE] = trk3_Incremental_Test;
    Trk4_ShapePerformance[Layer]["Incremental"][L_JE] = trk4_Incremental_Test;
    
    auto Printer =[&](TCanvas* can, TString Fitter, TString dir, VTH1F trk, VTH1F trk_T, TString Title, TString PDF)
    {
      PlotHists(trk, trk_T, Title, can); 
      can -> Print(dir);
      can -> Print(PDF);
      can -> Clear(); 
    }; 

    // Normal Fits
    TString Fitter = "Normal";
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_Normal_Test, ntrk1_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_Normal_Test, ntrk2_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_Normal_Test, ntrk3_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_Normal_Test, ntrk4_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    // Normal Fits
    Fitter = "ShiftNormal";  
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_ShiftNormal_Test, ntrk1_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_ShiftNormal_Test, ntrk2_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_ShiftNormal_Test, ntrk3_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_ShiftNormal_Test, ntrk4_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    // Normal Fits
    Fitter = "ShiftNormalFFT"; 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_ShiftNormalFFT_Test, ntrk1_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_ShiftNormalFFT_Test, ntrk2_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_ShiftNormalFFT_Test, ntrk3_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_ShiftNormalFFT_Test, ntrk4_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    // Normal Fits
    Fitter = "ShiftNormalWidthFFT"; 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_ShiftNormalWidthFFT_Test, ntrk1_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_ShiftNormalWidthFFT_Test, ntrk2_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_ShiftNormalWidthFFT_Test, ntrk3_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_ShiftNormalWidthFFT_Test, ntrk4_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    // Incremental 
    Fitter = "Incremental"; 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_Incremental_Test, ntrk1_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_Incremental_Test, ntrk2_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_Incremental_Test, ntrk3_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_Incremental_Test, ntrk4_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

  }
  can -> Print(PDF + "]"); 
  delete can; 
  
  TString PER = Mode + "-Performance.pdf"; 
  TString dir = "./ShapePerformance/";
  float min = 0.00001; 
  float max = 1;
  TCanvas* Lan = new TCanvas(); 
  Lan -> Print(PER + "["); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(4);
  
  std::vector<TGraph*> gr; 
  gr = GenerateMultiGraph(Trk1_ShapePerformance["IBL"], "1-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_IBL_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk1_ShapePerformance["Blayer"], "1-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk1_ShapePerformance["layer1"], "1-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk1_ShapePerformance["layer2"], "1-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 


  GenerateMultiGraph(Trk2_ShapePerformance["IBL"], "2-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_IBL_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 

  GenerateMultiGraph(Trk2_ShapePerformance["Blayer"], "2-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk2_ShapePerformance["layer1"], "2-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk2_ShapePerformance["layer2"], "2-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 


  GenerateMultiGraph(Trk3_ShapePerformance["IBL"], "3-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_IBL_Shape_ERROR.png");
  Lan -> Print(PER);
  Lan -> Clear(); 

  GenerateMultiGraph(Trk3_ShapePerformance["Blayer"], "3-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk3_ShapePerformance["layer1"], "3-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk3_ShapePerformance["layer2"], "3-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 


  GenerateMultiGraph(Trk4_ShapePerformance["IBL"], "4-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_IBL_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 

  GenerateMultiGraph(Trk4_ShapePerformance["Blayer"], "4-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk4_ShapePerformance["layer1"], "4-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  GenerateMultiGraph(Trk4_ShapePerformance["layer2"], "4-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 

  Lan -> Print(PER + "]");
  delete Lan;
}

MMVTH1F FLostPerformance(TString direct, TString Mode, bool Plotting)
{
  MMVTH1F F = ReadAlgorithmResults(direct);
  
  // Layer - Algo - Energy
  MMMF FLost2_Performance; 
  MMMF FLost3_Performance; 


  int start; 
  int end; 
  if (Mode == "Test")
  {
    start = 0; 
    end = 4; 
  }
  if (Mode == "Template")
  {
    start = 4; 
    end = 8; 
  }

  for (MMVTH1Fi x = F.begin(); x != F.end(); x++)
  {
   
    TString L_JE = x -> first; 
    TString Layer = "";  
    if (L_JE.Contains("IBL")){Layer = "IBL";}
    if (L_JE.Contains("Blayer")){Layer = "Blayer";}
    if (L_JE.Contains("layer1")){Layer = "layer1";}
    if (L_JE.Contains("layer2")){Layer = "layer2";}
    L_JE = L_JE.ReplaceAll(Layer + "_", ""); 
    
    if (!L_JE.Contains("GeV") || Layer == ""){continue;}

    MVTH1F M = F[x -> first];
    
    //Truth Histograms 
    VTH1F ntrk1_T = M["ntrk_1_Truth"];
    VTH1F ntrk2_T = M["ntrk_2_Truth"];
    VTH1F ntrk3_T = M["ntrk_3_Truth"];
    VTH1F ntrk4_T = M["ntrk_4_Truth"];

    if (ntrk1_T.size() == 0){continue;}
    int a = 0; 

    for (int k(0); k < Algos.size(); k++)
    {
      VTH1F np = GetSubVector(M[Algos[k]+"_ntrk_1"], start, end); 
      if (np.size() == 0){ continue;}
      a = k;
      break; 
    }

    // Normalization Histograms 
    VTH1F ntrk1_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_Normal_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 

    // ShiftNormal Histograms 
    a++; 
    VTH1F ntrk1_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_ShiftNormal_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 

    // ShiftNormalFFT Histograms
    a++;
    VTH1F ntrk1_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_ShiftNormalFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 

    // ShiftNormalWidthFFT Histograms
    a++;
    VTH1F ntrk1_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_ShiftNormalWidthFFT_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 

    // ShiftNormalWidthFFT Histograms
    a++;
    VTH1F ntrk1_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_Incremental_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
 
    
    // FLost2
    float FL2_Normal = Flost2({ntrk1_Normal_Test, ntrk2_Normal_Test}); 
    float FL3_Normal = Flost3({ntrk1_Normal_Test, ntrk2_Normal_Test, ntrk3_Normal_Test, ntrk4_Normal_Test}); 

    float FL2_ShiftNormal = Flost2({ntrk1_ShiftNormal_Test, ntrk2_ShiftNormal_Test}); 
    float FL3_ShiftNormal = Flost3({ntrk1_ShiftNormal_Test, ntrk2_ShiftNormal_Test, ntrk3_ShiftNormal_Test, ntrk4_ShiftNormal_Test}); 

    float FL2_ShiftNormalFFT = Flost2({ntrk1_ShiftNormalFFT_Test, ntrk2_ShiftNormalFFT_Test}); 
    float FL3_ShiftNormalFFT = Flost3({ntrk1_ShiftNormalFFT_Test, ntrk2_ShiftNormalFFT_Test, ntrk3_ShiftNormalFFT_Test, ntrk4_ShiftNormalFFT_Test}); 

    float FL2_ShiftNormalWidthFFT = Flost2({ntrk1_ShiftNormalWidthFFT_Test, ntrk2_ShiftNormalWidthFFT_Test}); 
    float FL3_ShiftNormalWidthFFT = Flost3({ntrk1_ShiftNormalWidthFFT_Test, ntrk2_ShiftNormalWidthFFT_Test, ntrk3_ShiftNormalWidthFFT_Test, ntrk4_ShiftNormalWidthFFT_Test}); 

    float FL2_Incremental = Flost2({ntrk1_Incremental_Test, ntrk2_Incremental_Test}); 
    float FL3_Incremental = Flost3({ntrk1_Incremental_Test, ntrk2_Incremental_Test, ntrk3_Incremental_Test, ntrk4_Incremental_Test}); 

    float FL2_Truth = Flost2({ntrk1_T, ntrk2_T}); 
    float FL3_Truth = Flost3({ntrk1_T, ntrk2_T, ntrk3_T, ntrk4_T}); 


    FLost2_Performance[Layer]["Normalization"][L_JE] = FL2_Normal; 
    FLost3_Performance[Layer]["Normalization"][L_JE] = FL3_Normal; 

    FLost2_Performance[Layer]["ShiftNormal"][L_JE] = FL2_ShiftNormal; 
    FLost3_Performance[Layer]["ShiftNormal"][L_JE] = FL3_ShiftNormal; 

    FLost2_Performance[Layer]["ShiftNormalFFT"][L_JE] = FL2_ShiftNormalFFT; 
    FLost3_Performance[Layer]["ShiftNormalFFT"][L_JE] = FL3_ShiftNormalFFT; 

    FLost2_Performance[Layer]["ShiftNormalWidthFFT"][L_JE] = FL2_ShiftNormalWidthFFT; 
    FLost3_Performance[Layer]["ShiftNormalWidthFFT"][L_JE] = FL3_ShiftNormalWidthFFT; 

    FLost2_Performance[Layer]["Incremental"][L_JE] = FL2_Incremental; 
    FLost3_Performance[Layer]["Incremental"][L_JE] = FL3_Incremental; 

    FLost2_Performance[Layer]["Truth"][L_JE] = FL2_Truth; 
    FLost3_Performance[Layer]["Truth"][L_JE] = FL3_Truth; 
  }

  TString dir = "./FLost/";
  
  MMVTH1F Output; 
  Output["IBL_FLost2"] = PrintFlost(FLost2_Performance, "IBL", dir + Mode + "_IBL_FLost2.txt"); 
  Output["Blayer_FLost2"] = PrintFlost(FLost2_Performance, "Blayer", dir + Mode + "_Blayer_FLost2.txt"); 
  Output["layer1_FLost2"] = PrintFlost(FLost2_Performance, "layer1", dir + Mode + "_layer1_FLost2.txt"); 
  Output["layer2_FLost2"] = PrintFlost(FLost2_Performance, "layer2", dir + Mode + "_layer2_FLost2.txt"); 

  Output["IBL_FLost3"] = PrintFlost(FLost3_Performance, "IBL", dir + Mode + "_IBL_FLost3.txt"); 
  Output["Blayer_FLost3"] = PrintFlost(FLost3_Performance, "Blayer", dir + Mode + "_Blayer_FLost3.txt"); 
  Output["layer1_FLost3"] = PrintFlost(FLost3_Performance, "layer1", dir + Mode + "_layer1_FLost3.txt"); 
  Output["layer2_FLost3"] = PrintFlost(FLost3_Performance, "layer2", dir + Mode + "_layer2_FLost3.txt"); 

  if (Plotting){ return Output; }

  TString PDF = Mode + "_FLost.pdf"; 
  TCanvas* Lan = new TCanvas(); 
  //Lan -> SetLogy(); 
  Lan -> Print(dir + PDF + "["); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(2);
  Lan -> SetRightMargin(0.01); 
  Lan -> SetTopMargin(0.1); 

  float min = 0.001; 
  float max = 2; 
  std::vector<TGraph*> gr; 
  gr = GenerateMultiGraph(FLost2_Performance["IBL"], "FLost2 Performance in the IBL", Lan, min, max, "FLost2"); 
  Lan -> Print(dir + Mode + "_IBL_FLOST2.png"); 
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  gr = GenerateMultiGraph(FLost2_Performance["Blayer"], "FLost2 Performance in the Blayer", Lan, min, max, "FLost2"); 
  Lan -> Print(dir + Mode + "_Blayer_FLOST2.png");
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  gr = GenerateMultiGraph(FLost2_Performance["layer1"], "FLost2 Performance in the layer1", Lan, min, max, "FLost2"); 
  Lan -> Print(dir + Mode + "_layer1_FLOST2.png");
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  gr = GenerateMultiGraph(FLost2_Performance["layer2"], "FLost2 Performance in the layer2", Lan, min, max, "FLost2"); 
  Lan -> Print(dir + Mode + "_layer2_FLOST2.png");
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 

  gr = GenerateMultiGraph(FLost3_Performance["IBL"], "FLost3 Performance in the IBL", Lan, min, max, "FLost3"); 
  Lan -> Print(dir + Mode + "_IBL_FLOST3.png");
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  gr = GenerateMultiGraph(FLost3_Performance["Blayer"], "FLost3 Performance in the Blayer", Lan, min, max, "FLost3"); 
  Lan -> Print(dir + Mode + "_Blayer_FLOST3.png"); 
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  gr = GenerateMultiGraph(FLost3_Performance["layer1"], "FLost3 Performance in the layer1", Lan, min, max, "FLost3"); 
  Lan -> Print(dir + Mode + "_layer1_FLOST3.png"); 
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  gr = GenerateMultiGraph(FLost3_Performance["layer2"], "FLost3 Performance in the layer2", Lan, min, max, "FLost3"); 
  Lan -> Print(dir + Mode + "_layer2_FLOST3.png");
  Lan -> Print(dir + PDF);
  Lan -> Clear(); 
  BulkDelete(gr); 
  
  Lan -> Print(dir + PDF + "]");
  delete Lan;
  
  return Output; 
}

void WriteFLostFile(TString Input, TString Output)
{ 

  MMVTH1F Out = FLostPerformance(Input, "Template"); 
  
  TFile* F = new TFile(Output + ".root", "RECREATE"); 
  for (MMVTH1Fi dir = Out.begin(); dir != Out.end(); dir++)
  {
    TString key = (dir -> first); 
    TString L; 
    TString FL; 
    
    MVTH1F Hist = dir -> second; 

    if (key.Contains("IBL")){ L = "IBL"; }
    if (key.Contains("Blayer")){ L = "Blayer"; }
    if (key.Contains("layer1")){ L = "layer1"; }
    if (key.Contains("layer2")){ L = "layer2"; }
    if (key.Contains("FLost2")){ FL = "FLost2"; }
    if (key.Contains("FLost3")){ FL = "FLost3"; }

    gDirectory -> mkdir(Output); 
    gDirectory -> mkdir(Output + "/" + L); 
    gDirectory -> mkdir(Output + "/" + L + "/Error_" + FL); 
    gDirectory -> mkdir(Output + "/" + L + "/Original_" + FL); 
    

    std::vector<TH1F*> Err = Hist["Error_" + FL]; 
    F -> cd(Output + "/" + L + "/Error_" + FL); 
    BulkWrite(Err); 
  
    std::vector<TH1F*> Ori = Hist["Original_" + FL]; 
    F -> cd(Output + "/" + L + "/Original_" + FL); 
    BulkWrite(Ori); 

    F -> cd("/"); 
  }
  F -> Close(); 
}

void ErrorFlostAlgo(TString dir)
{
  auto PlotMinimizer =[&] (MMMMF FL_Err, TString FL, TString pdf, TCanvas* can)
  {
    for (MMMMFi fl = FL_Err.begin(); fl != FL_Err.end(); fl++)
    {
      TString L = fl -> first; 
      for (MMMFi a_k = FL_Err[L].begin(); a_k != FL_Err[L].end(); a_k++)
      {
        MMF Mini = FL_Err[L][a_k -> first]; 
        std::vector<TGraph*> mini = GenerateMultiGraph(Mini,"Minimizer Comparison: " + L + " " + (a_k -> first), can, 0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
        can -> Print(pdf); 
        can -> Print("Minimizer/FLost" + FL + "_" + L + "_" + (a_k -> first) + ".png"); 
        
        for (TGraph* G : mini){ delete G; }
      }
    }
  }; 

  auto GetBestMinimizer =[&] (MMF FL_Err)
  {
    
    MF Best_Minimizer; 
    bool Pop = true; 
    for (TString JE : JetEnergy)
    {
      float err = -1; 
      TString Mini = "x"; 
      for (MMFi fl = FL_Err.begin(); fl != FL_Err.end(); fl++)
      {    
        float FL = FL_Err[fl -> first][JE]; 
        if (Pop == true){ Best_Minimizer[fl -> first] = 0; }
        if (float(FL) == float(0.)){continue;}
        if ( err < 0 || err > FL){ err = FL; Mini = fl -> first; }
        
      }
      Pop = false;
      Best_Minimizer[Mini]+=1;
    }
    return Best_Minimizer; 
  }; 

  auto CompileTable =[&] (MMF Minimizer, std::vector<TString> Ag, TString FL)
  {
    std::vector<TString> Out;  
    TString Title = FL; 
    int Margin = 24; 
    CoutText(&Title, 30 - Title.Sizeof(), " "); 
    Title += (" | "); 
    
    int l; 
    std::vector<TString> Mini = {};
    for (TString A : Ag)
    {
      for (TString L : Layer)
      {
        for (MFi mini = Minimizer[L + "+" + A].begin(); mini != Minimizer[L + "+" + A].end(); mini++)
        {
          TString Combi = mini -> first; 
          Title += Combi; 
          //CoutText(&Title, Margin - Combi.Sizeof(), " "); 
          Title += (" | ");

          Mini.push_back(Combi);
        }
        Out.push_back(Title); 
        l = Title.Sizeof(); 
        Title = ""; 
        CoutText(&Title, l, "="); 
        Out.push_back(Title); 
        break;
      }
      break;
    }
    
    MMF Best; 
    for (TString L : Layer)
    {
      for (TString A : Ag)
      {
        Title = L + "+" + A; 
        CoutText(&Title, 30 - Title.Sizeof(), " "); 
        Title += (" | "); 
        
        for (TString m : Mini)
        {
          TString f = ""; f += Minimizer[L + "+" + A][m]; 
          Title += f; 
          CoutText(&Title, m.Sizeof() - f.Sizeof(), " ");
          Title += (" | ");

          Best[L][m] += Minimizer[L + "+" + A][m]; 
        }
        Out.push_back(Title); 
      }

      Title = ""; 
      CoutText(&Title, l, "_"); 
      Out.push_back(Title);  
    
      Title = ""; 
      Title += ("Score"); 
      CoutText(&Title, 30 - Title.Sizeof(), " "); 
      Title += (" | "); 
      for (TString m : Mini)
      {
        TString Comb = ""; Comb += (Best[L][m]); 
        Title += Comb; 
        CoutText(&Title, m.Sizeof() - Comb.Sizeof(), " "); 
        Title += (" | "); 
      }
      Out.push_back(Title);  
      Out.push_back("");
    }
    MT Winner; 
    for (TString L : Layer)
    {
      int sc = 0; 
      TString B = ""; 
      for (TString M : Mini)
      {
        if (M == "x"){continue;}
        if (sc < Best[L][M]){ sc = Best[L][M]; B = M; }
      }
      Winner[L] = B;
    }
    std::ofstream print; 
    print.open(FL +  "_ScoreMatrix.txt"); 
    for (TString a : Out)
    {
      print << a << "\n"; 
    }
    print.close(); 
    return Winner; 
  };

  MMMMMF Obj = ReadMinimizers(dir); 
  MMMMF FL2_Error = Obj["Error_FLost2"]; 
  MMMMF FL3_Error = Obj["Error_FLost3"];
  
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(2);

  TCanvas* can = new TCanvas(); 
  
  TString pdf = "FLost2_Layer_Algo.pdf"; 
  can -> Print(pdf + "["); 
  PlotMinimizer(FL2_Error, "2", pdf, can); 
  can -> Print(pdf + "]"); 

  pdf = "FLost3_Layer_Algo.pdf"; 
  can -> Print(pdf + "["); 
  PlotMinimizer(FL3_Error, "3", pdf, can); 
  can -> Print(pdf + "]"); 
 
  MMF MinimizerMap_FL2; 
  MMF MinimizerMap_FL3; 

  std::vector<TString> Ag = {"Normalization", "ShiftNormal", "ShiftNormalFFT", "ShiftNormalWidthFFT", "Incremental"};
  for (TString L : Layer)
  {
    for (TString A : Ag)
    {
      MinimizerMap_FL2[L + "+" + A] = GetBestMinimizer(FL2_Error[L][A]); 
      MinimizerMap_FL3[L + "+" + A] = GetBestMinimizer(FL3_Error[L][A]); 
    }
  }

  MT FL2_Best = CompileTable(MinimizerMap_FL2, Ag, "FLost2");  
  MT FL3_Best = CompileTable(MinimizerMap_FL3, Ag, "FLost3");  
  
  MMMF FL2_Best_Mini; 
  MMMF FL3_Best_Mini; 
  for (TString L : Layer)
  {
    TString Mini_2 = FL2_Best[L]; 
    TString Mini_3 = FL3_Best[L]; 
    
    for (TString A : Ag)
    {
      FL2_Best_Mini[L][A] = FL2_Error[L][A][Mini_2];
      FL3_Best_Mini[L][A] = FL3_Error[L][A][Mini_3]; 
    }
  }
  
  std::vector<TGraph*> GR; 
  pdf = "Optimizer_Algo.pdf"; 
  can -> Print(pdf + "["); 
  TString Lay = "IBL"; 
  TString FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "Blayer"; 
  FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer1"; 
  FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer2"; 
  FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 


  Lay = "IBL"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "Blayer"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer1"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer2"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print("FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  can -> Print(pdf + "]"); 














}
