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

  can -> Print(PDF + "["); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(3);
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
    
    bool Skip = true;
    for (TString k : JetEnergy){ if (k == L_JE){Skip = false; break;} }
    if (!L_JE.Contains("GeV") || Layer == "" || Skip == true){continue;}

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

    std::cout << "-> " << L_JE << " -> " << Layer << " -> " << Algos[a] << std::endl;
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

    std::cout << "-> " << L_JE << " -> " << Layer << " -> " << Algos[a] << std::endl;
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

    std::cout << "-> " << L_JE << " -> " << Layer << " -> " << Algos[a] << std::endl;
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

    std::cout << "-> " << L_JE << " -> " << Layer << " -> " << Algos[a] << std::endl;
    // Experimental Histograms
    a++;
    VTH1F ntrk1_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], 0, 4); 
    VTH1F ntrk2_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], 0, 4); 
    VTH1F ntrk3_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], 0, 4); 
    VTH1F ntrk4_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], 0, 4); 
    float trk1_Experimental_Test = WeightedComparisonToTruth(ntrk1_Experimental_Test, ntrk1_T); 
    float trk2_Experimental_Test = WeightedComparisonToTruth(ntrk2_Experimental_Test, ntrk2_T); 
    float trk3_Experimental_Test = WeightedComparisonToTruth(ntrk3_Experimental_Test, ntrk3_T); 
    float trk4_Experimental_Test = WeightedComparisonToTruth(ntrk4_Experimental_Test, ntrk4_T); 
    
    if (trk1_Experimental_Test != -1 && Mode != "Test")
    {
      std::cout << "-> " << L_JE << " -> " << Layer << " -> " << Algos[a] << std::endl;
      Trk1_ShapePerformance[Layer]["Experimental"][L_JE] = trk1_Experimental_Test;
      Trk2_ShapePerformance[Layer]["Experimental"][L_JE] = trk2_Experimental_Test;
      Trk3_ShapePerformance[Layer]["Experimental"][L_JE] = trk3_Experimental_Test;
      Trk4_ShapePerformance[Layer]["Experimental"][L_JE] = trk4_Experimental_Test;
    }


    auto Printer =[&](TCanvas* can, TString Fitter, TString dir, VTH1F trk, VTH1F trk_T, TString Title, TString PDF)
    {
      if (trk.size() == 0){return;}
      PlotHists(trk, trk_T, Title, can); 
      can -> Print(dir);
      can -> Print(PDF);
      can -> Clear(); 
    }; 

    TString Fitter = "Normal";
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_Normal_Test, ntrk1_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_Normal_Test, ntrk2_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_Normal_Test, ntrk3_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_Normal_Test, ntrk4_T, "Normalization Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    Fitter = "ShiftNormal";  
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_ShiftNormal_Test, ntrk1_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_ShiftNormal_Test, ntrk2_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_ShiftNormal_Test, ntrk3_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_ShiftNormal_Test, ntrk4_T, "Normalization Shift Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    Fitter = "ShiftNormalFFT"; 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_ShiftNormalFFT_Test, ntrk1_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_ShiftNormalFFT_Test, ntrk2_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_ShiftNormalFFT_Test, ntrk3_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_ShiftNormalFFT_Test, ntrk4_T, "Normalization Shift FFT Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    Fitter = "ShiftNormalWidthFFT"; 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_ShiftNormalWidthFFT_Test, ntrk1_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_ShiftNormalWidthFFT_Test, ntrk2_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_ShiftNormalWidthFFT_Test, ntrk3_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_ShiftNormalWidthFFT_Test, ntrk4_T, "Normalization Shift Width FFT Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    Fitter = "Incremental"; 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_Incremental_Test, ntrk1_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_Incremental_Test, ntrk2_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_Incremental_Test, ntrk3_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_Incremental_Test, ntrk4_T, "Normalization Shift Width FFT Incremental Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

    Fitter = "Experimental"; 
    if (trk1_Experimental_Test == -1 && Mode == "Test") {continue;}
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk1/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk1_Experimental_Test, ntrk1_T, "Experimental Fit: " + Layer + " " + L_JE + " -> Track-1", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk2/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk2_Experimental_Test, ntrk2_T, "Experimental Fit: " + Layer + " " + L_JE + " -> Track-2", PDF);
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk3/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk3_Experimental_Test, ntrk3_T, "Experimental Fit: " + Layer + " " + L_JE + " -> Track-3", PDF); 
    Printer(can, Fitter, "./Histograms/" + Mode + "/trk4/" + Fitter + "/" + Layer + "_" + L_JE + ".png", ntrk4_Experimental_Test, ntrk4_T, "Experimental Fit: " + Layer + " " + L_JE + " -> Track-4", PDF);

  }
  can -> Print(PDF + "]"); 
  delete can; 
  
  TString PER = Mode + "-Performance.pdf"; 
  TString dir = "./ShapePerformance/";
  float min = 0.0001; 
  float max = 10;
  TCanvas* Lan = new TCanvas(); 
  Lan -> Print(PER + "["); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(4);
  
  std::vector<TGraph*> gr; 
  gr = GenerateMultiGraph(Trk1_ShapePerformance["IBL"], "1-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_IBL_Shape_ERROR.png"); 
  Lan -> Print(PER);
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk1_ShapePerformance["Blayer"], "1-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk1_ShapePerformance["layer1"], "1-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk1_ShapePerformance["layer2"], "1-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk1/" + Mode + "-trk1_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
 
  gr = GenerateMultiGraph(Trk2_ShapePerformance["IBL"], "2-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_IBL_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 

  gr = GenerateMultiGraph(Trk2_ShapePerformance["Blayer"], "2-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk2_ShapePerformance["layer1"], "2-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk2_ShapePerformance["layer2"], "2-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk2/" + Mode + "-trk2_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 


  gr = GenerateMultiGraph(Trk3_ShapePerformance["IBL"], "3-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_IBL_Shape_ERROR.png");
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 

  gr = GenerateMultiGraph(Trk3_ShapePerformance["Blayer"], "3-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk3_ShapePerformance["layer1"], "3-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk3_ShapePerformance["layer2"], "3-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk3/" + Mode + "-trk3_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 

  gr = GenerateMultiGraph(Trk4_ShapePerformance["IBL"], "4-Track " + Mode + " Shape Performance For Different Algorithms: IBL", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_IBL_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 

  gr = GenerateMultiGraph(Trk4_ShapePerformance["Blayer"], "4-Track " + Mode + " Shape Performance For Different Algorithms: Blayer", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_Blayer_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk4_ShapePerformance["layer1"], "4-Track " + Mode + " Shape Performance For Different Algorithms: layer1", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_layer1_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
  Lan -> Clear(); 
  
  gr = GenerateMultiGraph(Trk4_ShapePerformance["layer2"], "4-Track " + Mode + " Shape Performance For Different Algorithms: layer2", Lan, min, max, "Error"); 
  Lan -> Print(dir + "trk4/" + Mode + "-trk4_layer2_Shape_ERROR.png"); 
  Lan -> Print(PER);
  BulkDelete(gr); 
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
   
    bool Skip = true;
    for (TString k : JetEnergy){ if (k == L_JE){Skip = false; break;} }
    if (!L_JE.Contains("GeV") || Layer == "" || Skip == true){continue;}

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
    
    // Experimental 
    a++;
    VTH1F ntrk1_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_1"], start, end); 
    VTH1F ntrk2_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_2"], start, end); 
    VTH1F ntrk3_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_3"], start, end); 
    VTH1F ntrk4_Experimental_Test = GetSubVector(M[Algos[a]+"_ntrk_4"], start, end); 
    
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

    float FL2_Experimental = Flost2({ntrk1_Experimental_Test, ntrk2_Experimental_Test}); 
    float FL3_Experimental = Flost3({ntrk1_Experimental_Test, ntrk2_Experimental_Test, ntrk3_Experimental_Test, ntrk4_Experimental_Test}); 

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
    
    if (FL2_Experimental != -1)
    {
      FLost2_Performance[Layer]["Experimental"][L_JE] = FL2_Experimental; 
      FLost3_Performance[Layer]["Experimental"][L_JE] = FL3_Experimental; 
    }

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
    print.open(FL+"_Error/"+FL +  "_ScoreMatrix.txt"); 
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
  
  TCanvas* can = new TCanvas(); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(2);


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
  TString out = "FLost"+FL+"_Error/";
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "Blayer"; 
  FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer1"; 
  FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer2"; 
  FL = "2"; 
  GR = GenerateMultiGraph(FL2_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL2_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 


  Lay = "IBL"; 
  FL = "3"; 
  out = "FLost"+FL+"_Error/";
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "Blayer"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer1"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  Lay = "layer2"; 
  FL = "3"; 
  GR = GenerateMultiGraph(FL3_Best_Mini[Lay], "FLost" + FL + " Error Using the Best Minimizer: " + FL3_Best[Lay] + " - " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
  can -> Print(out + "FLost" + FL + "_" + Lay + "_Best.png");
  can -> Print(pdf);
  can -> Clear(); 
  BulkDelete(GR); 

  can -> Print(pdf + "]"); 
}

void BestMinimizerAlgo(TString dir)
{
  auto BestMiniForAlgo =[&] (MMMF FL)
  {
    MT BestCombination; 
    for (TString a : Ag)
    {
      MI CounterMatrix; 
      for (TString j : JetEnergy)
      {
        float err = -1; 
        TString Best = "x"; 
        for (MMFi M = FL[a].begin(); M != FL[a].end(); M++)
        {
          float e = FL[a][M -> first][j]; 
          if (e == 0){ continue; }
          if (err < 0 || err >  e){ Best = M -> first; err = e; }
        }
        CounterMatrix[Best]++; 
      }
      
      int counter = -1; 
      TString Best = "x"; 
      for (MIi c = CounterMatrix.begin(); c != CounterMatrix.end(); c++)
      {
        if (counter < 0 || counter < c -> second){ Best = c -> first; counter = c -> second; } 
      }
      BestCombination[a] = Best; 
    }

    return BestCombination; 
  }; 

  auto CleaningWriting =[&](MMMMF FL_Error, TCanvas* can, TString PDF, TString FL, TString out, TString Lay)
  {
    MMF FLost; 
    for (TString a : Ag)
    {
      TString Min = BestMiniForAlgo(FL_Error[Lay])[a];  
      FLost[Min + "+" + a] = FL_Error[Lay][a][Min];   
    }
    std::vector<TGraph*> GR = GenerateMultiGraph(FLost, "FLost" + FL + " Error of Algorithm with optimal Minimizer: " + Lay, can,  0.01, 10000, "#frac{#Delta F_{L" + FL + "}}{F_{T}} #times 100 (%)");
    can -> Print(out + "FLost" + FL + "_" + Lay + "_AlgoWithBestMin.png");
    can -> Print(PDF);
    can -> Clear(); 
    BulkDelete(GR); 
  
    MI CounterMatrix; 
    for (TString je : JetEnergy)
    {
      float err = -1; 
      TString Optimal = "x"; 
      for (MMFi B = FLost.begin(); B != FLost.end(); B++)
      {
        float e = FLost[B-> first][je]; 
        if (e == 0){continue;} 
        if (err < 0 || err > e){ Optimal = B -> first; err = e; }
      }
      CounterMatrix[Optimal]++; 
    }

    int counter = -1; 
    TString Best = "x"; 
    for (MIi c = CounterMatrix.begin(); c != CounterMatrix.end(); c++)
    {
      if (counter < 0 || counter < c -> second){ Best = c -> first; counter = c -> second; } 
    }
    return Best; 
  }; 



  MMMMMF Obj = ReadMinimizers(dir); 
  MMMMF FL2_Error = Obj["Error_FLost2"]; 
  MMMMF FL3_Error = Obj["Error_FLost3"];
  TString out = "OptimalMiniWithAlgo/";  
  TString pdf = "FLost2_AlgoWithOptimialMini.pdf"; 
  
  TCanvas* can = new TCanvas();
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(2);

  can -> Print(pdf + "["); 
  std::vector<TString> Out2; 
  Out2.push_back(CleaningWriting(FL2_Error, can, pdf, "2", out, "IBL")); 
  Out2.push_back(CleaningWriting(FL2_Error, can, pdf, "2", out, "Blayer")); 
  Out2.push_back(CleaningWriting(FL2_Error, can, pdf, "2", out, "layer1")); 
  Out2.push_back(CleaningWriting(FL2_Error, can, pdf, "2", out, "layer2")); 

  std::vector<TString> Out3; 
  Out3.push_back(CleaningWriting(FL3_Error, can, pdf, "3", out, "IBL")); 
  Out3.push_back(CleaningWriting(FL3_Error, can, pdf, "3", out, "Blayer")); 
  Out3.push_back(CleaningWriting(FL3_Error, can, pdf, "3", out, "layer1")); 
  Out3.push_back(CleaningWriting(FL3_Error, can, pdf, "3", out, "layer2")); 
  can -> Print(pdf + "]"); 

  std::ofstream print; 
  print.open(out + "FLost2_Verdict.txt"); 
  for (int i(0); i < Layer.size(); i++)
  {
    print << Layer[i] << ": " << Out2[i] << "\n"; 
  }
  print.close(); 

  print.open(out + "FLost3_Verdict.txt"); 
  for (int i(0); i < Layer.size(); i++)
  {
    print << Layer[i] << ": " << Out3[i] << "\n"; 
  }
  print.close(); 



}


void Clusters(TString dir)
{
  MMMTF F = ReadCTIDE(dir);  
    
  // Layer - trk - Energy - Float
  MMMF ClusterEntries; 
  for (MMMTFi f = F.begin(); f != F.end(); f++)
  {
    TString L = f -> first; 
    for (MMTFi en = F[L].begin(); en != F[L].end(); en++)
    {
      TString E = en -> first; 
      MTF Jet = F[L][E]; 
      std::vector<TH1F*> trk = Jet["Inside"]; 

      for (TH1F* H : trk)
      {
        TString Title = H -> GetTitle(); 

        if (Title.Contains("ntrk_1")){ ClusterEntries[L]["Track-1"][E] = H -> GetEntries(); }
        if (Title.Contains("ntrk_2")){ ClusterEntries[L]["Track-2"][E] = H -> GetEntries(); }
        if (Title.Contains("ntrk_3")){ ClusterEntries[L]["Track-3"][E] = H -> GetEntries(); }
        if (Title.Contains("ntrk_4")){ ClusterEntries[L]["Track-4"][E] = H -> GetEntries(); }
      }
    }
  }

  MMF IBL = ClusterEntries["IBL"]; 
  MMF Blayer = ClusterEntries["Blayer"]; 
  MMF layer1 = ClusterEntries["layer1"]; 
  MMF layer2 = ClusterEntries["layer2"]; 
  
  TCanvas* can = new TCanvas(); 
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(2);
  can -> Print("ClusterEntries.pdf[");  
  
  std::vector<TGraph*> GR;
  GR = GenerateMultiGraph(IBL, "Cluster Entries for Recorded Tracks in Detector Layer: IBL", can, 1, 1e9, "Entries"); 
  can -> Print("ClusterEntries.pdf");  
  can -> Print("IBL.png"); 
  BulkDelete(GR); 
  
  GR = GenerateMultiGraph(Blayer, "Cluster Entries for Recorded Tracks in Detector Layer: Blayer", can, 1, 1e9, "Entries"); 
  can -> Print("ClusterEntries.pdf");  
  can -> Print("Blayer.png"); 
  BulkDelete(GR); 
 
  GR = GenerateMultiGraph(layer1, "Cluster Entries for Recorded Tracks in Detector Layer: layer1", can, 1, 1e9, "Entries"); 
  can -> Print("ClusterEntries.pdf");  
  can -> Print("layer1.png"); 
  BulkDelete(GR); 
 
  GR = GenerateMultiGraph(layer2, "Cluster Entries for Recorded Tracks in Detector Layer: layer2", can, 1, 1e9, "Entries"); 
  can -> Print("ClusterEntries.pdf");  
  can -> Print("layer2.png"); 
  BulkDelete(GR); 
 
  can -> Print("ClusterEntries.pdf]");
}

void TemplateVariants(TString dir)
{
  TString Title = "Nominal"; 
  if (dir.Contains("Subtract")){Title = "Subtract";}
  if (dir.Contains("Smooth")){Title = "Smoothed"; }

  MMVTH1F F = ReadAlgorithmResults(dir);
  TH1F* Alg = F["Blayer_800_1000_GeV"]["Normalization_ntrk_1"][0]; 

  TCanvas* can = new TCanvas();
  can -> SetTopMargin(0.01);
  can -> SetRightMargin(0.02); 
  can -> SetLeftMargin(0.075); 
  can -> SetLogy(); 
  gStyle -> SetOptStat(0);
  gStyle -> SetImageScaling(4); 
  Alg -> GetXaxis() -> SetTitle("dE/dx [MeV g^{-1} cm^{2}]");
  Alg -> GetYaxis() -> SetTitle("Clusters");
  Alg -> GetXaxis() -> SetRangeUser(0, 12); 
  GeneratePlot(Alg, Title, can, kRed, kSolid, "HIST", 1); 
  Alg -> SetLineWidth(5); 
  can -> Update();
  can -> Print(Title + ".png");
}

void Debug_Cases(TString dir)
{
  MMMMTH1F debugs = Debug_File_Steps(dir);
 
  // Illustrate the cases we are trying to test and evaluate the fit predictions for 
  TString Layer_JE = "Blayer_1000_1200_GeV"; 
  MMMTH1F Algo_Cases = debugs[Layer_JE]; 

  // Truth Distributions 
  TH1F* ntrk1_ntru1_20 = Algo_Cases["NormalCase1"]["step_20"]["ntrk1-tru1"];
  TH1F* ntrk1_ntru2_20 = Algo_Cases["NormalCase1"]["step_20"]["ntrk1-tru2"];

  TH1F* ntrk1_ntru1_100 = Algo_Cases["NormalCase1"]["step_100"]["ntrk1-tru1"];
  TH1F* ntrk1_ntru2_100 = Algo_Cases["NormalCase1"]["step_100"]["ntrk1-tru2"];

  TH1F* ntrk1_ntru1_200 = Algo_Cases["NormalCase1"]["step_200"]["ntrk1-tru1"];
  TH1F* ntrk1_ntru2_200 = Algo_Cases["NormalCase1"]["step_200"]["ntrk1-tru2"];


  auto Printer =[&](TCanvas* can, TString dir, VTH1F trk, VTH1F trk_T, TString Title, TString PDF)
  {
    if (trk.size() == 0){return;}
    PlotHists(trk, trk_T, Title, can); 
    can -> Print(dir);
    can -> Print(PDF);
    can -> Clear(); 
  }; 

  // Case fits
  // 1. Use the 1-Truth as a starting template and build the 2-Truth template
  // 2. Use the nominal 1-Truth distribution and use the 2-Truth 
  // 3. Use the nominal 1 and 2-Truth templates
  TCanvas* can = new TCanvas();
  can -> SetTopMargin(0.01);
  can -> SetRightMargin(0.02); 
  //can -> SetLeftMargin(0.075); 
  can -> SetLogy(); 
  gStyle -> SetOptStat(0);
  gStyle -> SetImageScaling(4);
  TString pdf_n = "./IllustrationOfCases/Case"; 

  can -> Print(pdf_n + ".pdf["); 
  TH1F* NormalCase1_20_ntrk1_ntru1_Templ = Algo_Cases["NormalCase1"]["step_20"]["dEdx_ntrk_1_ntru_1_NormalCase11_step_20"]; 
  TH1F* NormalCase1_20_ntrk1_ntru2_Templ = Algo_Cases["NormalCase1"]["step_20"]["dEdx_ntrk_1_ntru_2_NormalCase12_step_20"]; 
  Printer(can, pdf_n + "1/NormalCase1.png", {NormalCase1_20_ntrk1_ntru1_Templ, NormalCase1_20_ntrk1_ntru2_Templ}, {ntrk1_ntru1_20, ntrk1_ntru2_20}, "Case-1: Normal", pdf_n + ".pdf");

  TH1F* NormalCase2_100_ntrk1_ntru1_Templ = Algo_Cases["NormalCase2"]["step_100"]["dEdx_ntrk_1_ntru_1_NormalCase21_step_100"]; 
  TH1F* NormalCase2_100_ntrk1_ntru2_Templ = Algo_Cases["NormalCase2"]["step_100"]["dEdx_ntrk_1_ntru_2_NormalCase22_step_100"]; 
  Printer(can, pdf_n + "2/NormalCase2.png", {NormalCase2_100_ntrk1_ntru1_Templ, NormalCase2_100_ntrk1_ntru2_Templ}, {ntrk1_ntru1_100, ntrk1_ntru2_100}, "Case-2: Normal", pdf_n + ".pdf");

  TH1F* NormalCase3_200_ntrk1_ntru1_Templ = Algo_Cases["NormalCase3"]["step_200"]["dEdx_ntrk_1_ntru_1_NormalCase31_step_200"]; 
  TH1F* NormalCase3_200_ntrk1_ntru2_Templ = Algo_Cases["NormalCase3"]["step_200"]["dEdx_ntrk_1_ntru_2_NormalCase32_step_200"]; 
  Printer(can, pdf_n + "3/NormalCase3.png", {NormalCase3_200_ntrk1_ntru1_Templ, NormalCase3_200_ntrk1_ntru2_Templ}, {ntrk1_ntru1_200, ntrk1_ntru2_200}, "Case-3: Normal", pdf_n + ".pdf");


  TH1F* ShiftNormalCase1_20_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalCase1"]["step_20"]["dEdx_ntrk_1_ntru_1_ShiftNormalCase11_step_20"]; 
  TH1F* ShiftNormalCase1_20_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalCase1"]["step_20"]["dEdx_ntrk_1_ntru_2_ShiftNormalCase12_step_20"]; 
  Printer(can, pdf_n + "1/ShiftNormalCase1.png", {ShiftNormalCase1_20_ntrk1_ntru1_Templ, ShiftNormalCase1_20_ntrk1_ntru2_Templ}, {ntrk1_ntru1_20, ntrk1_ntru2_20}, "Case-1: ShiftNormal", pdf_n + ".pdf");

  TH1F* ShiftNormalCase2_100_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalCase2"]["step_100"]["dEdx_ntrk_1_ntru_1_ShiftNormalCase21_step_100"]; 
  TH1F* ShiftNormalCase2_100_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalCase2"]["step_100"]["dEdx_ntrk_1_ntru_2_ShiftNormalCase22_step_100"]; 
  Printer(can, pdf_n + "2/ShiftNormalCase2.png", {ShiftNormalCase2_100_ntrk1_ntru1_Templ, ShiftNormalCase2_100_ntrk1_ntru2_Templ}, {ntrk1_ntru1_100, ntrk1_ntru2_100}, "Case-2: ShiftNormal", pdf_n + ".pdf");

  TH1F* ShiftNormalCase3_200_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalCase3"]["step_200"]["dEdx_ntrk_1_ntru_1_ShiftNormalCase31_step_200"]; 
  TH1F* ShiftNormalCase3_200_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalCase3"]["step_200"]["dEdx_ntrk_1_ntru_2_ShiftNormalCase32_step_200"]; 
  Printer(can, pdf_n + "3/ShiftNormalCase3.png", {ShiftNormalCase3_200_ntrk1_ntru1_Templ, ShiftNormalCase3_200_ntrk1_ntru2_Templ}, {ntrk1_ntru1_200, ntrk1_ntru2_200}, "Case-3: ShiftNormal", pdf_n + ".pdf");


  TH1F* ShiftNormalFFTCase1_20_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalFFTCase1"]["step_20"]["dEdx_ntrk_1_ntru_1_ShiftNormalFFTCase11_step_20"]; 
  TH1F* ShiftNormalFFTCase1_20_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalFFTCase1"]["step_20"]["dEdx_ntrk_1_ntru_2_ShiftNormalFFTCase12_step_20"]; 
  Printer(can, pdf_n + "1/ShiftNormalFFTCase1.png", {ShiftNormalFFTCase1_20_ntrk1_ntru1_Templ, ShiftNormalFFTCase1_20_ntrk1_ntru2_Templ}, {ntrk1_ntru1_20, ntrk1_ntru2_20}, "Case-1: ShiftNormalFFT", pdf_n + ".pdf");

  TH1F* ShiftNormalFFTCase2_100_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalFFTCase2"]["step_100"]["dEdx_ntrk_1_ntru_1_ShiftNormalFFTCase21_step_100"]; 
  TH1F* ShiftNormalFFTCase2_100_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalFFTCase2"]["step_100"]["dEdx_ntrk_1_ntru_2_ShiftNormalFFTCase22_step_100"]; 
  Printer(can, pdf_n + "2/ShiftNormalFFTCase2.png", {ShiftNormalFFTCase2_100_ntrk1_ntru1_Templ, ShiftNormalFFTCase2_100_ntrk1_ntru2_Templ}, {ntrk1_ntru1_100, ntrk1_ntru2_100}, "Case-2: ShiftNormalFFT", pdf_n + ".pdf");

  TH1F* ShiftNormalFFTCase3_200_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalFFTCase3"]["step_200"]["dEdx_ntrk_1_ntru_1_ShiftNormalFFTCase31_step_200"]; 
  TH1F* ShiftNormalFFTCase3_200_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalFFTCase3"]["step_200"]["dEdx_ntrk_1_ntru_2_ShiftNormalFFTCase32_step_200"]; 
  Printer(can, pdf_n + "3/ShiftNormalFFTCase3.png", {ShiftNormalFFTCase3_200_ntrk1_ntru1_Templ, ShiftNormalFFTCase3_200_ntrk1_ntru2_Templ}, {ntrk1_ntru1_200, ntrk1_ntru2_200}, "Case-3: ShiftNormalFFT", pdf_n + ".pdf");


  TH1F* ShiftNormalWidthFFTCase1_20_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalWidthFFTCase1"]["step_20"]["dEdx_ntrk_1_ntru_1_ShiftNormalWidthFFTCase11_step_20"]; 
  TH1F* ShiftNormalWidthFFTCase1_20_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalWidthFFTCase1"]["step_20"]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase12_step_20"]; 
  Printer(can, pdf_n + "1/ShiftNormalWidthFFTCase1.png", {ShiftNormalWidthFFTCase1_20_ntrk1_ntru1_Templ, ShiftNormalWidthFFTCase1_20_ntrk1_ntru2_Templ}, {ntrk1_ntru1_20, ntrk1_ntru2_20}, "Case-1: ShiftNormalWidthFFT", pdf_n + ".pdf");

  TH1F* ShiftNormalWidthFFTCase2_100_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalWidthFFTCase2"]["step_100"]["dEdx_ntrk_1_ntru_1_ShiftNormalWidthFFTCase21_step_100"]; 
  TH1F* ShiftNormalWidthFFTCase2_100_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalWidthFFTCase2"]["step_100"]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase22_step_100"]; 
  Printer(can, pdf_n + "2/ShiftNormalWidthFFTCase2.png", {ShiftNormalWidthFFTCase2_100_ntrk1_ntru1_Templ, ShiftNormalWidthFFTCase2_100_ntrk1_ntru2_Templ}, {ntrk1_ntru1_100, ntrk1_ntru2_100}, "Case-2: ShiftNormalWidthFFT", pdf_n + ".pdf");

  TH1F* ShiftNormalWidthFFTCase3_200_ntrk1_ntru1_Templ = Algo_Cases["ShiftNormalWidthFFTCase3"]["step_200"]["dEdx_ntrk_1_ntru_1_ShiftNormalWidthFFTCase31_step_200"]; 
  TH1F* ShiftNormalWidthFFTCase3_200_ntrk1_ntru2_Templ = Algo_Cases["ShiftNormalWidthFFTCase3"]["step_200"]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase32_step_200"]; 
  Printer(can, pdf_n + "3/ShiftNormalWidthFFTCase3.png", {ShiftNormalWidthFFTCase3_200_ntrk1_ntru1_Templ, ShiftNormalWidthFFTCase3_200_ntrk1_ntru2_Templ}, {ntrk1_ntru1_200, ntrk1_ntru2_200}, "Case-3: ShiftNormalWidthFFT", pdf_n + ".pdf");

  can -> Print(pdf_n + ".pdf]"); 

  // Demonstrate the excess in the 1-Track vs 1-Truth 
  TH1F* trk1_template = (TH1F*)Algo_Cases["NormalCase3"]["step_100"]["dEdx_ntrk_1_ntru_1_NormalCase31_step_100"] -> Clone("template1"); 
  TH1F* trk2_template = (TH1F*)Algo_Cases["NormalCase3"]["step_100"]["dEdx_ntrk_1_ntru_2_NormalCase32_step_100"] -> Clone("tempalte2");  
  TH1F* trk1_tru1 = (TH1F*)Algo_Cases["NormalCase3"]["step_100"]["ntrk1-tru1"] -> Clone("Truth"); 
  TH1F* trk1_excess = (TH1F*)trk1_template -> Clone("trk1_excess"); 
  trk2_template -> SetTitle("2-Track Template");
  trk1_excess -> Reset(); 
  trk1_excess -> SetTitle("Template-Truth Excess");
  for (int i(0); i < trk1_template -> GetNbinsX(); i++)
  {
    float diff = trk1_template -> GetBinContent(i+1) - trk1_tru1 -> GetBinContent(i+1); 
    if ( diff <= 0 ){ diff = 0; }
    trk1_excess -> SetBinContent(i+1, diff);
  }
  
  PlotHists({trk2_template}, {trk1_excess}, "Track-1 Template Excess", can); 
  trk2_template -> SetLineColor(kGreen); 
  trk2_template -> SetFillColorAlpha(kGreen, 0.2); 
  trk2_template -> SetLineWidth(2); 
  can -> Update(); 
  can -> Print("./Excess/Excess.png");
  
  delete trk1_template;
  delete trk2_template;
  delete trk1_tru1;
  delete trk1_excess;


  
  // Now show the integral difference between these cases from the truth as a function of the %
  int start = 20; 
  int end = 200; 
  int steps = 20; 
  int iter = (end - start) / steps; 
  TString pdf_steps = "./Steps/Case"; 
  float r; 

  can -> Print(pdf_steps + ".pdf["); 
  
  for (TString Lay : Layer)
  {
    MMF RatioMatrix;  
    MMMTH1F Algo_Case = debugs[Lay+ "_1000_1200_GeV"]; 
    for (int i(0); i <= iter; i++)
    {
      TString Steps = "step_"; Steps += (start + steps*i);   
      
      // Truth Increments
      TH1F* ntrk1_ntru1 = Algo_Case["NormalCase1"][Steps]["ntrk1-tru1"]; 
      TH1F* ntrk1_ntru2 = Algo_Case["NormalCase1"][Steps]["ntrk1-tru2"]; 

      TH1F* NormalCase11 = Algo_Case["NormalCase1"][Steps]["dEdx_ntrk_1_ntru_1_NormalCase11_"+Steps]; 
      TH1F* NormalCase12 = Algo_Case["NormalCase1"][Steps]["dEdx_ntrk_1_ntru_2_NormalCase12_"+Steps]; 
      
      TH1F* NormalCase21 = Algo_Case["NormalCase2"][Steps]["dEdx_ntrk_1_ntru_1_NormalCase21_"+Steps]; 
      TH1F* NormalCase22 = Algo_Case["NormalCase2"][Steps]["dEdx_ntrk_1_ntru_2_NormalCase22_"+Steps]; 

      TH1F* NormalCase31 = Algo_Case["NormalCase3"][Steps]["dEdx_ntrk_1_ntru_1_NormalCase31_"+Steps]; 
      TH1F* NormalCase32 = Algo_Case["NormalCase3"][Steps]["dEdx_ntrk_1_ntru_2_NormalCase32_"+Steps]; 

      TH1F* ShiftNormalWidthFFTCase11 = Algo_Case["ShiftNormalWidthFFTCase1"][Steps]["dEdx_ntrk_1_ntru_1_ShiftNormalWidthFFTCase11_"+Steps]; 
      TH1F* ShiftNormalWidthFFTCase12 = Algo_Case["ShiftNormalWidthFFTCase1"][Steps]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase12_"+Steps]; 
      
      TH1F* ShiftNormalWidthFFTCase21 = Algo_Case["ShiftNormalWidthFFTCase2"][Steps]["dEdx_ntrk_1_ntru_1_ShiftNormalWidthFFTCase21_"+Steps]; 
      TH1F* ShiftNormalWidthFFTCase22 = Algo_Case["ShiftNormalWidthFFTCase2"][Steps]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase22_"+Steps]; 
      
      TH1F* ShiftNormalWidthFFTCase31 = Algo_Case["ShiftNormalWidthFFTCase3"][Steps]["dEdx_ntrk_1_ntru_1_ShiftNormalWidthFFTCase31_"+Steps]; 
      TH1F* ShiftNormalWidthFFTCase32 = Algo_Case["ShiftNormalWidthFFTCase3"][Steps]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase32_"+Steps]; 
      
      can -> SetLogy();
      Printer(can, pdf_steps + "1/" + Lay + "_NormalCase1" + Steps + ".png", {NormalCase11, NormalCase12}, {ntrk1_ntru1, ntrk1_ntru2}, "Case 1: Normalization Fit at 20% Original 2-Truth", pdf_steps + ".pdf"); 
      Printer(can, pdf_steps + "2/" + Lay + "_NormalCase2" + Steps + ".png", {NormalCase21, NormalCase22}, {ntrk1_ntru1, ntrk1_ntru2}, "Case 2: Normalization Fit at 20% Original 2-Truth", pdf_steps + ".pdf"); 
      Printer(can, pdf_steps + "3/" + Lay + "_NormalCase3" + Steps + ".png", {NormalCase31, NormalCase32}, {ntrk1_ntru1, ntrk1_ntru2}, "Case 3: Normalization Fit at 20% Original 2-Truth", pdf_steps + ".pdf"); 
      
      Printer(can, pdf_steps + "1/" + Lay + "_ShiftNormalWidthFFTCase1" + Steps + ".png", {ShiftNormalWidthFFTCase11, ShiftNormalWidthFFTCase12}, {ntrk1_ntru1, ntrk1_ntru2}, "Case 1: ShiftNormalWidthFFT Fit at 20% Original 2-Truth", pdf_steps + ".pdf");
      Printer(can, pdf_steps + "2/" + Lay + "_ShiftNormalWidthFFTCase2" + Steps + ".png", {ShiftNormalWidthFFTCase21, ShiftNormalWidthFFTCase22}, {ntrk1_ntru1, ntrk1_ntru2}, "Case 2: ShiftNormalWidthFFT Fit at 20% Original 2-Truth", pdf_steps + ".pdf");   
      Printer(can, pdf_steps + "3/" + Lay + "_ShiftNormalWidthFFTCase3" + Steps + ".png", {ShiftNormalWidthFFTCase31, ShiftNormalWidthFFTCase32}, {ntrk1_ntru1, ntrk1_ntru2}, "Case 3: ShiftNormalWidthFFT Fit at 20% Original 2-Truth", pdf_steps + ".pdf"); 


      RatioMatrix["Normal-Case1"][Steps] = (NormalCase12 -> Integral()) / (ntrk1_ntru2 -> Integral());
      RatioMatrix["Normal-Case2"][Steps] = (NormalCase22 -> Integral()) / (ntrk1_ntru2 -> Integral());
      RatioMatrix["Normal-Case3"][Steps] = (NormalCase32 -> Integral()) / (ntrk1_ntru2 -> Integral());
     
      RatioMatrix["ShiftNormalWidthFFT-Case1"][Steps] = (ShiftNormalWidthFFTCase12 -> Integral()) / (ntrk1_ntru2 -> Integral());
      RatioMatrix["ShiftNormalWidthFFT-Case2"][Steps] = (ShiftNormalWidthFFTCase22 -> Integral()) / (ntrk1_ntru2 -> Integral());
      RatioMatrix["ShiftNormalWidthFFT-Case3"][Steps] = (ShiftNormalWidthFFTCase32 -> Integral()) / (ntrk1_ntru2 -> Integral());
    }
    
    can -> SetLogy(false);
    std::vector<TGraph*> GR = GenerateMultiGraph(RatioMatrix, "Integral Ratio Between Fit Prediction and Truth: " + Lay + " 1000-1200 GeV", can, 0, 2, "#frac{Pred}{Truth}"); 
    can -> Print("./RatioMatrix/" + Lay + "/RatioMatrix.png");
    can -> Print(pdf_steps + ".pdf");
    BulkDelete(GR);
  }
  can -> Print(pdf_steps + ".pdf]"); 
 
  for (TString Lay : Layer)
  {
    MMF RatioMatrixJet; 
    for (TString JE : JetEnergy)
    {
      // Truth Increments
      TH1F* ntrk1_ntru2_40 = debugs[Lay + "_"+JE]["ShiftNormalWidthFFTCase1"]["step_40"]["ntrk1-tru2"]; 
      TH1F* ntrk1_ntru2_60 = debugs[Lay + "_"+JE]["ShiftNormalWidthFFTCase1"]["step_60"]["ntrk1-tru2"]; 
      TH1F* ntrk1_ntru2_80 = debugs[Lay + "_"+JE]["ShiftNormalWidthFFTCase1"]["step_80"]["ntrk1-tru2"]; 

      TH1F* ShiftNormalWidthFFTCase12_40 = debugs[Lay + "_" + JE]["ShiftNormalWidthFFTCase1"]["step_40"]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase12_step_40"]; 
      TH1F* ShiftNormalWidthFFTCase22_60 = debugs[Lay + "_" + JE]["ShiftNormalWidthFFTCase2"]["step_60"]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase22_step_60"]; 
      TH1F* ShiftNormalWidthFFTCase32_80 = debugs[Lay + "_" + JE]["ShiftNormalWidthFFTCase3"]["step_80"]["dEdx_ntrk_1_ntru_2_ShiftNormalWidthFFTCase32_step_80"]; 

      TH1F* NormalCase12_40 = debugs[Lay + "_" + JE]["NormalCase1"]["step_40"]["dEdx_ntrk_1_ntru_2_NormalCase12_step_40"]; 
      TH1F* NormalCase22_60 = debugs[Lay + "_" + JE]["NormalCase2"]["step_60"]["dEdx_ntrk_1_ntru_2_NormalCase22_step_60"]; 
      TH1F* NormalCase32_80 = debugs[Lay + "_" + JE]["NormalCase3"]["step_80"]["dEdx_ntrk_1_ntru_2_NormalCase32_step_80"]; 

      if (ntrk1_ntru2_40 == 0){continue;}

      float r1 = (ShiftNormalWidthFFTCase12_40 -> Integral()) / (ntrk1_ntru2_40 -> Integral());
      float r2 = (ShiftNormalWidthFFTCase22_60 -> Integral()) / (ntrk1_ntru2_60 -> Integral());
      float r3 = (ShiftNormalWidthFFTCase32_80 -> Integral()) / (ntrk1_ntru2_80 -> Integral());
      
      RatioMatrixJet["ShiftNormalWidthFFT-Case1-40%"][JE] = r1;
      RatioMatrixJet["ShiftNormalWidthFFT-Case2-60%"][JE] = r2;
      RatioMatrixJet["ShiftNormalWidthFFT-Case3-80%"][JE] = r3;

      r1 = (NormalCase12_40 -> Integral()) / (ntrk1_ntru2_40 -> Integral());
      r2 = (NormalCase22_60 -> Integral()) / (ntrk1_ntru2_60 -> Integral());
      r3 = (NormalCase32_80 -> Integral()) / (ntrk1_ntru2_80 -> Integral());
      
      RatioMatrixJet["Normal-Case1-40%"][JE] = r1;
      RatioMatrixJet["Normal-Case2-60%"][JE] = r2;
      RatioMatrixJet["Normal-Case3-80%"][JE] = r3;

    }
    can -> SetLogy(false);
    std::vector<TGraph*> GR = GenerateMultiGraph(RatioMatrixJet, "Integral Ratio Between Prediction and Truth in the " + Lay, can, 0.1, 10, "#frac{Pred}{Truth}"); 
    can -> Print("./RatioMatrix/" + Lay + "/RatioMatrixJet.png");
    BulkDelete(GR);
    
    RatioMatrixJet.clear();
  }
  delete can;
}
