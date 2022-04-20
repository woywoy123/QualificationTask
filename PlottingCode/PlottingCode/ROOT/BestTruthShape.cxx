#include "../PlottingCode/BestTruthShape.h"

void GetWeightedError(_TS Template, 
                      std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST, 
                      std::map<_TS, std::vector<TH1F*>>* CTIDE, 
                      _TS trk, 
                      _TS FittingMode, 
                      std::map<_TS, float>* Error_Matrix
)
{
  for (_TS Layer : Layer_H)
  {
    for (_TS Energy : Energy_H)
    {
      std::vector<TH1F*> ntrk_T = (*CTIDE)[Layer + "/" + Energy + "/ntrk" + trk + "/InsideTruth"]; 
      for (_TS Alg : Algos_H)
      {
        std::vector<TH1F*> ntrk_F = POST[Layer + "/" + Energy + "/" + Alg][FittingMode + "_ntrk" + trk]; 
        _TS TMP_Key = Template + "/ntrk" + trk + "/" + Layer + "/" + Energy + "/" + Alg;
        if (ntrk_F.size() == 0){continue;}
        if (FittingMode == "Test" && Alg == "Experimental"){ (*Error_Matrix)[TMP_Key] = -1; continue;}
        
        (*Error_Matrix)[TMP_Key] = WeightedShapeError(ntrk_F, ntrk_T);
      }
    }
  }
}

void PrintEnergyAlgoTable(_TS root_key, _TS outdir, _TS mode, std::map<_TS, float> Matrix)
{
  Table T; 
  T.x_cellNames = Energy_H; 
  T.y_cellNames = Algos_H; 
  T.Round = 6; 
  T.Outdir = outdir; 
  T.Name = mode + ".txt";
  
  for (_TS Algo : Algos_H)
  {
    T.Newline();
    for (_TS Energy : Energy_H)
    {
      T.AddValue(Matrix[root_key + Energy + "/" + Algo]); 
    }
  }
  
  T.CompileTable(); 
}

void PrintSummaryOfTemplateVariation(_TS mode, _TS trk, _TS outdir, std::map<_TS, float> Matrix)
{
  MergeTables T_all; 
  T_all.Name = mode + ".txt"; 
  T_all.Outdir = outdir;
  T_all.x_cellNames = Energy_H; 
  T_all.y_cellNames = Algos_H; 
  T_all.z_cellNames = Layer_H;
  T_all.Round = 6;
  for (_TS Layer : Layer_H)
  {
    Table T; 
    T.x_cellNames = Energy_H; 
    T.y_cellNames = Algos_H; 
    T.Round = 6; 
    T.Outdir = outdir; 
    T.Name = "WeightedError.txt";
    
    for (_TS Algo : Algos_H)
    {
      T.Newline();
      for (_TS Energy : Energy_H)
      {
        T.AddValue(Matrix[mode + "/ntrk" + trk + "/" + Layer + "/" + Energy + "/" + Algo]); 
      }
    }
    T_all.InputTables.push_back(T); 
  } 
  T_all.CompileTable();
}

void PrintBestTemplateVariationForAlgorithms(_TS Algo, _TS trk, _TS Layer, std::vector<_TS> TemplatesVariation, _TS outdir, std::map<_TS, float> Matrix)
{

  Table T; 
  T.Name = Algo + ".txt"; 
  T.Outdir = outdir;
  T.x_cellNames = Energy_H;
  T.y_cellNames = TemplatesVariation; 
  T.Round = 6; 
  std::vector<TGraph*> GR;
  for (_TS mode : TemplatesVariation)
  {
    std::vector<float> tmp; 
    T.Newline();
    for (_TS Energy : Energy_H)
    {
      T.AddValue(Matrix[mode + "/ntrk" + trk + "/" + Layer + "/" + Energy + "/" + Algo]); 
      tmp.push_back(Matrix[mode + "/ntrk" + trk + "/" + Layer + "/" + Energy + "/" + Algo]); 
    } 
    TGraph* gr = MakeGraph(tmp, mode.ReplaceAll("_", "")); 
    GR.push_back(gr); 
  }
  T.CompileTable();
  PlotMultiGraph(&GR, "Minimizer Comparison: " + Layer + " " + Algo + " Track-" + trk, 0.00001, 1, outdir + Algo + ".png"); 
}

void PrintBestAlgorithms(_TS trk, _TS Layer, std::vector<_TS> TemplatesVariation, _TS outdir, std::map<_TS, float> Matrix)
{
  MergeTables T_all; 
  T_all.Name = "Tracks-" + trk + ".txt"; 
  T_all.Outdir = outdir;
  T_all.x_cellNames = Energy_H; 
  T_all.y_cellNames = Algos_H; 
  T_all.z_cellNames = TemplatesVariation;
  T_all.Round = 6;
  for (_TS mode : TemplatesVariation)
  {
    Table T; 
    T.Name = mode + ".txt"; 
    T.Outdir = outdir;
    T.x_cellNames = Energy_H;
    T.y_cellNames = Algos_H; 
    T.Round = 6; 
    for (_TS Algo : Algos_H)
    {
      T.Newline();
      for (_TS Energy : Energy_H)
      {
        T.AddValue(Matrix[mode + "/ntrk" + trk + "/" + Layer + "/" + Energy + "/" + Algo]); 
      } 
    }
    T_all.InputTables.push_back(T); 
  }
  T_all.CompileTable();
  
  T_all.y_i = 0; 
  std::vector<TGraph*> GR;
  for (_TS Alg : Algos_H)
  {
    T_all.x_i = 0; 
    std::vector<float> tmp; 
    for (_TS Energy : Energy_H)
    {
      tmp.push_back(T_all.values[T_all.y_i][T_all.x_i]); 
      T_all.x_i++; 
    }
    TGraph* gr = MakeGraph(tmp, Alg); 
    GR.push_back(gr); 
    T_all.y_i++; 
  }
  PlotMultiGraph(&GR, "Average Over Template Variations Using Different Fitting Algorithms: " + Layer + " Track-" + trk, 0.0001, 1, outdir + "Track-" + trk + ".png"); 
}




int main(int argc, char* argv[])
{

  _TS OutdirTruth = "FitToTruth"; 
  _TS OutdirData = "FitToData"; 
  _TS Truth = argv[1];
  
  TFile* C = new TFile(Truth); 
  std::map<_TS, std::vector<TH1F*>> CTIDE = ReadCTIDE(C); 
  std::map<_TS, float> Error_Matrix_Test; 
  std::map<_TS, float> Error_Matrix_Data; 
  std::vector<_TS> Modes; 
  std::vector<_TS> trks = {"1", "2", "3", "4"}; 

  std::vector<std::thread> th; 
  for (int i(2); i < argc; i++)
  {
    std::vector<_TS> f = Split(argv[i], "/");
    _TS mode = f[f.size()-2]; 
    Modes.push_back(mode); 
    
    TFile* F = new TFile(argv[i]); 
    std::map<_TS, std::map<_TS, std::vector<TH1F*>>> Read = ReadPostAnalysis(F);
    for (int k(0); k < trks.size(); k++)
    {
      th.push_back(std::thread(GetWeightedError, mode, Read, &CTIDE, trks[k], "Test", &Error_Matrix_Test));
      th.push_back(std::thread(GetWeightedError, mode, Read, &CTIDE, trks[k], "Template", &Error_Matrix_Data));
    }
    for (std::thread &t : th){ t.join(); }
    th.clear();
    Read.clear(); 
    F -> Close(); 
    delete F;
    std::cout << mode << std::endl;
  }

  for (_TS Layer : Layer_H)
  {
    for (_TS mode : Modes)
    {
      _TS trk1 = mode + "/ntrk1/" + Layer + "/";  
      _TS trk2 = mode + "/ntrk2/" + Layer + "/";
      _TS trk3 = mode + "/ntrk3/" + Layer + "/";
      _TS trk4 = mode + "/ntrk4/" + Layer + "/";

      _TS trk_out_T = OutdirTruth + "/" + Layer + "/Track-";   
      _TS trk_out_D = OutdirData + "/" + Layer + "/Track-";   

      th.push_back(std::thread(PrintEnergyAlgoTable, trk1, trk_out_T + "1/", mode, Error_Matrix_Test)); 
      th.push_back(std::thread(PrintEnergyAlgoTable, trk2, trk_out_T + "2/", mode, Error_Matrix_Test)); 
      th.push_back(std::thread(PrintEnergyAlgoTable, trk3, trk_out_T + "3/", mode, Error_Matrix_Test)); 
      th.push_back(std::thread(PrintEnergyAlgoTable, trk4, trk_out_T + "4/", mode, Error_Matrix_Test)); 

      th.push_back(std::thread(PrintEnergyAlgoTable, trk1, trk_out_D + "1/", mode, Error_Matrix_Data)); 
      th.push_back(std::thread(PrintEnergyAlgoTable, trk2, trk_out_D + "2/", mode, Error_Matrix_Data)); 
      th.push_back(std::thread(PrintEnergyAlgoTable, trk3, trk_out_D + "3/", mode, Error_Matrix_Data)); 
      th.push_back(std::thread(PrintEnergyAlgoTable, trk4, trk_out_D + "4/", mode, Error_Matrix_Data)); 
    }
  }

  for (_TS mode : Modes)
  {
    _TS trk_out_T = OutdirTruth + "/Summary/PerTemplateVariation/Track-"; 
    _TS trk_out_D = OutdirData + "/Summary/PerTemplateVariation/Track-"; 

    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "1", trk_out_T + "1", Error_Matrix_Test));
    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "2", trk_out_T + "2", Error_Matrix_Test));
    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "3", trk_out_T + "3", Error_Matrix_Test));
    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "4", trk_out_T + "4", Error_Matrix_Test));

    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "1", trk_out_D + "1", Error_Matrix_Data));
    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "2", trk_out_D + "2", Error_Matrix_Data));
    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "3", trk_out_D + "3", Error_Matrix_Data));
    th.push_back(std::thread(PrintSummaryOfTemplateVariation, mode, "4", trk_out_D + "4", Error_Matrix_Data));
  }
  
  for (_TS Layer : Layer_H)
  {
    for (_TS Algo : Algos_H)
    {
      _TS trk_out_T = OutdirTruth + "/Summary/PerAlgorithmVariation/" + Layer + "/Track-"; 
      _TS trk_out_D = OutdirData + "/Summary/PerAlgorithmVariation/" + Layer + "/Track-"; 
     
      PrintBestTemplateVariationForAlgorithms(Algo, "1", Layer, Modes, trk_out_T + "1/", Error_Matrix_Test); 
      PrintBestTemplateVariationForAlgorithms(Algo, "2", Layer, Modes, trk_out_T + "2/", Error_Matrix_Test); 
      PrintBestTemplateVariationForAlgorithms(Algo, "3", Layer, Modes, trk_out_T + "3/", Error_Matrix_Test); 
      PrintBestTemplateVariationForAlgorithms(Algo, "4", Layer, Modes, trk_out_T + "4/", Error_Matrix_Test); 

      PrintBestTemplateVariationForAlgorithms(Algo, "1", Layer, Modes, trk_out_D + "1/", Error_Matrix_Data); 
      PrintBestTemplateVariationForAlgorithms(Algo, "2", Layer, Modes, trk_out_D + "2/", Error_Matrix_Data); 
      PrintBestTemplateVariationForAlgorithms(Algo, "3", Layer, Modes, trk_out_D + "3/", Error_Matrix_Data); 
      PrintBestTemplateVariationForAlgorithms(Algo, "4", Layer, Modes, trk_out_D + "4/", Error_Matrix_Data); 
    }
  }

  for (_TS Layer : Layer_H)
  {
    _TS trk_out_T = OutdirTruth + "/Summary/OverallPerformance/" + Layer + "/"; 
    _TS trk_out_D = OutdirData + "/Summary/OverallPerformance/" + Layer + "/"; 

    PrintBestAlgorithms("1", Layer, Modes, trk_out_T, Error_Matrix_Test); 
    PrintBestAlgorithms("2", Layer, Modes, trk_out_T, Error_Matrix_Test); 
    PrintBestAlgorithms("3", Layer, Modes, trk_out_T, Error_Matrix_Test); 
    PrintBestAlgorithms("4", Layer, Modes, trk_out_T, Error_Matrix_Test); 

    PrintBestAlgorithms("1", Layer, Modes, trk_out_D, Error_Matrix_Data); 
    PrintBestAlgorithms("2", Layer, Modes, trk_out_D, Error_Matrix_Data); 
    PrintBestAlgorithms("3", Layer, Modes, trk_out_D, Error_Matrix_Data); 
    PrintBestAlgorithms("4", Layer, Modes, trk_out_D, Error_Matrix_Data); 
  }





  for (std::thread &t : th){ t.join(); }
  Error_Matrix_Test.clear(); 
  Error_Matrix_Data.clear(); 
  CTIDE.clear();
  C -> Close(); 


  return 0;
}
