#include "../PlottingCode/BestFLost.h"

void GetFLost(_TS Template, 
              std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST, 
              std::map<_TS, std::vector<TH1F*>>* CTIDE, 
              _TS FittingMode, 
              std::map<_TS, float>* Matrix
)
{
  for (_TS Layer : Layer_H)
  {
    for (_TS Energy : Energy_H)
    {
      std::vector<TH1F*> trk1_T = (*CTIDE)[Layer + "/" + Energy + "/ntrk1/InsideTruth"]; 
      std::vector<TH1F*> trk2_T = (*CTIDE)[Layer + "/" + Energy + "/ntrk2/InsideTruth"]; 
      std::vector<TH1F*> trk3_T = (*CTIDE)[Layer + "/" + Energy + "/ntrk3/InsideTruth"]; 
      std::vector<TH1F*> trk4_T = (*CTIDE)[Layer + "/" + Energy + "/ntrk4/InsideTruth"]; 

      for (_TS Alg : Algos_H)
      {
        std::vector<TH1F*> trk1_P = POST[Layer + "/" + Energy + "/" + Alg][FittingMode + "_ntrk1"]; 
        std::vector<TH1F*> trk2_P = POST[Layer + "/" + Energy + "/" + Alg][FittingMode + "_ntrk2"]; 
        std::vector<TH1F*> trk3_P = POST[Layer + "/" + Energy + "/" + Alg][FittingMode + "_ntrk3"]; 
        std::vector<TH1F*> trk4_P = POST[Layer + "/" + Energy + "/" + Alg][FittingMode + "_ntrk4"]; 

        _TS TMP_Key = Template + "/" + Layer + "/" + Energy + "/" + Alg;
        float FL2_T = Flost2(trk1_T, trk2_T); 
        float FL2_P = Flost2(trk1_P, trk2_P); 
        
        float FL3_T = Flost3(trk1_T, trk2_T, trk3_T); 
        float FL3_P = Flost3(trk1_P, trk2_P, trk3_P); 
     
        if (ErrorMode && Flost == "2")
        {
          (*Matrix)[TMP_Key] = 100*std::abs(FL2_T - FL2_P)/float(FL2_T); 
        }
        else if (ErrorMode && Flost == "3")
        {
          (*Matrix)[TMP_Key] = 100*std::abs(FL3_T - FL3_P)/float(FL3_T);
        }
        else if (!ErrorMode && Flost == "2")
        {
          (*Matrix)[TMP_Key] = FL2_P; 
        }
        else if (!ErrorMode && Flost == "3")
        {
          (*Matrix)[TMP_Key] = FL3_P;
        }
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

void PrintSummaryOfTemplateVariation(_TS mode, _TS outdir, std::map<_TS, float> Matrix)
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
    T.Name = mode + "_" + Layer + ".txt";
    
    for (_TS Algo : Algos_H)
    {
      T.Newline();
      for (_TS Energy : Energy_H)
      {
        T.AddValue(Matrix[mode + "/" + Layer + "/" + Energy + "/" + Algo]); 
      }
    }
    T_all.InputTables.push_back(T); 
  } 
  T_all.CompileTable();
}

void PrintBestTemplateVariationForAlgorithms(_TS Algo, _TS Layer, std::vector<_TS> TemplatesVariation, _TS outdir, std::map<_TS, float> Matrix)
{
  
  min = -1; 
  max = -1; 
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
      T.AddValue(Matrix[mode + "/" + Layer + "/" + Energy + "/" + Algo]); 
      tmp.push_back(Matrix[mode + "/" + Layer + "/" + Energy + "/" + Algo]); 
    } 
    TGraph* gr = MakeGraph(tmp, mode.ReplaceAll("_", "")); 
    GR.push_back(gr); 

    for (float f : tmp)
    {
      if (min <= 0){ min = f; }
      if (max <= 0){ max = f; }

      if (min > f && f > 0){ min = f; }
      if (max < f && f > 0){ max = f; }
    }
  }
  T.CompileTable();
  PlotMultiGraph(&GR, "Minimizer Comparison: " + Layer + " " + Algo + " FLost-" + Flost + "", outdir + Algo + ".png"); 
}

void PrintBestAlgorithms(_TS Layer, std::vector<_TS> TemplatesVariation, _TS outdir, std::map<_TS, float> Matrix)
{
  min = -1; 
  max = -1; 

  MergeTables T_all; 
  T_all.Name = "FLost.txt"; 
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
        T.AddValue(Matrix[mode + "/" + Layer + "/" + Energy + "/" + Algo]); 
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
    for (float f : tmp)
    {
      if (min <= 0){ min = f; }
      if (max <= 0){ max = f; }

      if (min > f && f > 0){ min = f; }
      if (max < f && f > 0){ max = f; }
    }
    
    TGraph* gr = MakeGraph(tmp, Alg); 
    GR.push_back(gr); 
    T_all.y_i++; 
  }
  PlotMultiGraph(&GR, "Average FLost-" + Flost + " Over Template Variations Using Different Fitting Algorithms: " + Layer, outdir + "FLost.png"); 
}




int main(int argc, char* argv[])
{

  _TS OutdirData = argv[1]; 
  _TS Truth = argv[2];

  if (OutdirData.Contains("FLost2")){ Flost = "2"; }
  if (OutdirData.Contains("FLost3")){ Flost = "3"; }

  if (OutdirData.Contains("Error")){ ErrorMode = true; }
  else { ErrorMode = false;}

  TFile* C = new TFile(Truth); 
  std::map<_TS, std::vector<TH1F*>> CTIDE = ReadCTIDE(C); 
  std::map<_TS, float> Matrix_Data; 
  std::vector<_TS> Modes; 
  std::vector<_TS> trks = {"1", "2", "3", "4"}; 

  std::vector<std::thread> th; 
  for (int i(3); i < argc; i++)
  {
    std::vector<_TS> f = Split(argv[i], "/");
    _TS mode = f[f.size()-2]; 
    Modes.push_back(mode); 
    
    TFile* F = new TFile(argv[i]); 
    std::map<_TS, std::map<_TS, std::vector<TH1F*>>> Read = ReadPostAnalysis(F);
    GetFLost(mode, Read, &CTIDE, "Template", &Matrix_Data);
    Read.clear(); 
    F -> Close(); 
    std::cout << "-> " << mode << std::endl;
    delete F;
  }

  for (_TS Layer : Layer_H)
  {
    for (_TS mode : Modes)
    {
      _TS flost = mode + "/" + Layer + "/";  
      _TS trk_out_D = OutdirData + "/" + Layer + "/";   
      PrintEnergyAlgoTable(flost, trk_out_D, mode, Matrix_Data); 
    }
  }

  for (_TS mode : Modes)
  {
    _TS trk_out_D = OutdirData + "/Summary/PerTemplateVariation/"; 

    PrintSummaryOfTemplateVariation( mode, trk_out_D, Matrix_Data);
  }
  for (std::thread &t : th){ t.join(); }
  
  for (_TS Layer : Layer_H)
  {
    for (_TS Algo : Algos_H)
    {
      _TS trk_out_D = OutdirData + "/Summary/PerAlgorithmVariation/" + Layer + "/"; 

      PrintBestTemplateVariationForAlgorithms(Algo, Layer, Modes, trk_out_D, Matrix_Data); 
    }
  }

  for (_TS Layer : Layer_H)
  {
    _TS trk_out_D = OutdirData + "/Summary/OverallPerformance/" + Layer + "/"; 
    PrintBestAlgorithms(Layer, Modes, trk_out_D, Matrix_Data); 
  }

  Matrix_Data.clear(); 
  CTIDE.clear();
  C -> Close(); 
  return 0;
}
