#include "../PlottingCode/BestTruthShape.h"

void GetWeightedError(_TS Template, 
                      std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST, 
                      std::map<_TS, std::vector<TH1F*>> CTIDE, 
                      _TS trk, 
                      _TS FittingMode, 
                      std::map<_TS, float>* Error_Matrix
)
{
  for (_TS Layer : Layer_H)
  {
    for (_TS Energy : Energy_H)
    {
      std::vector<TH1F*> ntrk_T = CTIDE[Layer + "/" + Energy + "/ntrk" + trk + "/InsideTruth"]; 
      for (_TS Alg : Algos_H)
      {
        std::vector<TH1F*> ntrk_F = POST[Layer + "/" + Energy + "/" + Alg][FittingMode + "_ntrk" + trk]; 
        if (ntrk_F.size() == 0){continue;}
        
        _TS TMP_Key = Template + "/ntrk" + trk + "/" + Layer + "/" + Energy + "/" + Alg; 
        (*Error_Matrix)[TMP_Key] = WeightedShapeError(ntrk_F, ntrk_T);
      }
    }
  }
}

void PrintEnergyAlgoTable(_TS root_key, _TS outdir, std::map<_TS, float> Matrix)
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
      T.AddValue(Matrix[root_key + Energy + "/" + Algo]); 
    }
  }
  
  T.CompileTable(); 
}

int main(int argc, char* argv[])
{
  _TS Truth = argv[1];
  std::map<_TS, std::vector<TH1F*>> CTIDE = ReadCTIDE(Truth); 
  std::map<_TS, std::map<_TS, std::vector<TH1F*>>> POST;  

  std::map<_TS, float> Error_Matrix_Test; 
  std::vector<_TS> Modes; 
  _TS FittingTo = "Test"; 
  for (int i(2); i < argc; i++)
  {
    std::vector<_TS> f = Split(argv[i], "/");
    _TS mode = f[f.size()-2]; 
    Modes.push_back(mode); 
    
    std::map<_TS, std::map<_TS, std::vector<TH1F*>>> Read = ReadPostAnalysis(argv[i]);
    for (MMVTFi x = Read.begin(); x != Read.end(); x++)
    {
      POST[mode + "/" + x -> first] = x -> second;
    }
    
    GetWeightedError(mode, Read, CTIDE, "1", FittingTo, &Error_Matrix_Test); 
    GetWeightedError(mode, Read, CTIDE, "2", FittingTo, &Error_Matrix_Test);
    GetWeightedError(mode, Read, CTIDE, "3", FittingTo, &Error_Matrix_Test);
    GetWeightedError(mode, Read, CTIDE, "4", FittingTo, &Error_Matrix_Test);
    break;
  }

  // Print out the error table /Layer/Track-n/Mode
  for (_TS Layer : Layer_H)
  {
    for (_TS mode : Modes)
    {
      _TS trk1 = mode + "/ntrk1/" + Layer + "/";  
      _TS trk2 = mode + "/ntrk2/" + Layer + "/";
      _TS trk3 = mode + "/ntrk3/" + Layer + "/";
      _TS trk4 = mode + "/ntrk4/" + Layer + "/";

      _TS trk1_out = "FittingAlgorithm/" + Layer + "/Track-1/" + mode;   
      _TS trk2_out = "FittingAlgorithm/" + Layer + "/Track-2/" + mode; 
      _TS trk3_out = "FittingAlgorithm/" + Layer + "/Track-3/" + mode; 
      _TS trk4_out = "FittingAlgorithm/" + Layer + "/Track-4/" + mode; 

      PrintEnergyAlgoTable(trk1, trk1_out, Error_Matrix_Test); 
      PrintEnergyAlgoTable(trk2, trk2_out, Error_Matrix_Test); 
      PrintEnergyAlgoTable(trk3, trk3_out, Error_Matrix_Test); 
      PrintEnergyAlgoTable(trk4, trk4_out, Error_Matrix_Test); 
    }
  }
  


  // Check the best performing fitting algorithm for each template variation per layer 
  














  return 0;
}
