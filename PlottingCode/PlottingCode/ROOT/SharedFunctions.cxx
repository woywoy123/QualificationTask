#include "../PlottingCode/SharedFunctions.h"

std::vector<_TS> Split(_TS Input, _TS Sub)
{
  std::vector<_TS> Output;
  if (!Input.Contains(Sub)){ return Output; }
  while (true)
  {
    int l = Input.First(Sub); 
    Output.push_back(Input(0, l)); 
    Input.Remove(0, l+Sub.Length());
    if (!Input.Contains(Sub)){ Output.push_back(Input); break;}
  }
  return Output; 
}

void BulkDelete(std::vector<TGraph*> gr)
{
  for (TGraph* g : gr){ delete g; }
}

std::vector<float> ReadTH1F(TH1F* TH)
{
  std::vector<float> out; 
  for (int i(0); i < TH -> GetNbinsX(); i++)
  {
    out.push_back(TH -> GetBinContent(i+1)); 
  }
  return out; 
}

void Normalize(std::vector<float>* v)
{
  float sum = 0; 
  for (float i : *v){ sum += i; }
  if (sum == 0){ return; }
  for (int i(0); i < v -> size(); i++){ (*v)[i] = (*v)[i] / sum; }
}

void Normalize(TH1F* Hist)
{
  float x = Hist -> Integral();
  if (x == 0){return;}
  Hist -> Scale(1/x); 
}

float WeightedShapeError(std::vector<TH1F*> Pred, std::vector<TH1F*> Truth)
{
  float n_int = 0; 
  for (TH1F* T : Truth){n_int += (T -> Integral());}
  
  float err = 0; 
  for (int i(0); i < Pred.size(); i++)
  {
    std::vector<float> P = ReadTH1F(Pred[i]); 
    std::vector<float> T = ReadTH1F(Truth[i]); 
    float m_int = Truth[i] -> Integral();
    Normalize(&P); 
    Normalize(&T); 
    
    float err_t = 0; 
    for (int j(0); j < P.size(); j++)
    {
      float sum = P[j] + T[j];
      if (sum == 0){continue;}
      float diff = P[j] - T[j]; 
      
      err_t += (std::pow(diff, 2) / sum);
    }
    err += 0.5*err_t*m_int;
  }
  
  err = err/n_int;
  return err; 
}

void Normalize(std::vector<TH1F*> Hists)
{
  for (TH1F* H : Hists){Normalize(H);}
}








void Table::CompileTable()
{
  Output.push_back("================= Raw Values =================="); 
  Output.push_back(Outdir);
  MakeHeader();
  for (y_i = 0; y_i < values.size(); y_i++){AppendLine(values);}
  BestValues();
  ResultList(); 
  OtherStats(); 
  WriteToFile();
  //Print(); 
}

void Table::MakeHeader()
{
  std::vector<_TS> All_X = x_cellNames; 
  for (y_i = 0; y_i < values.size(); y_i++)
  {
    for (x_i = 0; x_i < values[y_i].size(); x_i++)
    {
      All_X.push_back(Rounded(values[y_i][x_i])); 
    }
  }

  x_width = LargestElement(All_X); 
  y_width = LargestElement(y_cellNames);
  
  _TS TMP = ""; 
  RepeatString(&TMP, y_width, " "); 
  TMP += Spacer; 
  
  for (_TS i : x_cellNames)
  {
    TMP += i; 
    RepeatString(&TMP, x_width - i.Length(), " "); 
    TMP += Spacer; 
  }
  Output.push_back(TMP); 
}


void Table::AppendLine(std::vector<std::vector<float>> val)
{
  _TS TMP = ""; 
  TMP += y_cellNames[y_i]; 
  RepeatString(&TMP, y_width - y_cellNames[y_i].Length(), " ");
  TMP += Spacer; 
  
  for (x_i = 0; x_i < val[y_i].size(); x_i++)
  {
    _TS s = Rounded(val[y_i][x_i]); 
    if (val[y_i][x_i] == -1){ s = "x"; }
    TMP += s; 
    RepeatString(&TMP, x_width - s.Length(), " "); 
    TMP += Spacer; 
  }
  Output.push_back(TMP);
}

void Table::BestValues()
{
  std::vector<std::vector<float>> Score(y_cellNames.size(), std::vector<float>{}); 
  for (x_i = 0; x_i < x_cellNames.size(); x_i++)
  {
    float L = -1;
    for (y_i = 0; y_i < y_cellNames.size(); y_i++)
    {
      float v = values[y_i][x_i]; 

      if (L == -1){ L = v; }
      if ( L > v && v > 0){ L = v; } 
    }
    
    for (y_i = 0; y_i < y_cellNames.size(); y_i++)
    {
      if (values[y_i][x_i] == L)
      { 
        Score[y_i].push_back(1); 
        continue;
      }
      Score[y_i].push_back(0); 
    }
  }

  for (y_i = 0; y_i < y_cellNames.size(); y_i++)
  {
    int sum = 0; 
    for (x_i = 0; x_i < x_cellNames.size(); x_i++){sum += Score[y_i][x_i];}
    Score[y_i].push_back(sum); 
  }
  x_cellNames.push_back("Sum Scores"); 

  Output.push_back(""); 
  Output.push_back("________ Scores ________");

  MakeHeader(); 
  Integer = true; 
  for (y_i = 0; y_i < Score.size(); y_i++){AppendLine(Score);}
  Integer = false;
  x_cellNames.pop_back(); 
  TMP_Score = Score;
}

void Table::ResultList()
{
  Output.push_back(""); 
  Output.push_back("________ Best Values _______"); 
  for (x_i = 0; x_i < x_cellNames.size(); x_i++)
  {
    _TS TMP = ""; 
    TMP += x_cellNames[x_i]; 
    RepeatString(&TMP, x_width - TMP.Length(), " "); 
    TMP += Spacer; 
    
    for (y_i = 0; y_i < y_cellNames.size(); y_i++)
    {
      float v = TMP_Score[y_i][x_i]*values[y_i][x_i]; 
      if (v <= 0){continue;}

      TMP += y_cellNames[y_i]; 
      RepeatString(&TMP, LargestElement(y_cellNames) - y_cellNames[y_i].Length(), " ");
      TMP += Spacer; 

      TMP += "=> " + Rounded(v) + " ||| ";
    }
    Output.push_back(TMP);
  }
}

void Table::OtherStats() 
{
  std::vector<_TS> Temp;
  Output.push_back(""); 
  Output.push_back("________ Other Stats _______"); 
  Output.push_back("----> Absolute Lowest Error For Each Fitting Algorithm:"); 
  
  for (y_i = 0; y_i < y_cellNames.size(); y_i++)
  {
    _TS TMP = ""; 
    TMP += y_cellNames[y_i]; 
    RepeatString(&TMP, y_width - TMP.Length(), " "); 
    TMP += Spacer; 
   
    std::vector<int> ind = LowestValue(values[y_i]);
    float av = Sum(values[y_i]) / float(x_cellNames.size()); 
    if (av < 0){TMP += "Excluded!";}
    else { TMP += Rounded(av); }
    Temp.push_back(TMP); 
    
    TMP += Spacer; 
    for (int x_i : ind)
    {
      TMP += x_cellNames[x_i]; 
      RepeatString(&TMP, LargestElement(x_cellNames) - x_cellNames[x_i].Length(), " "); 
      TMP += Spacer;
      TMP += "=> " + Rounded( values[y_i][x_i] ) + " ||| "; 
    }

    Output.push_back(TMP);
  }
  Output.push_back(""); 
  Output.push_back("----> Average Error For Algorithm:"); 
  Output.insert(Output.end(), Temp.begin(), Temp.end()); 
}

int Table::LargestElement(std::vector<_TS> v)
{
  int max = 0; 
  for (_TS o : v)
  {
    if (max < o.Length()){ max = o.Length(); }
  }
  return max; 
}

_TS Table::Rounded(float l)
{
  _TS Out; 
  std::ostringstream p; 
  p.precision(Round); 
  if (Scientific){ p << std::scientific << l; }
  else if (Integer) { Out += (int(l)); return Out; }
  else { p << std::fixed << l; }
  Out += (p.str()); 
  return Out; 
}






void MergeTables::CompileTable()
{
  Output.push_back("================= Raw Values =================="); 
  Output.push_back(Outdir);

  TMP_Score = MakeEmptyVector(y_cellNames.size(), x_cellNames.size()); 
  values = MakeEmptyVector(y_cellNames.size(), x_cellNames.size()); 

  std::vector<std::vector<std::vector<float> > > Tensor; 
  for (int x(0); x < InputTables.size(); x++)
  {
    Table i = InputTables[x];
    _TS Title = z_cellNames[x]; 
    
    Output.push_back("-----> " + Title); 
    MakeHeader();
    for (y_i = 0; y_i < i.values.size(); y_i++){AppendLine(i.values);}

    i.BestValues(); 
    AddToTable(&TMP_Score, i.TMP_Score); 
    AddToTable(&values, i.values, 1.0/InputTables.size()); 
    Output.push_back(""); 
    Tensor.push_back(i.values); 
  }

  Output.push_back("------------ Best Values -------------"); 
  FindLowestErrorProjection(Tensor); 
  
  Output.push_back("");
  for (y_i = 0; y_i < y_cellNames.size(); y_i++)
  {
    int sum = 0; 
    for (x_i = 0; x_i < x_cellNames.size(); x_i++){sum += TMP_Score[y_i][x_i];}
    TMP_Score[y_i].push_back(sum); 
  }
  x_cellNames.push_back("Sum Scores"); 

  Output.push_back("--------- Accumuated Scores ---------"); 
  MakeHeader(); 
  Integer = true; 
  for (y_i = 0; y_i < TMP_Score.size(); y_i++){AppendLine(TMP_Score);}
  Integer = false;
  x_cellNames.pop_back();

  Output.push_back(""); 
  Output.push_back(""); 
  Output.push_back(""); 
  Output.push_back("========== Averaged Raw Values =========="); 
  MakeHeader();
  for (y_i = 0; y_i < values.size(); y_i++){AppendLine(values);}

  BestValues(); 
  ResultList(); 
  OtherStats();

  WriteToFile(); 


}

std::vector<_TS> MergeTables::FindLowestErrorProjection(std::vector<std::vector<std::vector<float>>> input)
{
  std::vector<_TS> TMP;  
  for (int x(0); x < x_cellNames.size(); x++)
  {
    std::vector<float> tmp_val; 
    std::vector<int> ind_z; 
    std::vector<int> ind_y; 

    for (int y(0); y < y_cellNames.size(); y++)
    {
      std::vector<float> tmp; 
      for (int z(0); z < z_cellNames.size(); z++){tmp.push_back(input[z][y][x]);}
      std::vector<int> z_index = LowestValue(tmp); 
      for (int i(0); i < z_index.size(); i++)
      {
        ind_z.push_back(z_index[i]); 
        ind_y.push_back(y); 
        tmp_val.push_back(input[z_index[i]][y][x]); 
      }
    }
    std::vector<int> zy_ind = LowestValue(tmp_val); 
    for (int k : zy_ind)
    {
      _TS out_str = x_cellNames[x]; 
      RepeatString(&out_str, x_width - out_str.Length(), " "); 
      out_str += Spacer; 
      
      out_str += y_cellNames[ind_y[k]];
      RepeatString(&out_str, LargestElement(y_cellNames) - y_cellNames[ind_y[k]].Length(), " "); 
      out_str += Spacer; 

      out_str += z_cellNames[ind_z[k]];
      RepeatString(&out_str, LargestElement(z_cellNames) - z_cellNames[ind_z[k]].Length(), " "); 
      out_str += Spacer; 

      out_str += "=> " + Rounded(tmp_val[k]); 
      TMP.push_back(out_str);
      Output.push_back(out_str); 
    }
  }
  return TMP; 
}
