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

float WeightedShapeError(std::vector<TH1F*> Pred, std::vector<TH1F*> Truth)
{
  float n_int = 0; 
  for (TH1F* T : Truth){n_int += (T -> Integral());}

  float err = 0; 
  for (int i(0); i < Pred.size(); i++)
  {
    TH1F* P = (TH1F*)Pred[i] -> Clone("p1"); 
    TH1F* T = (TH1F*)Truth[i] -> Clone("p1"); 
    float m_int = T -> Integral();
    Normalize(P); 
    Normalize(T); 
    
    float err_t = 0; 
    for (int j(0); j < P -> GetNbinsX(); j++)
    {
      float sum = P -> GetBinContent(j+1) + T -> GetBinContent(j+1);
      if (sum == 0){continue;}
      float diff = P -> GetBinContent(j+1) - T -> GetBinContent(j+1); 
      
      err_t += (std::pow(diff, 2) / sum);
      
    }
    err += 0.5*err_t*m_int;
    delete P, T;
  }
  
  err = err/n_int;
  return err; 
}

void Normalize(TH1F* Hist)
{
  float x = Hist -> Integral();
  if (x == 0){return;}
  Hist -> Scale(1/x); 
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
      if ( L > v ){ L = v; } 
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
      if (v == 0){continue;}
      _TS Alg = y_cellNames[y_i]; 
 
      TMP += Alg; TMP += " (" + Rounded(v) + ") ";
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
    
    float L = -1; 
    float sum = 0; 
    for (x_i = 0; x_i < x_cellNames.size(); x_i++)
    {
      if (L == -1){ L = values[y_i][x_i]; }
      if (L > values[y_i][x_i]){ L = values[y_i][x_i]; }
      sum += values[y_i][x_i]; 
    }
    sum = sum / float(x_cellNames.size()); 

    Temp.push_back(TMP + Rounded(sum)); 

    for (x_i = 0; x_i < x_cellNames.size(); x_i++)
    {
      if (L == values[y_i][x_i])
      { 
        TMP += x_cellNames[x_i] + "(" + Rounded(L) + ") "; 
      }
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
