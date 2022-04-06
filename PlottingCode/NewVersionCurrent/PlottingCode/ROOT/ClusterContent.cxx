#include "../PlottingCode/ClusterContent.h"

int main(int argc, char* argv[])
{
  TString Directory = argv[1];
  std::map<_TS, std::vector<TH1F*>> Content = ReadCTIDE(Directory);
  std::map<_TS, std::vector<float>> ClusterSize_IBL; 
  std::map<_TS, std::vector<float>> ClusterSize_Blayer; 
  std::map<_TS, std::vector<float>> ClusterSize_layer1; 
  std::map<_TS, std::vector<float>> ClusterSize_layer2; 
  
  float m_IBL = 0; 
  float m_Blayer = 0; 
  float m_layer1 = 0; 
  float m_layer2 = 0; 
  for (MVTFi it = Content.begin(); it != Content.end(); it++)
  {
    if (!(it -> first).Contains("InsideData")){continue;}
    std::vector<TH1F*> Tracks = it -> second; 
    
    std::vector<_TS> val = Split((it -> first), "/"); 
    _TS Layer = val[0]; 
    _TS Energy = val[1]; 
    _TS ntrk = val[2];

    for (TH1F* T : Tracks)
    {
      if (Layer == "IBL")
      { 
        ClusterSize_IBL[Energy].push_back( T -> Integral() ); 
        if (m_IBL < T -> Integral()){ m_IBL = T -> Integral(); }
      }

      if (Layer == "Blayer")
      { 
        ClusterSize_Blayer[Energy].push_back( T -> Integral() ); 
        if (m_Blayer < T -> Integral()){ m_Blayer = T -> Integral(); }
      }

      if (Layer == "layer1")
      { 
        ClusterSize_layer1[Energy].push_back( T -> Integral() ); 
        if (m_layer1 < T -> Integral()){ m_layer1 = T -> Integral(); }
      }

      if (Layer == "layer2")
      { 
        ClusterSize_layer2[Energy].push_back( T -> Integral() ); 
        if (m_layer2 < T -> Integral()){ m_layer2 = T -> Integral(); }
      }

    }
  }
  
  int mul = 100; 
  TCanvas* can = new TCanvas(); 
  can -> SetLogy();
  gStyle -> SetOptStat(0); 
  gStyle -> SetImageScaling(3); 
  can -> SetTopMargin(0.1); 
  
  std::vector<TGraph*> gIBL = PlotMultiGraph(ClusterSize_IBL, "Recorded Clusters in the IBL for n-Tracks", can, 4, 1, mul*m_IBL);
  can -> Print("IBL.png");
  can -> Clear();
  BulkDelete(gIBL);

  std::vector<TGraph*> gBlayer = PlotMultiGraph(ClusterSize_Blayer, "Recorded Clusters in the Blayer for n-Tracks", can, 4, 1, mul*m_Blayer);
  can -> Print("Blayer.png");
  can -> Clear();
  BulkDelete(gBlayer);
  
  std::vector<TGraph*> glayer1 = PlotMultiGraph(ClusterSize_layer1, "Recorded Clusters in the layer1 for n-Tracks", can, 4, 1, mul*m_layer1);
  can -> Print("layer1.png");
  can -> Clear();
  BulkDelete(glayer1);
  
  std::vector<TGraph*> glayer2 = PlotMultiGraph(ClusterSize_layer2, "Recorded Clusters in the layer2 for n-Tracks", can, 4, 1, mul*m_layer2);
  can -> Print("layer2.png");
  can -> Clear();
  BulkDelete(glayer2);
  return 0; 
}
