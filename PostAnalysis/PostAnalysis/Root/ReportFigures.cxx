#include<PostAnalysis/ReportFigures.h>
#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/AlgorithmFunctions.h>
#include<PostAnalysis/IO.h>

void Figure_Proxy()
{
  Figure_N_TrackTemplates();
  Figure_Trk1_Subtraction_Justification(); 

}






void Figure_N_TrackTemplates()
{

  std::map<TString, std::map<TString, std::vector<TH1F*>>> M = ReadCTIDE("Merged_MC.root");  
  TH1F* Blayer_trk1_outside = M[example]["ntrk_1_M_O"][0]; 
  TH1F* Blayer_trk2_outside = M[example]["ntrk_2_M_O"][0]; 
  TH1F* Blayer_trk3_outside = M[example]["ntrk_3_M_O"][0]; 
  TH1F* Blayer_trk4_outside = M[example]["ntrk_4_M_O"][0]; 
  
  std::vector<TH1F*> Blayer_trk1_O_M = {Blayer_trk1_outside, Blayer_trk2_outside, Blayer_trk3_outside, Blayer_trk4_outside};
  SubtractData(Blayer_trk1_O_M, Blayer_trk1_outside, 0, false); 
  std::vector<TH1F*> ntrk = BuildNtrkMtru(4, Blayer_trk1_outside, "T")[0];  

  for (int i(0); i < ntrk.size(); i++)
  {
    TString na; na += (i+1); na += ("-Track Template"); 
    ntrk[i] -> SetTitle(na);
  }
  
  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  GenerateNiceStacks(ntrk, "Example n-Track Template Distributions from Blayer at Jet PT 1200 - 1400 GeV", can, "dE/dx [MeV g^{-1} cm^{2}]", "", "nostack hist"); 
  can -> Print("ExampleTemplates.pdf"); 
}



void Figure_Trk1_Subtraction_Justification()
{
  std::map<TString, std::map<TString, std::vector<TH1F*>>> M = ReadCTIDE("Merged_MC.root");  
  TH1F* Blayer_trk1_outside = M[example]["ntrk_1_M_O"][0]; 
  TH1F* Blayer_trk2_outside = M[example]["ntrk_2_M_O"][0]; 
  TH1F* Blayer_trk3_outside = M[example]["ntrk_3_M_O"][0]; 
  TH1F* Blayer_trk4_outside = M[example]["ntrk_4_M_O"][0]; 

  std::vector<TH1F*> Blayer_trk1_O_Truth; 
  std::vector<TH1F*> Blayer_trk1_O_M; 
 
  std::vector<TH1F*> tmp_T = M[example]["ntrk_1_T_O"]; 
  std::vector<TH1F*> tmp = {Blayer_trk1_outside, Blayer_trk2_outside, Blayer_trk3_outside, Blayer_trk4_outside};
  for (int i(0); i < tmp_T.size(); i++)
  {
    TH1F* T = tmp_T[i]; 
    TH1F* M = tmp[i]; 
    if (M -> GetEntries() < 50){continue;}
    Blayer_trk1_O_Truth.push_back(T); 
    Blayer_trk1_O_M.push_back(M); 
    
    TString n_M = "dE/dx of "; n_M += (i+1); n_M += ("-Track Measured"); 
    M -> SetTitle(n_M); 

    n_M = "dE/dx of 1-Track,"; n_M += (i+1); n_M += ("-Truth"); 
    T -> SetTitle(n_M); 
  }

  TCanvas* can = new TCanvas(); 
  can -> SetLogy(); 
  Blayer_trk1_O_Truth[0] -> GetYaxis() -> SetTitle("Entries");
  PlotHists(Blayer_trk1_O_M, Blayer_trk1_O_Truth, "1-Track Truth Composition Superimposed with n-Track Cluster Measurements at \\DeltaR > 0.1", can);
  can -> Print("ExampleContamination.pdf"); 
  can -> Clear();
  //can -> Print("Output.pdf["); 
  
  std::map<TString, float> IBL_noSub; 
  std::map<TString, float> IBL_Sub;

  std::map<TString, float> Blayer_noSub; 
  std::map<TString, float> Blayer_Sub;

  std::map<TString, float> Layer1_noSub; 
  std::map<TString, float> Layer1_Sub;

  std::map<TString, float> Layer2_noSub; 
  std::map<TString, float> Layer2_Sub;

  // Normalization parameters
  std::map<TString, std::vector<float>> Params_N; 
  Params_N["Minimizer"] = {100000};

  for (MMVi x = M.begin(); x != M.end(); x++)
  {
    TString name = x -> first; 
    TH1F* trk1_inside = M[name]["ntrk_1_T_I"][0]; 
    TH1F* trk1_outside = M[name]["ntrk_1_M_O"][0]; 
    TH1F* trk2_outside = M[name]["ntrk_2_M_O"][0]; 
    TH1F* trk3_outside = M[name]["ntrk_3_M_O"][0]; 
    TH1F* trk4_outside = M[name]["ntrk_4_M_O"][0];
    Normalize(trk1_inside);
    TH1F* H1 = (TH1F*)trk1_outside -> Clone("H1"); 
    Normalize(H1);  

    std::vector<TH1F*> trk1_O_Truth = M[name]["ntrk_1_T_O"]; 
    std::vector<TH1F*> trk1_O_M = {trk1_outside, trk2_outside, trk3_outside, trk4_outside};
    if (trk1_outside -> GetEntries() < 20000){continue;}

    float err = ChiSquareError(trk1_inside, H1); 
    TString ibl = name; ibl = ibl.ReplaceAll("IBL_", ""); ibl = ibl.ReplaceAll("_", " ");
    TString bl = name; bl = bl.ReplaceAll("Blayer_", ""); bl = bl.ReplaceAll("_", " ");
    TString l1 = name; l1 = l1.ReplaceAll("Layer1_", ""); l1 = l1.ReplaceAll("_", " ");
    TString l2 = name; l2 = l2.ReplaceAll("Layer2_", ""); l2 = l2.ReplaceAll("_", " ");

    //PlotHists(trk1_O_M, trk1_O_Truth, name, can);
    //can -> Print("Output.pdf"); 
    //can -> Clear(); 

    if (name.Contains("IBL_")){  IBL_noSub[ibl] = err*100; }
    if (name.Contains("Blayer_")){ Blayer_noSub[bl] = err*100; }
    if (name.Contains("layer1_")){ Layer1_noSub[l1] = err*100; }
    if (name.Contains("layer2_")){ Layer2_noSub[l2] = err*100; }

    SubtractData(trk1_O_M, trk1_outside, 0, false); 
    //trk1_outside -> Add(trk2_outside, -1); 
    //trk1_outside -> Add(trk3_outside, -1); 
    //trk1_outside -> Add(trk4_outside, -1); 
    Normalize(trk1_outside); 
    
    float err_s = ChiSquareError(trk1_inside, trk1_outside); 
    
    if (name.Contains("IBL_")){ IBL_Sub[ibl] = err_s*100; }
    if (name.Contains("Blayer_")){ Blayer_Sub[bl] = err_s*100; }
    if (name.Contains("layer1_")){ Layer1_Sub[l1] = err_s*100; }
    if (name.Contains("layer2_")){ Layer2_Sub[l2] = err_s*100; }
      
    std::cout << err_s << " " << err << " " << name << " " << trk1_outside -> GetName() << std::endl;
    
    delete H1; 

  }
  
  TMultiGraph* Gr = new TMultiGraph(); 
  
  //can -> Print("Output.pdf]"); 
  can -> Clear(); 
  can -> SetLogy(false);
  
  can -> Print("ShapeToInsideTruth.pdf["); 

  TGraph* IBL = GenerateGraph(IBL_noSub, "IBL Track Jet-PT Ranges");
  TGraph* IBL_sub = GenerateGraph(IBL_Sub, "IBL Track Jet-PT Ranges (Subtract)"); 

  TGraph* Blayer = GenerateGraph(Blayer_noSub, "Blayer Track Jet-PT Ranges");
  TGraph* Blayer_sub = GenerateGraph(Blayer_Sub, "Blayer Track Jet-PT Ranges (Subtract)"); 
 
  TGraph* Layer1 = GenerateGraph(Layer1_noSub, "Layer1 Track Jet-PT Ranges");
  TGraph* Layer1_sub = GenerateGraph(Layer1_Sub, "Layer1 Track Jet-PT Ranges (Subtract)"); 
  
  TGraph* Layer2 = GenerateGraph(Layer2_noSub, "Layer2 Track Jet-PT Ranges");
  TGraph* Layer2_sub = GenerateGraph(Layer2_Sub, "Layer2 Track Jet-PT Ranges (Subtract)"); 
  
  IBL -> SetLineColor(kBlue); 
  IBL_sub -> SetLineColor(kRed); 

  Blayer -> SetLineColor(kBlue); 
  Blayer_sub -> SetLineColor(kRed); 

  Layer1 -> SetLineColor(kBlue); 
  Layer1_sub -> SetLineColor(kRed); 

  Layer2 -> SetLineColor(kBlue); 
  Layer2_sub -> SetLineColor(kRed); 

  std::vector<TGraph*> ng = {IBL, IBL_sub}; 
  IBL -> GetYaxis() -> SetRangeUser(0., 0.5); 
  IBL -> Draw("ALPSAME");
  IBL_sub -> Draw("SAME");
  GenerateLegend(ng, can); 
  can -> Print("ShapeToInsideTruth.pdf");
  can -> Clear(); 

  ng = {Blayer, Blayer_sub}; 
  Blayer -> GetYaxis() -> SetRangeUser(0., 0.5); 
  Blayer -> Draw("ALPSAME");
  Blayer_sub -> Draw("SAME");
  GenerateLegend(ng, can); 
  can -> Print("ShapeToInsideTruth.pdf");
  can -> Clear(); 

  ng = {Layer1, Layer1_sub}; 
  Layer1 -> GetYaxis() -> SetRangeUser(0., 0.5); 
  Layer1 -> Draw("ALPSAME");
  Layer1_sub -> Draw("SAME");
  GenerateLegend(ng, can); 
  can -> Print("ShapeToInsideTruth.pdf");
  can -> Clear(); 

  ng = {Layer2, Layer2_sub}; 
  Layer2 -> GetYaxis() -> SetRangeUser(0., 0.5); 
  Layer2 -> Draw("ALPSAME");
  Layer2_sub -> Draw("SAME");
  GenerateLegend(ng, can); 
  can -> Print("ShapeToInsideTruth.pdf");

  can -> Print("ShapeToInsideTruth.pdf]"); 










}
