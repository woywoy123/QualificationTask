#include<PostAnalysis/BaseFunctions.h>
#include<PostAnalysis/DerivedFunctions.h>
#include<PostAnalysis/Constants.h>
#include<PostAnalysis/Plotting.h>
#include<PostAnalysis/UnitTest.h>
#include<stdio.h>

int main(int argc, char *argv[])
{

  DerivedFunctions DF;

  int iter = 50;
  int cor_loop = 40; // Correction loop number 

  // ===== Good parameters that have been tested (out.root)
  //Gaussian Parameter used for deconvolution
	std::map<TString, std::vector<float>> Params;
  Params["Gaussian"] = {0, 6};
  Params["m_e"] = {1, 1, 1, 1, 1, 1};
  Params["m_s"] = {0, 0, 0, 0, 0, 0};        
  Params["s_s"] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  Params["s_e"] = {15, 15, 15, 15, 15, 15};
	int bins = 500; 
	int max = 20; 
	int min = 0; 
	float offset = 0.2; 

	TString Layer = argv[1]; 
	TString Energy = argv[2]; 

	std::cout << Layer << std::endl;

	auto Generate =[](std::map<TString, TH1F*> &Dic, std::vector<TString> Detector, std::vector<std::vector<TString>> Batch, std::vector<TString> Names, TString Dir, std::vector<TString> ntrk, int bins, float min, float max)
  {
    std::map<TString, int> Test; 
    TFile* File = new TFile(Dir); 
    for (int i(0); i < ntrk.size(); i++)
    {
      TString trk = ntrk[i] + "_All";
      Dic[trk] = new TH1F(trk, trk, bins, min, max); 
      Test[trk] = 0; 
      for (TString v : Names)
      {
        TString trk_all = trk + "_" + v;
        Dic[trk_all] = new TH1F(trk_all, trk_all, bins, min, max); 
        Test[trk_all] = 0; 
      }

      for (int v(0); v < Detector.size(); v++)
      {
        TString Lay = Detector[v]; 
        TString trk_L = ntrk[i] + "_" + Lay; 
        Dic[trk_L] = new TH1F(trk_L, trk_L, bins, min, max);
        Test[trk_L] = 0;  
        for (int x(0); x < Batch.size(); x++)
        {
          std::vector<TString> Bat = Batch[x];
          TString BN = Names[x]; 
          TString trk_bn = trk_L + "_" + BN;
          TString trk_all = trk + "_" + BN;
          Dic[trk_bn] = new TH1F(trk_bn, trk_bn, bins, min, max);  
          Test[trk_bn] = 0; 
          for (int b(0); b < Bat.size(); b++)
          {
            TString dir = Lay + Bat[b]; 
            File -> cd(dir);  
            TH1F* H = (TH1F*)gDirectory -> Get(ntrk[i]);
            Dic[trk_all] -> Add(H); 
            Dic[trk] -> Add(H); 
            Dic[trk_L] -> Add(H); 
            Dic[trk_bn] -> Add(H); 

            Test[trk_all]++;
            Test[trk]++;
            Test[trk_L]++;
            Test[trk_bn]++;

            File -> cd(); 
          }
        }
      }
    }
  };
 

	std::vector<TString> Detector_Layer = {"IBL", "Blayer", "layer1", "layer2"};
  std::vector<TString> E = Constants::energies;
  std::vector<std::vector<TString>> Batch = {{E[0]}, 
                                             {E[1], E[2]}, 
                                             {E[3], E[4], E[5]}, 
                                             {E[6], E[7], E[8], E[9], E[10], E[11], E[12], E[13], E[14], E[15]}};
  std::vector<TString> Batch_Names = {"200", "200-600", "600-1200", "1200+"};

  // Create the histograms for both the MC and the Data. The Data wont have truth but an empty set of histograms 
  std::map<TString, TH1F*> Data_MC;
  for (int i(0); i < Constants::All_MC.size(); i++)
  {
    std::vector<TString> ntrk = Constants::All_MC[i]; 
    Generate(Data_MC, Detector_Layer, Batch, Batch_Names, Constants::MC_dir, ntrk, bins, min, max);     
  }
  Generate(Data_MC, Detector_Layer, Batch, Batch_Names, Constants::MC_dir, Constants::Data_Names, bins, min, max);    

  Detector_Layer.push_back("All");

  // Iterate through the dictionary and build the data and truth sets used for the iteration 
  std::map<TString, std::vector<std::vector<TH1F*>>> Closure_Dic;
  std::map<TString, std::vector<TH1F*>> Data_Dic; 
  
  for (TString L : Detector_Layer)
  { 
    for (TString BN : Batch_Names)
    {
      for (std::vector<TString> tru : Constants::All_MC)
      {
        std::vector<TH1F*> truth;
        for (TString t : tru)
        { 
          TString condition = t + "_" + L + "_" + BN;
          std::map<TString, TH1F*>::iterator it;
          for (it = Data_MC.begin(); it != Data_MC.end(); it++)
          {
            TString name = it -> second -> GetName();
            if (name == condition)
            {
              truth.push_back(it -> second);
            }
          }
        }
        Closure_Dic[L + "_" + BN].push_back(truth); 
      }

      for (TString d : Constants::Data_Names)
      {
        TString condition = d + "_" + L + "_" + BN;
        std::map<TString, TH1F*>::iterator it;
        for (it = Data_MC.begin(); it != Data_MC.end(); it++)
        {
          TString name = it -> second -> GetName();
          if (name == condition)
          {
            Data_Dic[L + "_" + BN].push_back(it -> second);
          }
        }
      }
    }

    for (TString d : Constants::Data_Names)
    {
      TString condition = d + "_" + L;
      std::map<TString, TH1F*>::iterator it;
      for (it = Data_MC.begin(); it != Data_MC.end(); it++)
      {
        TString name = it -> second -> GetName();
        if (name == condition)
        {
          Data_Dic[L].push_back(it -> second); 
        }
      } 
    }

    for (std::vector<TString> d : Constants::All_MC)
    {
      std::vector<TH1F*> truth; 
      for (TString x : d)
      {
        TString condition = x + "_" + L;
        std::map<TString, TH1F*>::iterator it;
        for (it = Data_MC.begin(); it != Data_MC.end(); it++)
        {
          TString name = it -> second -> GetName();
          if (name == condition)
          {
            truth.push_back(it -> second); 
          }
        } 
      }
      Closure_Dic[L].push_back(truth); 
    }
  }

	TString Title = Layer + "_" + Energy + ".root"; 

	TString Execute = "";
	if (Energy == "")
	{
		Execute = Layer; 	
	}
	else 
	{
		Execute = Layer + "_" + Energy; 
	}

	std::vector<TH1F*> ntrk_Data = Data_Dic[Execute];
	std::vector<std::vector<TH1F*>> Truth_Sets = Closure_Dic[Execute];
	
  TFile Output(Title, "UPDATE"); 
  std::map<int, std::pair<TH1F*, std::vector<TH1F*>>> res =  DF.MainAlgorithm(ntrk_Data, Params, offset, iter, cor_loop, Truth_Sets); 
      
  Output.mkdir(Execute);
  for (int p(0); p < res.size(); p++)
  {
    Output.cd(Execute);   
    TString n = "Iteration_"; n+=(p); 
    gDirectory -> mkdir(n); 
    gDirectory -> cd(n); 
    TH1F* FLost = res[p].first;
    FLost -> Write(); 
    std::vector<TH1F*> Others = res[p].second; 
    for (TH1F* H_R : Others)
    {
      H_R -> Write(); 
    }
    gDirectory -> cd(); 
  }
  Output.Close(); 

	return 0; 
}
