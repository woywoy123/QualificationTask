#include<TStyle.h>

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants
{
  const std::vector<TString> Detector = {"IBL", "Blayer", "layer1", "layer2"};
  const std::vector<TString> Pure_Names = {"dEdx_ntrk_1_ntru_1", 
                                           "dEdx_ntrk_2_ntru_2", 
                                           "dEdx_ntrk_3_ntru_3", 
                                           "dEdx_ntrk_4_ntru_4",
                                           "dEdx_ntrk_5_ntru_5"};

  const std::vector<TString> trk_1 = {"dEdx_ntrk_1_ntru_1", 
                                      "dEdx_ntrk_1_ntru_2", 
                                      "dEdx_ntrk_1_ntru_3", 
                                      "dEdx_ntrk_1_ntru_4",
                                      "dEdx_ntrk_1_ntru_5",
                                      "dEdx_ntrk_1_ntru_6"};

  const std::vector<TString> trk_2 = {"dEdx_ntrk_2_ntru_1", 
                                      "dEdx_ntrk_2_ntru_2", 
                                      "dEdx_ntrk_2_ntru_3", 
                                      "dEdx_ntrk_2_ntru_4",
                                      "dEdx_ntrk_2_ntru_5",
                                      "dEdx_ntrk_2_ntru_6"};

  const std::vector<TString> trk_3 = {"dEdx_ntrk_3_ntru_1", 
                                      "dEdx_ntrk_3_ntru_2", 
                                      "dEdx_ntrk_3_ntru_3", 
                                      "dEdx_ntrk_3_ntru_4",
                                      "dEdx_ntrk_3_ntru_5",
                                      "dEdx_ntrk_3_ntru_6"};

  const std::vector<TString> trk_4 = {"dEdx_ntrk_4_ntru_1", 
                                      "dEdx_ntrk_4_ntru_2", 
                                      "dEdx_ntrk_4_ntru_3", 
                                      "dEdx_ntrk_4_ntru_4",
                                      "dEdx_ntrk_4_ntru_5",
                                      "dEdx_ntrk_4_ntru_6"};

  const std::vector<TString> Data_Names = {"dEdx_ntrk_1", //dEdx_out1_ntrk1_calib", 
                                           "dEdx_ntrk_2", 
                                           "dEdx_ntrk_3", 
                                           "dEdx_ntrk_4",
                                           "dEdx_ntrk_5"};

  const std::vector<TString> FitNames = {"ntrk_1", 
                                         "ntrk_2", 
                                         "ntrk_3", 
                                         "ntrk_4",
                                         "ntrk_5",
																				 "ntrk_6"};

  const std::vector<TString> energies = {"/200_up_GeV/", "/200_400_GeV/", 
                                         "/400_600_GeV/", "/600_800_GeV/", 
                                         "/800_1000_GeV/", "/1000_1200_GeV/", 
                                         "/1200_1400_GeV/", "/1400_1600_GeV/", 
                                         "/1600_1800_GeV/", "/1800_2000_GeV/", 
                                         "/2000_2200_GeV/", "/2200_2400_GeV/", 
                                         "/2400_2600_GeV/", "/2600_2800_GeV/", 
                                         "/2800_3000_GeV/", "/higher_GeV/"}; 

  const std::vector<std::vector<TString>> All_MC = {trk_1, trk_2, trk_3, trk_4}; 

  const std::vector<Color_t> Colors = {kRed, kGreen, kBlue, kCyan, kViolet, kOrange, kCoffee, kAurora, kMagenta, kViolet};
 
  const std::vector<float> COMP1 = {0.997, 0.001, 0.001, 0.001, 0.001};
  const std::vector<float> COMP2 = {0.25 , 0.60 , 0.05 , 0.05 , 0.05 };
  const std::vector<float> COMP3 = {0.05 , 0.15 , 0.6  , 0.15 , 0.05 };
  const std::vector<float> COMP4 = {0.05 , 0.05 , 0.2  , 0.5  , 0.2  };
  const std::vector<float> COMP5 = {0.005 , 0.005, 0.19 , 0.2  , 0.6  };

  const std::vector<float> LandauParameters = {1, 0.9, 0.1};
  const int GaussianToys = 500000;

  const TString MC_dir = "/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/MonteCarlo/Merged.root";
  const TString Data2016 = "/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/Data/2016.root";
  const TString Data2017 = "/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/Data/2017.root";
  const TString Data2018 = "/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/Data/2018.root";
}

#endif
