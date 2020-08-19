#include<TStyle.h>

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants
{
  const std::vector<TString> Detector = {"IBL", "Blayer", "layer1", "layer2"};
  const std::vector<TString> Pure_Names = {"dEdx_ntrk_1_ntru_1", 
                                           "dEdx_ntrk_2_ntru_2", 
                                           "dEdx_ntrk_3_ntru_3", 
                                           "dEdx_ntrk_4_ntru_4"};

  const std::vector<TString> trk_1 = {"dEdx_ntrk_1_ntru_1", 
                                      "dEdx_ntrk_1_ntru_2", 
                                      "dEdx_ntrk_1_ntru_3", 
                                      "dEdx_ntrk_1_ntru_4"};

  const std::vector<TString> trk_2 = {"dEdx_ntrk_2_ntru_1", 
                                      "dEdx_ntrk_2_ntru_2", 
                                      "dEdx_ntrk_2_ntru_3", 
                                      "dEdx_ntrk_2_ntru_4"};

  const std::vector<TString> trk_3 = {"dEdx_ntrk_3_ntru_1", 
                                      "dEdx_ntrk_3_ntru_2", 
                                      "dEdx_ntrk_3_ntru_3", 
                                      "dEdx_ntrk_3_ntru_4"};

  const std::vector<TString> trk_4 = {"dEdx_ntrk_4_ntru_1", 
                                      "dEdx_ntrk_4_ntru_2", 
                                      "dEdx_ntrk_4_ntru_3", 
                                      "dEdx_ntrk_4_ntru_4"};
  
  const std::vector<TString> Data_Names = {"dEdx_ntrk_1", 
                                           "dEdx_ntrk_2", 
                                           "dEdx_ntrk_3", 
                                           "dEdx_ntrk_4"};

  const std::vector<TString> FitNames = {"ntrk_1", 
                                         "ntrk_2", 
                                         "ntrk_3", 
                                         "ntrk_4"};



  const std::vector<Color_t> Colors = {kRed, kBlack, kBlue, kCyan, kTeal, kGreen, kSpring, kYellow, kGray, kOrange, kCoffee, kAurora, kMagenta, kViolet};
 
  const std::vector<float> COMP1 = {0.8, 0.1, 0.1, 0.1};
  const std::vector<float> COMP2 = {0.1, 0.89, 0.05, 0.05};
  const std::vector<float> COMP3 = {0.01, 0.2, 0.59, 0.2};
  const std::vector<float> COMP4 = {0.02, 0.2, 0.2, 0.58};
  const std::vector<float> LandauParameters = {1, 0.9, 0.1};

  const TString MC_dir = "/home/tnom6927/CTIDE/QualificationTask/PostAnalysisData/MonteCarlo/Merged.root";

}

#endif
