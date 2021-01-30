#include<PostAnalysis/BaseFunctionTest.h>
#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/PresentationFigures.h>
#include<PostAnalysis/AlgorithmTest.h>
#include<TFile.h>

int main(int argc, char** argv)
{
  bool rewrite = true; 
  int option = -5;

  TFile* F; 
  if (rewrite == true)
  {
    F = new TFile("output.root", "RECREATE");
    option = -2; 
  }
  else 
  {
    F = new TFile("output.root", "READ");
  }

  // Run Tests 
  if (option == -2)
  {
    //TestLandau(F); 
    //TestGaussian(F); 
    //TestGaussianXGaussian(F);
    //TestLandauXLandau(F); 
    //TestLandauXGaussian(F); 
    //TestDeconvGausXGaus(F); 
    //TestDeconvLandauXLandau(F); 
    //TestDeconvLandauXGaussian(F); 
    //TestGaussianDeconvolutionFit(F); 
    //TestLandauXGausFit(F); 
    //TestNLandauXNGausFit(F);
    //TestDeconvolutionFit(F);  
    //TestComparisonBinCenteringLandauXLandau(F); 
    //TestOscillationLucyRichardson(F); 
    //TestReadFile(F); 
    //TestReadFileTrackEnergy(F); 
    //TestMonteCarloMatchConvolution(F); 
    //TestMonteCarloFit(F); 

    TestAlgorithmMonteCarlo(); 
    DataAlgorithm(); 
  }

  // Run Experimental
  if (option == -3 && rewrite == false)
  {
    AlgorithmMonteCarlo();
  }
  
  if (option == -4)
  {
    // Compile the figures  
    FigureCompiler(F);   
  } 
  F -> Close();   
 
  if (option == -5)
  {
    //ProcessDataResults(); 
    ProcessMonteCarloResults(); 

  }
  std::cout << "fin" << std::endl;
  return 0; 
}
