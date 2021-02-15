#include<PostAnalysis/BaseFunctionTest.h>
#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/PresentationFigures.h>
#include<PostAnalysis/AlgorithmTest.h>
#include<TFile.h>
#include<PostAnalysis/Debugging.h>

int main(int argc, char** argv)
{
<<<<<<< Updated upstream
  bool rewrite = true; 
<<<<<<< Updated upstream
  int option = -3;
=======
  int option = -5;
=======
  bool rewrite = false; 
  int option = -6;
>>>>>>> Stashed changes
>>>>>>> Stashed changes

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
    //DataAlgorithm(); 
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
 
  if (option == -2)
  {
<<<<<<< Updated upstream
    //ProcessDataResults(); 
    ProcessMonteCarloResults(); 
=======
<<<<<<< Updated upstream
    ProcessDataResults(); 
    //ProcessMonteCarloResults(); 
>>>>>>> Stashed changes

=======
    //ProcessDataResults(); 
    ProcessMonteCarloResults(); 
>>>>>>> Stashed changes
  }
  
  if (option == -6)
  {
    Entry();  
  }
  
  
  
  
  
  std::cout << "fin" << std::endl;
  return 0; 
}
