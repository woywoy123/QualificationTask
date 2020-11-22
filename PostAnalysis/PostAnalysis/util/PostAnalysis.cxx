#include<PostAnalysis/BaseFunctionTest.h>
#include<PostAnalysis/Experimental.h>
#include<PostAnalysis/PresentationFigures.h>
#include<TFile.h>

int main(int argc, char** argv)
{
  bool rewrite = true; 
  int option = -2;

  TFile* F; 
  if (rewrite == true)
  {
    F = new TFile("output.root", "RECREATE"); 
  }
  else 
  {
    F = new TFile("output.root", "READ"); 
  }

  // Run Tests 
  if (option == -1)
  {
    //PlotLandau(); 
    //PlotGaussian(); 
    //PlotGaussianXGaussian();
    //PlotLandauXLandau(); 
    //PlotLandauXGaussian(); 
    //PlotDeconvGausXGaus(); 
    //PlotDeconvLandauXLandau(); 
    //PlotDeconvLandauXGaussian(); 
    //PlotGaussianDeconvolutionFit(); 
    //PlotLandauXGausFit(); 
    //PlotNLandauXNGausFit();
    //PlotDeconvolutionFit();  
  }

  // Produce Presentation Figures
  if (option == -2)
  {
    GaussianXGaussian(F); 
    LandauXLandau(); 
  }

  // Run Experimental
  if (option == -3)
  {
  }
 
  // Compile the figures  
  FigureCompiler(&(*F));   
 
 
 
 
 
 
 
  
  F -> Close();   
  
  std::cout << "fin" << std::endl;
  return 0; 

}
