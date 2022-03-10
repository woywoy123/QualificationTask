# Qualification Task - FLost2 and 3 Fit 
## Introduction:
This Qualification Task aims to improve the pre-existing approach of quantifying the fraction of 2 and 3-Tracks being lost during the Inner Detector reconstruction phase.
To begin with the Analysis, a sample needs to created using the CTIDE codebase, which outputs a ROOT file containing the 1, 2, 3 and 4-Track measurements for the IBL, B-Layer, Layer-1 and Layer-2 at different Jet energies. 
For the code to perform as expected, the CTIDE algorithm needs to be run over Monte Carlo samples to get the true 1, 2, 3 and 4-Track distributions, which are required for the code to perform a closure test.
Alternatively, the code can be modified to exclude truth information and simply run over track measurements, but without any closure. 

The aim of this project is to accurately generate dE/dx Track templates for different track multiplicities and subsequently fit them to a given measurement distribution. 
The technique employed in this code involves reading the 1-Track measurements and performing an n-Fold convolution of itself to generate initial 'approximate templates'. 
These are subsequently fitted using a variety of algorithms ranging in complexity. 
The most simplistic fit is a 'Normalization' fit, followed by an algorithm which takes template shifting into account, and finally more complex approaches involving Gaussian convolutions. 

The code is structured in two parts, within the 'PostAnalysis/PostAnalysis/Data' directory, the target sample is inserted and the bash script 'QuickSpawn.sh' initiates the fitting scripts. 
Within this script, the algorithm approach can be edited and executed in bulk. 
It should be noted, the spawning can be done in two way;
1. Running on CONDOR -> Once finished an 'example.submit' file is dumped.
2. Locally (Not recommended) -> This will create 'Spawn.sh' files which can be inidividually executed.

The second part of the code is concerned with the plotting of the results. 
Scripts to initialize the process can be found under the 'PlottingCode' directory.

## To Run:
### Creating Templates
- Naviate to the 'Data' directory within 'PostAnalysis' 
- Place the CTIDE output ROOT file into this directory
- Edit 'QuickSpawn.sh' accordingly: e.g. Filename, Algorithms, Jet energy, Layers etc.
- Execute 'QuickSpawn.sh' and wait.

### Analyse Template
- Execute 'MakeDir.sh' in 'CreatePlots' within the 'PlottingCode' directory.
- Place the output of the first part into the 'Input' folder.
- Execute 'Run.sh'



