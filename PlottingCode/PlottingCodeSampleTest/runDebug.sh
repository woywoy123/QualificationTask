#!/bin/bash 

source ~/.bashrc
cd ../
echo $PWD
source CreatePlots/Algorithms.sh
source CreatePlots/Functions.sh
BuildDirectory
CollectROOT
cd Output/
CollimateROOT
cp ./Output/_FitTo_Range/MultiTrackFit.root ProgramCode/
cd ProgramCode
CleanDebug
CompilePostPlot
MakeHistDir "."
./Plotter
exit
CompileBestAlgo

