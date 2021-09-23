#!/bin/bash
source ~/.bashrc
g++ PostPlot.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o Plotter 
#./Plotter

g++ BestAlgo.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o BestAlgo
#./BestAlgo _Subtract_Minimizer_Smooth
#./BestAlgo BestAlg.root

