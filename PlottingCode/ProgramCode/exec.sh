#!/bin/bash
source ~/.bashrc
g++ PostPlot.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o Plotter 
#./Plotter

g++ BestAlgo.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o BestAlgo
#./BestAlgo _Minimizer.root
#./BestAlgo BestAlgos.root
#./BestAlgo Blayer_800_1000_GeV_FitT_Normal_Subtract_Minimizer.root

g++ ClusterPlots.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o ClusterPlots
#./ClusterPlots Merged_MC_Negative.root

g++ Check.cxx `root-config --libs --cflags` -o Check
