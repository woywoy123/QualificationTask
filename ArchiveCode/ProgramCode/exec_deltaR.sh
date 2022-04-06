#!/bin/bash 
source ~/.bashrc
g++ deltaR.cxx BaseFunctions.cxx IO.cxx Plotting.cxx `root-config --libs --cflags` -o deltaR
mkdir Hists
./deltaR ../deltaR/Merger.root
