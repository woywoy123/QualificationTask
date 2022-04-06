#!/bin/bash
source ~/.bashrc
g++ Debug.cxx Evaluate.cxx BaseFunctions.cxx IO.cxx Plotting.cxx `root-config --libs --cflags` -o Debug
#./Debug _Debug.root
