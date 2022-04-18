#!/bin/bash 

cd $CODE_DIR/ShellCode
g++ ../ROOT/BestTruthShape.cxx ../ROOT/IO.cxx ../ROOT/SharedFunctions.cxx `root-config --libs --cflags` -o ../Binaries/BestTruthShape
cd $root_dir
