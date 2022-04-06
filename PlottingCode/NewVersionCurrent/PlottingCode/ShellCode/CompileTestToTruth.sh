#!/bin/bash 

cd $CODE_DIR/ShellCode
g++ ../ROOT/TestToTruthShape.cxx ../ROOT/IO.cxx ../ROOT/SharedFunctions.cxx `root-config --libs --cflags` -o ../Binaries/TestToTruthShape
cd $root_dir
