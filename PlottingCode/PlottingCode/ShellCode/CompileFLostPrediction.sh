#!/bin/bash 

cd $CODE_DIR/ShellCode
g++ ../ROOT/BestFLost.cxx ../ROOT/IO.cxx ../ROOT/SharedFunctions.cxx `root-config --libs --cflags` -o ../Binaries/BestFLost
cd $root_dir
