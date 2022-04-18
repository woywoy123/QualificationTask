#!/bin/bash 

cd $CODE_DIR/ShellCode
g++ ../ROOT/TestToDataShape.cxx ../ROOT/IO.cxx ../ROOT/SharedFunctions.cxx `root-config --libs --cflags` -o ../Binaries/TestToDataShape
cd $root_dir
