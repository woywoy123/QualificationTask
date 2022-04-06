#!/bin/bash 

cd $CODE_DIR/ShellCode
g++ ../ROOT/Debugging.cxx ../ROOT/IO.cxx ../ROOT/SharedFunctions.cxx `root-config --libs --cflags` -o ../Binaries/Debugging
cd $root_dir
