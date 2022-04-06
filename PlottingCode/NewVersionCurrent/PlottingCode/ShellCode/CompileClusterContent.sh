#!/bin/bash 

cd $CODE_DIR/ShellCode
g++ ../ROOT/ClusterContent.cxx ../ROOT/IO.cxx ../ROOT/SharedFunctions.cxx `root-config --libs --cflags` -o ../Binaries/ClusterContent
cd $root_dir
