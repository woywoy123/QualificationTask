#!/bin/bash 

root_dir=$PWD
cd ../
OUTPUT_DIR=$root_dir/Output
CODE_DIR=$PWD/PlottingCode
INPUT_DIR=$PWD/Input
ROOT_MERGE=/home/tnom6927/ROOT/bin/hadd
cd $root_dir 

source ~/.bashrc
source $CODE_DIR/ShellCode/Functions.sh
source $CODE_DIR/ShellCode/Algorithms.sh

#CleanWorkspace
mkdir $CODE_DIR/Binaries
mkdir $root_dir/Output

#source $CODE_DIR/ShellCode/CompileClusterContent.sh
#source $CODE_DIR/ShellCode/CompileDebugging.sh
source $CODE_DIR/ShellCode/CompileTestToTruth.sh

MakeOutputDir
#CollectCommonROOT
#MakeClusterContent
#MakeDebuggingContent
MakeTestToTruthShape

