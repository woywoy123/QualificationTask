#!/bin/bash
MergeROOT () {
  local str=$1
  local dir=$2
  cd $dir
  /home/tnom6927/ROOT/bin/hadd MultiTrackFit.root $str
  mv MultiTrackFit.root ../
}

Modes=( "_Debug" "_Debug_Smooth" "_Debug_Subtract_Smooth" "_Debug_Subtract")

root_dir=$PWD
rm -r Output_Debug
rm -r Final_Debug
mkdir Output_Debug
cd Output_Debug

for i in ${Modes[@]}
do
  echo "-> $i"
  mkdir $i
  mkdir "$i/ROOT"
done 

cd $root_dir
cd ProgramCode
bash exec.sh
cd ../

cd Input_Debug
for k in ${Modes[@]}
do 
  for j in *
  do
    if [[ "$j" == *"$k" ]]
    then 
      cp $root_dir/Input_Debug/$j/*.root $root_dir/Output_Debug/$k/ROOT/$j.root &
    fi
  done
done 
wait

cd $root_dir
cd Output_Debug

source ~/.bashrc
for i in *
do
  mkdir $i/IllustrationOfCases
  mkdir $i/IllustrationOfCases/Case1
  mkdir $i/IllustrationOfCases/Case2
  mkdir $i/IllustrationOfCases/Case3
 
  mkdir $i/Steps
  mkdir $i/Steps/Case1
  mkdir $i/Steps/Case2
  mkdir $i/Steps/Case3
  
  mkdir $i/RatioMatrix
  mkdir $i/RatioMatrix/IBL
  mkdir $i/RatioMatrix/Blayer
  mkdir $i/RatioMatrix/layer1
  mkdir $i/RatioMatrix/layer2

  mkdir $i/Excess

  cd $i/ROOT/
  str=""
  for j in *
  do
    str="$str $j"
  done
  
  MergeROOT "$str" "$PWD"
  cd $root_dir/Output_Debug
done
cd $root_dir

cd ProgramCode
bash exec_debug.sh
cd ../

cd Output_Debug
for i in *
do
  cd $i 
  cp ../../ProgramCode/Debug ./
  chmod +x Debug 
  ./Debug MultiTrackFit.root &
  cd ../
done 
wait

