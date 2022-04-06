#!/bin/bash 

source ~/.bashrc
cd ../
source CreatePlots/Algorithms.sh
source CreatePlots/Functions.sh
#BuildDirectory
#CollectROOT
#cd Output/
#CollimateROOT
#cp ./Output/_FitTo/MultiTrackFit.root ProgramCode/
cd ProgramCode
#CompilePostPlot
CompileBestAlgo
cd ../
cd Output
#
#it=0
#for i in ${Modes[@]}
#do
#  cd $i
#  cp ../../ProgramCode/Plotter .
#  ./Plotter &
#
#  if [[ $it -gt 8 ]]
#  then 
#    it=0
#    wait
#  fi
#  it=$((it+1))
#  cd ../
#done
#wait

for i in ${Modes[@]}
do
  cd $i
  rm $i.root
  cp ../../ProgramCode/BestAlgo .
  ./BestAlgo $i &
  cd ../
done
wait
cd ../

BuildFinalDirectory "Blayer_800_1000_GeV" "_FitTo_Range"
CollectFinalPlots "Blayer_800_1000_GeV" "_FitTo_Range"



