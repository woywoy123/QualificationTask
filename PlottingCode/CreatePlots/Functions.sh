#!/bin/bash

MakeHistDir () {
  local i=$1
  echo "-> $i"
  mkdir $i/ROOT
  trks=("trk1" "trk2" "trk3" "trk4")
  mkdir $i/Histograms
  mkdir $i/Histograms/Test
  mkdir $i/Histograms/Template
  mkdir $i/ShapePerformance
  mkdir $i/FLost
  for j in ${trks[@]}
  do
    mkdir $i/Histograms/Test/$j
    mkdir $i/Histograms/Template/$j
    mkdir $i/ShapePerformance/$j
    for l in ${Algs[@]}
    do
      mkdir $i/Histograms/Test/$j/$l 
      mkdir $i/Histograms/Template/$j/$l
    done
  done
}

BuildDirectory () {
  mkdir Input
  rm -rf Output
  mkdir Output
  cd Output
  
  for i in ${Modes[@]}
  do
    mkdir $i
    MakeHistDir $i
  done 
  cd ../
}

CollectROOT () {
  cd Input
  for i in "${Algs[@]}"
  do
    for k in ${Modes[@]}
    do
      for j in *
      do
        str="FitT_"$i$k
        if [[ "$j" == *"$str" || "$j" == *"Truth" ]]
        then 
          cp $j/Fit_Tracks.root ../Output/$k/ROOT/$j.root
          echo "---> $j"
        fi
      done
    done
  done
  cd ../
}

MergeROOT () {
  local str=$1
  /home/tnom6927/ROOT/bin/hadd MultiTrackFit.root $str
  mv MultiTrackFit.root ../
}

CollimateROOT () {
  for i in *
  do
    cd $i/ROOT
    local str=""
    for j in *
    do
      str="$str $j"
    done
    MergeROOT "$str"
    cd ../../
  done
  cd ../
}

CompilePostPlot () {
  g++ PostPlot.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o Plotter
}

CompileBestAlgo () {
  g++ BestAlgo.cxx BaseFunctions.cxx IO.cxx Evaluate.cxx Plotting.cxx `root-config --libs --cflags` -o BestAlgo
}

CleanDebug () {
  rm -r FLost 
  rm -r ROOT
  rm -r Histograms
  rm -r ShapePerformance
  rm -r *.pdf
}
