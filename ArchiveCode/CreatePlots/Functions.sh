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
  for i in ${Algs[@]}
  do
    for k in ${Modes[@]}
    do
      for j in *
      do
        str="FitT_"$i$k"_ntrk"
        
        if [[ "$j" == "Truth" ]]
        then 
          cp $j/Fit_Tracks.root ../Output/$k/ROOT/$j.root
          echo "---> $j $str"
        fi
        
        if [[ "$j" == *"FitT_$i$k" && "$i" == "Experimental" ]]
        then
          cp $j/Fit_Tracks.root ../Output/$k/ROOT/$j.root
          echo "---> $j"
        fi

        if [[ "$j" != *"$str"* ]]
        then 
          continue
        fi

        cp $j/Fit_Tracks.root ../Output/$k/ROOT/$j.root
        echo "---> $j"
      done
    done
  done
  cd ../
}

MergeROOT () {
  local col=$1
  local name=$2

  /home/tnom6927/ROOT/bin/hadd $name.root $col
}

CollimateROOT () {
  for i in *
  do
    cd $i/ROOT

    for l in ${Layer[@]}
    do
      local str=""
      for j in *
      do
        if [[ "$j" != *"$l"* ]]
        then 
          continue
        fi
        str="$str $j"
      done
      MergeROOT "$str" "$l"
    done
    local str=""
    for l in ${Layer[@]}
    do
      str="$str $l.root"
    done
    str="$str Truth.root"
    MergeROOT "$str" "MultiTrackFit"
    mv MultiTrackFit.root ../

    cd ../../
  done
  cd ../
  wait
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

BuildFinalDirectory() {
  local energy=$1
  local mini=$2
  root_dir=$PWD

  rm -rf Final
  mkdir Final 
  mkdir Final/IllustrationsOfAlgorithmFits
  mkdir Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth
  mkdir Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth/Nominal-NoSubtraction-NoSmoothing
  ill=Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth/Nominal-NoSubtraction-NoSmoothing
  
  for i in "${TRK[@]}"
  do
    mkdir Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth/Nominal-NoSubtraction-NoSmoothing/Track-$i
  done

  mkdir Final/BestAlgos
  mkdir Final/BestAlgos/ROOT

  cd $root_dir/Final/BestAlgos
  mkdir Minimizer
  mkdir OptimalMiniWithAlgo
  mkdir FLost2_Error
  mkdir FLost3_Error
  cd $root_dir
}

CollectFinalPlots () {
  local energy=$1
  local mini=$2
  
  root_dir=$PWD
  ill=Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth/Nominal-NoSubtraction-NoSmoothing
  for i in "${Algs[@]}"
  do
    for j in "${TRK[@]}"
    do
      cp Output/$mini/Histograms/Test/trk$j/$i/$energy.png "$ill/Track-$j/$energy""_""$i.png"
    done
  done

  for i in "${TRK[@]}"
  do
    cp Output/$mini/ShapePerformance/trk$i/* $ill/Track-$i
  done

  for i in "${Modes[@]}"
  do
    cp Output/$i/$i.root Final/BestAlgos/ROOT
  done

  cd $root_dir/Final/BestAlgos/ROOT
  str=""
  for i in *
  do
    str="$str $i"
  done
  
  MergeROOT "$str" "BestAlgos"
  mv BestAlgos.root ../
  cp $root_dir/ProgramCode/BestAlgo ../
  cd ../
  chmod +x BestAlgo
  ./BestAlgo BestAlgos.root


  



}
