#!/bin/bash

CleanWorkspace(){
  cd $root_dir
  rm -rf ../PlottingCode/Binaries
  rm -rf $root_dir/Output
}

MakeClusterContent(){
  mkdir $OUTPUT_DIR/Appendix/ClusterContent
  mkdir $OUTPUT_DIR/Appendix/ClusterContent/ROOT
  cp $INPUT_DIR/Merged.root $OUTPUT_DIR/Appendix/ClusterContent/ROOT/Merged.root
  cp -f $CODE_DIR/Binaries/ClusterContent $OUTPUT_DIR/Appendix/ClusterContent/ClusterContent
  cd $OUTPUT_DIR/Appendix/ClusterContent
  chmod +x ClusterContent
  ./ClusterContent $PWD/ROOT/Merged.root
  rm -rf ROOT
  rm -rf ClusterContent
  cd $root_dir
}

MergeROOT(){
  local r=$PWD
  cd $1
  mkdir ../TMP
  x=0
  it=0
  mkdir ../TMP/$x
  for i in *
  do
    mv $i ../TMP/$x
    if [[ $it -gt 100 ]]
    then 
      x=$((x+1))
      it=0
      mkdir ../TMP/$x
    fi
    it=$((it+1))
  done
  cd ../TMP
  
  str_tmp=""
  for i in *
  do
    cd $i
    str=""
    for j in *
    do
      str="$str $j"
    done
    $ROOT_MERGE $i.root $str > /dev/null
    mv $i.root ../
    str_out="$str_out $i.root"
    cd ../
  done
  $ROOT_MERGE Output.root $str_out > /dev/null
  mv Output.root ../
}

CollectDebuggingROOT(){
  cd $INPUT_DIR
  for i in *
  do 
    for j in ${!TRKTEMPLATE_D[@]}
    do
      mod=${TRKTEMPLATE_D[j]}
      dir=${TRKTEMPLATE_F[j]}
      
      if [[ "$mod" == "_" ]]
      then 
        n="Debug"
      else
        n="Debug$mod"
      fi
      
      if [[ "$i" == *"$n" ]]
      then 
        cp $PWD/$i/*.root $OUTPUT_DIR/Appendix/Debugging/$dir/ROOT/$i.root
      fi
    done
  done 
  
  cd $OUTPUT_DIR

  for i in ${TRKTEMPLATE_F[@]}
  do
    MergeROOT $OUTPUT_DIR/Appendix/Debugging/$i/ROOT &
  done
  wait
}

MakeDebuggingContent(){
  mkdir $OUTPUT_DIR/Appendix/Debugging
  
  for i in ${TRKTEMPLATE_F[@]}
  do
    mkdir $OUTPUT_DIR/Appendix/Debugging/$i
    mkdir $OUTPUT_DIR/Appendix/Debugging/$i/ROOT
    cp -f $CODE_DIR/Binaries/Debugging $OUTPUT_DIR/Appendix/Debugging/$i/Debugging
   
    mkdir $OUTPUT_DIR/Appendix/Debugging/$i/Ratio_2Truth
    mkdir $OUTPUT_DIR/Appendix/Debugging/$i/Ratio_2Truth_Layer

    for j in ${Layer[@]}
    do
      mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j
      for k in ${ENERGY[@]}
      do
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/Normal
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormal
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormalFFT

        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/Normal/Case1
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormal/Case1
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormalFFT/Case1

        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/Normal/Case2
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormal/Case2
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormalFFT/Case2
       
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/Normal/Case3
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormal/Case3
        mkdir $OUTPUT_DIR/Appendix/Debugging/$i/$j/$k/ShiftNormalFFT/Case3
      done
    done
  done
  
  CollectDebuggingROOT

  for i in ${TRKTEMPLATE_F[@]}
  do
    cd $OUTPUT_DIR/Appendix/Debugging/$i
    chmod +x Debugging
    ./Debugging Output.root &
  done
  wait
}

CollectCommonROOT()
{
  mkdir $OUTPUT_DIR/Appendix/Common
  
  cd $INPUT_DIR
  for i in ${Modes[@]}
  do
    mkdir $OUTPUT_DIR/Appendix/Common/$i
    mkdir $OUTPUT_DIR/Appendix/Common/$i/ROOT
    
    for k in *
    do
      if [[ $k == *"TRUTH"* ]]
      then
        continue
      fi

      if [[ $k == *"$i""_ntrk"* ]]
      then 
        cp $k/Fit_Tracks.root $OUTPUT_DIR/Appendix/Common/$i/ROOT/$k.root
        continue
      elif [[ $k == *"$i" ]]
      then 
        cp $k/Fit_Tracks.root $OUTPUT_DIR/Appendix/Common/$i/ROOT/$k.root
      fi
    done
  done
  
  cd $OUTPUT_DIR/Appendix/Common
  for i in *
  do
    echo $i
    MergeROOT "$OUTPUT_DIR/Appendix/Common/$i/ROOT" &
  done
  wait
  cp $INPUT_DIR/Merged.root $OUTPUT_DIR/Appendix/Common/Merged.root
  cd $root_dir
}



MakeTestToTruthShape(){
  mkdir $OUTPUT_DIR/Appendix/TestToTruthShape
  for i in ${Modes[@]}
  do
    mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i
    mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/Distributions
    
    cd $OUTPUT_DIR/Appendix/TestToTruthShape/$i/
    cp -f $CODE_DIR/Binaries/TestToTruthShape $OUTPUT_DIR/Appendix/TestToTruthShape/$i/TestToTruthShape 
    chmod +x TestToTruthShape     
   
    for j in ${Algs_CPP[@]}
    do
      mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/Distributions/$j 
      mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/Distributions/$j/Track-1
      mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/Distributions/$j/Track-2
      mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/Distributions/$j/Track-3
      mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/Distributions/$j/Track-4
    done
    
    if [[ "$i" != "_Minimizer" ]]
    then 
      continue
    fi
    ./TestToTruthShape $OUTPUT_DIR/Appendix/Common/$i/Output.root $OUTPUT_DIR/Appendix/Common/Merged.root
    
    exit

  done







}



MakeOutputDir(){
  mkdir $root_dir/Output/Appendix
  #mkdir $root_dir/Output/Appendix
  #mkdir $root_dir/Output/ClusterContent
  #MakeTrack $root_dir/Output/Appendix
  #mkdir $root_dir/Output/Debugging
  #mkdir $root_dir/Output/ErrorDirectTruthFits_FitTo
  #mkdir $root_dir/Output/ErrorOfDataTemplateFits
  #mkdir $root_dir/Output/ExampleOfExperimental
  #mkdir $root_dir/Output/FLost

}

