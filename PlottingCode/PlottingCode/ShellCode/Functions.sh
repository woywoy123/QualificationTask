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
  local x=0
  local it=0
  mkdir ../TMP/$x
  for i in *
  do
    cp $i ../TMP/$x
    if [[ $it -gt 100 ]]
    then 
      x=$((x+1))
      it=0
      mkdir ../TMP/$x
    fi
    it=$((it+1))
  done
  cd ../TMP
  
  local str_tmp=""
  for i in *
  do
    cd $i
    local str=""
    for j in *
    do
      str="$str $j"
    done
    $ROOT_MERGE $i.root $str > /dev/null
    mv $i.root ../
    str_tmp="$str_tmp $i.root"
    cd ../
  done
  $ROOT_MERGE Output.root $str_tmp > /dev/null
  mv Output.root ../
  cd ../
  rm -rf TMP
}

CollectDebuggingROOT(){
  mkdir $OUTPUT_DIR/Appendix/Debugging
  for i in ${TRKTEMPLATE_F[@]}
  do
    mkdir $OUTPUT_DIR/Appendix/Debugging/$i
    mkdir $OUTPUT_DIR/Appendix/Debugging/$i/ROOT
  done

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
  #mkdir $OUTPUT_DIR/Appendix/Debugging
  
  for i in ${TRKTEMPLATE_F[@]}
  do
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
    
    mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/ShapePerformance
    for j in ${Layer[@]}
    do
      mkdir $OUTPUT_DIR/Appendix/TestToTruthShape/$i/ShapePerformance/$j
    done
    
    ./TestToTruthShape $OUTPUT_DIR/Appendix/Common/$i/Output.root $OUTPUT_DIR/Appendix/Common/Merged.root &
  done
  wait
  cd $root_dir
}

MakeTestToDataShape(){
  mkdir $OUTPUT_DIR/Appendix/TestToDataShape
  for i in ${Modes[@]}
  do
    mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i
    mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/Distributions
    
    cd $OUTPUT_DIR/Appendix/TestToDataShape/$i/
    cp -f $CODE_DIR/Binaries/TestToDataShape $OUTPUT_DIR/Appendix/TestToDataShape/$i/TestToDataShape 
    chmod +x TestToDataShape     
   
    for j in ${Algs_CPP[@]}
    do
      mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/Distributions/$j 
      mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/Distributions/$j/Track-1
      mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/Distributions/$j/Track-2
      mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/Distributions/$j/Track-3
      mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/Distributions/$j/Track-4
    done
    
    mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/ShapePerformance
    for j in ${Layer[@]}
    do
      mkdir $OUTPUT_DIR/Appendix/TestToDataShape/$i/ShapePerformance/$j
    done
    
    ./TestToDataShape $OUTPUT_DIR/Appendix/Common/$i/Output.root $OUTPUT_DIR/Appendix/Common/Merged.root &
  done
  wait
  cd $root_dir
}

MakeFindBestTruthShape(){
  mkdir $OUTPUT_DIR/Appendix/BestTruthShape
  DIR_NAME=$OUTPUT_DIR/Appendix/BestTruthShape
  cd $DIR_NAME

  local str=""
  mkdir $DIR_NAME/FitToTruth 
  mkdir $DIR_NAME/FitToData 
  for l in ${Layer[@]}
  do 
    mkdir $DIR_NAME/FitToTruth/$l
    mkdir $DIR_NAME/FitToTruth/$l/Track-1
    mkdir $DIR_NAME/FitToTruth/$l/Track-2
    mkdir $DIR_NAME/FitToTruth/$l/Track-3
    mkdir $DIR_NAME/FitToTruth/$l/Track-4

    mkdir $DIR_NAME/FitToData/$l
    mkdir $DIR_NAME/FitToData/$l/Track-1
    mkdir $DIR_NAME/FitToData/$l/Track-2
    mkdir $DIR_NAME/FitToData/$l/Track-3
    mkdir $DIR_NAME/FitToData/$l/Track-4

  done
  mkdir $DIR_NAME/FitToTruth/Summary/
  mkdir $DIR_NAME/FitToTruth/Summary/Track-1
  mkdir $DIR_NAME/FitToTruth/Summary/Track-2
  mkdir $DIR_NAME/FitToTruth/Summary/Track-3
  mkdir $DIR_NAME/FitToTruth/Summary/Track-4

  mkdir $DIR_NAME/FitToData/Summary/
  mkdir $DIR_NAME/FitToData/Summary/Track-1
  mkdir $DIR_NAME/FitToData/Summary/Track-2
  mkdir $DIR_NAME/FitToData/Summary/Track-3
  mkdir $DIR_NAME/FitToData/Summary/Track-4


  for i in ${Modes[@]}
  do
    for l in ${Layer[@]}
    do 
      mkdir $DIR_NAME/FitToTruth/$l/Track-1/$i
      mkdir $DIR_NAME/FitToTruth/$l/Track-2/$i
      mkdir $DIR_NAME/FitToTruth/$l/Track-3/$i
      mkdir $DIR_NAME/FitToTruth/$l/Track-4/$i

      mkdir $DIR_NAME/FitToData/$l/Track-1/$i
      mkdir $DIR_NAME/FitToData/$l/Track-2/$i
      mkdir $DIR_NAME/FitToData/$l/Track-3/$i
      mkdir $DIR_NAME/FitToData/$l/Track-4/$i
    done
  done
  
  for i in ${Modes[@]}
  do
    str="$OUTPUT_DIR/Appendix/Common/$i/Output.root $str"
  done
  
  cp -f $CODE_DIR/Binaries/BestTruthShape $DIR_NAME/BestTruthShape 
  chmod +x BestTruthShape 
  ./BestTruthShape $OUTPUT_DIR/Appendix/Common/Merged.root $str

  cd $root_dir
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

