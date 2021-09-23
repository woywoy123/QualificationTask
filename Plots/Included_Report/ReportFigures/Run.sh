#!/bin/bash
MergeROOT () {
  local str=$1
  local dir=$2
  cd $dir
  /home/tnom6927/ROOT/bin/hadd MultiTrackFit.root $str
  mv MultiTrackFit.root ../
}

Algs=("Normal" "ShiftNormal" "ShiftNormalFFT" "ShiftNormalWidthFFT" "Incremental" "Simultaneous" "Experimental" "Truth")
Modes=("_Subtract_Smooth_FitTo" "_Subtract_Minimizer_Smooth" "_Subtract_Minimizer" "_Subtract_FitTo" "_Smooth_FitTo" "_Minimizer_Smooth" "_Minimizer" "_FitTo" "_Smooth_FitTo_Range" "_Minimizer_Smooth_Range" "_Minimizer_Range" "_FitTo_Range" "_TRUTH_Minimizer" "_TRUTH_FitTo")

root_dir=$PWD
rm -r Output
mkdir Output
cd Output

for i in ${Modes[@]}
do
  echo "-> $i"
  mkdir $i
  mkdir "$i/ROOT"
done 

cd $root_dir
cd Input

for i in ${Algs[@]}
do
  for k in ${Modes[@]}
  do 
    str="FitT_"$i$k
    for j in *
    do
      copy=false
      if [[ "$j" == *"$str" ]]
      then 
        copy=true 
        #echo $j " ----> " $str
      elif [[ "$j" == *"Truth" ]]
      then 
        copy=true
      fi
      
      if [[ $copy == true ]]
      then 
        f=$(wc -c $PWD/$j/Fit_Tracks.root)
        size=($f)
        size=${size[0]}
        if [[ $size != 63 ]]
        then 
          cp $PWD/$j/Fit_Tracks.root $root_dir/Output/$k/ROOT/$j.root
        fi
      fi

    done
  done 
done 
cd $root_dir
cd Output

source ~/.bashrc
for i in *
do
  cd $i/ROOT/
  str=""
  for j in *
  do
    str="$str $j"
  done
  
  MergeROOT "$str" "$PWD"
  cd $root_dir/Output
done
cd $root_dir

cd ProgramCode
bash exec.sh
prog=$PWD/Plotter
alg=$PWD/BestAlgo
cd ../
chmod +x $prog

cd Output
for i in ${Modes[@]}
do
  echo $i
  cd $i
  mkdir ShapePerformance
  mkdir IntegralPerformance
  mkdir Histograms

  mkdir ./Histograms/Test
  mkdir ./Histograms/Test/trk1
  mkdir ./Histograms/Test/trk1/Normal
  mkdir ./Histograms/Test/trk1/ShiftNormal
  mkdir ./Histograms/Test/trk1/ShiftNormalFFT
  mkdir ./Histograms/Test/trk1/ShiftNormalWidthFFT
  mkdir ./Histograms/Test/trk1/Incremental

  mkdir ./Histograms/Test/trk2
  mkdir ./Histograms/Test/trk2/Normal
  mkdir ./Histograms/Test/trk2/ShiftNormal
  mkdir ./Histograms/Test/trk2/ShiftNormalFFT
  mkdir ./Histograms/Test/trk2/ShiftNormalWidthFFT
  mkdir ./Histograms/Test/trk2/Incremental

  mkdir ./Histograms/Test/trk3
  mkdir ./Histograms/Test/trk3/Normal
  mkdir ./Histograms/Test/trk3/ShiftNormal
  mkdir ./Histograms/Test/trk3/ShiftNormalFFT
  mkdir ./Histograms/Test/trk3/ShiftNormalWidthFFT
  mkdir ./Histograms/Test/trk3/Incremental

  mkdir ./Histograms/Test/trk4
  mkdir ./Histograms/Test/trk4/Normal
  mkdir ./Histograms/Test/trk4/ShiftNormal
  mkdir ./Histograms/Test/trk4/ShiftNormalFFT
  mkdir ./Histograms/Test/trk4/ShiftNormalWidthFFT
  mkdir ./Histograms/Test/trk4/Incremental

  mkdir ./Histograms/Template
  mkdir ./Histograms/Template/trk1
  mkdir ./Histograms/Template/trk1/Normal
  mkdir ./Histograms/Template/trk1/ShiftNormal
  mkdir ./Histograms/Template/trk1/ShiftNormalFFT
  mkdir ./Histograms/Template/trk1/ShiftNormalWidthFFT
  mkdir ./Histograms/Template/trk1/Incremental
 
  mkdir ./Histograms/Template/trk2
  mkdir ./Histograms/Template/trk2/Normal
  mkdir ./Histograms/Template/trk2/ShiftNormal
  mkdir ./Histograms/Template/trk2/ShiftNormalFFT
  mkdir ./Histograms/Template/trk2/ShiftNormalWidthFFT
  mkdir ./Histograms/Template/trk2/Incremental
                     
  mkdir ./Histograms/Template/trk3
  mkdir ./Histograms/Template/trk3/Normal
  mkdir ./Histograms/Template/trk3/ShiftNormal
  mkdir ./Histograms/Template/trk3/ShiftNormalFFT
  mkdir ./Histograms/Template/trk3/ShiftNormalWidthFFT
  mkdir ./Histograms/Template/trk3/Incremental
                    
  mkdir ./Histograms/Template/trk4
  mkdir ./Histograms/Template/trk4/Normal
  mkdir ./Histograms/Template/trk4/ShiftNormal
  mkdir ./Histograms/Template/trk4/ShiftNormalFFT
  mkdir ./Histograms/Template/trk4/ShiftNormalWidthFFT
  mkdir ./Histograms/Template/trk4/Incremental


  mkdir ./ShapePerformance/trk1
  mkdir ./ShapePerformance/trk2
  mkdir ./ShapePerformance/trk3
  mkdir ./ShapePerformance/trk4

  mkdir FLost

  cp $prog ./
  cp $alg ./
  chmod +rwx Plotter
  chmod +rwx BestAlgo
  ./Plotter &
  ./BestAlgo $i &
  cd ../
done
wait 

cd $root_dir
rm -r Final
mkdir Final 
mkdir Final/IllustrationsOfAlgorithmFits
mkdir Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth
mkdir Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth/Nominal-NoSubtraction-NoSmoothing
cd Final/IllustrationsOfAlgorithmFits/DirectFitsToTruth/Nominal-NoSubtraction-NoSmoothing
IllFit=$PWD
mkdir Track-1
mkdir Track-2
mkdir Track-3
mkdir Track-4
cd $root_dir

cp ./Output/_FitTo_Range/Histograms/Test/trk1/Incremental/Blayer_800_1000_GeV.png $IllFit/Track-1/Blayer_800_1000_GeV_Incremental.png
cp ./Output/_FitTo_Range/Histograms/Test/trk2/Incremental/Blayer_800_1000_GeV.png $IllFit/Track-2/Blayer_800_1000_GeV_Incremental.png
cp ./Output/_FitTo_Range/Histograms/Test/trk3/Incremental/Blayer_800_1000_GeV.png $IllFit/Track-3/Blayer_800_1000_GeV_Incremental.png
cp ./Output/_FitTo_Range/Histograms/Test/trk4/Incremental/Blayer_800_1000_GeV.png $IllFit/Track-4/Blayer_800_1000_GeV_Incremental.png

cp ./Output/_FitTo_Range/Histograms/Test/trk1/Normal/Blayer_800_1000_GeV.png $IllFit/Track-1/Blayer_800_1000_GeV_Normal.png
cp ./Output/_FitTo_Range/Histograms/Test/trk2/Normal/Blayer_800_1000_GeV.png $IllFit/Track-2/Blayer_800_1000_GeV_Normal.png
cp ./Output/_FitTo_Range/Histograms/Test/trk3/Normal/Blayer_800_1000_GeV.png $IllFit/Track-3/Blayer_800_1000_GeV_Normal.png
cp ./Output/_FitTo_Range/Histograms/Test/trk4/Normal/Blayer_800_1000_GeV.png $IllFit/Track-4/Blayer_800_1000_GeV_Normal.png

cp ./Output/_FitTo_Range/Histograms/Test/trk1/ShiftNormal/Blayer_800_1000_GeV.png $IllFit/Track-1/Blayer_800_1000_GeV_ShiftNormal.png
cp ./Output/_FitTo_Range/Histograms/Test/trk2/ShiftNormal/Blayer_800_1000_GeV.png $IllFit/Track-2/Blayer_800_1000_GeV_ShiftNormal.png
cp ./Output/_FitTo_Range/Histograms/Test/trk3/ShiftNormal/Blayer_800_1000_GeV.png $IllFit/Track-3/Blayer_800_1000_GeV_ShiftNormal.png
cp ./Output/_FitTo_Range/Histograms/Test/trk4/ShiftNormal/Blayer_800_1000_GeV.png $IllFit/Track-4/Blayer_800_1000_GeV_ShiftNormal.png

cp ./Output/_FitTo_Range/Histograms/Test/trk1/ShiftNormalFFT/Blayer_800_1000_GeV.png $IllFit/Track-1/Blayer_800_1000_GeV_ShiftNormalFFT.png
cp ./Output/_FitTo_Range/Histograms/Test/trk2/ShiftNormalFFT/Blayer_800_1000_GeV.png $IllFit/Track-2/Blayer_800_1000_GeV_ShiftNormalFFT.png
cp ./Output/_FitTo_Range/Histograms/Test/trk3/ShiftNormalFFT/Blayer_800_1000_GeV.png $IllFit/Track-3/Blayer_800_1000_GeV_ShiftNormalFFT.png
cp ./Output/_FitTo_Range/Histograms/Test/trk4/ShiftNormalFFT/Blayer_800_1000_GeV.png $IllFit/Track-4/Blayer_800_1000_GeV_ShiftNormalFFT.png

cp ./Output/_FitTo_Range/Histograms/Test/trk1/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $IllFit/Track-1/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png
cp ./Output/_FitTo_Range/Histograms/Test/trk2/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $IllFit/Track-2/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png
cp ./Output/_FitTo_Range/Histograms/Test/trk3/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $IllFit/Track-3/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png
cp ./Output/_FitTo_Range/Histograms/Test/trk4/ShiftNormalWidthFFT/Blayer_800_1000_GeV.png $IllFit/Track-4/Blayer_800_1000_GeV_ShiftNormalWidthFFT.png

cp ./Output/_FitTo_Range/ShapePerformance/trk1/* $IllFit/Track-1
cp ./Output/_FitTo_Range/ShapePerformance/trk2/* $IllFit/Track-2
cp ./Output/_FitTo_Range/ShapePerformance/trk3/* $IllFit/Track-3
cp ./Output/_FitTo_Range/ShapePerformance/trk4/* $IllFit/Track-4

cd $root_dir 
mkdir Final/BestAlgos
mkdir Final/BestAlgos/ROOT

for i in ${Modes[@]}
do
  
  if [[ $i != *"TRUTH"* ]]
  then 
    cp Output/$i/$i.root Final/BestAlgos/ROOT
  fi
done 
cd Final/BestAlgos/ROOT

str=""
for i in *
do
  str="$str $i"
done

echo $str
MergeROOT "$str" ./
cd ../

mkdir Minimizer
mv MultiTrackFit.root BestAlgos.root
cp $root_dir/ProgramCode/BestAlgo .
chmod +x BestAlgo
./BestAlgo BestAlgos.root



