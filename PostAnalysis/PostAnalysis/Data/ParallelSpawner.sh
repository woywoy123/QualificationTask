#!/bin/bash


function CreateBatches_Local
{
  
  echo "#!/bin/bash" >> Spawn.sh
  echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> Spawn.sh
  echo "source $""{ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" >> Spawn.sh
  echo "cd $4" >> Spawn.sh
  echo "cur=$""PWD" >> Spawn.sh
  echo "cd ../PostAnalysis && asetup --restore" >> Spawn.sh
  echo "cd ../build/" >> Spawn.sh
  echo "source x86_64-centos7-gcc62-opt/setup.sh" >> Spawn.sh
  echo "cd $""cur" >> Spawn.sh
  echo "PostAnalysis $1 $2 $3/$5" >> Spawn.sh
  #echo "mkdir /eos/home-t/tnommens/Analysis/" >> Spawn.sh
  #echo "cp ./* /eos/home-t/tnommens/Analysis/" >> Spawn.sh
}

function CondorBuild
{
  echo "executable = Spawn.sh" >> example.submit
  #echo "output = ./results.output.$""(ClusterID)"  >> example.submit
  echo "error =  ./results.error.$""(ClusterID)"  >> example.submit
  echo "log =  ./results.log.$""(ClusterID)"  >> example.submit
  echo "Request_Cpus = 1"  >> example.submit
  echo "Request_Memory = 1024MB" >> example.submit
  echo "+RequestRunTime= 172800"  >> example.submit
  echo "queue 1"  >> example.submit
}

#Constants that we need to generate the names 
Condor_active=true
compiler="PostAnalysisCompiler_Negative"
filename="Merged_MC_Negative.root"
Layer=("IBL" "Blayer" "layer1" "layer2") 
JetEnergy=("200_400_GeV" "400_600_GeV" "600_800_GeV" "800_1000_GeV" "1000_1200_GeV" "1200_1400_GeV" "1400_1600_GeV" "1600_1800_GeV" "1800_2000_GeV" "2000_2200_GeV" "2200_2400_GeV" "2400_2600_GeV" "2600_2800_GeV" "2800_3000_GeV" "higher_GeV")
Mode=("Truth")

# Default no subtract
Mode1=("FitT_Normal_Minimizer" "FitT_ShiftNormal_Minimizer" "FitT_ShiftNormalFFT_Minimizer" "FitT_ShiftNormalWidthFFT_Minimizer" "FitT_Incremental_Minimizer")
Mode2=("FitT_Normal_Minimizer_Smooth" "FitT_ShiftNormal_Minimizer_Smooth" "FitT_ShiftNormalFFT_Minimizer_Smooth" "FitT_ShiftNormalWidthFFT_Minimizer_Smooth" "FitT_Incremental_Minimizer_Smooth")
Mode3=("FitT_Normal_Smooth_FitTo" "FitT_ShiftNormal_Smooth_FitTo" "FitT_ShiftNormalFFT_Smooth_FitTo" "FitT_ShiftNormalWidthFFT_Smooth_FitTo" "FitT_Incremental_Smooth_FitTo")
Mode4=("FitT_Normal_FitTo" "FitT_ShiftNormal_FitTo" "FitT_ShiftNormalFFT_FitTo" "FitT_ShiftNormalWidthFFT_FitTo" "FitT_Incremental_FitTo")
#
## Subtraction
Mode5=("FitT_Normal_Subtract_Minimizer" "FitT_ShiftNormal_Subtract_Minimizer" "FitT_ShiftNormalFFT_Subtract_Minimizer" "FitT_ShiftNormalWidthFFT_Subtract_Minimizer" "FitT_Incremental_Subtract_Minimizer")
Mode6=("FitT_Normal_Subtract_Minimizer_Smooth" "FitT_ShiftNormal_Subtract_Minimizer_Smooth" "FitT_ShiftNormalFFT_Subtract_Minimizer_Smooth" "FitT_ShiftNormalWidthFFT_Subtract_Minimizer_Smooth" "FitT_Incremental_Subtract_Minimizer_Smooth")
Mode7=("FitT_Normal_Subtract_Smooth_FitTo" "FitT_ShiftNormal_Subtract_Smooth_FitTo" "FitT_ShiftNormalFFT_Subtract_Smooth_FitTo" "FitT_ShiftNormalWidthFFT_Subtract_Smooth_FitTo" "FitT_Incremental_Subtract_Smooth_FitTo")
Mode8=("FitT_Normal_Subtract_FitTo" "FitT_ShiftNormal_Subtract_FitTo" "FitT_ShiftNormalFFT_Subtract_FitTo" "FitT_ShiftNormalWidthFFT_Subtract_FitTo" "FitT_Incremental_Subtract_FitTo" "FitT_Simultaneous_Subtract_FitTo")
#
## Truth Fits 
Mode9=("FitT_Normal_TRUTH_Minimizer" "FitT_ShiftNormal_TRUTH_Minimizer" "FitT_ShiftNormalFFT_TRUTH_Minimizer" "FitT_ShiftNormalWidthFFT_TRUTH_Minimizer" "FitT_Incremental_TRUTH_Minimizer")
Mode10=("FitT_Normal_TRUTH_FitTo" "FitT_ShiftNormal_TRUTH_FitTo" "FitT_ShiftNormalFFT_TRUTH_FitTo" "FitT_ShiftNormalWidthFFT_TRUTH_FitTo" "FitT_Incremental_TRUTH_FitTo")

# Default Range Fits
Mode11=("FitT_Normal_Minimizer_Range" "FitT_ShiftNormal_Minimizer_Range" "FitT_ShiftNormalFFT_Minimizer_Range" "FitT_ShiftNormalWidthFFT_Minimizer_Range" "FitT_Incremental_Minimizer_Range")
Mode12=("FitT_Normal_Minimizer_Smooth_Range" "FitT_ShiftNormal_Minimizer_Smooth_Range" "FitT_ShiftNormalFFT_Minimizer_Smooth_Range" "FitT_ShiftNormalWidthFFT_Minimizer_Smooth_Range" "FitT_Incremental_Minimizer_Smooth_Range")
Mode13=("FitT_Normal_Smooth_FitTo_Range" "FitT_ShiftNormal_Smooth_FitTo_Range" "FitT_ShiftNormalFFT_Smooth_FitTo_Range" "FitT_ShiftNormalWidthFFT_Smooth_FitTo_Range" "FitT_Incremental_Smooth_FitTo_Range")
Mode14=("FitT_Normal_FitTo" "FitT_ShiftNormal_FitTo_Range" "FitT_ShiftNormalFFT_FitTo_Range" "FitT_ShiftNormalWidthFFT_FitTo_Range" "FitT_Incremental_FitTo_Range")

Mode+=(${Mode1[@]})
Mode+=(${Mode2[@]})
Mode+=(${Mode3[@]})
Mode+=(${Mode4[@]})
Mode+=(${Mode5[@]})
Mode+=(${Mode6[@]})
Mode+=(${Mode7[@]})
Mode+=(${Mode8[@]})
Mode+=(${Mode9[@]})
Mode+=(${Mode10[@]})

root_dir=$PWD
echo $root_dir

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

cd ../../
PostAnalysis_root_dir=$PWD
echo $PostAnalysis_root_dir

echo "You have 30 seconds to cancel the deleting of $compiler !"
sleep 10

cd $HOME
rm -r $compiler
mkdir $compiler

cd $compiler
cp -r $PostAnalysis_root_dir . 

cd PostAnalysis && asetup AnalysisBase,21.2.58,here
cd ../
mkdir build 
cd build 
cmake ../PostAnalysis 
make -j12 
cd ../
cp $root_dir/$filename ./
File=$PWD
for L in ${Layer[@]}
do

  for E in ${JetEnergy[@]}
  do
    for M in ${Mode[@]}
    do
      Line=$L"_"$E"_"$M
     
      echo $Line
      mkdir $Line
      cd $Line 
      LJE=$L"_"$E
      
      CreateBatches_Local $LJE $M $File $PWD $filename
      CondorBuild
      chmod +x Spawn.sh
      
      if [[ $Condor_active == true && $M != "Truth" ]]
      then 
        condor_submit example.submit 
      else
        bash Spawn.sh
      fi 
      
      cd ../
    done 
  done 


  for M in ${Mode[@]}
  do
    Line=$L"_"$M
   
    echo $Line
    mkdir $Line
    cd $Line 
    LJE=$L
    
    CreateBatches_Local $LJE $M $File $PWD $filename
    CondorBuild
    chmod +x Spawn.sh
    
    if [[ $Condor_active == true && $M != "Truth" ]]
    then 
      condor_submit example.submit 
    else
      bash Spawn.sh
    fi 
    
    cd ../
  done 

done

JetEnergy+=("All")
for E in ${JetEnergy[@]}
do
  for M in ${Mode[@]}
  do
    Line=$E"_"$M
   
    echo $Line
    mkdir $Line
    cd $Line 
    LJE=$E
    
    CreateBatches_Local $LJE $M $File $PWD $filename
    CondorBuild
    chmod +x Spawn.sh
    
    if [[ $Condor_active == true && $M != "Truth" ]]
    then 
      condor_submit example.submit 
    else
      bash Spawn.sh
    fi 
    
    cd ../
  done 
done 


