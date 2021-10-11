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
Condor_active=false
compiler="PostAnalysisCompiler_Negative"
filename="Merged_MC_Negative.root"
Layer=("IBL" "Blayer" "layer1" "layer2") 
JetEnergy=("200_400_GeV" "400_600_GeV" "600_800_GeV" "800_1000_GeV" "1000_1200_GeV" "1200_1400_GeV" "1400_1600_GeV" "1600_1800_GeV" "1800_2000_GeV" "2000_2200_GeV" "2200_2400_GeV" "2400_2600_GeV" "2600_2800_GeV" "2800_3000_GeV" "higher_GeV")
Mode=("Truth")
Algos=("Experimental" "Normal" "ShiftNormal" "ShiftNormalFFT" "ShiftNormalWidthFFT" "Incremental")

# Default templates - No Subtract
Nominal=("Minimizer" "FitTo" "Minimizer_Smooth" "FitTo_Smooth")
Subtract=("Minimizer_Subtract" "FitTo_Subtract" "Minimizer_Subtract_Smooth" "FitTo_Subtract_Smooth")
TRU=("TRUTH_Minimizer" "TRUTH_FitTo")
Range=("Minimizer_Range" "FitTo_Range" "Minimizer_Range_Smooth" "FitTo_Range_Smooth")



for m in ${Algos[@]}
do
  for i in ${Nominal[@]}
  do
      str="FitT_""$m""_$i"
      Mode+=("$str")
  done

  for i in ${Subtract[@]}
  do
      str="FitT_""$m""_$i"
      Mode+=("$str")
  done

  for i in ${TRU[@]}
  do
      str="FitT_""$m""_$i"
      Mode+=("$str")
  done

  for i in ${Range[@]}
  do
      str="FitT_""$m""_$i"
      Mode+=("$str")
  done
done

#Mode=("Debug" "Debug_Subtract" "Debug_Subtract_Smooth" "Debug_Smooth")

root_dir=$PWD
echo $root_dir

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

cd ../../
PostAnalysis_root_dir=$PWD
echo $PostAnalysis_root_dir
cd $HOME
#rm -r $compiler
#mkdir $compiler

cd $compiler
rm -r PostAnalysis
rm -r build
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
      
      rm Spawn.sh
      CreateBatches_Local $LJE $M $File $PWD $filename
      CondorBuild
      chmod +x Spawn.sh
      
      if [[ $Condor_active == true && $M != "Truth" ]]
      then 
        condor_submit example.submit 
      else
        : #bash Spawn.sh
      fi 
      
      cd ../
    done 
  done 
done

