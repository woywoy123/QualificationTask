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
  echo "PostAnalysis $1 $2 $3/Merged_MC.root" >> Spawn.sh
}

function CondorBuild
{
  echo "executable = Spawn.sh" >> example.submit
  #echo "output = ./results.output.$""(ClusterID)"  >> example.submit
  echo "error =  ./results.error.$""(ClusterID)"  >> example.submit
  echo "log =  ./results.log.$""(ClusterID)"  >> example.submit
  echo "Request_Cpus = 1"  >> example.submit
  echo "Request_Memory = 100MB" >> example.submit
  echo "+RequestRunTime= 172800"  >> example.submit
  echo "queue 1"  >> example.submit
}


#Constants that we need to generate the names 
Condor_active=true
compiler="PostAnalysisCompiler"
Layer=("IBL" "Blayer" "layer1" "layer2") 
JetEnergy=("200_400_GeV" "400_600_GeV" "600_800_GeV" "800_1000_GeV" "1000_1200_GeV" "1200_1400_GeV" "1400_1600_GeV" "1600_1800_GeV" "1800_2000_GeV" "2000_2200_GeV" "2200_2400_GeV" "2400_2600_GeV" "2600_2800_GeV" "2800_3000_GeV" "higher_GeV")

Mode=("Truth" "Normal" "ShiftNormal" "ShiftNormalFFT" "ShiftNormalWidthFFT" "Incremental" "Simultaneous" "Experimental" "Debug")
root_dir=$PWD
echo $root_dir

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

cd ../../
PostAnalysis_root_dir=$PWD
echo $PostAnalysis_root_dir

echo "You have 30 seconds to cancel the deleting of $compiler !"
sleep 30

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
cp $root_dir/Merged_MC.root ./
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
      
      CreateBatches_Local $LJE $M $File $PWD
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
    
    CreateBatches_Local $LJE $M $File $PWD
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
    
    CreateBatches_Local $LJE $M $File $PWD
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


