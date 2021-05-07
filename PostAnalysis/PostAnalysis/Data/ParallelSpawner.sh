#!/bin/bash


function CreateBatches_Local
{
  
  echo "#!/bin/bash" >> Spawn.sh
  echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> Spawn.sh
  echo "source $""{ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" >> Spawn.sh
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
  #echo "log =  ./results.log.$""(ClusterID)"  >> example.submit
  echo "Request_Cpus = 4"  >> example.submit
  echo "Request_Memory = 1GB" >> example.submit
  echo "+RequestRunTime= 43200"  >> example.submit
  echo "queue 1"  >> example.submit
}


#Constants that we need to generate the names 
Condor_active=true
Layer=("IBL" "Blayer" "layer1" "layer2") 
JetEnergy=("200_up_GeV" "200_400_GeV" "400_600_GeV" "600_800_GeV" "800_1000_GeV" "1000_1200_GeV" "1200_1400_GeV" "1400_1600_GeV" "1600_1800_GeV" "1800_2000_GeV" "2000_2200_GeV" "2200_2400_GeV" "2400_2600_GeV" "2600_2800_GeV" "2800_3000_GeV" "higher_GeV" "All")
Mode=("ShiftNormal" "ShiftNormalFFT" "ShiftNormalWidthFFT" "Experimental" "Normal" "Truth")
#Mode=("Normal" "Truth"); 
root_dir=$PWD
echo $root_dir

cd ../../
PostAnalysis_root_dir=$PWD
echo $PostAnalysis_root_dir

cd $HOME
rm -r PostAnalysisCompiler
mkdir PostAnalysisCompiler

cd PostAnalysisCompiler
cp -r $PostAnalysis_root_dir . 
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
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
      
      CreateBatches_Local $LJE $M $File
      CondorBuild
      chmod +x Spawn.sh
      
      if [[ $Condor_active == true ]]
      then 
        condor_submit example.submit 
      else
        bash Spawn.sh
      fi 
      
      cd ../
    done 
  done 
done

for L in ${Layer[@]}
do
  for M in ${Mode[@]}
  do
    Line=$L"_"$M
   
    echo $Line
    mkdir $Line
    cd $Line 

    CreateBatches_Local $L $M $File
    CondorBuild
    chmod +x Spawn.sh

    if [[ $Condor_active == true ]]
    then 
      condor_submit example.submit 
    else
      bash Spawn.sh
    fi 

    
    cd ../
  done 
done 

for L in ${JetEnergy[@]}
do
  for M in ${Mode[@]}
  do
    Line=$L"_"$M
   
    echo $Line
    mkdir $Line
    cd $Line 

    CreateBatches_Local $L $M $File
    CondorBuild
    chmod +x Spawn.sh

    if [[ $Condor_active == true ]]
    then 
      condor_submit example.submit 
    else
      bash Spawn.sh
    fi 


    cd ../
  done 
done 


