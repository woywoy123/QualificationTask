#!/bin/bash


function CreateBatches
{
  
  echo "#!/bin/bash" >> Spawn.sh
  echo "ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> Spawn.sh
  echo "source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" >> Spawn.sh
  echo "mkdir build" >> Spawn.sh
  echo "cd build" >> Spawn.sh
  echo "cp $2/Merged_MC.root ." >> Spawn.sh
  echo "cd $1 && asetup --restore" >> Spawn.sh
  echo "cd $PWD/build" >> Spawn.sh
  echo "cmake $1" >> Spawn.sh
  echo "make -j 12" >> Spawn.sh
  echo "source ./x86_64-centos7-gcc62-opt/setup.sh" >> Spawn.sh
  echo "PostAnalysis $Line" >> Spawn.sh
}

#Constants that we need to generate the names 
Layer=("IBL" "Blayer" "layer1" "layer2") 
JetEnergy=("200_up_GeV" "200_400_GeV" "400_600_GeV" "600_800_GeV" "800_1000_GeV" "1000_1200_GeV" "1200_1400_GeV" "1400_1600_GeV" "1600_1800_GeV" "1800_2000_GeV" "2000_2200_GeV" "2200_2400_GeV" "2400_2600_GeV" "2600_2800_GeV" "2800_3000_GeV" "higher_GeV")

oper="x86_64-centos7-gcc62-opt"
#oper="x86_64-slc6-gcc62-opt"

root_dir=$PWD
echo $root_dir

cd ../../
PostAnalysis_root_dir=$PWD
echo $PostAnalysis_root_dir

cd $HOME
mkdir PostAnalysisCompiler

cd PostAnalysisCompiler

for L in ${Layer[@]}
do

  for E in ${JetEnergy[@]}
  do

    Line=$L"_"$E
    
    rm -rf $Line
    
    echo $Line
    mkdir $Line

    cd $Line 
    
    CreateBatches $PostAnalysis_root_dir $root_dir $Line $oper
    bash Spawn.sh
    sleep 15

    cd ../
  done 
done 
