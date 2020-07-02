#!/bin/bash

NumberOfCores=8
Files_Per_Core=3
SampleDir="/CERN/QT/mc16_13TeV.361026.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ6W.recon.DAOD_IDTIDE.e3569_s3126_r10504"
CollectionDir="Output/data-PixhistOutput"

mkdir Processed
mkdir ToProcess

Output=$PWD/Processed

cd ../
RootDir=$PWD
cd src/

# Set the current environment for the code to run
ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"

# Source the code for the CTIDE build 
echo $PWD
asetup --restore 
source $RootDir/build/x86_64-centos7-gcc62-opt/setup.sh

# execute the analysis code
cd $SampleDir/
Files=()
for item in *
do 
  Files+=("$item")
  echo "$item"
done

cd $RootDir/Parallel_Compute

for ((i=1; i<=$NumberOfCores; i++))
do
  mkdir "ToProcess/$i"
done 

Completed=()
RunTime=true
Analysis_index=0
while [ $RunTime != 0 ]
do 

  # Run the analysis
  PIDs=()
  for ((core=1; core<=$NumberOfCores; core++))
  do
    index=1
    for i in ${Files[@]}
    do
      if [[ "${Completed[@]}" =~ "${i}" ]]
      then
        :
      else 
        echo "Copying $i to core $core folder."
        cp $SampleDir/$i ToProcess/$core/
        Completed+=("$i") 
        index=$(($index+1))
      fi 
  
      if (( $index > $Files_Per_Core ))
      then
        break
      fi
    done

    mkdir $Output/$Analysis_index 
    runAnalysis -inDS $PWD/ToProcess/$core -submitDir $Output/$Analysis_index/Output -driver local -analysis dEdx >> $Output/$Analysis_index/log.txt &    
    PIDs+=("$!") 
    Analysis_index=$(($Analysis_index+1)) 

  done


  while true
  do
    alive=()
    for i in ${PIDs[@]}
    do

      if ps -p $i > /dev/null
      then
        echo "Running: $i"
        alive+=("$i")
      fi

    done
 
    echo "${#alive[@]}" 
    if (( ${#alive[@]} == 0 ))
    then
      break
    fi
    sleep 30
  
  done
  
  # remove the files 
  for ((core=1; core<=$NumberOfCores; core++))
  do
    rm -r $PWD/ToProcess/$core/*
    sleep 2
    
    # Run the analysis
    PIDs=()
  done

  cl="${#Completed[@]}"
  fl="${#Files[@]}"
  if [[ "$cl" == "$fl" ]];
  then 
    RunTime=0
  fi
done

cd $Output

mkdir ../Merger
for item in *
do
  echo $item
  cp $item/$CollectionDir/* ../Merger/$item.root
done

cd ../Merger
echo $PWD
string=""
for item in *
do
  string="$item $string"
done 

echo $PWD
$ROOTSYS/bin/hadd Merged.root $string


