#!/bin/bash

Data=./CalibrationData

mkdir Merger
cd Merger
Merge=$PWD
cd ..

cd $Data
data=$PWD
for item in *
do
  cd $item
  for file in *
  do
    files=$file
  done
  cp $data/$item/$files $Merge/

  cd ../
done

cd $Merge
string=""
for item in *
do
  string="$item $string"
done

ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "root 6.20.06-x86_64-centos7-gcc-opt"

$ROOTSYS/bin/hadd Merger.root $string


