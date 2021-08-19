#!/bin/bash 

compiler="gfortran"
options="-Wall -Wno-unused-dummy-argument -Wno-integer-division -Wunused-variable -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero -g -pedantic -Wno-uninitialized"
source="ModIntVec.f90 ModVec.f90 ModMtrx.f90 MerTwist.f90 Main.f90"
object="ModIntVec.o ModVec.o ModMtrx.o MerTwist.o Main.o"
program="Main_gnu.exe"

mkdir -p obj_gnu;
for file in $source; do
  $compiler -J./obj_gnu -I./obj_gnu $options -c src/$file
  if [[ $? != 0 ]]; then
    exit
  fi
done
$compiler -I./obj_gnu -o $program $object -llapack -lblas
mv $object obj_gnu
./$program
