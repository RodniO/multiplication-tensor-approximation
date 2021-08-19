# Multiplication Tensor Approximation

Hello!
This repository contains a Fortran code for the approximation of small n^2 x m^2 x k^2 multiplication tensors. For 9x9x9 tensors it should take on average 30 minutes to find a low-rank approximation, containing only 0, 1 and -1.

Algorithm was checked to work on all sizes, such that nmk <= 27, without changing any parameters.

Algorithm code is in incfiles/ExampleA.f90 and incfiles/ExampleTEN.f90.

The main program is Main.f90. By default, it will do a single try for 9x9x9 tensor.

To compile and run it, use either

make gnu_run

(Requires gfortran, BLAS and LAPACK)

or

make run

(Requires ifort and mkl)

The code uses matrix modules from "projective-volume-low-rank" repository.

Structure:

The code is organized to enable compilation of multiple programs with a single makefile. This is achieved through 'src/incfiles' directory, where one can write functions and then include and call them in 'Main.f90'. Makefile pays attention only to files with '.f90' extension, so remember to name them appropriately.
Module and object files are saved in 'obj_gnu' and 'obj_intel' folders (created during compilation), enabling to have gfortran and ifort versions compiled and saved simultaneously.

Large sizes:

Code uses 'allocatable' arrays, which are usually put in stack. Therefore, for large amounts of data, you may require to use 'ulimit -s unlimited' to allow unlimited stack size.

Debug:

debugrun.sh script uses gfortran with debug options to recompile and run everything. Remember to run 'make clean' after debug before running 'make gnu'.
