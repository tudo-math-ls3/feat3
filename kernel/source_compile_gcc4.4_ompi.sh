#!/bin/tcsh

# kick out all modules
module purge

# load required modules
module load gcc/4.4.0
module load openmpi/gcc/1.3.2

#compile prototype
#mpic++ -O0 -Wall -I../ universe_test.cpp -o universe_test
