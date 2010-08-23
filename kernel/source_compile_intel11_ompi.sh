#!/bin/tcsh

# kick out all modules
module purge

# load required modules
module load intel/cc/11.1.072
module load intel/fc/11.1.072
module load openmpi/intel11

#compile prototype
#mpic++ -O0 -Wall -I../ universe_test.cpp -o universe_test
