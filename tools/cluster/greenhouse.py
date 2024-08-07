#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

import subprocess
import sys
import os
import time

sys.dont_write_bytecode = True

for i in range(1, 5):
  with open("temp", "w") as f:
    f.write("#!/bin/bash" + os.linesep)
    f.write("#PBS -q batch" + os.linesep)
    f.write("#PBS -l vmem=100gb,nodes=" + str(i) + ":ppn=16:bttf-cpu" + os.linesep)
    #f.write("#PBS -l vmem=100gb,nodes=" + str(i) + ":ppn=4:gpus=1:shared:h5-cpu-ib-gpu" + os.linesep)
    f.write("#PBS -l walltime=01:00:00" + os.linesep)
    f.write("#PBS -N feat" + os.linesep)
    f.write("#PBS -j oe" + os.linesep)

    f.write("export OMP_NUM_THREADS=1" + os.linesep)
    f.write("source ~/.bashrc &>/dev/null" + os.linesep)
    f.write("mpirun -np " + str(i*16) + " --map-by node ~/feat/trunk/applications/poisson_dirichlet_2d  --level " +  str(9) + " 3 --part_min_elems 500" + os.linesep)

  #sbatch
  subprocess.call(["qsub", "temp"])
  time.sleep(1)

os.remove("temp")
