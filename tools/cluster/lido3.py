#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

import subprocess
import sys
import os
import time

sys.dont_write_bytecode = True

for i in range(0, 5):
  with open("temp", "w") as f:
    f.write("#!/bin/bash -x" + os.linesep)

    f.write("#SBATCH --nodes=" + str(4**i) + os.linesep)
    f.write("#SBATCH --ntasks-per-node=10" + os.linesep)
    f.write("#SBATCH --cpus-per-task=1" + os.linesep)
    f.write("#SBATCH --time=00:59:00" + os.linesep)
    f.write("#SBATCH --partition=short" + os.linesep)
    f.write("#SBATCH --exclusive" + os.linesep)

   # Possible 'constraint' values:
   #               all, public
   #               cgpu01  OR  gpu     OR  tesla_k40
   #               cquad01 OR  xeon_e54640v4
   #               cquad02 OR  xeon_e54640v4
   #               cstd01  OR  xeon_e52640v4   OR  ib_1to3
   #               cstd02  OR  xeon_e52640v4   OR  ib_1to1 OR  nonblocking_comm
   #SBATCH --constraint="cstd01"
   # Maximum 'mem' values depending on constraint (values in MB):
   #               cstd01/xeon_e52640v4/ib_1to3 AND
   #               cstd02/xeon_e52640v4/ib_1to1/nonblocking_comm: 62264
   #               cquad01: 255800
   #               cquad02: 1029944
   #               gpu: 62261
    f.write("#SBATCH --constraint=cstd01" + os.linesep)
    #f.write("#SBATCH --constraint=all" + os.linesep)
    #f.write("#SBATCH --constraint=gpu" + os.linesep)
    #f.write("#SBATCH --gres='gpu:1'" + os.linesep)
    f.write("#SBATCH --mem=62264" + os.linesep) #cstd01
    #f.write("#SBATCH --mem=62261" + os.linesep) #gpu / all
    f.write("#SBATCH --mem_bind=verbose,local --hint=memory_bound" + os.linesep)
    f.write("#SBATCH --job-name=feat" + os.linesep)

   #SBATCH --mail-user=<IHRE.E-MAILADRESSE>@tu-dortmund.de
   # Possible 'mail-type' values: NONE, BEGIN, END, FAIL, ALL (= BEGIN,END,FAIL)
   #SBATCH --mail-type=NONE

    f.write("export OMP_NUM_THREADS=1" + os.linesep)
    f.write("source ~/.bashrc &>/dev/null" + os.linesep)
    f.write("mpirun --map-by socket ~/nobackup/feat-build/applications/poisson_dirichlet  --level " +  str(11+i) + " 7 0 --mesh ~/feat3/data/meshes/unit-square-quad.xml" + os.linesep)
    #f.write("mpirun --map-by socket ~/nobackup/feat-build/applications/poisson_solver_factory  --level " +  str(10+i) + " 5 --mesh ~/feat3/data/meshes/unit-square-quad.xml --solver-ini ~/feat3/data/ini/solver_example.ini" + os.linesep)

  #sbatch
  subprocess.call(["sbatch", "temp"])
  time.sleep(1)

os.remove("temp")
