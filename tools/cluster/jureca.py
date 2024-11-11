#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

import subprocess
import sys
import os
import time

sys.dont_write_bytecode = True

for i in range(1, 5):
  with open("temp", "w") as f:
    f.write("#!/bin/bash -x" + os.linesep)

    f.write("#SBATCH --nodes=" + str(4**i) + os.linesep)

    f.write("#SBATCH --ntasks-per-node=16" + os.linesep)
    f.write("#SBATCH --time=00:20:00" + os.linesep)
    f.write("#SBATCH --partition=batch" + os.linesep)
    # http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JURECA/UserInfo/QuickIntroduction.html?nn=1803700#ReqGres
    #f.write("#SBATCH --partition=gpus" + os.linesep)
    # http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JURECA/UserInfo/GpuNodes.html?nn=1803700
    #f.write("#SBATCH --gres=gpu:1" + os.linesep)
    f.write("#SBATCH --job-name=feat" + os.linesep)

    f.write("export PSP_ONDEMAND=1" + os.linesep)
    f.write("export OMP_NUM_THREADS=1" + os.linesep)
    #f.write("export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" + os.linesep)
    f.write("source ~/.bashrc &>/dev/null" + os.linesep)
    f.write("export SCOREP_TOTAL_MEMORY=3900MB" + os.linesep)
    f.write("export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=40" + os.linesep)
    #f.write("export SCOREP_ENABLE_PROFILING=true" + os.linesep)
    f.write("export SCOREP_ENABLE_TRACING=true" + os.linesep)
    f.write("export SCOREP_FILTERING_FILE=" + os.getcwd() + "/scorep.filt" + os.linesep)
    f.write("srun --hint=nomultithread --cpu_bind=map_cpu:0,1,2,3,4,5,6,7,12,13,14,15,16,17,18,19 ~/feat/applications/poisson_dirichlet_2d  --level " +  str(11+i) + " 4 --part_min_elems 1000" + os.linesep)

  #sbatch
  subprocess.call(["sbatch", "temp"])
  time.sleep(1)

os.remove("temp")
