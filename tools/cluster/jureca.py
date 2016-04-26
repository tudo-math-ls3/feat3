#!/usr/bin/env python
# vim: set filetype=python sw=2 sts=2 et nofoldenable :
# This is a python script.
# If you encounter problemes when executing it on its own, start it with a python interpreter
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
    f.write("#SBATCH --job-name=feast" + os.linesep)

    f.write("export OMP_NUM_THREADS=1" + os.linesep)
    #f.write("export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}" + os.linesep)
    f.write("source ~/.bashrc &>/dev/null" + os.linesep)
    f.write("export SCOREP_TOTAL_MEMORY=7g" + os.linesep)
    f.write("export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=40" + os.linesep)
    f.write("export SCOREP_ENABLE_TRACING=true" + os.linesep)
    f.write("export SCOREP_FILTERING_FILE=stokes_poiseuille_2d.filt" + os.linesep)
    f.write("cd ~/feast/applications/" + os.linesep)
    f.write("srun parti_poisson_dirichlet_2d  --level " +  str(9+i) + " 4 --part_min_elems 500" + os.linesep)

  #sbatch
  subprocess.call(["sbatch", "temp"])
  time.sleep(1)

os.remove("temp")
