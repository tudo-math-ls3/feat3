# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

import subprocess

def detect_cuda_arch():
  cuda_arch = "native"
  cuda_arch_int = 1000000
  call_args = ("nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader")
  popen = subprocess.Popen(call_args, stdout=subprocess.PIPE, text=True)
  popen.wait()
  output = popen.stdout.read()
  if("command not found..." in output or "error" in output.lower()):
    cuda_arch_int = 60
  else:
    output_split_glob = output.splitlines()
    for output_str in output_split_glob:
      output_split = output_str.split(".")
      if(len(output_split) != 2):
        print("ERROR: Could not parse output of nvidia-smi querry" + output)
        sys.exit(1)
      cuda_arch_int = min(cuda_arch_int, int(output_split[0])*10 + int(output_split[1]))
  if(cuda_arch_int < 60):
    print("Warning: Lowest device compute capability at " + str(cuda_arch_int) + "\n feat cuda will not work as expected")
    cuda_arch_int = 60

  return cuda_arch, cuda_arch_int
