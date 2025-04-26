# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import platform
import sys
from build_system.feat_util import get_output

def configure_icx(cpu, buildid, compiler, system_host_compiler, restrict_errors):
  version = get_output(compiler + " --version")
  version = version[0].split()[4]
  major = int(version.split('.')[0])
  minor = int(version.split('.')[1])
  minor2 = int(version.split('.')[2])
  print ("Detected icx version: " + str(major) + " " + str(minor) + " " + str(minor2))

  cxxflags = "-std=c++17 -g -Wall -Wextra -Wdeprecated -Wnon-virtual-dtor -Wpointer-arith -Wreturn-type -Wshadow -Wshorten-64-to-32 -Wuninitialized"

  if restrict_errors:
    cxxflags += " -diag-error-limit1"

  if system_host_compiler:
    cxxflags += " -gcc-name=" + system_host_compiler

  if "debug" in buildid or "nopt" in buildid:
    cxxflags += "  -O0 -debug all -fp-model=strict"
    #TODO missing ftrapuv, not know by icpx: Initializes stack local variables to an unusual value to aid error detection.

  elif "opt" in buildid:
    cxxflags += " -O3 -fp-model=precise"

    if "lto" in buildid :
      cxxflags += " -ipo "

    if cpu == "unknown":
      # generate code for every simd unit, existing so far
      cxxflags += " -xhost"
      print ("Warning: cpu type not detected, using -xhost instead.")

    # INTEL
    elif cpu == "haswell":
      cxxflags += " -xhaswell -m64"
    elif cpu == "broadwell":
      cxxflags += " -xbroadwell -m64"
    elif cpu == "skylake":
      cxxflags += " -xskylake -m64"
    elif cpu == "skylake-sp":
      cxxflags += " -xskylake-avx512 -m64"
    elif cpu == "kaby-lake":
      cxxflags += " -xskylake -m64" #no special kabylake support, yet
    elif cpu == "ice-lake":
      cxxflags += " -xicelake-server -m64"
    elif cpu == "coffee-lake":
      cxxflags += " -xskylake -m64"
    elif cpu == "cascadelake":
      cxxflags += " -xcascadelake -m64"
    elif cpu == "sapphirerapids":
      cxxflags += " -xsapphirerapids -m64"
    elif cpu == "alder-lake":
      cxxflags += " -xalderlake -m64"

    # AMD
    elif cpu == "zen":
      cxxflags += " -mavx2 -mfma -axCORE-AVX2,CORE-AVX512" # not really supported by intel
    elif cpu == "zen2":
      cxxflags += " -mavx2 -mfma -axCORE-AVX2,CORE-AVX512" # not really supported by intel
    elif cpu == "zen3":
      cxxflags += " -mavx2 -mfma -axCORE-AVX2,CORE-AVX512" # not really supported by intel
    elif cpu == "zen4":
      cxxflags += " -mavx2 -mfma -axCORE-AVX2,CORE-AVX512" # not really supported by intel
    elif cpu == "zen5":
      cxxflags += " -mavx2 -mfma -axCORE-AVX2,CORE-AVX512" # not really supported by intel
    else:
      print ("Warning: Detected cpu type not supported by configure_icx.py, using generic vectorisation support instead")
      # generate code for every simd unit, existing so far
      cxxflags += " -mavx2 -mfma -axCORE-AVX2,CORE-AVX512"

  return cxxflags
