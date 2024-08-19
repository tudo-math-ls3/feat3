# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import platform
from build_system.feat_util import get_output

def configure_icc(cpu, buildid, compiler, system_host_compiler, restrict_errors):
  version = get_output(compiler + " -dM -E - ")
  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__INTEL_COMPILER"][0:2])
  minor = int(version["__INTEL_COMPILER"][-2:])
  minor2 = int(version["__INTEL_COMPILER_UPDATE"])
  print ("Detected icc version: " + str(major) + " " + str(minor) + " " + str(minor2))

  if major < 15:
    print ("Error: Intel Compiler version less then 15 is not supported, please update your compiler!")
    sys.exit(1)

  cxxflags = "-std=c++17 -g -Wall -Wextra -Wcheck -Wdeprecated -Wnon-virtual-dtor -Wpointer-arith -Wreturn-type -Wshadow -Wp64 -Wshorten-64-to-32 -Wuninitialized -diag-disable 2304,2305"

  if restrict_errors:
    cxxflags += " -diag-error-limit1"

  if system_host_compiler:
    cxxflags += " -gcc-name=" + system_host_compiler

  # do not use stone old clang libc++ from apple with icc14, hopefully any gcc headers are present
  if platform.system() == "Darwin" and major == 14:
    cxxflags += " -no-use-clang-env"

  if "debug" in buildid or "nopt" in buildid:
    cxxflags += "  -O0 -debug all -ftrapuv"

  elif "opt" in buildid:
    cxxflags += " -no-prec-div -diag-disable 11074,11076,25463,25464 -diag-disable openmp"
    if "lto" in buildid and major >= 18:
      cxxflags += " -ipo -diag-disable 11000,11001,11006"
    # we don't set -ip because it causes horrific compile times during nightly regression tests
    #else:
    #  cxxflags += " -ip"

    if "opt" in buildid:
      cxxflags += " -O3"

    if cpu == "unknown":
      # generate code for every simd unit, existing so far
      cxxflags += " -axCORE-AVX512"
      print ("Warning: cpu type not detected, using generic vectorisation support instead.")

    # INTEL
    elif cpu == "i486":
      pass
    elif cpu == "pentium":
      cxxflags += " -mia32"
    elif cpu == "pentiumpro":
      cxxflags += " -mia32"
    elif cpu == "pentium2":
      cxxflags += " -mia32"
    elif cpu == "pentium3":
      cxxflags += " -mia32"
    elif cpu == "pentiumm":
      cxxflags += " -xsse2"
    elif cpu == "pentiu4m":
      cxxflags += " -xsse2"
    elif cpu == "coresolo":
      cxxflags += " -xsse2"
    elif cpu == "coreduo":
      cxxflags += " -xsse3"
    elif cpu == "penryn":
      cxxflags += " -xsse4.1"
    elif cpu == "nehalem":
      cxxflags += " -xsse4.2"
    elif cpu == "westmere":
      cxxflags += " -xsse4.2"
    elif cpu == "sandybridge":
      cxxflags += " -xAVX"
    elif cpu == "ivybridge":
      cxxflags += " -xCORE-AVX-I"
    elif cpu == "haswell":
      cxxflags += " -xCORE-AVX2"
    elif cpu == "broadwell":
      cxxflags += " -xCORE-AVX2"
    elif cpu == "knightslanding":
      cxxflags += " -xMIC-AVX512"
    elif cpu == "skylake":
      cxxflags += " -xCORE-AVX2"
    elif cpu == "skylake-sp":
      cxxflags += " -xCORE-AVX512"
    elif cpu == "cascadelake":
      cxxflags += " -xCORE-AVX512"
    elif cpu == "coffee-lake":
      cxxflags += " -xCORE-AVX2"
    elif cpu == "kaby-lake":
      cxxflags += " -xCORE-AVX2"
    elif cpu == "ice-lake":
      cxxflags += " -xCORE-AVX2"
    elif cpu == "sapphirerapids":
      cxxflags += " -xCORE-AVX512"
    elif cpu == "alder-lake":
      cxxflags += " -xCORE-AVX512"
    elif cpu == "itanium":
      # no setting necessary, the itanium version of the intel compiler
      # sets everything automatically
      pass
    elif cpu == "pentium4":
      cxxflags += " -xsse2"
    elif cpu == "nocona":
      cxxflags += " -xsse3"
    elif cpu == "itanium":
      # no setting necessary, the itanium version of the intel compiler
      # sets everything automatically
      pass

    # AMD
    elif cpu == "amd486":
      cxxflags += " -mia32"
    elif cpu == "k5":
      cxxflags += " -mia32"
    elif cpu == "k6":
      cxxflags += " -mia32"
    elif cpu == "athlon":
      cxxflags += " -mia32"
    elif cpu == "athlonxp":
      cxxflags += " -mia32"
    elif cpu == "opteron":
      cxxflags += " -msse2"
    elif cpu == "athlon64":
      cxxflags += " -msse2"
    elif cpu == "opteronx2":
      cxxflags += " -msse3"
    elif cpu == "turionx2":
      cxxflags += " -msse3"
    elif cpu == "barcelona":
      cxxflags += " -msse3"
    elif cpu == "shanghai":
      cxxflags += " -msse3"
    elif cpu == "istanbul":
      cxxflags += " -msse3"
    elif cpu == "magnycours":
      cxxflags += " -msse3"
    elif cpu == "zen":
      cxxflags += " -march=core-avx2" # not really supported by intel
    elif cpu == "zen2":
      cxxflags += " -march=core-avx2" # not really supported by intel
    elif cpu == "zen3":
      cxxflags += " -march=core-avx2" # not really supported by intel
    elif cpu == "zen4":
      cxxflags += " -march=core-avx2" # not really supported by intel
    else:
      print ("Warning: Detected cpu type not supported by configure_icc.py, using generic vectorisation support instead")
      # generate code for every simd unit, existing so far
      cxxflags += " -axCORE-AVX512"

  return cxxflags
