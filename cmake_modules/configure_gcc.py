import platform
import sys
from cmake_modules.feast_util import get_output

def configure_gcc(cpu, buildid):
  version = get_output("g++ -dM -E - ")
  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__GNUC__"])
  minor = int(version["__GNUC_MINOR__"])
  minor2 = int(version["__GNUC_PATCHLEVEL__"])
  print ("Detected gcc version: " + str(major) + " " + str(minor) + " " + str(minor2))

  if major < 4  or (major == 4 and minor < 7):
    print ("GNU Compiler version less then 4.7 is not supported, please update your compiler!")
    sys.exit(1)

  cxxflags = "-pipe -std=c++11 -ggdb -gdwarf-4"

  if major >= 4  and minor >= 9:
    cxxflags += " -fdiagnostics-color=always"

  if "coverage" in buildid:
    cxxflags += " -fprofile-arcs -ftest-coverage"

  # For 'quadmath', we need to enable extended numeric literals, as otherwise the g++ will
  # not recognise the 'q' suffix for __float128 constants when compiling with --std=c++11 starting from gcc version 4.8.
  if major >= 4 and minor >= 8 and "quadmath" in buildid:
    cxxflags += " -fext-numeric-literals"

  if "debug" in buildid:
    cxxflags += " -O0 -Wall -Wextra -Wundef -Wshadow -Woverloaded-virtual -Wuninitialized -fdiagnostics-show-option -fno-omit-frame-pointer"
    #do not use stl debug libs under darwin, as these are as buggy as everything else in macos
    if platform.system() != "Darwin":
      cxxflags += " -D_GLIBCXX_DEBUG"
    if major >= 4  and minor >= 8 and "serial" in buildid and not "cuda" in buildid:
      cxxflags += " -fsanitize=address"
    if major >= 4  and minor >= 9:
      cxxflags += " -fsanitize=undefined"
  elif "opt" in buildid:
    cxxflags += " -O3"
    if cpu == "unknown":
      cxxflags += " -march=native"
    # INTEL
    elif cpu == "i486":
      cxxflags += " -march=i486 -m32"
    elif cpu == "pentium":
      cxxflags += " -march=pentium -m32"
    elif cpu == "pentiumpro":
      cxxflags += " -march=pentiumpro -m32"
    elif cpu == "pentium2":
      cxxflags += " -march=pentium2 -m32"
    elif cpu == "pentium3":
      cxxflags += " -march=pentium3 -m32"
    elif cpu == "pentiumm":
      cxxflags += " -march=pentium-m -m32"
    elif cpu == "pentium4m":
      cxxflags += " -march=pentium4m -m32"
    elif cpu == "coresolo":
      cxxflags += " -march=prescott -msse2"
    elif cpu == "coreduo":
      cxxflags += " -march=core2 -m64"
    elif cpu == "penryn":
      cxxflags += " -march=core2 -msse4.1 -m64"
    elif cpu == "nehalem":
      if major == 4 and minor < 6:
        cxxflags += " -march=core2 -msse4.2 -m64"
      else:
        cxxflags += " -march=corei7 -m64"
    elif cpu == "westmere":
      if major == 4 and minor < 6:
        cxxflags += " -march=core2 -msse4.2 -m64"
      else:
        cxxflags += " -march=corei7 -msse4.2 -m64"
    elif cpu == "sandybridge":
      if major == 4 and minor < 6:
        cxxflags += " -march=core2 -msse4.2 -m64"
      else:
        if platform.system() == "Darwin":
          cxxflags += " -march=corei7 -msse4 -msse4.1 -msse4.2 -m64"
        else:
          cxxflags += " -march=corei7-avx -msse4 -msse4.1 -msse4.2 -m64"
    elif cpu == "ivybridge":
      if major == 4 and minor < 6:
        cxxflags += " -march=core2 -msse4.2 -m64"
      else:
        if platform.system() == "Darwin":
          cxxflags += " -march=corei7 -msse4 -msse4.1 -msse4.2 -m64"
        else:
          cxxflags += " -march=corei7-avx -msse4 -msse4.1 -msse4.2 -m64"
    elif cpu == "haswell":
      if major == 4 and minor < 6:
        cxxflags += " -march=core2 -msse4.2 -m64"
      else:
        if platform.system() == "Darwin":
          cxxflags += " -march=corei7 -msse4 -msse4.1 -msse4.2 -m64"
        else:
          cxxflags += " -march=core-avx2 -msse4 -msse4.1 -msse4.2 -m64"
    elif cpu == "itanium":
      cxxflags += " -march=itanium"
    elif cpu == "pentium4":
      cxxflags += " -march=pentium4m -m64"
    elif cpu == "nocona":
      cxxflags += " -march=nocona -m64"
    elif cpu == "itanium2":
      cxxflags += " -march=itanium2"

    # AMD
    elif cpu == "amd486":
      cxxflags += " -m32"
    elif cpu == "k5":
      cxxflags += " -m32"
    elif cpu == "k6":
      cxxflags += " -m32 -march=k6"
    elif cpu == "athlon":
      cxxflags += " -m32 -march=athlon"
    elif cpu == "athlonxp":
      cxxflags += " -m32 -march=athlonxp"
    elif cpu == "opteron":
      cxxflags += " -m64 -march=k8"
    elif cpu == "athlon64":
      cxxflags += " -m64 -march=k8"
    elif cpu == "opteronx2":
      cxxflags += " -m64 -march=k8-sse3"
    elif cpu == "turionx2":
      cxxflags += " -m64 -march=k8-sse3"
    elif cpu == "turionx2":
      cxxflags += " -m64 -march=barcelona"
    elif cpu == "barcelona":
      cxxflags += " -m64 -march=barcelona"
    elif cpu == "shanghai":
      cxxflags += " -m64 -march=barcelona"
    elif cpu == "istanbul":
      cxxflags += " -m64 -march=barcelona"
    elif cpu == "magnycours":
      cxxflags += " -m64 -march=barcelona"

    else:
      cxxflags += " -march=native"
      print ("Detected cpu type not supported by configure_gcc.py, using -march=native instead.")

  return cxxflags
