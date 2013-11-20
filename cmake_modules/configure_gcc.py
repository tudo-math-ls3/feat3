import subprocess
import platform

def configure_gcc(cpu, buildmode):
  if "check_output" not in dir( subprocess ): #deactivated as its not available bevor python 2.7
    pipe = subprocess.Popen("g++ -dM -E - ".split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pipe.stdin.close()
    version = pipe.stdout.read().splitlines()
  else:
    version = subprocess.check_output("echo | g++ -dM -E - ", shell=True).splitlines()

  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__GNUC__"])
  minor = int(version["__GNUC_MINOR__"])
  minor2 = int(version["__GNUC_PATCHLEVEL__"])

  if major <= 4 and minor <= 4:
    print ("GNU Compiler version less then 4.4 is not supported, please update your compiler")

  cxxflags = "-pipe -std=c++0x -ggdb"

  if "coverage" in buildmode:
    cxxflags += " -fprofile-arcs -ftest-coverage"
  if "debug" in buildmode:
    cxxflags += " -O0 -Wall -Wextra -Wundef -Wno-unused-parameter"
  elif "opt" in buildmode:
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
      cxxflags += " -march=native -m64"
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
          cxxflags += " -march=corei7-avx -msse4 -msse4.1 -msse4.2 -m64"
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
      print "Detected cpu type not supported by configure_gcc.py"

  return cxxflags
