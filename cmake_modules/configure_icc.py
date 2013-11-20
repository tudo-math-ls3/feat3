import subprocess
import platform

def configure_icc(cpu, buildmode):
  if "check_output" not in dir( subprocess ): #deactivated as its not available bevor python 2.7
    pipe = subprocess.Popen("icpc -dM -E - ".split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    pipe.stdin.close()
    version = pipe.stdout.read().splitlines()
  else:
    version = subprocess.check_output("echo | icpc -dM -E - ", shell=True).splitlines()

  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__INTEL_COMPILER"][0:2])
  minor = int(version["__INTEL_COMPILER"][-2:])
  minor2 = int(version["__INTEL_COMPILER_UPDATE"])

  cxxflags = "-std=c++11 -g"
  if buildmode == "debug":
    cxxflags += "  -O0 -Wall -Wp64 -Wshorten-64-to-32 -debug all -ftrapuv"
  elif "opt" in buildmode:
    cxxflags += " -O3"
    if cpu == "unknown":
      pass

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
      cxxflags += " -xsse3"
    elif cpu == "coreduo":
      cxxflags += " -xsse3"
    elif cpu == "pentryn":
      cxxflags += " -xsse4.1"
    elif cpu == "nehalem":
      cxxflags += " -xsse4.2"
    elif cpu == "westmere":
      cxxflags += " -xsse4.2"
    elif cpu == "sandybridge":
      cxxflags += " -xsse4.2 -axAVX"
    elif cpu == "ivybridge":
      cxxflags += " -xsse4.2 -axAVX"
    elif cpu == "haswell":
      cxxflags += " -xsse4.2 -axCORE-AVX2"
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
      cxxflags += " -msse4.1"
    elif cpu == "shanghai":
      cxxflags += " -msse4.1"
    elif cpu == "istanbul":
      cxxflags += " -msse4.1"
    elif cpu == "magnycours":
      cxxflags += " -msse4.1"
    else:
      print ("Detected cpu type not supported by configure_gcc.py")

  return cxxflags
