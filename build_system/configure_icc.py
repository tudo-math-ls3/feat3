import platform
from build_system.feat_util import get_output

def configure_icc(cpu, buildid, compiler, system_host_compiler):
  version = get_output(compiler + " -dM -E - ")
  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__INTEL_COMPILER"][0:2])
  minor = int(version["__INTEL_COMPILER"][-2:])
  minor2 = int(version["__INTEL_COMPILER_UPDATE"])
  print ("Detected icc version: " + str(major) + " " + str(minor) + " " + str(minor2))

  if major < 15:
    print ("Error: Intel Compiler version less then 15 is not supported, please update your compiler!")
    sys.exit(1)

  cxxflags = "-std=c++11 -g -Wall -Wextra -Wcheck -Wdeprecated -Wnon-virtual-dtor -Wpointer-arith -Wreturn-type -Wshadow -Wp64 -Wshorten-64-to-32 -Wuninitialized -diag-disable 2304,2305"
  if system_host_compiler:
    cxxflags += " -gcc-name=" + system_host_compiler

  # do not use stone old clang libc++ from apple with icc14, hopefully any gcc headers are present
  if platform.system() == "Darwin" and major == 14:
    cxxflags += " -no-use-clang-env"

  if "debug" in buildid or "noop" in buildid:
    cxxflags += "  -O0 -debug all -ftrapuv"

  elif "opt" in buildid or "fast" in buildid:
    cxxflags += " -no-prec-div"
    if "lto" in buildid:
      cxxflags += " -ipo -diag-disable 11074,11076,11000,11001,11006"
    else:
      cxxflags += " -ip -diag-disable 11074,11076"

    if "opt" in buildid:
      cxxflags += " -O3"
    elif "fast" in buildid:
      cxxflags += " -Ofast"

    if cpu == "unknown":
      # generate code for every simd unit, existing so far
      cxxflags += " -axCORE-AVX2"
      print ("Warning: cpu type not detected, using generic sse support instead.")

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
    # Haswell, Broadwell and Skylake do support -xcore-AVX2, but all
    # Intel compilers up to 16.0.4 produce code with significantly
    # lower precision, most notably in kernel/solver/optimiser-test.
    # This is the case for both Haswell and Skylake CPUs.
    # At this point it is not clear if there is a compiler bug or this
    # is normal behaviour and AVX2 reduces the precision, so deactivate
    # this for the time being.
    elif cpu == "haswell":
      cxxflags += " -xCORE-AVX-I"
    elif cpu == "broadwell":
      cxxflags += " -xCORE-AVX-I"
    elif cpu == "knightslanding":
      cxxflags += " -xMIC-AVX512"
    elif cpu == "skylake":
      cxxflags += " -xCORE-AVX-I"
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
    else:
      print ("Warning: Detected cpu type not supported by configure_icc.py, using generic sse support instead")
      # generate code for every simd unit, existing so far
      cxxflags += " -axCORE-AVX2"

  return cxxflags
