import platform
from build_system.feast_util import get_output

def configure_clang(cpu, buildid, compiler):
  if not compiler:
    compiler = "clang++"
  version = get_output(compiler + " -dM -E - ")
  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__clang_major__"])
  minor = int(version["__clang_minor__"])
  minor2 = int(version["__clang_patchlevel__"])
  print ("Detected clang version: " + str(major) + " " + str(minor) + " " + str(minor2))

  if major < 3  or (major == 3 and minor < 3):
    print ("Error: Clang Compiler version less then 3.3 is not supported, please update your compiler!")
    sys.exit(1)

  cxxflags = "-pipe  -std=c++11 -ggdb -fcolor-diagnostics -m64"
  if "ccache" in buildid:
    cxxflags += " -Qunused-arguments"
  if "debug" in buildid:
    cxxflags += " -O0 -Wall -Wextra -Wshadow -Wundef -Wshorten-64-to-32 -Wconversion -Wstrict-aliasing=2 -Wunknown-pragmas -Wundef -Wno-unused-value -Wno-unused-parameter -Wuninitialized -fdiagnostics-show-template-tree -fdiagnostics-show-category=name -fno-omit-frame-pointer -fno-optimize-sibling-calls -Wswitch -Wunused-label -Woverloaded-shift-op-parentheses -Wempty-body -Wheader-guard -Wimplicit-fallthrough -Wloop-analysis -Wheader-hygiene"
    if platform.system() != "Darwin":
      cxxflags += " -fsanitize=undefined" # darwin clang does not like sanitize=undefined
    if "mpi" not in buildid and "cuda" not in buildid and "valgrind" not in buildid and "xcode" not in buildid:
      cxxflags += " -fsanitize=address" # -fsanitize=memory" -fsanitize=address-full

  elif "opt" in buildid:
    cxxflags += " -O3"
    if cpu == "unknown":
      cxxflags += " -mtune=generic"
      print ("Warning: cpu type not detected, using -mtune=generic instead.")

    # INTEL
    elif cpu == "i486":
      cxxflags += " -march=i486"
    elif cpu == "pentium":
      cxxflags += " -march=pentium"
    elif cpu == "pentiumpro":
      cxxflags += " -march=pentiumpro"
    elif cpu == "pentium2":
      cxxflags += " -march=pentium2"
    elif cpu == "pentium3":
      cxxflags += " -march=pentium3"
    elif cpu == "pentiumm":
      cxxflags += " -march=pentium-m"
    elif cpu == "pentium4m":
      cxxflags += " -march=pentium4m"
    elif cpu == "coresolo":
      cxxflags += " -march=prescott"
    elif cpu == "coreduo":
      cxxflags += " -march=core2"
    elif cpu == "penryn":
      cxxflags += " -march=core2"
    elif cpu == "nehalem":
      cxxflags += " -march=corei7"
    elif cpu == "westmere":
      cxxflags += " -march=corei7"
    elif cpu == "sandybridge":
      cxxflags += " -march=corei7-avx"
    elif cpu == "ivybridge":
      cxxflags += " -march=corei7-avx"
    elif cpu == "haswell":
      cxxflags += " -march=core-avx2"
    elif cpu == "itanium":
      cxxflags += " -march=itanium"
    elif cpu == "pentium4":
      cxxflags += " -march=pentium4m"
    elif cpu == "nocona":
      cxxflags += " -march=nocona"
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
      print ("Warning: Detected cpu type not supported by configure_clang.py, using -march=native instead.")

  return cxxflags