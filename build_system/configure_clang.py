# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import platform
import sys
from build_system.feat_util import get_output

def configure_clang(cpu, buildid, compiler, system_host_compiler, restrict_errors):
  version = get_output(compiler + " -dM -E - ")
  version = dict(map(lambda x : (x[1], " ".join(x[2:])), [line.split() for line in version]))
  major = int(version["__clang_major__"])
  minor = int(version["__clang_minor__"])
  minor2 = int(version["__clang_patchlevel__"])
  print ("Detected clang version: " + str(major) + " " + str(minor) + " " + str(minor2))

  if major < 3  or (major == 3 and minor < 3):
    print ("Error: Clang Compiler version less then 3.3 is not supported, please update your compiler!")
    sys.exit(1)

  standard = "-std=c++17"
  if(platform.system() == "Windows"):
    standard = "-std=c++14"
  cxxflags = standard + " -ggdb -fcolor-diagnostics -m64 -Wall -Wextra -Wshadow -Wundef -Wshorten-64-to-32 -Wconversion -Wstrict-aliasing=2 -Wunknown-pragmas -Wundef -Wuninitialized -Wswitch -Wunused-label -Woverloaded-shift-op-parentheses -Wempty-body -Wheader-guard -Wimplicit-fallthrough -Wloop-analysis -Wheader-hygiene -Wpedantic"

  if(platform.system() != "Windows"):
    cxxflags += " -pthread"


  if restrict_errors:
    cxxflags += " -Wfatal-errors"

  if major > 3 or (major == 3 and minor > 6):
    cxxflags += " -Wrange-loop-analysis -Wobjc-circular-container"

  if major > 3 or (major == 3 and minor >= 9):
    cxxflags += " -Wcomma"

  if major >= 5:
    cxxflags += " -Wcast-qual -Wunused-lambda-capture -Wstrict-prototypes"

  if major >= 6:
    cxxflags += " -Wtautological-compare"

  if major >= 7:
    cxxflags += " -Woverloaded-virtual -Wnon-virtual-dtor"

  if major >= 8:
    cxxflags += " -Wextra-semi-stmt -Wextra-semi"

  if major >= 10:
    cxxflags += " -Wtautological-overlap-compare -Wtautological-bitwise-compare -Wmisleading-indentation -Wsizeof-array-div -Wsizeof-pointer-div -Wxor-used-as-pow -Wfinal-dtor-non-final-class"

  if major >= 11:
    cxxflags += " -Wpointer-to-int-cast"

  if major >= 13:
    cxxflags += " -Wunused-but-set-parameter -Wunused-but-set-variable -Wnull-pointer-subtraction"

  if major >= 16:
    cxxflags += " -Wno-gnu-line-marker"

  if system_host_compiler:
    cxxflags += " --gcc-toolchain=" + system_host_compiler

  if "ccache" in buildid:
    cxxflags += " -Qunused-arguments"

  if "deathhandler" in buildid:
    cxxflags += " -ldl"

  if "debug" in buildid or "nopt" in buildid:

    if major < 4:
      cxxflags += " -O0"
    else:
      cxxflags += " -Og"

    if major >= 5:
      if "sanitizer" in buildid:
        cxxflags += " -fsanitize=pointer-overflow -fsanitize=nullability"

    if major >= 7:
      if "sanitizer" in buildid:
        cxxflags += " -fsanitize=implicit-conversion"

    cxxflags += " -ftemplate-backtrace-limit=0 -fdiagnostics-show-template-tree -fdiagnostics-show-category=name -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    if platform.system() != "Darwin":
      if "sanitizer" in buildid:
        cxxflags += " -fsanitize=undefined" # darwin clang does not like sanitize=undefined
    if "mpi" not in buildid and "cuda" not in buildid and "valgrind" not in buildid and "xcode" not in buildid:
      if "sanitizer" in buildid:
        cxxflags += " -fsanitize=address" #" -fsanitize=memory" #-fsanitize=address-full  #Problem with clang LLVM on windows with VS22..
      #see https://github.com/llvm/llvm-project/issues/56300

  elif "opt" in buildid or "fast" in buildid:
    if major >= 7:
      cxxflags += " -funroll-loops"
      if "lto" in buildid:
        cxxflags += " -flto"

    if "opt" in buildid:
      cxxflags += " -O3"
    elif "fast" in buildid:
      if major == 3 and minor < 7:
        cxxflags += " -O3"
      else:
        # work around for https://llvm.org/bugs/show_bug.cgi?id=13745 - math.h/math-finite-h broken in gcc 4.6 (host compiler of ubuntu 12.04)
        cxxflags += " -Ofast -D__extern_always_inline='extern __always_inline' -ffp-model=fast"

    if cpu == "unknown":
      cxxflags += " -march=native"
      print ("Warning: cpu type not detected, using -march=native instead.")

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
    elif cpu == "broadwell":
      cxxflags += " -march=core-avx2"
    elif cpu == "skylake":
      cxxflags += " -march=core-avx2"
    elif cpu == "skylake-sp":
      cxxflags += " -march=core-avx2"
    elif cpu == "cascadelake":
      cxxflags += " -march=core-avx2"
    elif cpu == "coffee-lake":
      cxxflags += " -march=core-avx2"
    elif cpu == "kaby-lake":
      cxxflags += " -march=core-avx2"
    elif cpu == "ice-lake":
      cxxflags += " -march=core-avx2"
    elif cpu == "sapphirerapdis":
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
      cxxflags += " -m32 -march=athlon-xp"
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
    elif cpu == "zen":
      cxxflags += " -m64 -march=znver1"
    elif cpu == "zen2":
      if major >= 9:
        cxxflags += " -m64 -march=znver2"
      else:
        cxxflags += " -m64 -march=znver1"
    elif cpu == "zen3":
      if major >= 12:
        cxxflags += " -m64 -march=znver3"
      elif major >= 9:
        cxxflags += " -m64 -march=znver2"
      else:
        cxxflags += " -m64 -march=znver1"
    elif cpu == "zen4":
      cxxflags += " -m64 -march=znver4"
    #ARM
    elif cpu == "cortexa53":
      cxxflags += " -march=a53"
    elif cpu == "armv8":
      cxxflags += " -m64 -march=v8.2a"
    elif cpu == "a64fx":
      cxxflags += " -mcpu=a64fx"
    elif cpu == "thunderx2":
      cxxflags += " -mcpu=thunderx2t99"
    else:
      cxxflags += " -march=native"
      print ("Warning: Detected cpu type not supported by configure_clang.py, using -march=native instead.")

  return cxxflags
