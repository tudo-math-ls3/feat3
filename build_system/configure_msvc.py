# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import platform
from build_system.feat_util import get_output

def configure_msvc(cpu, buildid, compiler, restrict_errors):
  ldflags = []

  cxxflags = [
    "/Zi",
    "/nologo",
    "/W4",
    "/WX-",
    "/diagnostics:column",
    "/Gm-",
    "/EHsc",
    "/fp:except",
    #"/Za",
    "/Zc:wchar_t",
    "/Zc:forScope",
    "/Zc:inline",
    "/GR",
    "/std:c++17",
    "/permissive-",
    "/external:W0",
    "/Gd",
    "/TP",
    "/FC",
    "/errorReport:queue",
    "/bigobj"
  ]

  if "debug" in buildid:
    cxxflags += [
      "/JMC",
      "/Od",
      "/Ob1",
      "/Oy-",
      "/RTC1",
      "/MDd",
      "/GS",
      "/Gy"
    ]

  if "opt" in buildid:
    cxxflags += [
      "/MP",
      "/Ox",
      "/Ob2",
      "/Oi",
      "/Ot",
      "/Oy",
      #"/GL",
      "/MD",
      "/GS-",
      "/Gy-"
    ]

  if "lto" in buildid:
    cxxflags.append("/GL")
    ldflags.append("/LTCG")

  if "omp" in buildid:
    cxxflags.append("/openmp:llvm")
  else:
    cxxflags.append("/openmp-")

  cudaflags = []
  if "cuda" in buildid:
    ldflags.append("/NODEFAULTLIB:libcmt")
    cudaflags += [
      "--std=c++17",
      "--ptxas-options",
      "-suppress-stack-size-warning"
    ]

  return cxxflags, ldflags, cudaflags
