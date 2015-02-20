import platform
import sys
from build_system.feast_util import get_output

def configure_msc(cpu, buildid):

  # set default flags
  cxxflags = "/DWIN32 /D_WINDOWS /Wall /GR /EHsc /wd4350"
  cxxflags += " /Zi /bigobj /errorReport:none /nologo"

  if "debug" in buildid:
    cxxflags += " /Od /Gm"
  elif "opt" in buildid:
    cxxflafs += " /Ox /GL /MP"

  if "omp" in buildid:
    cxxflags += " /openmp"

  return cxxflags
