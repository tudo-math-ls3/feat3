# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import platform
from build_system.feat_util import get_output

def configure_cray(cpu, buildid, compiler, restrict_errors):

  cxxflags = "-h std=c++11"

  if "debug" in buildid:
    cxxflags += "  -O0 -g"

  elif "opt" in buildid:
    cxxflags += " -O3"

  return cxxflags
