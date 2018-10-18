import platform
from build_system.feat_util import get_output

def configure_pgi(cpu, buildid, compiler, restrict_errors):
  # TODO min versiong 16.4

  cxxflags = "--c++11"

  if restrict_errors:
    cxxflags += " -e1"

  if "debug" in buildid:
    cxxflags += "  -O0 -g"

  elif "opt" in buildid:
    cxxflags += " -O3 -mcmodel=medium -Minline -fastsse -Mipa=fast,inline -Msmartalloc -Mprefetch -Msafeptr" #-gopt

  return cxxflags
