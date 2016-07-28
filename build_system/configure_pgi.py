import platform
from build_system.feat_util import get_output

def configure_pgi(cpu, buildid, compiler):
  # TODO min versiong 16.4

  cxxflags = "--c++11"
  if "debug" in buildid:
    cxxflags += " -g -O0"

  elif "opt" in buildid:
    cxxflags += " -O3 -gopt -mcmodel=medium -Minline -fastsse -Mipa=fast,inline -Msmartalloc -Mprefetch -Msafeptr"

  return cxxflags
