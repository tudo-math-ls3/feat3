# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Peter Zajac"
__date__   = "February 2023"
from build_system.thirdparty_package import ThirdpartyPackage
import os
import platform

class SuperLU(ThirdpartyPackage):
  def __init__(self,trunk_dirname):
    self.version = "8.1.1"
    self.names = ["superlu"]
    self.name = "SuperLU"
    self.dirname = "superlu_dist-" + self.version
    self.realdirname = trunk_dirname + "\\" + self.dirname
    self.filename = "superlu_dist-" + self.version + ".tar.gz"
    self.url = "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v" + self.version + ".tar.gz"
    self.cmake_flags = " -DFEAT_HAVE_SUPERLU_DIST:BOOL=ON -DFEAT_DIRNAME_SUPERLU:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
    self.cmake_find_package_avail = False

  def patch(self):
    if platform.system() != "Windows":
      return
    print("\nPatching " + self.name + " " + self.version + " sources...")
    # don't use ParMETIS
    self.patch_file(os.path.join("SRC", "superlu_dist_config.h.in"), [
      [7, '#cmakedefine HAVE_PARMETIS @HAVE_PARMETIS@', '/* #undef HAVE_PARMETIS */'],
    ])
    # don't include unistd.h header
    self.patch_file(os.path.join("SRC", "util.c"), [
      [25, '#include <unistd.h>', '/*#include <unistd.h>*/'],
    ])