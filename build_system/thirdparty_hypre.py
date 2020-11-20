# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock, Peter Zajac"
__date__   = "July 2018"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class HYPRE(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "2.20.0"
    self.names = ["hypre"]
    self.dirname = "hypre-" + self.version
    self.filename = "hypre-" + self.version + ".tar.gz"
    self.url = "https://github.com/hypre-space/hypre/archive/v" + self.version + ".tar.gz"
    self.cmake_flags = " -DFEAT_HAVE_HYPRE:BOOL=ON -DFEAT_DIRNAME_HYPRE:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname

  def patch(self):
    print("Patching HYPRE sources...")
    # there is a shell script named 'version' in "thirdparty/hypre/hypre-2.20.0/src/utilities"
    # which has to be renamed because it collides with clang's internal <version> header
    x = os.path.join(self.target_dirname, self.dirname, "src", "utilities", "version")
    if os.path.isfile(x):
      print("Renaming file '%s'" % x)
      os.rename(x, x + ".sh")
