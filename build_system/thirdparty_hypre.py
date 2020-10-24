# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "July 2018"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class HYPRE(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["hypre"]
    self.dirname = "hypre"
    self.filename = "hypre-2.11.2.tar.gz"
    self.url = "https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_HYPRE:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = os.path.join(trunk_dirname, self.dirname)

  def patch(self):
    print("Patching HYPRE sources...")
    # there is a shell script named 'version' in "thirdparty/hypre/hypre-2.11.2/src/utilities"
    # which has to be renamed because it collides with clang's internal <version> header
    x = os.path.join(self.target_dirname, "hypre-2.11.2", "src", "utilities", "version")
    if os.path.isfile(x):
      print("Renaming file '%s'" % x)
      os.rename(x, x + ".sh")
