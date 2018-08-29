# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Peter Zajac"
__date__   = "October 2021"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class SuperLU(ThirdpartyPackage):
  def __init__(self,trunk_dirname):
    self.version = "7.1.0"
    self.names = ["superlu"]
    self.dirname = "superlu_dist-" + self.version
    self.filename = "superlu_dist-" + self.version + ".tar.gz"
    self.url = "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v" + self.version + ".tar.gz"
    self.cmake_flags = " -DFEAT_HAVE_SUPERLU_DIST:BOOL=ON -DFEAT_DIRNAME_SUPERLU:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
