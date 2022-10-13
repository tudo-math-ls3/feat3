# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2022 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Peter Zajac"
__date__   = "October 2022"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class Zoltan(ThirdpartyPackage):
  def __init__(self,trunk_dirname):
    self.version = "3.901"
    self.names = ["zoltan"]
    self.dirname = "Zoltan-" + self.version
    self.filename = "Zoltan-" + self.version + ".tar.gz"
    self.url = "https://github.com/sandialabs/Zoltan/archive/refs/tags/v" + self.version + ".tar.gz"
    self.cmake_flags = " -DFEAT_HAVE_ZOLTAN:BOOL=ON -DFEAT_DIRNAME_ZOLTAN:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
