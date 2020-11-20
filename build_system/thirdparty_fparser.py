# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Peter Zajac"
__date__   = "2020-11-20"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class FParser(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "4.5.2"
    self.names = ["fparser"]
    self.dirname = "fparser-" + self.version
    self.filename = "fparser" + self.version + ".zip"
    self.url = "http://warp.povusers.org/FunctionParser/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_FPARSER:BOOL=ON -DFEAT_DIRNAME_FPARSER:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = os.path.join(trunk_dirname, self.dirname)
