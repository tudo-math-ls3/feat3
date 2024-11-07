# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Maximilian Esser"
__date__   = "October 2022"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ZFP(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "1.0.0"
    self.names = ["zfp"]
    self.dirname = "zfp-" + self.version
    self.filename ="zfp-" + self.version + ".zip"
    self.url ="https://github.com/LLNL/zfp/archive/" + self.version + ".zip"
    self.cmake_flags = " -DFEAT_HAVE_ZFP:BOOL=ON -DFEAT_DIRNAME_ZFP:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
    self.cmake_find_package_avail = False
