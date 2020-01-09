# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "January 2020"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ZFP(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["zfp"]
    self.dirname = "zfp"
    self.filename = "zfp-0.5.5.tar.gz"
    self.url = "https://computing.llnl.gov/projects/floating-point-compression/download/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_ZFP:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
