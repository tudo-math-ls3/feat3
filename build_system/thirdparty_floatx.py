# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Peter Zajac"
__date__   = "2018-03-19"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class FParser(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["floatx"]
    self.dirname = "floatx"
    self.filename = "FloatX-develop.zip"
    self.url = "https://codeload.github.com/oprecomp/FloatX/zip/develop"
    self.cmake_flags = " -DFEAT_HAVE_FLOATX:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
