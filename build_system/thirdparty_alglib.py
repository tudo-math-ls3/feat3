# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Jordi Paul"
__date__   = "April 2014"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ALGLIB(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["alglib"]
    self.dirname = "alglib"
    self.filename = "alglib-3.13.0.cpp.gpl.zip"
    self.url = "http://www.alglib.net/translator/re/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_ALGLIB:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
