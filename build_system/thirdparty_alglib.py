# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Jordi Paul"
__date__   = "April 2014"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ALGLIB(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "3.16.0"
    self.names = ["alglib"]
    self.dirname = "alglib-" + self.version
    self.filename = "alglib-" + self.version + ".cpp.gpl.zip"
    self.url = "http://www.alglib.net/translator/re/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_ALGLIB:BOOL=ON -DFEAT_DIRNAME_ALGLIB:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = os.path.join(trunk_dirname, self.dirname)
