# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock, Peter Zajac"
__date__   = "April 2017"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ZLIB(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "1.2.13"
    self.names = ["zlib"]
    self.dirname = "zlib-" + self.version
    self.filename = "zlib-" + self.version + ".tar.gz"
    self.url = "http://www.zlib.net/fossils/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_ZLIB:BOOL=ON -DFEAT_DIRNAME_ZLIB:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
    self.cmake_find_package_avail = False
