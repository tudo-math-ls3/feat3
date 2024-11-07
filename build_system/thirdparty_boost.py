# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "Oct 2022"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class BOOST(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "1.80.0"
    self.names = ["boost"]
    self.dirname = "boost_1_80_0"
    self.filename = "boost_1_80_0.tar.gz"
    self.url = "https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_BOOST:BOOL=ON -DFEAT_DIRNAME_BOOST:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
    self.cmake_find_package_avail = True