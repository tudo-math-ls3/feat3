# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "July 2015"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class PARMETIS(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["parmetis", "metis"]
    self.dirname = "parmetis"
    self.filename = "parmetis-4.0.3.tar.gz"
    self.url = "http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_PARMETIS:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
