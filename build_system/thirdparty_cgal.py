# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "Nov 2018"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class CGAL(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "4.14"
    self.names = ["cgal"]
    self.dirname = "CGAL-" + self.version
    self.filename = "CGAL-" + self.version + ".zip"
    self.url = "https://github.com/CGAL/cgal/releases/download/releases/CGAL-" + self.version + "/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_CGAL:BOOL=ON -DFEAT_DIRNAME_CGAL:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
