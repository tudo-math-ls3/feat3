# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Peter Zajac"
__date__   = "April 2018"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class Triangle(ThirdpartyPackage):
  def __init__(self,trunk_dirname):
    self.names = ["triangle"]
    self.dirname = "triangle"
    self.filename = "triangle.zip"
    self.url = "http://www.netlib.org/voronoi/triangle.zip"
    self.cmake_flags = " -DFEAT_HAVE_TRIANGLE:BOOL=ON -DFEAT_DIRNAME_TRIANGLE:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = os.path.join(trunk_dirname, self.dirname)
