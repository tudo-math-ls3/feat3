# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Jordi Paul, Peter Zajac"
__date__   = "April 2014"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class SuiteSparse(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["umfpack"]
    self.dirname = "SuiteSparse"
    self.filename = "SuiteSparse-5.2.0.tar.gz"
    self.url = "http://faculty.cse.tamu.edu/davis/SuiteSparse/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_UMFPACK:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
