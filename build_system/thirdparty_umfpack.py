# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Jordi Paul, Peter Zajac"
__date__   = "April 2014"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class SuiteSparse(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "5.8.1"
    self.names = ["umfpack","suitesparse"]
    self.dirname = "SuiteSparse-" + self.version
    self.filename = "SuiteSparse-" + self.version + ".tar.gz"
    self.url = "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v" + self.version + ".tar.gz"
    self.cmake_flags = " -DFEAT_HAVE_UMFPACK:BOOL=ON -DFEAT_DIRNAME_SUITESPARSE:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
