# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Jordi Paul"
__date__   = "April 2014"
from cmake_modules.thirdparty_package import ThirdpartyPackage
import os

class UMFPACK(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "umfpack"
    self.dirname = "UMFPACK"
    self.filename = "UMFPACK-5.6.2.tar.gz"
    self.url = "http://www.cise.ufl.edu/research/sparse/umfpack/" + self.filename
    self.cmake_flags = " -DFEAST_HAVE_UMFPACK:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname

# Overwrite add function to take care of the dependencies AMD and SuiteSparse_config
  def add(self):
    ThirdpartyPackage.add(self)
    amd = AMD(self.trunk_dirname)
    amd.add()
    suitesparseconfig = SuiteSparse_config(self.trunk_dirname)
    suitesparseconfig.add()

class AMD(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "amd"
    self.dirname = "AMD"
    self.filename = "AMD-2.3.1.tar.gz"
    self.url = "http://www.cise.ufl.edu/research/sparse/amd/" + self.filename
    self.cmake_flags = ""
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname

class SuiteSparse_config(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "suitesparse_config"
    self.dirname = "SuiteSparse_config"
    self.filename = "SuiteSparse_config-4.2.1.tar.gz"
    self.url = "http://www.cise.ufl.edu/research/sparse/SuiteSparse_config/" + self.filename
    self.cmake_flags = ""
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
