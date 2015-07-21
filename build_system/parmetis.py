# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Dirk Ribbrock"
__date__   = "July 2015"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class PARMETIS(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "parmetis"
    self.dirname = "parmetis"
    self.filename = "parmetis-4.0.3.tar.gz"
    self.url = "http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/" + self.filename
    self.cmake_flags = " -DFEAST_HAVE_PARMETIS:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
