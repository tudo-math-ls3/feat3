# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Dirk Ribbrock"
__date__   = "April 2017"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ZLIB(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "zlib"
    self.dirname = "zlib"
    self.filename = "zlib-1.2.11.tar.gz"
    self.url = "http://www.zlib.net/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_ZLIB:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
