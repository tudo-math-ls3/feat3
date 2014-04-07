from thirdparty.thirdparty_package import *
import os

class ALGLIB(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "alglib"
    self.dirname = "alglib"
    self.filename = "alglib-3.8.2.cpp.zip"
    self.url = "http://www.alglib.net/translator/re/" + self.filename
    self.cmake_flags = " -DFEAST_HAVE_ALGLIB:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.unpacker = "unzip -q "
    self.unpack_target = " -d " + self.trunk_dirname + os.sep + self.dirname
