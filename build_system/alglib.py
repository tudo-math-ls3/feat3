# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Jordi Paul"
__date__   = "April 2014"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ALGLIB(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "alglib"
    self.dirname = "alglib"
    self.filename = "alglib-3.10.0.cpp.gpl.zip"
    self.url = "http://www.alglib.net/translator/re/" + self.filename
    self.cmake_flags = " -DFEAST_HAVE_ALGLIB:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
