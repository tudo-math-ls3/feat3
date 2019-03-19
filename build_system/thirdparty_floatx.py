# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Peter Zajac"
__date__   = "2018-03-19"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class FParser(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["floatx"]
    self.dirname = "floatx"
    self.filename = "FloatX-master.zip"
    self.url = "https://codeload.github.com/oprecomp/FloatX/zip/master"
    self.cmake_flags = " -DFEAT_HAVE_FLOATX:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
