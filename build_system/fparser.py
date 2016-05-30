# vim: set filetype=python sw=2 sts=2 et nofoldenable :
from build_system.thirdparty_package import ThirdpartyPackage
import os

class FParser(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.name = "fparser"
    self.dirname = "fparser"
    self.filename = "fparser4.5.2.zip"
    self.url = "http://warp.povusers.org/FunctionParser/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_FPARSER:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
