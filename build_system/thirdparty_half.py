# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Dirk Ribbrock"
__date__   = "Feb 2016"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class HALF(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["half"]
    self.dirname = "half"
    self.filename = "half-1.12.0.zip"
    self.url = "http://downloads.sourceforge.net/project/half/half/1.12.0/half-1.12.0.zip?r=http%3A%2F%2Fhalf.sourceforge.net%2Findex.html&ts=1455282582&use_mirror=netcologne"
    self.cmake_flags = " -DFEAT_HAVE_HALFMATH:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
