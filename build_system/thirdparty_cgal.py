# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Dirk Ribbrock"
__date__   = "Nov 2018"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class CGAL(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["cgal"]
    self.dirname = "cgal"
    self.filename = "CGAL-4.13.zip"
    self.url = "https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_CGAL:BOOL=ON"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname+os.sep+self.dirname
