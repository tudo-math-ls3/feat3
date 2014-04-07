# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Jordi Paul"
__date__   = "April 2014"

import os
import imp
import inspect
import glob

# Baseclass for third party packages
class ThirdpartyPackage(object):
  def __init__(self,trunk_dirname):
    self.name = "NullPackage"
    self.dirname = "NullDirName"
    self.url = "NullUrl"
    self.filename = "NullFilename"
    self.trunk_dirname = trunk_dirname
    self.cmake_flags = "NullCmakeFlags"
    self.unpacker = "NullUnpacker"

# How to add a third party package to the build process
  def add(self):
    # Directory where the package is expected
    target_dirname = self.trunk_dirname+os.sep+self.dirname
    target_filename = self.trunk_dirname+os.sep+self.filename
    if not os.path.isdir(target_dirname):
      print (self.name +" enabled, but could not find a directory with name " + target_dirname+ ", checking for file...")
      if not os.path.isfile(target_filename):
	print(target_filename + " not found, attempting to automatically download it from " + self.url + "...")
        dlcmd = "wget " + self.url + " -O " + target_filename
        os.system(dlcmd)

      unpkcmd = self.unpacker + target_filename + self.unpack_target
      os.system(unpkcmd)
    return

def available_packages(path,name=''):
  available_files = glob.glob(path+os.sep+'*.py')
  available_packages = []
  for name in available_files:
    my_path, filename = os.path.split(name)
    (f,p,d) = imp.find_module(filename[:-3],[my_path])
    mod = imp.load_module(filename[:-3],f,p,d)
    members = inspect.getmembers(mod)
    classes = []
    for x in members:
      if type(x) == tuple:
        if len(x) == 2:
          if inspect.isclass(x[1]):
	    if issubclass(x[1],ThirdpartyPackage) and (x[1] != ThirdpartyPackage):
               classes.append(x[1])
    for x in classes:
      result = x(path)
      available_packages.append(result)
  return available_packages
