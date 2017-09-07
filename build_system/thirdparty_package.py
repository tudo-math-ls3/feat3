# vim: set filetype=python sw=2 sts=2 et nofoldenable :
__author__ = "Jordi Paul"
__date__   = "April 2014"

import os
import subprocess
import imp
import inspect
import glob

# Baseclass for third party packages
class ThirdpartyPackage(object):
  def __init__(self,trunk_dirname):
    # Package name as given in the build id.
    self.name = "NullPackage"
    # Name of folder where the package is expected to be found.
    self.dirname = "NullDirName"
    # Url from where to download the file, not including the filename.
    self.url = "NullUrl"
    # Filename
    self.filename = "NullFilename"
    # Absolute path of the FEAT2 directory. Crude, but works.
    self.trunk_dirname = trunk_dirname
    # Where to extract the archive. Needed because sometimes the archive's files are not in an appropriately named
    # subfolder
    self.target_dirname = trunk_dirname
    # This is added to the cmake_flags
    self.cmake_flags = "NullCmakeFlags"

# How to add a third party package to the build process
  def add(self):
    # Directory where the package is expected
    target_filename = self.trunk_dirname+os.sep+self.filename
    # This is the folder where the third party package is assumed to be. Not to be confused with self.target_dirname
    expected_dirname = self.trunk_dirname+os.sep+self.dirname
    if not os.path.isdir(expected_dirname):
      print (self.name +" enabled, but could not find a directory with name " + expected_dirname+ ", checking for file...")
      if not os.path.isfile(target_filename):
        print(target_filename + " not found, attempting to automatically download it from " + self.url + "...")
        download(self.url,target_filename)

      self.unpack()
    return

  def unpack(self):
    target_filename = self.trunk_dirname+os.sep+self.filename

    if(self.filename.endswith(".zip")):
      import zipfile
      archive = zipfile.ZipFile(target_filename, "r")
      for f in archive.namelist():
        if f.startswith('..') or f.startswith(os.sep):
          print("Error: File "+self.filename+" contains absolute path or one starting with '..' ")
          print("This is an unsafe operation, aborting.")
          sys.exit(1)
    elif(self.filename.endswith("tar.gz")):
      import tarfile
      archive = tarfile.open(target_filename, "r")
      for f in archive.getnames():
        if f.startswith('..') or f.startswith(os.sep):
          print("Error: File "+self.filename+" contains absolute path or one starting with '..' ")
          print("This is an unsafe operation, aborting.")
          sys.exit(1)

    archive.extractall(self.target_dirname)

    return

# Downloads a file from url and saves it to filename
def download(url,filename):
  import sys
  if sys.version_info[0] < 3:
    from urllib import FancyURLopener
    class MyOpener(FancyURLopener):
      version = "Wget/1.13.4 (linux-gnu)"
    urlretrieve = MyOpener().retrieve
    urlretrieve(url, filename)
  else:
    import urllib.request
    import shutil
    with urllib.request.urlopen(url) as response:
      with open(filename, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
  return

# Find available third party packages by parsing all .py files in the build_system folder
def available_packages(files_path,target_path,name=''):
  available_files = glob.glob(files_path+os.sep+'*.py')
  available_packages = []
  for name in available_files:
    my_path, filename = os.path.split(name)
    (f,p,d) = imp.find_module(filename[:-3],[my_path])
    mod = imp.load_module(filename[:-3],f,p,d)
    members = inspect.getmembers(mod)
    classes = []
    # For all files, find all the classes defined in them...
    for x in members:
      if type(x) == tuple:
        if len(x) == 2:
          if inspect.isclass(x[1]):
            # ... that are subclasses of ThirdpartyPackage but not ThirdpartyPackage itself
            if issubclass(x[1],ThirdpartyPackage) and (x[1] != ThirdpartyPackage):
               classes.append(x[1])
    for x in classes:
      result = x(target_path)
      available_packages.append(result)
  return available_packages
