# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Jordi Paul, Peter Zajac"
__date__   = "April 2014"

import os
import subprocess
import importlib
import importlib.util
import inspect
import glob
import ssl
import sys

# Baseclass for third party packages
class ThirdpartyPackage(object):
  def __init__(self,trunk_dirname):
    # Package name as given in the build id.
    self.name = ["NullPackage"]
    # Name of folder where the package is expected to be found.
    self.dirname = "NullDirName"
    self.realdirname = self.dirname
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
    target_filename = os.path.join(self.trunk_dirname, self.filename)
    # This is the folder where the third party package is assumed to be. Not to be confused with self.target_dirname
    expected_dirname = os.path.join(self.trunk_dirname, self.dirname)
    if not os.path.isdir(expected_dirname):
      print (self.names[0] +" enabled, but could not find a directory with name " + expected_dirname + ", checking for file...")
      if not os.path.isfile(target_filename):
        print(target_filename + " not found, attempting to automatically download it from " + self.url + "...")
        download(self.url,target_filename)
        print("Download finished, so unpacking the sources now...")
      else:
        print(target_filename + " found, so unpacking the sources now...")
      self.unpack()
    return

  # unpacks the TAR/GZ/ZIP source archive
  def unpack(self):
    target_filename = os.path.join(self.trunk_dirname, self.filename)

    is_zip = False
    is_tar = False
    is_zip = is_zip or self.filename.endswith(".zip")
    is_tar = is_tar or self.filename.endswith(".tar.gz")
    is_tar = is_tar or self.filename.endswith(".tar.bz2")
    is_tar = is_tar or self.filename.endswith(".tar.xz")

    if(is_zip):
      import zipfile
      archive = zipfile.ZipFile(target_filename, "r")
      for f in archive.namelist():
        if f.startswith('..') or f.startswith(os.sep):
          print("Error: File "+self.filename+" contains absolute path or one starting with '..' ")
          print("This is an unsafe operation, aborting.")
          sys.exit(1)
    elif(is_tar):
      import tarfile
      archive = tarfile.open(target_filename, "r")
      for f in archive.getnames():
        if f.startswith('..') or f.startswith(os.sep):
          print("Error: File "+self.filename+" contains absolute path or one starting with '..' ")
          print("This is an unsafe operation, aborting.")
          sys.exit(1)
    else:
      print("Error: Unknown archive extension: '" + self.filename + "'")
      print("Don't know what to do, aborting.")
      sys.exit(1)

    archive.extractall(self.target_dirname)

    # maybe we have to patch the sources manually
    self.patch()

    return

    # patch_lst entry: [line-no, original, patched]
  def patch_file(self, filename, patch_lst):
    filename = os.path.join(self.realdirname, filename)
    # create backup file if it doesn't exist
    filename_b = filename + ".backup"
    if not os.path.isfile(filename_b):
      os.rename(filename, filename_b)
    print("Patching '%s'..." % filename)
    # open backup file for reading
    fi = open(filename_b, "rt")
    # open file for writing
    fo = open(filename, "wt")
    # loop over all input file lines
    lno = 0
    for line in fi:
      lno = lno + 1
      if (len(patch_lst) > 0) and (patch_lst[0][0] == lno):
        # this line is to be patched
        if line.strip() != patch_lst[0][1]:
          print("ERROR: when processing file '%s': in line %i" % (filename, lno))
          print("expected : '%s'" % patch_lst[0][1])
          print("but found: '%s'" % line.strip())
          print("Patch aborted!")
          sys.exit(1)
        # okay replace line
        fo.write(patch_lst[0][2] + "\n")
        # remove patch line
        patch_lst = patch_lst[1:]
      else:
        fo.write(line)
    # ensure that all patches were applied
    if len(patch_lst) > 0:
      print("ERROR: when processing file '%s': end of file found, but there are still patches left")
      print("Patch aborted!")
      sys.exit(1)
    # okay, that's it
    fo.close()
    fi.close()

  # patches the source code after unpacking;
  # this may be overriden by a derived class
  def patch(self):
    # nothing to patch by default
    pass

# Downloads a file from url and saves it to filename
def download(url,filename):
  import sys
  if sys.version_info[0] < 3:
    from urllib import FancyURLopener
    class MyOpener(FancyURLopener):
      version = "Wget/1.13.4 (linux-gnu)"
    urlretrieve = MyOpener().retrieve
    context = ssl.SSLContext(ssl.PROTOCOL_SSLv23)
    context.verify_mode=ssl.CERT_NONE
    ssl._create_default_https_context = ssl._create_unverified_context
    urlretrieve(url, filename)
  else:
    import urllib.request
    import shutil
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36'})
    with urllib.request.urlopen(req) as response:
      with open(filename, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
  return

# Find available third party packages by parsing all .py files in the build_system folder
def available_packages(files_path,target_path,name=''):
  available_files = glob.glob(files_path+os.sep+'thirdparty_*.py')
  available_packages = []
  for name in available_files:
    filename = os.path.basename(name)
    loader = importlib.machinery.SourceFileLoader(filename, name)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    members = inspect.getmembers(mod, inspect.isclass)
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
