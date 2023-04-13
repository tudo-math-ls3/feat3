#!/usr/bin/env python
########################################################################################################################
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
########################################################################################################################
# Third-Party Library build script for Visual Studio 2022 (VC17)
# ----------------------------------------------------------------------------------------------------------------------
#
# \author Peter Zajac
########################################################################################################################
import os
import sys
import re
import shutil
import glob
import subprocess
import tkinter.filedialog

sys.dont_write_bytecode = True

########################################################################################################################

# skip Visual Studio preparation?
skip_prepare = ("--skip-prepare" in sys.argv[1:])

# MS-MPI installed?
path_mpi = None
path_vc17 = None

# map of available packages
packages = {}

########################################################################################################################
########################################################################################################################
# Auxiliary Functions
########################################################################################################################
########################################################################################################################

def ask_yes_no_cancel():
  print("> ", end="", flush=True)
  for line in sys.stdin:
    if line.lower().startswith("y"):
      return True
    if line.lower().startswith("n"):
      return False
    if line.lower().startswith("q"):
      return None
    print("Please enter [y]es, [n]o or [q]uit\n> ", end="")

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

print("")
#      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
print("Welcome to the Visual Studio 2022 third-party library assistant script for FEAT 3")
print("")

# check whether we have
if os.environ.get("VS17PATH") == None:
  # check the default installation path
  def_path = os.path.join(os.environ.get("ProgramFiles"), "Microsoft Visual Studio", "2022", "Community")
  if os.path.isdir(def_path):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The 'VS17PATH' environment variable has not been set yet, however, it seems that")
    print("Visual Studio 2022 is installed in its default location:")
    print(def_path)
    print("\nDo you want to set the environment variable permanently now? ([y]es, [n]o, [q]uit)")
    rtn = ask_yes_no_cancel()
    if(rtn == None):
      print("User cancelled...")
      sys.exit(1)
    elif(rtn == True):
      print("Trying to set 'VS17PATH' environment variable permanently...")
      os.system("setx VS17PATH \"%s\"" % def_path)
      os.environ["VS17PATH"] = def_path
      path_vc17 = def_path
    else:
      print("Setting 'VS17PATH' environment variable for this session only")
      os.environ["VS17PATH"] = def_path
      path_vc17 = def_path
  else:
    # Visual Studio not installed in default location
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The 'VS17PATH' environment variable has not been set yet and it seems that")
    print("Visual Studio 2022 is not installed in its default location, so you will need")
    print("to specify the directory of your Visual Studio 2022 installation manually.\n")
    print("You need to locate the directory that contains (among others) the 'Common7' and")
    print("'VC' directories, which is usually named 'Community' or whatever edition you")
    print("have installed.\n")
    print("Do you want to set the environment variable permanently once you have selected")
    print("the installation directory?  ([y]es, [n]o, [q]uit)")
    rtn = ask_yes_no_cancel()
    if(rtn == None):
      print("User cancelled...")
      sys.exit(1)
    path_vc17 = tkinter.filedialog.askdirectory()
    if(len(path_vc17) == 0):
      print("User cancelled...")
      sys.exit(1)
    if(rtn == True):
      print("Trying to set 'VS17PATH' environment variable permanently...")
      os.system("setx VS17PATH \"%s\"" % path_vc17)
      os.environ["VS17PATH"] = path_vc17
    elif(rtn == False):
      print("Setting 'VS17PATH' environment variable for this session only")
      os.environ["VS17PATH"] = path_vc17
else:
  # 'VS17PATH' environment var is set
  path_vc17 = os.environ.get("VS17PATH")

# check whether we have MPI installed
if os.environ.get("MSMPI_INC") != None:
  path_mpi = os.path.abspath(os.path.join(os.environ.get("MSMPI_INC"), ".."))

# print current settings
print("The following paths have been found:")
print("Visual Studio 2022: " + (path_vc17 if path_vc17 != None else "- N/A -"))
print("Microsoft MPI.....: " + (path_mpi if path_mpi != None else "- N/A -"))

########################################################################################################################

# skip preparation?
if not skip_prepare:
  # try to load compiler
  print("\nPreparing Visual Studio 2022 build environment...")

  # call configuration
  try:
    config_call = "cmd.exe /C @call \"" + os.path.join(path_vc17, "VC", "Auxiliary", "Build", "vcvars64.bat") +"\" > NUL && set"
    vcvars_out = subprocess.check_output(config_call, universal_newlines=True)
  except:
    print("ERROR: Failed to configure Visual Studio 2022 for build!")
    sys.exit(1)

  # split and compile environment
  vcvars_list = vcvars_out.split("\n")
  vcenv = {}
  for v in vcvars_list:
    if len(v) < 1:
      continue
    k = v.split("=")
    if len(k) < 1:
      continue
    vcenv[k[0]] = k[1]

  # update our environment to include VC stuff
  os.environ.update(vcenv)

# create object directory
if not os.path.isdir(os.path.join("..", "lib")):
  os.makedirs(os.path.join("..", "lib"), exist_ok=True)
# create library directory
if not os.path.isdir(os.path.join("..", "lib")):
  os.makedirs(os.path.join("..", "lib"), exist_ok=True)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyPackage(object):
  def __init__(self):
    self.name = None
    self.version = None
    self.date = None
    self.dir = None
    self.file = None
    self.trunk = None
    self.url = None
    self.page = None
    self.prefix = ""
    self.suffix = ""
    self.need_mpi = False # MPI mandatory?
    # path to license file in extracted directory
    self.license_files = []
    # addition license info
    self.license_info = None
    # linker flags
    self.libflags = "/MACHINE:X64 /NOLOGO"
    # compiler flags that are always set
    self.baseflags = " /c /GS /Gd /W3 /WX- /Zc:wchar_t /Zi /Zc:forScope /errorReport:none /EHsc /nologo"
    # compiler flags to use for debug builds
    self.dbgflags = " /Od /RTC1 /MDd"
    # compiler flags to use for optimized builds
    self.optflags = " /MP /Gy /O2 /Ob2 /Oi /MD"
    # extended flags (set by derived class)
    self.extflags = ""
    # output file flags (set by prepare function)
    self.outflags = None
    self.objpath = None
    self.libpath = None

  def info(self):
    print("No info available for package '" + self.name + "'")

  def license(self):
    # print additional license info if available
    if self.license_info != None:
      print(self.license_info)
    # no license file available?
    if len(self.license_files) < 1:
      if self.license_info != None:
        print("No license info available for package '" + self.name + "'")
      return
    # download source archive if necessary
    if not os.path.isfile(self.file):
      print("Archive '" + self.file + "' not found, so downloading it now...")
      self.download()
    # unpack source archive if necessary
    if not os.path.isdir(self.dir):
      print("Directory '" + self.dir + "' not found, so unpacking archive now...")
      self.unpack()
    # loop over all license files
    for lic_file in self.license_files:
      # try to open license file
      lic_path = os.path.join(self.dir, lic_file)
      if not os.path.isfile(lic_path):
        print("ERROR: license file '" + lic_path + "' not found")
        return
      # print license to console
      with open(lic_path) as fin:
        print(fin.read())
      print("-"*100)

  def ready(self):
    return os.path.isfile(os.path.join("..", "lib", self.name + self.prefix + "-" + self.version + ".vc17-opt-x64.lib"))

  def remove(self):
    print("Removing package '%s'..." % self.name)
    for fn in glob.glob(os.path.join("..", "lib", self.name + "-" + self.version + ".vc17-*.lib")):
      os.remove(fn)
    for fn in glob.glob(os.path.join("..", "obj", self.name + "-" + self.version + ".vc17-*")):
      shutil.rmtree(fn)
    if os.path.isdir(self.dir):
      shutil.rmtree(self.dir)

  def install(self):
    if (path_vc17 == None):
      print("ERROR: Visual Studio 2022 installation not found, aborting...")
      sys.exit(1)
    if self.need_mpi and (path_mpi == None):
      print("WARNING: Cannot install " + self.name + " because MPI was not found, skipping...")
      return
    # download source archive if necessary
    if not os.path.isfile(self.file):
      print("Archive '" + self.file + "' not found, so downloading it now...")
      self.download()
    # unpack source archive if necessary
    if not os.path.isdir(self.dir):
      print("Directory '" + self.dir + "' not found, so unpacking archive now...")
      self.unpack()
    # patch sources if necessary
    self.patch()
    # build debug library
    self.suffix = "vc17-dbg-x64"
    self.modeflags = self.dbgflags
    self.build()
    # build opt library
    self.suffix = "vc17-opt-x64"
    self.modeflags = self.optflags
    self.build()

  # download the source archive from website
  def download(self):
    import urllib.request
    with urllib.request.urlopen(self.url) as response:
      with open(self.file, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    return

  # unpacks the TAR/GZ/ZIP source archive
  def unpack(self):
    is_zip = False
    is_tar = False
    is_zip = is_zip or self.file.endswith(".zip")
    is_tar = is_tar or self.file.endswith(".tar.gz")
    is_tar = is_tar or self.file.endswith(".tar.bz2")
    is_tar = is_tar or self.file.endswith(".tar.xz")

    if(is_zip):
      import zipfile
      archive = zipfile.ZipFile(self.file, "r")
      for f in archive.namelist():
        if f.startswith('..') or f.startswith(os.sep):
          print("ERROR: File '" + self.file + "' contains absolute path or one starting with '..' ")
          print("This is an unsafe operation, aborting.")
          sys.exit(1)
    elif(is_tar):
      import tarfile
      archive = tarfile.open(self.file, "r")
      for f in archive.getnames():
        if f.startswith('..') or f.startswith(os.sep):
          print("ERROR: File '" + self.file + "' contains absolute path or one starting with '..' ")
          print("This is an unsafe operation, aborting.")
          sys.exit(1)
    else:
      print("ERROR: Unknown archive extension: '" + self.file + "'")
      print("Don't know what to do, aborting.")
      sys.exit(1)

    archive.extractall(self.trunk)
    return

  # patch_lst entry: [line-no, original, patched]
  def patch_file(self, filename, patch_lst):
    filename = os.path.join(self.dir, filename)
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

  # patches the sources
  def patch(self):
    # nothing to patch by default
    pass

  # build the actual sources
  def prepare(self):
    # set output name
    self.outname = self.name + self.prefix + "-" + self.version + "." + self.suffix
    # set output paths
    self.objpath = os.path.join("..", "obj", self.outname)
    self.libpath = os.path.join("..", "lib", self.outname + ".lib")
    # set output file flags
    self.outflags  = ' /Fd"' + self.objpath + '/' + self.outname + '.pdb"'
    self.outflags += ' /Fp"' + self.objpath + '/' + self.outname + '.pch"'
    # set compiler flags
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    # create object directory
    if not os.path.isdir(self.objpath):
      os.makedirs(self.objpath, exist_ok=True)
    # check whether the output file already exists
    if os.path.isfile(self.libpath):
      print("\nLibrary '%s' already exists, so skipping build..." % self.libpath)
      return False
    # ok, let's get crackin'
    return True

  # compiles a single source file
  def compile(self, srcfile, objfile, extflags = ""):
    #print("cl.exe " + self.cxxflags + " " + extflags + " " + os.path.join(self.dir, srcfile) + ' /Fo"' + os.path.join(self.objpath, objfile) + '"')
    subprocess.run("cl.exe " + self.cxxflags + " " + extflags + " " + os.path.join(self.dir, srcfile) + ' /Fo"' + os.path.join(self.objpath, objfile) + '"')

  # compiles a single C source file with default object filename
  def compile_c(self, filename, extflags = ""):
    self.compile(filename + ".c", filename.replace("/","_") + ".obj", extflags + " /TC")

  # compiles a single CC source file with default object filename
  def compile_cc(self, filename, extflags = ""):
    self.compile(filename + ".cc", filename.replace("/","_") + ".obj", extflags + " /TP")

  # compiles a single CPP source file with default object filename
  def compile_cpp(self, filename, extflags = ""):
    self.compile(filename + ".cpp", filename.replace("/","_") + ".obj", extflags + " /TP")

  # links all object files to a library
  def link(self):
    print("\nLinking '" + self.outname + ".lib'...")
    subprocess.run("lib.exe " + self.libflags + ' /OUT:"' + self.libpath + '" ' + os.path.join(self.objpath,"*.obj"))


########################################################################################################################
########################################################################################################################
########################################################################################################################
# B O O S T  -  B O O S T  -  B O O S T  -  B O O S T  -  B O O S T  -  B O O S T  -  B O O S T  -  B O O S T  -  B O O
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyBoost(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "boost"
    self.version = "1.80.0"
    self.date = "2022-08-10"
    self.file = self.name + "_" + self.version.replace(".","_") + ".zip"
    self.dir = self.name + "_" + self.version.replace(".","_")
    self.trunk = "."
    self.license_files = ["LICENSE_1_0.txt"]
    self.url = "https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/" + self.file
    self.page = "https://www.boost.org"

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("boost is a large set of libraries providing various features. As of now, FEAT does not use the")
    print("boost libraries by itself, but some other third-party libraries (e.g. CGAL) may require boost")
    print("as a mandatory prerequisite.")
    print("")
    print("Prerequisites/Constraints: none")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_BOOST' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def ready(self):
    # check for version header
    return os.path.isfile(os.path.join(self.dir, "boost", "version.hpp"))

  def remove(self):
    if(os.path.isdir(self.dir)):
      shutil.rmtree(self.dir)

  def build(self):
    # boost is a header-only library
    pass

packages["boost"] = ThirdPartyBoost()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# C G A L  -  C G A L  -  C G A L  -  C G A L  -  C G A L  -  C G A L  -  C G A L  -  C G A L  -  C G A L  -  C G A L  -
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyCGAL(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "CGAL"
    self.version = "5.5.1"
    self.date = "2022-10-12"
    self.file = self.name + "-" + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.license_files = ["LICENSE"]
    self.url = "https://github.com/CGAL/cgal/releases/download/v" + self.version + "/" + self.file
    self.page = "https://www.cgal.org/"

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The Computational Geometry Algorithms Library (CGAL) is a library that offers a large variety of")
    print("of geometry-related algorithms. Installing this package enables the use of the Geometry::CGALWrapper")
    print("class that can be used for inside/outside tests of geometries represented by a surface triangulation.")
    print("")
    print("Prerequisites/Constraints: This library requires the 'boost' third-party library.")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_CGAL' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def ready(self):
    # check for versiion header
    return os.path.isfile(os.path.join(self.dir, "include", "CGAL", "version.h"))

  def remove(self):
    if(os.path.isdir(self.dir)):
      shutil.rmtree(self.dir)

  def build(self):
    # CGAL is a header-only library
    pass

packages["cgal"] = ThirdPartyCGAL()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# F P A R S E R  -  F P A R S E R  -  F P A R S E R  -  F P A R S E R  -  F P A R S E R  -  F P A R S E R  -  F P A R S
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyFParser(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "fparser"
    self.version = "4.5.2"
    self.date = "2015-06-07"
    self.file = self.name + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = self.dir
    self.license_files = [os.path.join("docs", "lgpl.txt")]
    self.url = "http://warp.povusers.org/FunctionParser/" + self.file
    self.page = "http://warp.povusers.org/FunctionParser"
    self.baseflags += ' /wd"4068" /wd"4244" /wd"4838"'
    self.baseflags += ' /I"./' + self.dir + '"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The FParser library offers a light-weight function parser that can be used to evaluate floating")
    print("point expressions at runtime. Installing this package enables the use of the")
    print("Analytic::ParsedScalarFunction and Analytic::ParsedVectorFunction classes, which are very commonly")
    print("used by many tutorials (e.g. Tutorial 04) and many applications.")
    print("")
    print("Prerequisites/Constraints: none")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_FPARSER' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " sources...")
    self.compile_cc("fparser")
    self.compile_cc("fpoptimizer")
    self.link()

packages["fparser"] = ThirdPartyFParser()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# H Y P R E  -  H Y P R E  -  H Y P R E  -  H Y P R E  -  H Y P R E  -  H Y P R E  -  H Y P R E  -  H Y P R E  -  H Y P
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyHypre(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "hypre"
    self.version = "2.26.0"
    self.date = "2022-10-14"
    self.file = self.name + "-" + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.license_files = ["LICENSE-MIT", "LICENSE-APACHE"]
    self.license_info = \
      "Note: HYPRE is available under both the MIT and the Apache 2.0 licenses, which are both shown here.\n" +\
      "You can pick whichever of these two licenses suits you best.\n\n"
    self.url = "https://github.com/hypre-space/hypre/archive/refs/tags/v" + self.version + ".zip"
    self.page = "http://www.llnl.gov/casc/hypre"
    self.baseflags += ' /wd"4028" /wd"4244" /wd"4293" /wd"4305" /wd"4334" /wd"4700" /wd"4996"'
    self.baseflags += ' /D "HAVE_CONFIG_H" /D "WIN32"'
    self.dbgflags += ' /D "HYPRE_DEBUG"'
    self.baseflags += ' /I"./' + self.dir + '/src"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The HYPRE library is a large collection of high performance serial/parallel preconditioners, which")
    print("includes (amongst others) the MPI-parallel algebraic multigrid preconditioner named 'BoomerAMG'.")
    print("Installing this package enables the use of the HYPRE-wrapper solver classes, which currently include")
    print("the Solver::BoomerAMG, Solver::EuclidPrecond and Solver::ParaSailsPrecond classes.")
    print("FEAT includes backends for both the sequential as well as the MPI-parallel builds of HYPRE and the")
    print("appropriate backend is chosen depending on whether FEAT is build with MPI enabled or disabled.")
    print("")
    print("Prerequisites/Constraints: parallel version requires MPI to be installed (obviously)")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_HYPRE' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def patch(self):
    print("\nPatching " + self.name + " " + self.version + " sources...")
    # write configure file
    cfg_file = os.path.join(self.dir, "src", "HYPRE_config.h")
    if not os.path.isfile(cfg_file):
      print("Writing Config Header '%s'..." % cfg_file)
      fo = open(cfg_file, "wt")
      fo.write("/* HYPRE config header generated by FEAT3 'vc17-thirdparty.py' script */\n")
      fo.write("/* Note: HYPRE_HAVE_MPI and HYPRE_SEQUENTIAL are defined by the */\n")
      fo.write("/*       build system and are therefore not defined here. */\n");
      fo.write("#define HAVE_INTTYPES_H\n")
      fo.write("#define HAVE_MEMORY_H\n")
      fo.write("#define HAVE_STDINT_H\n")
      fo.write("#define HAVE_STDLIB_H\n")
      fo.write("#define HAVE_STRING_H\n")
      fo.write("#define HYPRE_BIGINT 1\n")
      fo.write("#define HYPRE_MAXDIM 3\n")
      fo.write("#define HYPRE_RELEASE_BUGS \"https://github.com/hypre-space/hypre/issues\"\n")
      fo.write("#define HYPRE_RELEASE_DATE \"2022/10/14\"\n")
      fo.write("#define HYPRE_RELEASE_NAME \"hypre\"\n")
      fo.write("#define HYPRE_RELEASE_NUMBER 22600\n")
      fo.write("#define HYPRE_RELEASE_TIME \"00:00:00\"\n")
      fo.write("#define HYPRE_RELEASE_VERSION \"2.26.0\"\n")
      fo.write("#define HYPRE_USING_HOST_MEMORY 1\n")
      fo.write("#define HYPRE_USING_HYPRE_BLAS 1\n")
      fo.write("#define HYPRE_USING_HYPRE_LAPACK 1\n")
      fo.write("#define PACKAGE_BUGREPORT \"\"\n")
      fo.write("#define PACKAGE_NAME \"hypre\"\n")
      fo.write("#define PACKAGE_STRING \"hypre 2.26.0\"\n")
      fo.write("#define PACKAGE_TARNAME \"hypre\"\n")
      fo.write("#define PACKAGE_URL \"\"\n")
      fo.write("#define PACKAGE_VERSION \"2.26.0\"\n")
      fo.close()
    # patch out non-existing header includes
    self.patch_file(os.path.join("src", "matrix_matrix", "HYPRE_ConvertParCSRMatrixToDistributedMatrix.c"), [
      [15, '#include <gmalloc.h>', '/*#include <gmalloc.h>*/']
    ])
    self.patch_file(os.path.join("src", "matrix_matrix", "HYPRE_ConvertPETScMatrixToDistributedMatrix.c"), [
      [15, '#include <gmalloc.h>', '/*#include <gmalloc.h>*/']
    ])
    # patch 32 vs 64 bit int issues in utilities
    self.patch_file(os.path.join("src", "utilities", "_hypre_utilities.h"), [
      [2171, 'hypre_ulongint h64 = HYPRE_XXH_PRIME64_5 + sizeof(input);', '  hypre_ulonglongint h64 = HYPRE_XXH_PRIME64_5 + sizeof(input);'],
      [2173, 'hypre_ulongint k1 = input;', '  hypre_ulonglongint k1 = input;']
    ])

  def ready(self):
    return os.path.isfile(os.path.join("..", "lib", "hypre-seq-" + self.version + ".vc17-opt-x64.lib"))

  def remove(self):
    for fn in glob.glob(os.path.join("..", "lib","hypre-*-" + self.version + ".vc17-*-x64.lib")):
      os.remove(fn)
    for fn in glob.glob(os.path.join("..", "obj", "hypre-*-" + self.version + ".vc17-*-x64")):
      shutil.rmtree(fn)
    if os.path.isdir(self.dir):
      shutil.rmtree(self.dir)

  def build(self):
    self.build_serial()
    if(path_mpi != None):
      self.build_parallel()
    else:
      print("WARNING: Can't install MPI-parallel version of HYPRE because MPI was not found")

  def build_serial(self):
    # build serial version
    print("Building sequential HYPRE library...")
    self.prefix = "-seq"
    self.extflags = ' /DHYPRE_SEQUENTIAL'
    self.build_sources()

  def build_parallel(self):
    print("Building MPI-parallel HYPRE library...")
    self.prefix = "-mpi"
    self.extflags = ' /DHYPRE_HAVE_MPI'
    self.extflags += ' /I"' + os.path.join(path_mpi, "Include") + '"'
    self.build_sources()

  def build_sources(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " blas sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.compile_c("src/blas/dasum")
    self.compile_c("src/blas/daxpy")
    self.compile_c("src/blas/dcopy")
    self.compile_c("src/blas/ddot")
    self.compile_c("src/blas/dgemm")
    self.compile_c("src/blas/dgemv")
    self.compile_c("src/blas/dger")
    self.compile_c("src/blas/dnrm2")
    self.compile_c("src/blas/drot")
    self.compile_c("src/blas/dscal")
    self.compile_c("src/blas/dswap")
    self.compile_c("src/blas/dsymm")
    self.compile_c("src/blas/dsymv")
    self.compile_c("src/blas/dsyr2")
    self.compile_c("src/blas/dsyr2k")
    self.compile_c("src/blas/dsyrk")
    self.compile_c("src/blas/dtrmm")
    self.compile_c("src/blas/dtrmv")
    self.compile_c("src/blas/dtrsm")
    self.compile_c("src/blas/dtrsv")
    self.compile_c("src/blas/f2c")
    self.compile_c("src/blas/idamax")
    self.compile_c("src/blas/lsame")
    self.compile_c("src/blas/xerbla")

    print("\nBuilding " + self.name + " " + self.version + " lapack sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/blas"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.compile_c("src/lapack/dbdsqr")
    self.compile_c("src/lapack/dgebd2")
    self.compile_c("src/lapack/dgebrd")
    self.compile_c("src/lapack/dgelq2")
    self.compile_c("src/lapack/dgelqf")
    self.compile_c("src/lapack/dgels")
    self.compile_c("src/lapack/dgeqr2")
    self.compile_c("src/lapack/dgeqrf")
    self.compile_c("src/lapack/dgesvd")
    self.compile_c("src/lapack/dgetrf")
    self.compile_c("src/lapack/dgetri")
    self.compile_c("src/lapack/dgetrs")
    self.compile_c("src/lapack/dgetf2")
    self.compile_c("src/lapack/dlabad")
    self.compile_c("src/lapack/dlabrd")
    self.compile_c("src/lapack/dlacpy")
    self.compile_c("src/lapack/dlae2")
    self.compile_c("src/lapack/dlaev2")
    self.compile_c("src/lapack/dlange")
    self.compile_c("src/lapack/dlanst")
    self.compile_c("src/lapack/dlansy")
    self.compile_c("src/lapack/dlapy2")
    self.compile_c("src/lapack/dlarfb")
    self.compile_c("src/lapack/dlarf")
    self.compile_c("src/lapack/dlarfg")
    self.compile_c("src/lapack/dlarft")
    self.compile_c("src/lapack/dlartg")
    self.compile_c("src/lapack/dlas2")
    self.compile_c("src/lapack/dlascl")
    self.compile_c("src/lapack/dlaset")
    self.compile_c("src/lapack/dlasq1")
    self.compile_c("src/lapack/dlasq2")
    self.compile_c("src/lapack/dlasq3")
    self.compile_c("src/lapack/dlasq4")
    self.compile_c("src/lapack/dlasq5")
    self.compile_c("src/lapack/dlasq6")
    self.compile_c("src/lapack/dlasr")
    self.compile_c("src/lapack/dlasrt")
    self.compile_c("src/lapack/dlassq")
    self.compile_c("src/lapack/dlaswp")
    self.compile_c("src/lapack/dlasv2")
    self.compile_c("src/lapack/dlatrd")
    self.compile_c("src/lapack/dorg2l")
    self.compile_c("src/lapack/dorg2r")
    self.compile_c("src/lapack/dorgbr")
    self.compile_c("src/lapack/dorgl2")
    self.compile_c("src/lapack/dorglq")
    self.compile_c("src/lapack/dorgql")
    self.compile_c("src/lapack/dorgqr")
    self.compile_c("src/lapack/dorgtr")
    self.compile_c("src/lapack/dorm2r")
    self.compile_c("src/lapack/dormbr")
    self.compile_c("src/lapack/dorml2")
    self.compile_c("src/lapack/dormlq")
    self.compile_c("src/lapack/dormqr")
    self.compile_c("src/lapack/dpotf2")
    self.compile_c("src/lapack/dpotrf")
    self.compile_c("src/lapack/dpotrs")
    self.compile_c("src/lapack/dsteqr")
    self.compile_c("src/lapack/dsterf")
    self.compile_c("src/lapack/dsyev")
    self.compile_c("src/lapack/dsygs2")
    self.compile_c("src/lapack/dsygst")
    self.compile_c("src/lapack/dsygv")
    self.compile_c("src/lapack/dsytd2")
    self.compile_c("src/lapack/dsytrd")
    self.compile_c("src/lapack/dtrtri")
    self.compile_c("src/lapack/dtrti2")
    self.compile_c("src/lapack/ieeeck")
    self.compile_c("src/lapack/ilaenv")
    self.compile_c("src/lapack/lsame")
    self.compile_c("src/lapack/xerbla")
    self.compile_c("src/lapack/dlamch")

    print("\nBuilding " + self.name + " " + self.version + " struct_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_mv"'
    self.compile_c("src/utilities/HYPRE_handle")
    self.compile_c("src/utilities/HYPRE_version")
    self.compile_c("src/utilities/amg_linklist")
    self.compile_c("src/utilities/binsearch")
    self.compile_c("src/utilities/exchange_data")
    self.compile_c("src/utilities/fortran_matrix")
    self.compile_c("src/utilities/ap")
    self.compile_c("src/utilities/log")
    self.compile_c("src/utilities/complex")
    self.compile_c("src/utilities/error")
    self.compile_c("src/utilities/hopscotch_hash")
    self.compile_c("src/utilities/memory_tracker")
    self.compile_c("src/utilities/merge_sort")
    self.compile_c("src/utilities/mmio")
    self.compile_c("src/utilities/mpi_comm_f2c")
    self.compile_c("src/utilities/prefix_sum")
    self.compile_c("src/utilities/printf")
    self.compile_c("src/utilities/qsort")
    self.compile_c("src/utilities/utilities")
    self.compile_c("src/utilities/mpistubs")
    self.compile_c("src/utilities/qsplit")
    self.compile_c("src/utilities/random")
    self.compile_c("src/utilities/threading")
    self.compile_c("src/utilities/timer")
    self.compile_c("src/utilities/timing")
    self.compile_c("src/utilities/device_utils")
    self.compile_c("src/utilities/general")
    self.compile_c("src/utilities/handle")
    self.compile_c("src/utilities/int_array")
    self.compile_c("src/utilities/memory")
    self.compile_c("src/utilities/omp_device")
    self.compile_c("src/utilities/nvtx")

    print("\nBuilding " + self.name + " " + self.version + " multivector sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/multivector"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.compile_c("src/multivector/multivector")
    self.compile_c("src/multivector/temp_multivector")

    print("\nBuilding " + self.name + " " + self.version + " krylov sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/krylov"'
    self.cxxflags += ' /I"./' + self.dir + '/src/blas"'
    self.cxxflags += ' /I"./' + self.dir + '/src/lapack"'
    self.cxxflags += ' /I"./' + self.dir + '/src/multivector"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.compile_c("src/krylov/bicgstab")
    self.compile_c("src/krylov/cgnr")
    self.compile_c("src/krylov/gmres")
    self.compile_c("src/krylov/cogmres")
    self.compile_c("src/krylov/flexgmres")
    self.compile_c("src/krylov/lgmres")
    self.compile_c("src/krylov/HYPRE_bicgstab")
    self.compile_c("src/krylov/HYPRE_cgnr")
    self.compile_c("src/krylov/HYPRE_gmres")
    self.compile_c("src/krylov/HYPRE_cogmres")
    self.compile_c("src/krylov/HYPRE_lgmres")
    self.compile_c("src/krylov/HYPRE_flexgmres")
    self.compile_c("src/krylov/HYPRE_pcg")
    self.compile_c("src/krylov/pcg")
    self.compile_c("src/krylov/HYPRE_lobpcg")
    self.compile_c("src/krylov/lobpcg")

    print("\nBuilding " + self.name + " " + self.version + " seq_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.compile_c("src/seq_mv/csr_matop")
    self.compile_c("src/seq_mv/csr_matrix")
    self.compile_c("src/seq_mv/csr_matvec")
    self.compile_c("src/seq_mv/genpart")
    self.compile_c("src/seq_mv/HYPRE_csr_matrix")
    self.compile_c("src/seq_mv/HYPRE_mapped_matrix")
    self.compile_c("src/seq_mv/HYPRE_multiblock_matrix")
    self.compile_c("src/seq_mv/HYPRE_vector")
    self.compile_c("src/seq_mv/mapped_matrix")
    self.compile_c("src/seq_mv/multiblock_matrix")
    self.compile_c("src/seq_mv/vector")
    self.compile_c("src/seq_mv/vector_batched")
    #self.compile_c("src/seq_mv/csr_matop_device")
    #self.compile_c("src/seq_mv/csr_matrix_cuda_utils")
    #self.compile_c("src/seq_mv/csr_matvec_device")
    self.compile_c("src/seq_mv/csr_matvec_oomp")
    #self.compile_c("src/seq_mv/csr_spadd_device")
    #self.compile_c("src/seq_mv/csr_spgemm_device")
    #self.compile_c("src/seq_mv/csr_spgemm_device_cusparse")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numblocks")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer1")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer2")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer3")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer4")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer5")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer6")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer7")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer8")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer9")
    #self.compile_c("src/seq_mv/csr_spgemm_device_numer10")
    #self.compile_c("src/seq_mv/csr_spgemm_device_onemklsparse")
    #self.compile_c("src/seq_mv/csr_spgemm_device_rocsparse")
    #self.compile_c("src/seq_mv/csr_spgemm_device_rowest")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl1")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl2")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl3")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl4")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl5")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl6")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl7")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl8")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl9")
    #self.compile_c("src/seq_mv/csr_spgemm_device_symbl10")
    #self.compile_c("src/seq_mv/csr_spgemm_device_util")
    #self.compile_c("src/seq_mv/csr_spmv_device")
    #self.compile_c("src/seq_mv/csr_sptrans_device")
    #self.compile_c("src/seq_mv/vector_device")

    print("\nBuilding " + self.name + " " + self.version + " parcsr_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/blas"'
    self.cxxflags += ' /I"./' + self.dir + '/src/lapack"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.compile_c("src/parcsr_mv/communicationT")
    self.compile_c("src/parcsr_mv/HYPRE_parcsr_matrix")
    self.compile_c("src/parcsr_mv/HYPRE_parcsr_vector")
    self.compile_c("src/parcsr_mv/gen_fffc")
    self.compile_c("src/parcsr_mv/new_commpkg")
    self.compile_c("src/parcsr_mv/numbers")
    self.compile_c("src/parcsr_mv/par_csr_aat")
    self.compile_c("src/parcsr_mv/par_csr_assumed_part")
    self.compile_c("src/parcsr_mv/par_csr_bool_matop")
    self.compile_c("src/parcsr_mv/par_csr_bool_matrix")
    self.compile_c("src/parcsr_mv/par_csr_communication")
    self.compile_c("src/parcsr_mv/par_csr_matop")
    self.compile_c("src/parcsr_mv/par_csr_matrix")
    self.compile_c("src/parcsr_mv/par_csr_matvec")
    self.compile_c("src/parcsr_mv/par_csr_matop_marked")
    self.compile_c("src/parcsr_mv/par_csr_triplemat")
    self.compile_c("src/parcsr_mv/par_make_system")
    self.compile_c("src/parcsr_mv/par_vector")
    self.compile_c("src/parcsr_mv/par_vector_batched")
    #self.compile_c("src/parcsr_mv/par_csr_matvec_device")
    #self.compile_c("src/parcsr_mv/par_csr_fffc_device")
    #self.compile_c("src/parcsr_mv/par_csr_matop_device")
    #self.compile_c("src/parcsr_mv/par_csr_triplemat_device")
    #self.compile_c("src/parcsr_mv/par_vector_device")

    print("\nBuilding " + self.name + " " + self.version + " parcsr_block_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_block_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/multivector"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_ls"'
    self.cxxflags += ' /I"./' + self.dir + '/src/IJ_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/krylov"'
    self.compile_c("src/parcsr_block_mv/csr_block_matrix")
    self.compile_c("src/parcsr_block_mv/csr_block_matvec")
    self.compile_c("src/parcsr_block_mv/par_csr_block_matrix")
    self.compile_c("src/parcsr_block_mv/par_csr_block_matvec")
    self.compile_c("src/parcsr_block_mv/par_csr_block_comm")
    self.compile_c("src/parcsr_block_mv/par_csr_block_rap")
    self.compile_c("src/parcsr_block_mv/par_csr_block_rap_communication")
    self.compile_c("src/parcsr_block_mv/par_csr_block_interp")
    self.compile_c("src/parcsr_block_mv/par_csr_block_relax")
    self.compile_c("src/parcsr_block_mv/par_block_nodal_systems")

    print("\nBuilding " + self.name + " " + self.version + " distributed_matrix sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_matrix"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.compile_c("src/distributed_matrix/distributed_matrix")
    self.compile_c("src/distributed_matrix/HYPRE_distributed_matrix")
    self.compile_c("src/distributed_matrix/distributed_matrix_ISIS")
    self.compile_c("src/distributed_matrix/distributed_matrix_PETSc")
    self.compile_c("src/distributed_matrix/distributed_matrix_parcsr")

    print("\nBuilding " + self.name + " " + self.version + " IJ_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/IJ_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_ls"'
    self.compile_c("src/IJ_mv/aux_parcsr_matrix")
    self.compile_c("src/IJ_mv/aux_par_vector")
    self.compile_c("src/IJ_mv/HYPRE_IJMatrix")
    self.compile_c("src/IJ_mv/HYPRE_IJVector")
    self.compile_c("src/IJ_mv/IJ_assumed_part")
    self.compile_c("src/IJ_mv/IJMatrix")
    self.compile_c("src/IJ_mv/IJMatrix_parcsr")
    self.compile_c("src/IJ_mv/IJVector")
    self.compile_c("src/IJ_mv/IJVector_parcsr")
    #self.compile_c("src/IJ_mv/IJMatrix_parcsr_device")
    #self.compile_c("src/IJ_mv/IJVector_parcsr_device")

    print("\nBuilding " + self.name + " " + self.version + " matrix_matrix sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/matrix_matrix"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_matrix"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.compile_c("src/matrix_matrix/HYPRE_ConvertParCSRMatrixToDistributedMatrix")
    self.compile_c("src/matrix_matrix/HYPRE_ConvertPETScMatrixToDistributedMatrix")

    print("\nBuilding " + self.name + " " + self.version + " parcsr_ls sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_ls"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_ls/Euclid"'
    self.cxxflags += ' /I"./' + self.dir + '/src/blas"'
    self.cxxflags += ' /I"./' + self.dir + '/src/lapack"'
    self.cxxflags += ' /I"./' + self.dir + '/src/multivector"'
    self.cxxflags += ' /I"./' + self.dir + '/src/krylov"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_matrix"'
    self.cxxflags += ' /I"./' + self.dir + '/src/matrix_matrix"'
    self.cxxflags += ' /I"./' + self.dir + '/src/IJ_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_block_mv"'
    self.compile_c("src/parcsr_ls/amg_hybrid")
    self.compile_c("src/parcsr_ls/aux_interp")
    self.compile_c("src/parcsr_ls/gen_redcs_mat")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_amg")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_amgdd")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_bicgstab")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_block")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_cgnr")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_Euclid")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_gmres")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_cogmres")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_flexgmres")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_lgmres")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_hybrid")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_int")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_mgr")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_ilu")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_fsai")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_ParaSails")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_pcg")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_pilut")
    self.compile_c("src/parcsr_ls/HYPRE_parcsr_schwarz")
    self.compile_c("src/parcsr_ls/HYPRE_ams")
    self.compile_c("src/parcsr_ls/HYPRE_ads")
    self.compile_c("src/parcsr_ls/HYPRE_ame")
    self.compile_c("src/parcsr_ls/par_2s_interp")
    self.compile_c("src/parcsr_ls/par_amg")
    self.compile_c("src/parcsr_ls/par_amgdd")
    self.compile_c("src/parcsr_ls/par_amgdd_comp_grid")
    self.compile_c("src/parcsr_ls/par_amgdd_setup")
    self.compile_c("src/parcsr_ls/par_amgdd_solve")
    self.compile_c("src/parcsr_ls/par_amgdd_fac_cycle")
    self.compile_c("src/parcsr_ls/par_amgdd_helpers")
    self.compile_c("src/parcsr_ls/par_amg_solve")
    self.compile_c("src/parcsr_ls/par_amg_solveT")
    self.compile_c("src/parcsr_ls/par_fsai")
    self.compile_c("src/parcsr_ls/par_fsai_setup")
    self.compile_c("src/parcsr_ls/par_fsai_solve")
    self.compile_c("src/parcsr_ls/par_cg_relax_wt")
    self.compile_c("src/parcsr_ls/par_coarsen")
    self.compile_c("src/parcsr_ls/par_cgc_coarsen")
    self.compile_c("src/parcsr_ls/par_cheby")
    self.compile_c("src/parcsr_ls/par_coarse_parms")
    self.compile_c("src/parcsr_ls/par_coordinates")
    self.compile_c("src/parcsr_ls/par_cr")
    self.compile_c("src/parcsr_ls/par_cycle")
    self.compile_c("src/parcsr_ls/par_add_cycle")
    self.compile_c("src/parcsr_ls/par_difconv")
    self.compile_c("src/parcsr_ls/par_gauss_elim")
    self.compile_c("src/parcsr_ls/par_gsmg")
    self.compile_c("src/parcsr_ls/par_indepset")
    self.compile_c("src/parcsr_ls/par_interp")
    self.compile_c("src/parcsr_ls/par_jacobi_interp")
    self.compile_c("src/parcsr_ls/par_krylov_func")
    self.compile_c("src/parcsr_ls/par_mod_lr_interp")
    self.compile_c("src/parcsr_ls/par_multi_interp")
    self.compile_c("src/parcsr_ls/par_mod_multi_interp")
    self.compile_c("src/parcsr_ls/par_laplace")
    self.compile_c("src/parcsr_ls/par_laplace_27pt")
    self.compile_c("src/parcsr_ls/par_laplace_9pt")
    self.compile_c("src/parcsr_ls/par_lr_interp")
    self.compile_c("src/parcsr_ls/par_mgr")
    self.compile_c("src/parcsr_ls/par_mgr_setup")
    self.compile_c("src/parcsr_ls/par_mgr_solve")
    self.compile_c("src/parcsr_ls/par_nongalerkin")
    self.compile_c("src/parcsr_ls/par_nodal_systems")
    self.compile_c("src/parcsr_ls/par_rap")
    self.compile_c("src/parcsr_ls/par_rap_communication")
    self.compile_c("src/parcsr_ls/par_rotate_7pt")
    self.compile_c("src/parcsr_ls/par_relax")
    self.compile_c("src/parcsr_ls/par_relax_more")
    self.compile_c("src/parcsr_ls/par_relax_interface")
    self.compile_c("src/parcsr_ls/par_scaled_matnorm")
    self.compile_c("src/parcsr_ls/par_schwarz")
    self.compile_c("src/parcsr_ls/par_stats")
    self.compile_c("src/parcsr_ls/par_strength")
    self.compile_c("src/parcsr_ls/par_sv_interp")
    self.compile_c("src/parcsr_ls/par_sv_interp_ln")
    self.compile_c("src/parcsr_ls/par_vardifconv")
    self.compile_c("src/parcsr_ls/par_vardifconv_rs")
    self.compile_c("src/parcsr_ls/partial")
    self.compile_c("src/parcsr_ls/schwarz")
    self.compile_c("src/parcsr_ls/block_tridiag")
    self.compile_c("src/parcsr_ls/par_restr")
    self.compile_c("src/parcsr_ls/par_lr_restr")
    self.compile_c("src/parcsr_ls/dsuperlu")
    self.compile_c("src/parcsr_ls/ads")
    self.compile_c("src/parcsr_ls/ams")
    self.compile_c("src/parcsr_ls/ame")
    self.compile_c("src/parcsr_ls/par_amg_setup")
    self.compile_c("src/parcsr_ls/par_ilu")
    self.compile_c("src/parcsr_ls/par_ilu_setup")
    self.compile_c("src/parcsr_ls/par_ilu_solve")
    #self.compile_c("src/parcsr_ls/par_cheby_device")
    #self.compile_c("src/parcsr_ls/par_relax_more_device")
    #self.compile_c("src/parcsr_ls/par_coarsen_device")
    #self.compile_c("src/parcsr_ls/par_coarse_parms_device")
    #self.compile_c("src/parcsr_ls/par_indepset_device")
    #self.compile_c("src/parcsr_ls/par_interp_device")
    #self.compile_c("src/parcsr_ls/par_lr_restr_device")
    #self.compile_c("src/parcsr_ls/par_interp_trunc_device")
    #self.compile_c("src/parcsr_ls/par_lr_interp_device")
    #self.compile_c("src/parcsr_ls/par_strength_device")
    #self.compile_c("src/parcsr_ls/par_strength2nd_device")
    #self.compile_c("src/parcsr_ls/par_amgdd_fac_cycle_device")
    #self.compile_c("src/parcsr_ls/par_2s_interp_device")
    #self.compile_c("src/parcsr_ls/par_relax_device")
    #self.compile_c("src/parcsr_ls/par_mod_multi_interp_device")
    #self.compile_c("src/parcsr_ls/par_mgr_device")

    print("\nBuilding " + self.name + " " + self.version + " struct_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.compile_c("src/struct_mv/assumed_part")
    self.compile_c("src/struct_mv/box_algebra")
    self.compile_c("src/struct_mv/box_boundary")
    self.compile_c("src/struct_mv/box")
    self.compile_c("src/struct_mv/box_manager")
    self.compile_c("src/struct_mv/communication_info")
    self.compile_c("src/struct_mv/computation")
    self.compile_c("src/struct_mv/HYPRE_struct_grid")
    self.compile_c("src/struct_mv/HYPRE_struct_matrix")
    self.compile_c("src/struct_mv/HYPRE_struct_stencil")
    self.compile_c("src/struct_mv/HYPRE_struct_vector")
    self.compile_c("src/struct_mv/project")
    self.compile_c("src/struct_mv/struct_grid")
    self.compile_c("src/struct_mv/struct_io")
    self.compile_c("src/struct_mv/struct_matrix_mask")
    self.compile_c("src/struct_mv/struct_stencil")
    self.compile_c("src/struct_mv/struct_axpy")
    self.compile_c("src/struct_mv/struct_communication")
    self.compile_c("src/struct_mv/struct_copy")
    self.compile_c("src/struct_mv/struct_innerprod")
    self.compile_c("src/struct_mv/struct_matrix")
    self.compile_c("src/struct_mv/struct_matvec")
    self.compile_c("src/struct_mv/struct_scale")
    self.compile_c("src/struct_mv/struct_vector")

    print("\nBuilding " + self.name + " " + self.version + " struct_ls sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_ls"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/multivector"'
    self.cxxflags += ' /I"./' + self.dir + '/src/krylov"'
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_mv"'
    self.compile_c("src/struct_ls/coarsen")
    self.compile_c("src/struct_ls/hybrid")
    self.compile_c("src/struct_ls/HYPRE_struct_bicgstab")
    self.compile_c("src/struct_ls/HYPRE_struct_cycred")
    self.compile_c("src/struct_ls/HYPRE_struct_flexgmres")
    self.compile_c("src/struct_ls/HYPRE_struct_gmres")
    self.compile_c("src/struct_ls/HYPRE_struct_hybrid")
    self.compile_c("src/struct_ls/HYPRE_struct_jacobi")
    self.compile_c("src/struct_ls/HYPRE_struct_lgmres")
    self.compile_c("src/struct_ls/HYPRE_struct_pfmg")
    self.compile_c("src/struct_ls/HYPRE_struct_smg")
    self.compile_c("src/struct_ls/HYPRE_struct_sparse_msg")
    self.compile_c("src/struct_ls/jacobi")
    self.compile_c("src/struct_ls/pcg_struct")
    self.compile_c("src/struct_ls/pfmg")
    self.compile_c("src/struct_ls/pfmg_relax")
    self.compile_c("src/struct_ls/pfmg_setup_rap")
    self.compile_c("src/struct_ls/pfmg_solve")
    self.compile_c("src/struct_ls/semi")
    self.compile_c("src/struct_ls/smg_relax")
    self.compile_c("src/struct_ls/smg_setup")
    self.compile_c("src/struct_ls/smg_setup_rap")
    self.compile_c("src/struct_ls/smg_setup_restrict")
    self.compile_c("src/struct_ls/smg_solve")
    self.compile_c("src/struct_ls/sparse_msg")
    self.compile_c("src/struct_ls/sparse_msg_setup")
    self.compile_c("src/struct_ls/sparse_msg_setup_rap")
    self.compile_c("src/struct_ls/sparse_msg_solve")
    self.compile_c("src/struct_ls/cyclic_reduction")
    self.compile_c("src/struct_ls/HYPRE_struct_int")
    self.compile_c("src/struct_ls/HYPRE_struct_pcg")
    self.compile_c("src/struct_ls/pfmg2_setup_rap")
    self.compile_c("src/struct_ls/pfmg3_setup_rap")
    self.compile_c("src/struct_ls/pfmg_setup")
    self.compile_c("src/struct_ls/pfmg_setup_interp")
    self.compile_c("src/struct_ls/pfmg_setup_rap5")
    self.compile_c("src/struct_ls/pfmg_setup_rap7")
    self.compile_c("src/struct_ls/point_relax")
    self.compile_c("src/struct_ls/red_black_constantcoef_gs")
    self.compile_c("src/struct_ls/red_black_gs")
    self.compile_c("src/struct_ls/semi_interp")
    self.compile_c("src/struct_ls/semi_restrict")
    self.compile_c("src/struct_ls/semi_setup_rap")
    self.compile_c("src/struct_ls/smg2_setup_rap")
    self.compile_c("src/struct_ls/smg3_setup_rap")
    self.compile_c("src/struct_ls/smg")
    self.compile_c("src/struct_ls/smg_axpy")
    self.compile_c("src/struct_ls/smg_residual")
    self.compile_c("src/struct_ls/smg_setup_interp")
    self.compile_c("src/struct_ls/sparse_msg2_setup_rap")
    self.compile_c("src/struct_ls/sparse_msg3_setup_rap")
    self.compile_c("src/struct_ls/sparse_msg_filter")
    self.compile_c("src/struct_ls/sparse_msg_interp")
    self.compile_c("src/struct_ls/sparse_msg_restrict")

    print("\nBuilding " + self.name + " " + self.version + " sstruct_mv sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/sstruct_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/IJ_mv"'
    self.compile_c("src/sstruct_mv/HYPRE_sstruct_graph")
    self.compile_c("src/sstruct_mv/HYPRE_sstruct_grid")
    self.compile_c("src/sstruct_mv/HYPRE_sstruct_matrix")
    self.compile_c("src/sstruct_mv/HYPRE_sstruct_stencil")
    self.compile_c("src/sstruct_mv/HYPRE_sstruct_vector")
    self.compile_c("src/sstruct_mv/sstruct_axpy")
    self.compile_c("src/sstruct_mv/sstruct_copy")
    self.compile_c("src/sstruct_mv/sstruct_graph")
    self.compile_c("src/sstruct_mv/sstruct_grid")
    self.compile_c("src/sstruct_mv/sstruct_innerprod")
    self.compile_c("src/sstruct_mv/sstruct_matvec")
    self.compile_c("src/sstruct_mv/sstruct_scale")
    self.compile_c("src/sstruct_mv/sstruct_stencil")
    self.compile_c("src/sstruct_mv/sstruct_matrix")
    self.compile_c("src/sstruct_mv/sstruct_vector")

    print("\nBuilding " + self.name + " " + self.version + " sstruct_ls sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/sstruct_ls"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/multivector"'
    self.cxxflags += ' /I"./' + self.dir + '/src/krylov"'
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/IJ_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_block_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/sstruct_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/struct_ls"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_ls"'
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_bicgstab")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_gmres")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_flexgmres")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_lgmres")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_InterFAC")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_int")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_maxwell")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_pcg")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_split")
    self.compile_c("src/sstruct_ls/HYPRE_sstruct_sys_pfmg")
    self.compile_c("src/sstruct_ls/bsearch")
    self.compile_c("src/sstruct_ls/fac")
    self.compile_c("src/sstruct_ls/fac_amr_zero_data")
    self.compile_c("src/sstruct_ls/fac_cf_coarsen")
    self.compile_c("src/sstruct_ls/fac_cfstencil_box")
    self.compile_c("src/sstruct_ls/fac_CFInterfaceExtents")
    self.compile_c("src/sstruct_ls/fac_interp2")
    self.compile_c("src/sstruct_ls/fac_relax")
    self.compile_c("src/sstruct_ls/fac_solve3")
    self.compile_c("src/sstruct_ls/fac_zero_cdata")
    self.compile_c("src/sstruct_ls/krylov")
    self.compile_c("src/sstruct_ls/krylov_sstruct")
    self.compile_c("src/sstruct_ls/eliminate_rowscols")
    self.compile_c("src/sstruct_ls/maxwell_grad")
    self.compile_c("src/sstruct_ls/maxwell_physbdy")
    self.compile_c("src/sstruct_ls/maxwell_PNedelec")
    self.compile_c("src/sstruct_ls/maxwell_PNedelec_bdy")
    self.compile_c("src/sstruct_ls/maxwell_semi_interp")
    self.compile_c("src/sstruct_ls/maxwell_solve")
    self.compile_c("src/sstruct_ls/maxwell_solve2")
    self.compile_c("src/sstruct_ls/maxwell_TV")
    self.compile_c("src/sstruct_ls/maxwell_TV_setup")
    self.compile_c("src/sstruct_ls/maxwell_zeroBC")
    self.compile_c("src/sstruct_ls/nd1_amge_interpolation")
    self.compile_c("src/sstruct_ls/sstruct_amr_intercommunication")
    self.compile_c("src/sstruct_ls/sstruct_owninfo")
    self.compile_c("src/sstruct_ls/sstruct_recvinfo")
    self.compile_c("src/sstruct_ls/sstruct_sendinfo")
    self.compile_c("src/sstruct_ls/sstruct_sharedDOFComm")
    self.compile_c("src/sstruct_ls/sys_pfmg")
    self.compile_c("src/sstruct_ls/sys_pfmg_relax")
    self.compile_c("src/sstruct_ls/sys_pfmg_setup")
    self.compile_c("src/sstruct_ls/sys_pfmg_setup_interp")
    self.compile_c("src/sstruct_ls/sys_pfmg_setup_rap")
    self.compile_c("src/sstruct_ls/sys_pfmg_solve")
    self.compile_c("src/sstruct_ls/sys_semi_interp")
    self.compile_c("src/sstruct_ls/sys_semi_restrict")
    self.compile_c("src/sstruct_ls/fac_amr_fcoarsen")
    self.compile_c("src/sstruct_ls/fac_amr_rap")
    self.compile_c("src/sstruct_ls/fac_restrict2")
    self.compile_c("src/sstruct_ls/fac_setup2")
    self.compile_c("src/sstruct_ls/fac_zero_stencilcoef")
    self.compile_c("src/sstruct_ls/node_relax")

    print("\nBuilding " + self.name + " " + self.version + " distributed_ls/pilut sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_ls/pilut"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/blas"'
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_matrix"'
    self.compile_c("src/distributed_ls/pilut/comm")
    self.compile_c("src/distributed_ls/pilut/debug")
    self.compile_c("src/distributed_ls/pilut/distributed_qsort")
    self.compile_c("src/distributed_ls/pilut/distributed_qsort_si")
    self.compile_c("src/distributed_ls/pilut/HYPRE_DistributedMatrixPilutSolver")
    self.compile_c("src/distributed_ls/pilut/ilut")
    self.compile_c("src/distributed_ls/pilut/parilut")
    self.compile_c("src/distributed_ls/pilut/parutil")
    self.compile_c("src/distributed_ls/pilut/pblas1")
    self.compile_c("src/distributed_ls/pilut/serilut")
    self.compile_c("src/distributed_ls/pilut/trifactor")
    self.compile_c("src/distributed_ls/pilut/util")

    print("\nBuilding " + self.name + " " + self.version + " distributed_ls/ParaSails sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_ls/ParaSails"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/blas"'
    self.cxxflags += ' /I"./' + self.dir + '/src/lapack"'
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_matrix"'
    self.compile_c("src/distributed_ls/ParaSails/ConjGrad")
    self.compile_c("src/distributed_ls/ParaSails/DiagScale")
    self.compile_c("src/distributed_ls/ParaSails/FGmres")
    self.compile_c("src/distributed_ls/ParaSails/Hash")
    self.compile_c("src/distributed_ls/ParaSails/hypre_ParaSails")
    self.compile_c("src/distributed_ls/ParaSails/LoadBal")
    self.compile_c("src/distributed_ls/ParaSails/Matrix")
    self.compile_c("src/distributed_ls/ParaSails/Mem")
    self.compile_c("src/distributed_ls/ParaSails/Numbering")
    self.compile_c("src/distributed_ls/ParaSails/OrderStat")
    self.compile_c("src/distributed_ls/ParaSails/ParaSails")
    self.compile_c("src/distributed_ls/ParaSails/PrunedRows")
    self.compile_c("src/distributed_ls/ParaSails/RowPatt")
    self.compile_c("src/distributed_ls/ParaSails/StoredRows")

    print("\nBuilding " + self.name + " " + self.version + " distributed_ls/Euclid sources...")
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_ls/Euclid"'
    self.cxxflags += ' /I"./' + self.dir + '/src/utilities"'
    self.cxxflags += ' /I"./' + self.dir + '/src/parcsr_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/seq_mv"'
    self.cxxflags += ' /I"./' + self.dir + '/src/distributed_matrix"'
    self.compile_c("src/distributed_ls/Euclid/blas_dh")
    self.compile_c("src/distributed_ls/Euclid/Euclid_apply")
    self.compile_c("src/distributed_ls/Euclid/Euclid_dh")
    self.compile_c("src/distributed_ls/Euclid/ExternalRows_dh")
    self.compile_c("src/distributed_ls/Euclid/Factor_dh")
    self.compile_c("src/distributed_ls/Euclid/getRow_dh")
    self.compile_c("src/distributed_ls/Euclid/globalObjects")
    self.compile_c("src/distributed_ls/Euclid/Hash_dh")
    self.compile_c("src/distributed_ls/Euclid/Hash_i_dh")
    self.compile_c("src/distributed_ls/Euclid/ilu_mpi_bj")
    self.compile_c("src/distributed_ls/Euclid/ilu_mpi_pilu")
    self.compile_c("src/distributed_ls/Euclid/ilu_seq")
    self.compile_c("src/distributed_ls/Euclid/io_dh")
    self.compile_c("src/distributed_ls/Euclid/krylov_dh")
    self.compile_c("src/distributed_ls/Euclid/Mat_dh")
    self.compile_c("src/distributed_ls/Euclid/mat_dh_private")
    self.compile_c("src/distributed_ls/Euclid/MatGenFD")
    self.compile_c("src/distributed_ls/Euclid/Mem_dh")
    self.compile_c("src/distributed_ls/Euclid/Numbering_dh")
    self.compile_c("src/distributed_ls/Euclid/Parser_dh")
    self.compile_c("src/distributed_ls/Euclid/shellSort_dh")
    self.compile_c("src/distributed_ls/Euclid/sig_dh")
    self.compile_c("src/distributed_ls/Euclid/SortedList_dh")
    self.compile_c("src/distributed_ls/Euclid/SortedSet_dh")
    self.compile_c("src/distributed_ls/Euclid/SubdomainGraph_dh")
    self.compile_c("src/distributed_ls/Euclid/TimeLog_dh")
    self.compile_c("src/distributed_ls/Euclid/Timer_dh")
    self.compile_c("src/distributed_ls/Euclid/Vec_dh")

    self.link()

packages["hypre"] = ThirdPartyHypre()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# P A R M E T I S  -  P A R M E T I S  -  P A R M E T I S  -  P A R M E T I S  -  P A R M E T I S  -  P A R M E T I S  -
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyParMETIS(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "ParMETIS"
    self.version = "4.0.3"
    self.date = "2013-03-30"
    self.file = "parmetis-" + self.version + ".tar.gz"
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.need_mpi = True
    self.license_info = \
      "ParMETIS has a rather restrictive license, which is the first one shown here. Additionally,\n" +\
      "ParMETIS relies on METIS. which has its own license, which is the second license shown below\n\n"
    self.license_files = ["LICENSE.txt", os.path.join("metis", "LICENSE.txt")]
    self.url = "http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/" + self.file
    self.page = "http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview"
    self.baseflags += " /DUSE_GKREGEX"
    self.baseflags += ' /wd"4028" /wd"4018" /wd"4101" /wd"4102" /wd"4244" /wd"4267" /wd"4305" /wd"4996"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The ParMETIS library offers a set of MPI-parallel graph and mesh partitioning algorithms.")
    print("Installing this package enables the use of the geometric 'parmetis' partitioner as a valid choice")
    print("for the Control::PartiDomainControl class.")
    print("")
    print("Important: Please note that ParMETIS has a rather restrictive license model.")
    print("")
    print("Prerequisites/Constraints: requires MPI to be installed (obviously)")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_PARMETIS' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def patch(self):
    print("\nPatching " + self.name + " " + self.version + " sources...")
    # remove some outdated MS-specific preprocessor mumbo-jumbo
    self.patch_file(os.path.join("metis", "GKlib", "gk_arch.h"), [
      [35, '#include "ms_stdint.h"', '  #include <stdint.h>'],
      [36, '#include "ms_inttypes.h"', '  #include <inttypes.h>'],
      [48, '', '#if defined(_WIN32) && !defined(WIN32)\n  #define WIN32\n#endif'],
      [61, '#ifdef __MSC__', '#if 0'],
      [70, '', '#ifdef __MSC__\n#define __thread __declspec(thread)\n#endif']
    ])
    # patch some broken pointer casts
    self.patch_file(os.path.join("metis", "GKlib", "gkregex.c"), [
      [5089, 'postorder (elem, mark_opt_subexp, (void *) (long) elem->token.opr.idx);',
         '    postorder (elem, mark_opt_subexp, (void *) (intptr_t) elem->token.opr.idx);'],
      [6301, 'int idx = (int) (long) extra;', '  int idx = (int) (intptr_t) extra;']
    ])
    # patch more outdated MS-specific preprocessor mumbo-jumbo
    self.patch_file(os.path.join("metis", "include", "metis.h"), [
      [66, '#ifdef COMPILER_MSC', '#if 0 /*def COMPILER_MSC*/']
    ])

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " GKlib sources...")
    self.extflags = ' /I"./' + os.path.join(self.dir, "metis", "GKlib") + '"'
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.compile_c("metis/GKlib/b64")
    self.compile_c("metis/GKlib/blas")
    self.compile_c("metis/GKlib/csr")
    self.compile_c("metis/GKlib/error")
    self.compile_c("metis/GKlib/evaluate")
    self.compile_c("metis/GKlib/fkvkselect")
    self.compile_c("metis/GKlib/fs")
    self.compile_c("metis/GKlib/getopt")
    self.compile_c("metis/GKlib/gkregex")
    self.compile_c("metis/GKlib/graph")
    self.compile_c("metis/GKlib/htable")
    self.compile_c("metis/GKlib/io")
    self.compile_c("metis/GKlib/itemsets")
    self.compile_c("metis/GKlib/mcore")
    self.compile_c("metis/GKlib/memory")
    self.compile_c("metis/GKlib/omp")
    self.compile_c("metis/GKlib/pdb")
    self.compile_c("metis/GKlib/pqueue")
    self.compile_c("metis/GKlib/random")
    self.compile_c("metis/GKlib/rw")
    self.compile_c("metis/GKlib/seq")
    self.compile_c("metis/GKlib/sort")
    self.compile_c("metis/GKlib/string")
    self.compile_c("metis/GKlib/timers")
    self.compile_c("metis/GKlib/tokenizer")
    self.compile_c("metis/GKlib/util")
    print("\nBuilding " + self.name + " " + self.version + " METIS sources...")
    self.extflags  = ' /I"./' + os.path.join(self.dir, "metis", "include") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "metis", "libmetis") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "metis", "GKlib") + '"'
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.compile_c("metis/libmetis/auxapi")
    self.compile_c("metis/libmetis/balance")
    self.compile_c("metis/libmetis/bucketsort")
    self.compile_c("metis/libmetis/checkgraph")
    self.compile_c("metis/libmetis/coarsen")
    self.compile_c("metis/libmetis/compress")
    self.compile_c("metis/libmetis/contig")
    self.compile_c("metis/libmetis/debug")
    self.compile_c("metis/libmetis/fm")
    self.compile_c("metis/libmetis/fortran")
    self.compile_c("metis/libmetis/frename")
    self.compile_c("metis/libmetis/gklib")
    self.compile_c("metis/libmetis/graph")
    self.compile_c("metis/libmetis/initpart")
    self.compile_c("metis/libmetis/kmetis")
    self.compile_c("metis/libmetis/kwayfm")
    self.compile_c("metis/libmetis/kwayrefine")
    self.compile_c("metis/libmetis/mcutil")
    self.compile_c("metis/libmetis/mesh")
    self.compile_c("metis/libmetis/meshpart")
    self.compile_c("metis/libmetis/minconn")
    self.compile_c("metis/libmetis/mincover")
    self.compile_c("metis/libmetis/mmd")
    self.compile_c("metis/libmetis/ometis")
    self.compile_c("metis/libmetis/options")
    self.compile_c("metis/libmetis/parmetis")
    self.compile_c("metis/libmetis/pmetis")
    self.compile_c("metis/libmetis/refine")
    self.compile_c("metis/libmetis/separator")
    self.compile_c("metis/libmetis/sfm")
    self.compile_c("metis/libmetis/srefine")
    self.compile_c("metis/libmetis/stat")
    self.compile_c("metis/libmetis/timing")
    self.compile_c("metis/libmetis/util")
    self.compile_c("metis/libmetis/wspace")
    print("\nBuilding " + self.name + " " + self.version + " ParMETIS sources...")
    self.extflags  = ' /I"' + os.path.join(path_mpi, "Include") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "include") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "libparmetis") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "metis", "include") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "metis", "libmetis") + '"'
    self.extflags += ' /I"./' + os.path.join(self.dir, "metis", "GKlib") + '"'
    self.cxxflags = self.baseflags + self.modeflags + self.extflags + self.outflags
    self.compile_c("libparmetis/akwayfm")
    self.compile_c("libparmetis/ametis")
    self.compile_c("libparmetis/balancemylink")
    self.compile_c("libparmetis/comm")
    self.compile_c("libparmetis/csrmatch")
    self.compile_c("libparmetis/ctrl")
    self.compile_c("libparmetis/debug")
    self.compile_c("libparmetis/diffutil")
    self.compile_c("libparmetis/frename")
    self.compile_c("libparmetis/gkmetis")
    self.compile_c("libparmetis/gkmpi")
    self.compile_c("libparmetis/graph")
    self.compile_c("libparmetis/initbalance")
    self.compile_c("libparmetis/initmsection")
    self.compile_c("libparmetis/initpart")
    self.compile_c("libparmetis/kmetis")
    self.compile_c("libparmetis/kwayrefine")
    self.compile_c("libparmetis/match")
    self.compile_c("libparmetis/mdiffusion")
    self.compile_c("libparmetis/mesh")
    self.compile_c("libparmetis/mmetis")
    self.compile_c("libparmetis/move")
    self.compile_c("libparmetis/msetup")
    self.compile_c("libparmetis/node_refine")
    self.compile_c("libparmetis/ometis")
    self.compile_c("libparmetis/pspases")
    self.compile_c("libparmetis/redomylink")
    self.compile_c("libparmetis/remap")
    self.compile_c("libparmetis/renumber")
    self.compile_c("libparmetis/rmetis")
    self.compile_c("libparmetis/selectq")
    self.compile_c("libparmetis/serial")
    self.compile_c("libparmetis/stat")
    self.compile_c("libparmetis/timer")
    self.compile_c("libparmetis/util")
    self.compile_c("libparmetis/wave")
    self.compile_c("libparmetis/weird")
    self.compile_c("libparmetis/wspace")
    self.compile_c("libparmetis/xyzpart")

    self.link()

packages["parmetis"] = ThirdPartyParMETIS()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# S U I T E S P A R S E  -  S U I T E S P A R S E  -  S U I T E S P A R S E  -  S U I T E S P A R S E  -  S U I T E S P
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartySuiteSparse(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "SuiteSparse"
    self.version = "5.13.0"
    self.date = "2022-08-25"
    self.file = self.name + "-" + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.license_info = \
      "FEAT only uses the AMD and UMFPACK libraries from SuiteSparse, so only the licenses for these\n" +\
      "two packages are shown here:\n\n"
    self.license_files = [os.path.join("AMD", "Doc", "License.txt"), os.path.join("UMFPACK", "Doc", "License.txt")]
    self.url = "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v" + self.version + ".zip"
    self.page = "https://people.engr.tamu.edu/davis/suitesparse.html"
    self.baseflags += ' /DNBLAS /DNCHOLMOD /D_MBCS'
    self.baseflags += ' /wd"4068" /wd"4101" /wd"4244" /wd"4267" /wd"4996"'
    self.baseflags += ' /I"./' + self.dir + '/UMFPACK/Include"'
    self.baseflags += ' /I"./' + self.dir + '/AMD/Include"'
    self.baseflags += ' /I"./' + self.dir + '/SuiteSparse_config"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The SuiteSparse library is a large collection of various linear solver related algorithms, most")
    print("notably the UMFPACK sparse direct solver, which is currently the only part of the SuiteSparse")
    print("library collection that is used by FEAT, although this may change in future. Installing this package")
    print("enables the use of the Solver::Umfpack class as well as its related solver classes.")
    print("")
    print("Prerequisites/Constraints: none")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_UMFPACK' as well as the 'FEAT_HAVE_SUITESPARSE'")
    print("preprocessor macros.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  # compile 'DINT' and 'DLONG' version of source file
  def compile_c_il(self, filename, extflags = "", extname = ""):
    self.compile(filename + ".c", filename.replace("/","_") + extname + "_i.obj", extflags + ' /TC /DDINT')
    self.compile(filename + ".c", filename.replace("/","_") + extname + "_l.obj", extflags + ' /TC /DDLONG')

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " SuiteSparse_config sources...")
    self.compile_c("SuiteSparse_config/SuiteSparse_config")
    print("\nBuilding " + self.name + " " + self.version + " AMD sources...")
    self.compile_c("AMD/Source/amd_global")
    self.compile_c_il("AMD/Source/amd_aat")
    self.compile_c_il("AMD/Source/amd_1")
    self.compile_c_il("AMD/Source/amd_2")
    self.compile_c_il("AMD/Source/amd_dump")
    self.compile_c_il("AMD/Source/amd_postorder")
    self.compile_c_il("AMD/Source/amd_post_tree")
    self.compile_c_il("AMD/Source/amd_defaults")
    self.compile_c_il("AMD/Source/amd_order")
    self.compile_c_il("AMD/Source/amd_control")
    self.compile_c_il("AMD/Source/amd_info")
    self.compile_c_il("AMD/Source/amd_valid")
    self.compile_c_il("AMD/Source/amd_preprocess")
    print("\nBuilding " + self.name + " " + self.version + " UMFPACK sources...")
    self.compile_c("UMFPACK/Source/umfpack_timer")
    self.compile_c("UMFPACK/Source/umfpack_tictoc")
    self.compile_c_il("UMFPACK/Source/umf_analyze")
    self.compile_c_il("UMFPACK/Source/umf_apply_order")
    self.compile_c_il("UMFPACK/Source/umf_colamd")
    self.compile_c_il("UMFPACK/Source/umf_cholmod")
    self.compile_c_il("UMFPACK/Source/umf_free")
    self.compile_c_il("UMFPACK/Source/umf_fsize")
    self.compile_c_il("UMFPACK/Source/umf_is_permutation")
    self.compile_c_il("UMFPACK/Source/umf_malloc")
    self.compile_c_il("UMFPACK/Source/umf_realloc")
    self.compile_c_il("UMFPACK/Source/umf_report_perm")
    self.compile_c_il("UMFPACK/Source/umf_singletons")
    self.compile_c_il("UMFPACK/Source/umf_assemble")
    self.compile_c_il("UMFPACK/Source/umf_blas3_update")
    self.compile_c_il("UMFPACK/Source/umf_build_tuples")
    self.compile_c_il("UMFPACK/Source/umf_create_element")
    self.compile_c_il("UMFPACK/Source/umf_dump")
    self.compile_c_il("UMFPACK/Source/umf_extend_front")
    self.compile_c_il("UMFPACK/Source/umf_garbage_collection")
    self.compile_c_il("UMFPACK/Source/umf_get_memory")
    self.compile_c_il("UMFPACK/Source/umf_init_front")
    self.compile_c_il("UMFPACK/Source/umf_kernel")
    self.compile_c_il("UMFPACK/Source/umf_kernel_init")
    self.compile_c_il("UMFPACK/Source/umf_kernel_wrapup")
    self.compile_c_il("UMFPACK/Source/umf_local_search")
    self.compile_c_il("UMFPACK/Source/umf_lsolve")
    self.compile_c_il("UMFPACK/Source/umf_ltsolve")
    self.compile_c_il("UMFPACK/Source/umf_mem_alloc_element")
    self.compile_c_il("UMFPACK/Source/umf_mem_alloc_head_block")
    self.compile_c_il("UMFPACK/Source/umf_mem_alloc_tail_block")
    self.compile_c_il("UMFPACK/Source/umf_mem_free_tail_block")
    self.compile_c_il("UMFPACK/Source/umf_mem_init_memoryspace")
    self.compile_c_il("UMFPACK/Source/umf_report_vector")
    self.compile_c_il("UMFPACK/Source/umf_row_search")
    self.compile_c_il("UMFPACK/Source/umf_scale_column")
    self.compile_c_il("UMFPACK/Source/umf_set_stats")
    self.compile_c_il("UMFPACK/Source/umf_solve")
    self.compile_c_il("UMFPACK/Source/umf_symbolic_usage")
    self.compile_c_il("UMFPACK/Source/umf_transpose")
    self.compile_c_il("UMFPACK/Source/umf_tuple_lengths")
    self.compile_c_il("UMFPACK/Source/umf_usolve")
    self.compile_c_il("UMFPACK/Source/umf_utsolve")
    self.compile_c_il("UMFPACK/Source/umf_valid_numeric")
    self.compile_c_il("UMFPACK/Source/umf_valid_symbolic")
    self.compile_c_il("UMFPACK/Source/umf_grow_front")
    self.compile_c_il("UMFPACK/Source/umf_start_front")
    self.compile_c_il("UMFPACK/Source/umf_store_lu")
    self.compile_c_il("UMFPACK/Source/umf_scale")
    self.compile_c_il("UMFPACK/Source/umfpack_col_to_triplet")
    self.compile_c_il("UMFPACK/Source/umfpack_defaults")
    self.compile_c_il("UMFPACK/Source/umfpack_free_numeric")
    self.compile_c_il("UMFPACK/Source/umfpack_free_symbolic")
    self.compile_c_il("UMFPACK/Source/umfpack_get_numeric")
    self.compile_c_il("UMFPACK/Source/umfpack_get_lunz")
    self.compile_c_il("UMFPACK/Source/umfpack_get_symbolic")
    self.compile_c_il("UMFPACK/Source/umfpack_get_determinant")
    self.compile_c_il("UMFPACK/Source/umfpack_numeric")
    self.compile_c_il("UMFPACK/Source/umfpack_qsymbolic")
    self.compile_c_il("UMFPACK/Source/umfpack_report_control")
    self.compile_c_il("UMFPACK/Source/umfpack_report_info")
    self.compile_c_il("UMFPACK/Source/umfpack_report_matrix")
    self.compile_c_il("UMFPACK/Source/umfpack_report_numeric")
    self.compile_c_il("UMFPACK/Source/umfpack_report_perm")
    self.compile_c_il("UMFPACK/Source/umfpack_report_status")
    self.compile_c_il("UMFPACK/Source/umfpack_report_symbolic")
    self.compile_c_il("UMFPACK/Source/umfpack_report_triplet")
    self.compile_c_il("UMFPACK/Source/umfpack_report_vector")
    self.compile_c_il("UMFPACK/Source/umfpack_solve")
    self.compile_c_il("UMFPACK/Source/umfpack_symbolic")
    self.compile_c_il("UMFPACK/Source/umfpack_transpose")
    self.compile_c_il("UMFPACK/Source/umfpack_triplet_to_col")
    self.compile_c_il("UMFPACK/Source/umfpack_scale")
    self.compile_c_il("UMFPACK/Source/umfpack_load_numeric")
    self.compile_c_il("UMFPACK/Source/umfpack_save_numeric")
    self.compile_c_il("UMFPACK/Source/umfpack_load_symbolic")
    self.compile_c_il("UMFPACK/Source/umfpack_save_symbolic")
    # preprocessor abuse at its finest...
    self.compile_c_il("UMFPACK/Source/umf_ltsolve")
    self.compile_c_il("UMFPACK/Source/umf_ltsolve", " /DCONJUGATE_SOLVE", '_cs')
    self.compile_c_il("UMFPACK/Source/umf_utsolve")
    self.compile_c_il("UMFPACK/Source/umf_utsolve", " /DCONJUGATE_SOLVE", '_cs')
    self.compile_c_il("UMFPACK/Source/umf_assemble")
    self.compile_c_il("UMFPACK/Source/umf_assemble", '/DFIXQ', '_fixq')
    self.compile_c_il("UMFPACK/Source/umf_store_lu")
    self.compile_c_il("UMFPACK/Source/umf_store_lu", '/DDROP', '_drop')
    self.compile_c_il("UMFPACK/Source/umfpack_solve")
    self.compile_c_il("UMFPACK/Source/umfpack_solve", " /DWSOLVE", '_w')
    self.compile_c_il("UMFPACK/Source/umf_triplet")
    self.compile_c_il("UMFPACK/Source/umf_triplet", " /DDO_MAP", "_map_novalues")
    self.compile_c_il("UMFPACK/Source/umf_triplet", " /DDO_VALUES", "_nomap_values")
    self.compile_c_il("UMFPACK/Source/umf_triplet", " /DDO_MAP /DDO_VALUES", "_map_values")
    self.link()

packages["suitesparse"] = ThirdPartySuiteSparse()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# S U P E R L U  -  S U P E R L U  -  S U P E R L U  -  S U P E R L U  -  S U P E R L U  -  S U P E R L U  -  S U P E R
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartySuperLU(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "SuperLU"
    self.version = "8.1.1"
    self.date = "2022-10-01"
    self.prefix = "_dist"
    self.file = self.name + self.prefix + "-" + self.version + ".zip"
    self.dir = self.name.lower() + self.prefix + "-" + self.version
    self.trunk = "."
    self.need_mpi = True
    self.license_files = ["License.txt"]
    self.url = "https://github.com/xiaoyeli/superlu_dist/archive/refs/tags/v" + self.version + ".zip"
    self.page = "https://github.com/xiaoyeli/superlu_dist"
    self.baseflags += ' /wd"4996" /wd"4005" /wd"4013" /wd"4068" /wd"4101" /wd"4146" /wd"4244" /wd"4267" /wd"4305" /wd"4334" /wd"4477" /wd"4700" /wd"4715"'
    self.baseflags += ' /I"' + os.path.join(path_mpi, "Include") + '"'
    self.baseflags += ' /I"./' + self.dir + '/CBLAS"'
    self.baseflags += ' /I"./' + self.dir + '/SRC"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The SuperLU library (or more precisely: the SuperLU_dist library) offers an MPI-parallel")
    print("implementation of LU-decomposition based direct solver. Note that although a sequential version")
    print("of SuperLU exists, FEAT only uses the MPI-parallel version; if you require a sequential direct")
    print("solver, please refer to the superior UMPFACK third-party library instead. Installing this package")
    print("enables the use of the Solver::SuperLU solver class.")
    print("")
    print("Prerequisites/Constraints: requires MPI to be installed (obviously)")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_SUPERLU' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def patch(self):
    print("\nPatching " + self.name + " " + self.version + " sources...")
    # don't use ParMETIS
    self.patch_file(os.path.join("SRC", "superlu_dist_config.h"), [
      [10, '#define HAVE_PARMETIS TRUE', '/* #undef HAVE_PARMETIS */'],
    ])
    # don't include unistd.h header
    self.patch_file(os.path.join("SRC", "util.c"), [
      [25, '#include <unistd.h>', '/*#include <unistd.h>*/'],
    ])

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " CBLAS sources...")
    self.compile_c("CBLAS/dgemm")
    self.compile_c("CBLAS/dgemv")
    self.compile_c("CBLAS/dger")
    self.compile_c("CBLAS/dtrsm")
    self.compile_c("CBLAS/dtrsv")
    self.compile_c("CBLAS/daxpy")
    self.compile_c("CBLAS/dscal")
    self.compile_c("CBLAS/input_error_dist")

    print("\nBuilding " + self.name + " " + self.version + " sources...")
    self.compile_c("SRC/sp_ienv")
    self.compile_c("SRC/etree")
    self.compile_c("SRC/sp_colorder")
    self.compile_c("SRC/get_perm_c")
    self.compile_c("SRC/mmd")
    self.compile_c("SRC/comm")
    self.compile_c("SRC/memory")
    self.compile_c("SRC/util")
    self.compile_c("SRC/gpu_api_utils")
    self.compile_c("SRC/superlu_grid")
    self.compile_c("SRC/pxerr_dist")
    self.compile_c("SRC/superlu_timer")
    self.compile_c("SRC/symbfact")
    self.compile_c("SRC/psymbfact")
    self.compile_c("SRC/psymbfact_util")
    self.compile_c("SRC/get_perm_c_parmetis")
    self.compile_c("SRC/mc64ad_dist")
    self.compile_c("SRC/xerr_dist")
    self.compile_c("SRC/smach_dist")
    self.compile_c("SRC/dmach_dist")
    self.compile_c("SRC/superlu_dist_version")
    self.compile_c("SRC/comm_tree")
    self.compile_c("SRC/superlu_grid3d")
    self.compile_c("SRC/supernodal_etree")
    self.compile_c("SRC/supernodalForest")
    self.compile_c("SRC/trfAux")
    self.compile_c("SRC/communication_aux")
    self.compile_c("SRC/treeFactorization")
    self.compile_c("SRC/sec_structs")
    self.compile_c("SRC/wingetopt")
    self.compile_c("SRC/dlangs_dist")
    self.compile_c("SRC/dgsequ_dist")
    self.compile_c("SRC/dlaqgs_dist")
    self.compile_c("SRC/dutil_dist")
    self.compile_c("SRC/dmemory_dist")
    self.compile_c("SRC/dmyblas2_dist")
    self.compile_c("SRC/dsp_blas2_dist")
    self.compile_c("SRC/dsp_blas3_dist")
    self.compile_c("SRC/pdgssvx")
    self.compile_c("SRC/pdgssvx_ABglobal")
    self.compile_c("SRC/dreadhb")
    self.compile_c("SRC/dreadrb")
    self.compile_c("SRC/dreadtriple")
    self.compile_c("SRC/dreadtriple_noheader")
    self.compile_c("SRC/dbinary_io")
    self.compile_c("SRC/dreadMM")
    self.compile_c("SRC/pdgsequ")
    self.compile_c("SRC/pdlaqgs")
    self.compile_c("SRC/dldperm_dist")
    self.compile_c("SRC/pdlangs")
    self.compile_c("SRC/pdutil")
    self.compile_c("SRC/pdsymbfact_distdata")
    self.compile_c("SRC/ddistribute")
    self.compile_c("SRC/pddistribute")
    self.compile_c("SRC/pdgstrf")
    self.compile_c("SRC/dstatic_schedule")
    self.compile_c("SRC/pdgstrf2")
    self.compile_c("SRC/pdgstrs")
    self.compile_c("SRC/pdgstrs1")
    self.compile_c("SRC/pdgstrs_lsum")
    self.compile_c("SRC/pdgstrs_Bglobal")
    self.compile_c("SRC/pdgsrfs")
    self.compile_c("SRC/pdgsmv")
    self.compile_c("SRC/pdgsrfs_ABXglobal")
    self.compile_c("SRC/pdgsmv_AXglobal")
    self.compile_c("SRC/pdGetDiagU")
    self.compile_c("SRC/pdgssvx3d")
    self.compile_c("SRC/dnrformat_loc3d")
    self.compile_c("SRC/pdgstrf3d")
    self.compile_c("SRC/dtreeFactorization")
    self.compile_c("SRC/dtreeFactorizationGPU")
    self.compile_c("SRC/dgather")
    self.compile_c("SRC/dscatter3d")
    self.compile_c("SRC/pd3dcomm")
    self.compile_c("SRC/dtrfAux")
    self.compile_c("SRC/dcommunication_aux")
    self.compile_c("SRC/dtrfCommWrapper")
    self.compile_c("SRC/dsuperlu_blas")
    self.link()

packages["superlu"] = ThirdPartySuperLU()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# T R I A N G L E  -  T R I A N G L E  -  T R I A N G L E  -  T R I A N G L E  -  T R I A N G L E  -  T R I A N G L E  -
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyTriangle(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "triangle"
    self.version = "1.6"
    self.date = "2005-07-28"
    self.file = self.name + "-" + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = self.dir
    self.license_info = "Note: triangle does not come with a dedicated license file, so here comes the entire README instead.\n"
    self.license_files = ["README"] # no dedicated license file available
    self.url = "http://www.netlib.org/voronoi/triangle.zip"
    self.page = "https://www.cs.cmu.edu/~quake/triangle.html"
    self.baseflags += ' /wd"4996" /wd"4267" /wd"4244" /wd"4311" /wd"4312"'
    self.baseflags += ' /DANSI_DECLARATORS /DTRILIBRARY /DNO_TIMER /DCPU86'
    self.baseflags += ' /I"./' + self.dir + '"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The triangle library is a lightweight collection of 2D triangulation algorithms, which can be used")
    print("to generate triangle meshes. Installing this package currently does not enable anything in FEAT,")
    print("but this will hopefully change in future once some student workers have been put on this job.")
    print("")
    print("Prerequisites/Constraints: none")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_TRIANGLE' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " sources...")
    self.compile_c("triangle")
    self.link()

packages["triangle"] = ThirdPartyTriangle()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -  Z F P  -
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyZfp(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "zfp"
    self.version = "1.0.0"
    self.date = "2022-08-01"
    self.file = self.name + "-" + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.license_files = ["LICENSE"]
    self.url = "https://github.com/LLNL/zfp/releases/download/" + self.version + "/" + self.file
    self.page = "https://computing.llnl.gov/projects/zfp"
    self.baseflags += ' /DZFP_INT64="long long" /DZFP_INT64_SUFFIX="ll"'
    self.baseflags += ' /DZFP_UINT64="unsigned long long" /DZFP_UINT64_SUFFIX="ull"'
    self.baseflags += ' /wd"4996" /wd"4267" /wd"4146"'
    self.baseflags += ' /I"./' + self.dir + '/include"'
    self.baseflags += ' /I"./' + self.dir + '/src"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The ZFP library is small library offering various lossy floating point compression algorithms, which")
    print("can be used to compress large binary output files. Installing this package enables FEAT to use ZFP")
    print("compression algorithms for the binary I/O of LAFEM containers.")
    print("")
    print("Prerequisites/Constraints: none")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_ZFP' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " sources...")
    self.compile_c("src/bitstream")
    self.compile_c("src/decode1d")
    self.compile_c("src/decode1f")
    self.compile_c("src/decode1i")
    self.compile_c("src/decode1l")
    self.compile_c("src/decode2d")
    self.compile_c("src/decode2f")
    self.compile_c("src/decode2i")
    self.compile_c("src/decode2l")
    self.compile_c("src/decode3d")
    self.compile_c("src/decode3f")
    self.compile_c("src/decode3i")
    self.compile_c("src/decode3l")
    self.compile_c("src/decode4d")
    self.compile_c("src/decode4f")
    self.compile_c("src/decode4i")
    self.compile_c("src/decode4l")
    self.compile_c("src/encode1d")
    self.compile_c("src/encode1f")
    self.compile_c("src/encode1i")
    self.compile_c("src/encode1l")
    self.compile_c("src/encode2d")
    self.compile_c("src/encode2f")
    self.compile_c("src/encode2i")
    self.compile_c("src/encode2l")
    self.compile_c("src/encode3d")
    self.compile_c("src/encode3f")
    self.compile_c("src/encode3i")
    self.compile_c("src/encode3l")
    self.compile_c("src/encode4d")
    self.compile_c("src/encode4f")
    self.compile_c("src/encode4i")
    self.compile_c("src/encode4l")
    self.compile_c("src/zfp")
    self.link()

packages["zfp"] = ThirdPartyZfp()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -  Z L I B  -
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyZlib(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "zlib"
    self.version = "1.2.13"
    self.date = "2022-10-13"
    self.file = self.name + self.version.replace(".", "") + ".zip" # filename without version dots
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.license_files = ["LICENSE"]
    self.url = "https://zlib.net/fossils/" + self.file
    self.page = "https://zlib.net"
    self.baseflags += ' /wd"4996" /wd"4244" /wd"4267"'
    self.baseflags += ' /I"./' + self.dir + '"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The Zlib library is small library offering lossless compression algorithms based on the 'deflate'")
    print("method, which is also the backbone of the famous ZIP file format. Installing this package enables")
    print("FEAT to use ZLIB compression algorithms for the binary I/O of LAFEM containers.")
    print("")
    print("Prerequisites/Constraints: none")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_ZLIB' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " sources...")
    self.compile_c("adler32")
    self.compile_c("compress")
    self.compile_c("crc32")
    self.compile_c("deflate")
    self.compile_c("gzclose")
    self.compile_c("gzlib")
    self.compile_c("gzread")
    self.compile_c("gzwrite")
    self.compile_c("infback")
    self.compile_c("inffast")
    self.compile_c("inflate")
    self.compile_c("inftrees")
    self.compile_c("trees")
    self.compile_c("uncompr")
    self.compile_c("zutil")
    self.link()

packages["zlib"] = ThirdPartyZlib()

########################################################################################################################
########################################################################################################################
########################################################################################################################
# Z O L T A N  -  Z O L T A N  -  Z O L T A N  -  Z O L T A N  -  Z O L T A N  -  Z O L T A N  -  Z O L T A N  -  Z O L
########################################################################################################################
########################################################################################################################
########################################################################################################################

class ThirdPartyZoltan(ThirdPartyPackage):
  def __init__(self):
    super().__init__()
    self.name = "Zoltan"
    self.version = "3.901"
    self.date = "2022-02-01"
    self.file = self.name + "-" + self.version + ".zip"
    self.dir = self.name + "-" + self.version
    self.trunk = "."
    self.need_mpi = True
    self.license_files = ["COPYRIGHT_AND_LICENSE"]
    self.url = "https://github.com/sandialabs/Zoltan/archive/refs/tags/v" + self.version + ".zip"
    self.page = "https://sandialabs.github.io/Zoltan"
    self.baseflags += ' /wd"4028" /wd"4018" /wd"4090" /wd"4101" /wd"4102" /wd"4113" /wd"4244" /wd"4267" /wd"4305" /wd"4311" /wd"4477" /wd"4716" /wd"4996"'
    self.baseflags += ' /I"' + os.path.join(path_mpi, "Include") + '"'
    self.baseflags += ' /I"./' + self.dir + '/src"'
    self.baseflags += ' /I"./' + self.dir + '/src/all"'
    self.baseflags += ' /I"./' + self.dir + '/src/coloring"'
    self.baseflags += ' /I"./' + self.dir + '/src/graph"'
    self.baseflags += ' /I"./' + self.dir + '/src/ha"'
    self.baseflags += ' /I"./' + self.dir + '/src/hier"'
    self.baseflags += ' /I"./' + self.dir + '/src/hsfc"'
    self.baseflags += ' /I"./' + self.dir + '/src/include"'
    self.baseflags += ' /I"./' + self.dir + '/src/lb"'
    self.baseflags += ' /I"./' + self.dir + '/src/matrix"'
    self.baseflags += ' /I"./' + self.dir + '/src/order"'
    self.baseflags += ' /I"./' + self.dir + '/src/par"'
    self.baseflags += ' /I"./' + self.dir + '/src/params"'
    self.baseflags += ' /I"./' + self.dir + '/src/phg"'
    self.baseflags += ' /I"./' + self.dir + '/src/rcb"'
    self.baseflags += ' /I"./' + self.dir + '/src/reftree"'
    self.baseflags += ' /I"./' + self.dir + '/src/simple"'
    self.baseflags += ' /I"./' + self.dir + '/src/timer"'
    self.baseflags += ' /I"./' + self.dir + '/src/tpls"'
    self.baseflags += ' /I"./' + self.dir + '/src/Utilities/shared"'
    self.baseflags += ' /I"./' + self.dir + '/src/Utilities/Timer"'
    self.baseflags += ' /I"./' + self.dir + '/src/zz"'

  def info(self):
    #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    print("The Zoltan library offers a set of MPI-parallel hypergraph partitioning algorithms.")
    print("Installing this package enables the use of the 'zoltan' partitioner as a valid choice")
    print("for the Control::PartiDomainControl class.")
    print("")
    print("Prerequisites/Constraints: requires MPI to be installed (obviously)")
    print("")
    print("Installing this package defines the 'FEAT_HAVE_ZOLTAN' preprocessor macro.")
    print("")
    print("Visit the " + self.name + " homepage for more information:")
    print(self.page)

  def patch(self):
    print("\nPatching " + self.name + " " + self.version + " sources...")
    # write configure file
    cfg_file = os.path.join(self.dir, "src", "include", "Zoltan_config.h")
    if not os.path.isfile(cfg_file):
      print("Writing Config Header '%s'..." % cfg_file)
      fo = open(cfg_file, "wt")
      fo.write("/* Zoltan config header generated by FEAT3 'vc17-thirdparty.py' script */\n")
      fo.write("#define HAVE_MPI 1\n")
      fo.write("#define UNSIGNED_LONG_LONG_GLOBAL_IDS 1\n")
      fo.write("#define strcasecmp _stricmp\n")
      fo.write("#define strncasecmp _strnicmp\n")
      fo.close()

    # patch "phg/phg_util.c"
    self.patch_file(os.path.join("src", "phg", "phg_util.c"), [
      [75, "{", "{/*"],
      [85, "}", "*/}"]
    ])
    # patch "rcb/inertial2d.c"
    self.patch_file(os.path.join("src", "rcb", "inertial2d.c"), [
      [66, "#define max(a, b) ((a) < (b) ? (b) : (a))", "/*#define max(a, b) ((a) < (b) ? (b) : (a))*/"],
      [67, "#define min(a, b) ((a) > (b) ? (b) : (a))", "/*#define min(a, b) ((a) > (b) ? (b) : (a))*/"]
    ])
    # patch "rcb/inertial3d.c"
    self.patch_file(os.path.join("src", "rcb", "inertial3d.c"), [
      [54, "#define max(a, b) ((a) < (b) ? (b) : (a))", "/*#define max(a, b) ((a) < (b) ? (b) : (a))*/"],
      [55, "#define min(a, b) ((a) > (b) ? (b) : (a))", "/*#define min(a, b) ((a) > (b) ? (b) : (a))*/"]
    ])
    # patch "zz/zz_const.c"
    self.patch_file(os.path.join("src", "zz", "zz_const.h"), [
      [55, "#include <strings.h>", "/*#include <strings.h>*/"]
    ])
    # patch "zz/zz_util.c"
    self.patch_file(os.path.join("src", "zz", "zz_util.c"), [
      [55, "#include <unistd.h>", "/*#include <unistd.h>*/"]
    ])

  def build(self):
    if not self.prepare():
      return
    print("\nBuilding " + self.name + " " + self.version + " sources...")
    self.compile_c("src/all/all_allo")
    self.compile_c("src/coloring/coloring")
    self.compile_c("src/coloring/color_test")
    self.compile_c("src/coloring/bucket")
    self.compile_c("src/coloring/g2l_hash")
    self.compile_c("src/graph/graph")
    self.compile_c("src/ha/divide_machine")
    self.compile_c("src/ha/get_processor_name")
    self.compile_c("src/ha/ha_ovis")
    self.compile_c("src/hier/hier")
    self.compile_c("src/hier/hier_free_struct")
    self.compile_c("src/hsfc/hsfc_box_assign")
    self.compile_c("src/hsfc/hsfc")
    self.compile_c("src/hsfc/hsfc_hilbert")
    self.compile_c("src/hsfc/hsfc_point_assign")
    self.compile_c("src/lb/lb_balance")
    self.compile_c("src/lb/lb_box_assign")
    self.compile_c("src/lb/lb_copy")
    self.compile_c("src/lb/lb_eval")
    self.compile_c("src/lb/lb_free")
    self.compile_c("src/lb/lb_init")
    self.compile_c("src/lb/lb_invert")
    self.compile_c("src/lb/lb_migrate")
    self.compile_c("src/lb/lb_part2proc")
    self.compile_c("src/lb/lb_point_assign")
    self.compile_c("src/lb/lb_remap")
    self.compile_c("src/lb/lb_set_fn")
    self.compile_c("src/lb/lb_set_method")
    self.compile_c("src/lb/lb_set_part_sizes")
    self.compile_c("src/matrix/matrix_build")
    self.compile_c("src/matrix/matrix_distribute")
    self.compile_c("src/matrix/matrix_operations")
    self.compile_c("src/matrix/matrix_sym")
    self.compile_c("src/matrix/matrix_utils")
    self.compile_c("src/order/order")
    self.compile_c("src/order/order_struct")
    self.compile_c("src/order/order_tools")
    self.compile_c("src/order/hsfcOrder")
    self.compile_c("src/order/perm")
    self.compile_c("src/par/par_average")
    self.compile_c("src/par/par_bisect")
    self.compile_c("src/par/par_median")
    self.compile_c("src/par/par_median_randomized")
    self.compile_c("src/par/par_stats")
    self.compile_c("src/par/par_sync")
    self.compile_c("src/par/par_tflops_special")
    self.compile_c("src/params/assign_param_vals")
    self.compile_c("src/params/bind_param")
    self.compile_c("src/params/check_param")
    self.compile_c("src/params/free_params")
    self.compile_c("src/params/key_params")
    self.compile_c("src/params/print_params")
    self.compile_c("src/params/set_param")
    self.compile_c("src/tpls/build_graph")
    self.compile_c("src/tpls/postprocessing")
    self.compile_c("src/tpls/preprocessing")
    self.compile_c("src/tpls/scatter_graph")
    self.compile_c("src/tpls/third_library")
    self.compile_c("src/tpls/verify_graph")
    self.compile_c("src/phg/phg_build")
    self.compile_c("src/phg/phg_build_calls")
    self.compile_c("src/phg/phg")
    self.compile_c("src/phg/phg_lookup")
    self.compile_c("src/phg/phg_verbose")
    self.compile_c("src/phg/phg_coarse")
    self.compile_c("src/phg/phg_comm")
    self.compile_c("src/phg/phg_distrib")
    self.compile_c("src/phg/phg_gather")
    self.compile_c("src/phg/phg_hypergraph")
    self.compile_c("src/phg/phg_match")
    self.compile_c("src/phg/phg_order")
    self.compile_c("src/phg/phg_parkway")
    self.compile_c("src/phg/phg_patoh")
    self.compile_c("src/phg/phg_plot")
    self.compile_c("src/phg/phg_rdivide")
    self.compile_c("src/phg/phg_refinement")
    self.compile_c("src/phg/phg_scale")
    self.compile_c("src/phg/phg_serialpartition")
    self.compile_c("src/phg/phg_util")
    self.compile_c("src/phg/phg_tree")
    self.compile_c("src/phg/phg_Vcycle")
    self.compile_c("src/rcb/box_assign")
    self.compile_c("src/rcb/create_proc_list")
    self.compile_c("src/rcb/inertial1d")
    self.compile_c("src/rcb/inertial2d")
    self.compile_c("src/rcb/inertial3d")
    self.compile_c("src/rcb/point_assign")
    self.compile_c("src/rcb/rcb_box")
    self.compile_c("src/rcb/rcb")
    self.compile_c("src/rcb/rcb_util")
    self.compile_c("src/rcb/rib")
    self.compile_c("src/rcb/rib_util")
    self.compile_c("src/rcb/shared")
    self.compile_c("src/reftree/reftree_build")
    self.compile_c("src/reftree/reftree_coarse_path")
    self.compile_c("src/reftree/reftree_hash")
    self.compile_c("src/reftree/reftree_part")
    self.compile_c("src/simple/block")
    self.compile_c("src/simple/cyclic")
    self.compile_c("src/simple/random")
    self.compile_c("src/timer/timer_params")
    self.compile_c("src/Utilities/Communication/comm_exchange_sizes")
    self.compile_c("src/Utilities/Communication/comm_invert_map")
    self.compile_c("src/Utilities/Communication/comm_do")
    self.compile_c("src/Utilities/Communication/comm_do_reverse")
    self.compile_c("src/Utilities/Communication/comm_info")
    self.compile_c("src/Utilities/Communication/comm_create")
    self.compile_c("src/Utilities/Communication/comm_resize")
    self.compile_c("src/Utilities/Communication/comm_sort_ints")
    self.compile_c("src/Utilities/Communication/comm_destroy")
    self.compile_c("src/Utilities/Communication/comm_invert_plan")
    self.compile_c("src/Utilities/Timer/zoltan_timer")
    self.compile_c("src/Utilities/Timer/timer")
    self.compile_c("src/Utilities/DDirectory/DD_Memory")
    self.compile_c("src/Utilities/DDirectory/DD_Find")
    self.compile_c("src/Utilities/DDirectory/DD_Destroy")
    self.compile_c("src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn3")
    self.compile_c("src/Utilities/DDirectory/DD_Remove")
    self.compile_c("src/Utilities/DDirectory/DD_Create")
    self.compile_c("src/Utilities/DDirectory/DD_Update")
    self.compile_c("src/Utilities/DDirectory/DD_Stats")
    self.compile_c("src/Utilities/DDirectory/DD_Hash2")
    self.compile_c("src/Utilities/DDirectory/DD_Print")
    self.compile_c("src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn2")
    self.compile_c("src/Utilities/DDirectory/DD_Set_Hash_Fn")
    self.compile_c("src/Utilities/DDirectory/DD_Set_Neighbor_Hash_Fn1")
    self.compile_c("src/Utilities/Memory/mem")
    self.compile_c("src/Utilities/shared/zoltan_align")
    self.compile_c("src/Utilities/shared/zoltan_id")
    self.compile_c("src/zz/zz_coord")
    self.compile_c("src/zz/zz_gen_files")
    self.compile_c("src/zz/zz_hash")
    self.compile_c("src/zz/murmur3")
    self.compile_c("src/zz/zz_map")
    self.compile_c("src/zz/zz_heap")
    self.compile_c("src/zz/zz_init")
    self.compile_c("src/zz/zz_obj_list")
    self.compile_c("src/zz/zz_rand")
    self.compile_c("src/zz/zz_set_fn")
    self.compile_c("src/zz/zz_sort")
    self.compile_c("src/zz/zz_struct")
    self.compile_c("src/zz/zz_back_trace")
    self.compile_c("src/zz/zz_util")
    self.link()

packages["zoltan"] = ThirdPartyZoltan()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

def split_cmd(line):
  ll = line.strip().lower().split()
  if len(ll) < 1:
    return None, None
  cmd = ll[0]
  # all packages?
  if "all" in ll[1:]:
    return cmd, list(packages.keys())
  # check pacakage names
  for pn in ll[1:]:
    if not pn in packages:
      print("ERROR: unknown library name '%s'" % pn)
      return None, None
  # okay
  return cmd, ll[1:]

while(True):
  print("="*100)
  print("\nAvailable/Installed Packages:")
  n = 0
  for pn in packages:
    pack = packages[pn]
    n += 1
    b = 'X' if pack.ready() else ' '
    print("%2i [%s] %s [v%s]" % (n, b, pack.name, pack.version))

  print("\nType 'help' for a list of all available commands, or 'quit' or 'exit' or 'stop' to quit/exit/stop.")
  print("By installing a third-party library, you explicitly agree to its license, which can be displayed")
  print("by typing 'lic <name>' for a specific library or by typing 'lic all' for all available libraries.")
  print("\n> ", end="", flush=True)

  # read command lines
  for line in sys.stdin:
    cmd, args = split_cmd(line)
    if cmd == None:
      print("\n> ", end="", flush=True)
      continue
    ### quit/exit/stop
    elif (cmd == "quit"):
      print("Quitting...")
      sys.exit(0)
    elif (cmd == "exit"):
      print("Exiting...")
      sys.exit(0)
    elif (cmd == "stop"):
      print("Stopping...")
      sys.exit(0)
    ### help
    elif cmd == "help":
      print("Available commands:")
      #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
      print("help                          Displays this message.")
      print("info <libs...>                Prints information about a set of libraries.")
      print("lic <libs...>                 Prints the licenses of a set of libraries.")
      print("install <libs...>             Installs a set of libraries.")
      print("remove <libs...>              Uninstalls a set of libraries.")
      print("quit                          Quits this program.")
      print("exit                          Exits this program.")
      print("stop                          Stops this program.")
      print("\nYou can also specify 'all' as an argument to execute a command for all available libraries.")
      input("\nPress Enter to continue...")
      break
    ### info
    elif cmd == 'info':
      if len(args) > 0:
        for pn in args:
          print("-"*100)
          packages[pn].info()
        input("\nPress Enter to continue...")
        break
      print("Missing arguments: 'info all' or 'info <libs...>' for specific libraries")
    ### lic / license
    elif cmd.startswith('lic'): # license
      if len(args) > 0:
        for pn in args:
          print("-"*100)
          packages[pn].license()
        input("\nPress Enter to continue...")
        break
      print("Missing arguments: 'lic all' or 'lic <libs...>' for specific libraries")
    ### install
    elif cmd == "install":
      if len(args) > 0:
        for pn in args:
          print("-"*100)
          packages[pn].install()
        break
      print("Missing arguments: 'install all' or 'install <libs...>' for specific libraries")
    ### install
    elif cmd == "remove":
      if len(args) > 0:
        for pn in args:
          print("-"*100)
          packages[pn].remove()
        break
      print("Missing arguments: 'remove all' or 'remove <libs...>' for specific libraries")
    # other fallbacks
    elif cmd == "make":
      print("The command to install packages is 'install', not 'make'")
    elif cmd == "rm":
      print("The command to uninstall packages is 'remove', not 'rm'")
    # unknown command
    else:
      print("ERROR: unknown command '%s'" % cmd)
    print("\n> ", end="", flush=True)
