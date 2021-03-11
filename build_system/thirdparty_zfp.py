# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock, Peter Zajac"
__date__   = "January 2020"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class ZFP(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    # Note: the ZFP authors are extremely lazy with new releases, so we have download a specific
    # commit and not a release, because 0.5.5 does not compile and 0.5.6 is not released yet.
    # Change this to version  0.5.6 once it has been released.
    self.version = "08adb27dbc7ed91373ba76312aa7e48504a02190"
    self.names = ["zfp"]
    self.dirname = "zfp-" + self.version
    self.filename ="zfp-" + self.version + ".zip"
    self.url ="https://github.com/LLNL/zfp/archive/" + self.version + ".zip"
    self.cmake_flags = " -DFEAT_HAVE_ZFP:BOOL=ON -DFEAT_DIRNAME_ZFP:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = trunk_dirname
