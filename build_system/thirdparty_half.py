# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "Feb 2016"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class HALF(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.version = "2.1.0"
    self.names = ["half"]
    self.dirname = "half-" + self.version
    self.filename = "half-" + self.version + ".zip"
    self.url = "http://downloads.sourceforge.net/project/half/half/" + self.version + "/" + self.filename
    self.cmake_flags = " -DFEAT_HAVE_HALFMATH:BOOL=ON -DFEAT_DIRNAME_HALF:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = os.path.join(trunk_dirname, self.dirname)
