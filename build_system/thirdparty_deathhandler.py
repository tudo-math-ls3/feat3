# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
__author__ = "Dirk Ribbrock"
__date__   = "2019-11-21"
from build_system.thirdparty_package import ThirdpartyPackage
import os

class FParser(ThirdpartyPackage):

  def __init__(self,trunk_dirname):
    self.names = ["deathhandler"]
    self.dirname = "deathhandler"
    self.filename = "DeathHandler.zip"
    self.url = "https://github.com/vmarkovtsev/DeathHandler/archive/master.zip"
    self.cmake_flags = " -DFEAT_HAVE_DEATH_HANDLER:BOOL=ON -DFEAT_DIRNAME_DEATHHANDLER:STRING='" + self.dirname + "'"
    self.trunk_dirname = trunk_dirname
    self.target_dirname = os.path.join(trunk_dirname, self.dirname)
