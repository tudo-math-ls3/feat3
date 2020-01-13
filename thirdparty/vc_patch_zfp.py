#!/usr/bin/env python
################################################################################
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
################################################################################
# ZFP 0.5.5 patch for Visual Studio
# ------------------------------------------------------------------------------
# This script patches a header file in the source code of the ZFP library
# to enable static library builds.
#
# Note:
# As of 2020-01-13, the latest stable version of zfp is 0.5.5, which needs this
# patch to successfully link the library. However, the source code that is
# hosted on GitHub is already fixed, so I expect that this patch will be
# obsolete from version 0.5.6 onwards.
#
# \author Peter Zajac
################################################################################
import os
import sys

# set header paths
path_o = os.path.join(".","zfp","include","zfp","system.h")
path_i = os.path.join(".","zfp","include","zfp","system.h.backup")

# check whether ZFP is unpacked
if not os.path.isfile(os.path.join(".","zfp","utils","zfp.c")):
  print("ERROR: ZFP source not found; nothing to patch...")
  sys.exit(1)

# create backup files unless they already exist
if not os.path.isfile(path_i):
  os.rename(path_o, path_i)

################################################################################
##### patch 'system.h' #####
print("Patching '%s'..." % path_o)

# open backup file
fi = open(path_i, "rt")

# open output file
fo = open(path_o, "wt")

# loop over all lines
lno = 0
for line in fi:
  lno = lno + 1
  # replace line 11
  if (lno == 11) and (line == "#ifdef _MSC_VER\n"):
    fo.write("#if 0 //#ifdef _MSC_VER\n")
    continue
  # okay, write line
  fo.write(line)

fo.close()
fi.close()

# okay
print("Patch applied successfully")
