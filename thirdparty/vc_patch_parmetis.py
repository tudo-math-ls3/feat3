#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
################################################################################
# ParMETIS 4.0.3 patch for Visual Studio
# ------------------------------------------------------------------------------
# This script patches two of the GKlib files to ensure that the library
# can be compiled under Visual Studio 14.
#
# \author Peter Zajac
################################################################################
import os
import sys

# set GKlib paths
gk_path = os.path.join(".","parmetis","metis","GKlib")
gk_path1o = os.path.join(gk_path,"gk_arch.h")
gk_path1i = os.path.join(gk_path,"gk_arch.h.backup")
gk_path2o = os.path.join(gk_path,"gkregex.c")
gk_path2i = os.path.join(gk_path,"gkregex.c.backup")

# check whether ParMETIS is undpacked
if not os.path.isfile(os.path.join(".","parmetis","include","parmetis.h")):
  print("ERROR: ParMETIS source not found; nothing to patch...")
  sys.exit(1)

# create backup files unless they already exist
if not os.path.isfile(gk_path1i):
  os.rename(gk_path1o, gk_path1i)
if not os.path.isfile(gk_path2i):
  os.rename(gk_path2o, gk_path2i)

################################################################################
##### patch 'gk_arch.h' #####
print("Patching '%s'..." % gk_path1o)

# open backup file
fi = open(gk_path1i, "rt")

# open output file
fo = open(gk_path1o, "wt")

# loop over all lines
lno = 0
for line in fi:
  lno = lno + 1
  # insert macro before line 48
  if (lno == 48):
    fo.write("#if defined(_WIN32) && !defined(WIN32)\n")
    fo.write("  #define WIN32\n")
    fo.write("#endif\n")
  # remove lines 62-70
  if (lno > 61) and (lno < 69):
    continue
  # insert line
  if (lno == 69):
    fo.write("#define __thread __declspec(thread)\n")
  # okay, write line
  fo.write(line)

fo.close()
fi.close()

################################################################################
##### patch 'gkregex.c' #####
print("Patching '%s'..." % gk_path2o)

# open backup file
fi = open(gk_path2i, "rt")

# open output file
fo = open(gk_path2o, "wt")

# loop over all lines
lno = 0
for line in fi:
  lno = lno + 1
  # replace line 5089
  if (lno == 5089):
    fo.write("    postorder (elem, mark_opt_subexp, (void *) (intptr_t) elem->token.opr.idx);\n")
    continue
  # replace line 6301
  if (lno == 6301):
    fo.write("  int idx = (int) (intptr_t) extra;\n")
    continue
  # okay, write line
  fo.write(line)

fo.close()
fi.close()

# okay
print("Patch applied successfully")
