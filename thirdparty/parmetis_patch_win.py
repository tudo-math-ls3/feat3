################################################################################
# ParMETIS 4.0.3 patch for Visual Studio
# ------------------------------------------------------------------------------
# This script patches one of the GKlib header files to ensure that the library
# can be compiled under Visual Studio 14.
#
# \author Peter Zajac
################################################################################
import os
import sys

# set GKlib paths
gk_path = os.path.join(".","parmetis","metis","GKlib")
gk_path1 = os.path.join(gk_path,"gk_arch.h")
gk_path2 = os.path.join(gk_path,"gk_arch.h.backup")

# check whether 'gk_arch.h' exists
if not os.path.isfile(gk_path1):
  print("ERROR: ParMETIS source not found; nothing to patch...")
  sys.exit(1)

# create backup file unless it already exists
if not os.path.isfile(gk_path2):
  os.rename(gk_path1, gk_path2)

# open backup file
fi = open(gk_path2, "rt")

# open output file
fo = open(gk_path1, "wt")

# loop over all lines
lno = 0
for line in fi:
  lno = lno + 1
  # insert macro before line 48
  if (lno == 48):
    fo.write("#ifdef _WIN32\n  #define WIN32\n#endif\n")
  # remove lines 62-70
  if (lno > 61) and (lno < 71):
    continue
  # insert line
  if (lno == 71):
    fo.write("#define __thread __declspec(thread)\n")
  # okay, write line
  fo.write(line)

fo.close()
fi.close()

# okay
print("Patch applied successfully")
