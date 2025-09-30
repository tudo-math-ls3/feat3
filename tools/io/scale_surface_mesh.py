#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
################################################################################
# This tool scales a 3D surface mesh, represented for now by an .off file,
# by a given factor and saves it into a new file.
#
# USAGE: scale_surface_mesh.py filename scale [outfile]
#
# Options:
# filename   The filename of the surface mesh. Checked if it ends with .off .
# scale      The scaling factor to be applied. Should be positive.
# outfile    Optional filename for the output mesh. If not provided
#            use the filename_scaled.off .
#
# Note: This script is compatible with both Python 2.x and Python 3.x
#
# \author Maximilian Esser
#
import os.path
import sys


# get the script filename
myself = os.path.basename(sys.argv[0])

# do we have enough arguments?
if len(sys.argv) < 3:
  print("")
  print("USAGE: " + myself + " filename scale [outfile]")
  print("")
  print("Options:")
  print("filename  The (relative) filepath to the surface mesh")
  print("scale     The scaling factor to be applied. Should be greater 0")
  print("outfile   Optional: The (relative) filepath to the output file")
  print("")
  print("Information:")
  print("This tool scales a 3D surface mesh, represented for now by an .off file")
  print("by a given factor and saves it into a new file.")
  sys.exit(0)

# parse input arguments
filename = sys.argv[1]
scale = float(sys.argv[2])
outfile = ""
if(len(sys.argv) > 3):
  outfile = sys.argv[3]

# some basic sanity checks
if(scale <= 0.):
  print("ERROR: invalid scaling factor. Should be larger 0")
  sys.exit(1)

split_filename = filename.split(".")
if(len(split_filename) != 2 or split_filename[-1] != "off"):
  print("ERROR: invalid file")
  sys.exit(1)
if(not outfile):
  outfile = split_filename[0] + "_scaled.off"

with open(outfile, "wt") as out_f:

  with open(filename, "rt") as in_f:
    tmp = in_f.readline()
    if(tmp != "OFF\n"):
      print("ERROR: " + filename + " is not an actual off file.")
      sys.exit()
    tmp = in_f.readline()
    out_f.write("OFF\n")
    out_f.write(tmp)
    for line in in_f:
      split_line = line.split()
      if(len(split_line) != 3):
        out_f.write(line)
        break
      x = scale * float(split_line[0])
      y = scale * float(split_line[1])
      z = scale * float(split_line[2])
      out_f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

    for line in in_f:
      out_f.write(line)


print("Done!")
