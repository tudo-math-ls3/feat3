#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
########################################################################################################################
# This tool generates a 2D circle mesh.
# The generated meshfile also contains a chart as well as a meshpart named 'bnd:outer' representing the outer boundary.
#
# The mesh that is created by this script consists of a single square quadilateral that is rotated by 45 degrees in
# the center of the domain and which is enclosed by N layers of 4 quadrilaterals, where N+1 is the number of radii
# given as arguments to the script. The radius R0 is the radius of the inner square element, whereas each additional
# 4-element layer #i has the given radius Ri.
#
#           .
#        ./'|'\.
#      ./'  |  '\.
#    ./'  ./'\   '\.
#  ./'  ./'   '\.  '\.
# +----+         +----+
#  '\.  '\.   ./'  ./'
#    '\.  '\./'  ./'
#      '\   |  ./'
#        '\.|./'
#           '
#
# USAGE: 2d_quad_circle cx cy R0 R1 ... RN filename
#
# Options:
# cx cy        The coordinates of the circle barycentre
# R0 R1... RN  A list of at least 2 increasing radii of the circle segments,
#              where R0 is the radius of the inner square element
# filename     The name of the meshfile to be created
#
# Note: This script is compatible with both Python 2.x and Python 3.x
#
# \author Peter Zajc
#
import os.path
import sys
import math

########################################################################################################################
########################################################################################################################
########################################################################################################################

# get the script filename
myself = os.path.basename(sys.argv[0])

# do we have enough arguments?
if len(sys.argv) < 5:
  print("USAGE: " + myself + " cx cy R0 R1... RN filename")
  print("")
  sys.exit(0)

# parse center coords
cx = float(sys.argv[1])
cy = float(sys.argv[2])

# radii array
radii = []

# parse input arguments
last_r = 0.0
for i in range(3, len(sys.argv)-1):
  r = float(sys.argv[i])
  if not r > last_r:
    print("ERROR: invalid radius: %f is not bigger than previous radius %f" % (r, last_r))
    sys.exit(1)
  radii += [r]
filename = sys.argv[-1]

# print radii
n = len(radii)-1
if n < 1:
  print("ERROR: you must specify at least 2 radii")
  sys.exit(1)

# compute vertex, edge and quad counts
nv = 4 + 4*n
ne = 4 + 8*n
nq = 1 + 4*n

# print some basic information
print("Radii: " + (", ".join(["%g" % r for r in radii])))
print("Verts: %i" % nv)
print("Edges: %i" % ne)
print("Quads: %i" % nq)
print("")

# try to create the output file
print("Creating '" + filename + "'...")
f = open(filename, "wt")
# write header
f.write("<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:2:2\">\n")
# write outer circle chart
f.write("<Chart name=\"outer\">\n")
f.write("  <Circle radius=\"%g\" midpoint=\"%g %g\" domain=\"0 4\" />\n" % (radii[-1], cx, cy))
f.write("</Chart>\n")
# write mesh
f.write("  <Mesh type=\"conformal:hypercube:2:2\" size=\"%i %i %i\">\n" % (nv, ne, nq))
f.write("    <Vertices>\n")
# write vertices in counter-clock-wise order
for i in range(0,n+1):
  r = radii[i]
  f.write("      %g %g\n" % (cx+r, cy  ))
  f.write("      %g %g\n" % (cx  , cy+r))
  f.write("      %g %g\n" % (cx-r, cy  ))
  f.write("      %g %g\n" % (cx  , cy-r))
f.write("    </Vertices>\n")
# write edges
f.write("    <Topology dim=\"1\">\n")
f.write("      0 1\n")
f.write("      1 2\n")
f.write("      2 3\n")
f.write("      3 0\n")
for i in range(0,n):
  f.write("      %i %i\n" % (4*i+0, 4*i+4))
  f.write("      %i %i\n" % (4*i+1, 4*i+5))
  f.write("      %i %i\n" % (4*i+2, 4*i+6))
  f.write("      %i %i\n" % (4*i+3, 4*i+7))
  f.write("      %i %i\n" % (4*i+4, 4*i+5))
  f.write("      %i %i\n" % (4*i+5, 4*i+6))
  f.write("      %i %i\n" % (4*i+6, 4*i+7))
  f.write("      %i %i\n" % (4*i+7, 4*i+4))
f.write("    </Topology>\n")
f.write("    <Topology dim=\"2\">\n")
f.write("      0 1 3 2\n")
# write onion layer quads
for i in range(0,n):
  f.write("      %i %i %i %i\n" % (4*i+0, 4*i+1, 4*i+4, 4*i+5))
  f.write("      %i %i %i %i\n" % (4*i+1, 4*i+2, 4*i+5, 4*i+6))
  f.write("      %i %i %i %i\n" % (4*i+2, 4*i+3, 4*i+6, 4*i+7))
  f.write("      %i %i %i %i\n" % (4*i+3, 4*i+0, 4*i+7, 4*i+4))
# write inner quad
f.write("    </Topology>\n")
f.write("  </Mesh>\n")
# write outer mesh part
f.write("  <MeshPart name=\"bnd:outer\" parent=\"root\" chart=\"outer\" topology=\"full\" size=\"5 4\">\n")
f.write("    <Mapping dim=\"0\">\n")
for i in range(0,5):
  f.write("      %i\n" % (4*n+ i%4))
f.write("    </Mapping>\n")
f.write("    <Mapping dim=\"1\">\n")
for i in range(0,4):
  f.write("      %i\n" % (8*n + i))
f.write("    </Mapping>\n")
f.write("    <Topology dim=\"1\">\n")
f.write("      0 1\n")
f.write("      1 2\n")
f.write("      2 3\n")
f.write("      3 4\n")
f.write("    </Topology>\n")
f.write("    <Attribute name=\"param\" dim=\"1\">\n")
f.write("      0\n")
f.write("      1\n")
f.write("      2\n")
f.write("      3\n")
f.write("      4\n")
f.write("    </Attribute>\n")
f.write("  </MeshPart>\n")
f.write("</FeatMeshFile>\n")
print("Done!")
