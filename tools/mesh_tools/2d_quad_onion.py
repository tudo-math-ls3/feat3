#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
################################################################################
# This tool generates a 2D rectangular mesh representing the domain
# [x0,x1]x[y0,y1] consisting of one rectangle and 4*n surrounding trapezoidal
# quads forming an "onion":
#
# For n=2:   +-------------------+
#            |'\.             ./'|
#            |   +-----------+   |
#            |   |'\.     ./'|   |
#            |   |   +---+   |   |
#            |   |   |   |   |   |
#            |   |   +---+   |   |
#            |   |./'     '\.|   |
#            |   +-----------+   |
#            |./'             '\.|
#            +-------------------+
#
# The resulting mesh has a aspect ratio depending on the scaling factor 's',
# but it is independent of 'n'.
#
# Note: This script is compatible with both Python 2.x and Python 3.x
#
# \author Peter Zajc
#
import os.path
import sys
import math

################################################################################
################################################################################
################################################################################

# get the script filename
myself = os.path.basename(sys.argv[0])

# do we have enough arguments?
if len(sys.argv) < 8:
  print("USAGE: " + myself + " x0 x1 y0 y1 s n filename")
  print("")
  print("Options:")
  print("x0 x1     The X-range of the rectangular domain")
  print("y0 y1     The Y-range of the rectangular domain")
  print("s         The scaling factor for onion layers")
  print("n         The number of onion layers around the inner quad")
  print("filename  The name of the meshfile to be created")
  print("")
  print("Information:")
  print("This tool generates a 2D rectangular mesh representing the domain")
  print("[x0,x1]x[y0,y1] consisting of  one rectangle and 4*n surrounding")
  print("trapezoidal quads forming an \"onion\".")
  print("The generated meshfile also contains a chart for the outer boundary")
  print("named 'bnd' as well as four meshparts representing the four outer")
  print("edges of the domain named:")
  print("bnd:l     The left boundary edge meshpart")
  print("bnd:r     The right boundary edge meshpart")
  print("bnd:t     The top boundary edge meshpart")
  print("bnd:b     The bottom boundary edge meshpart")
  sys.exit(0)

# parse input arguments
x0 = float(sys.argv[1])
x1 = float(sys.argv[2])
y0 = float(sys.argv[3])
y1 = float(sys.argv[4])
s  = float(sys.argv[5])
n  = int(sys.argv[6])
filename = sys.argv[7]

# some basic sanity checks
if(x0 >= x1):
  print("ERROR: invalid X-range: x0 must be < x1")
  sys.exit(1)
if(y0 >= y1):
  print("ERROR: invalid Y-range: y0 must be < y1")
  sys.exit(1)
if (s < 0.0001) or (s > 0.9999):
  print("ERROR: 's' must be in range (0, 1)")
  sys.exit(1)
if(n < 0):
  print("ERROR: 'n' must be non-negative")
  sys.exit(1)

# compute vertex, edge and quad counts
nv = 4 + 4*n
ne = 4 + 8*n
nq = 1 + 4*n

# print some basic information
print("Domain: [%g , %g] x [%g , %g]" % (x0, x1, y0, y1))
print("Layers: %i" % (n))
print("Scale.: %g" % (s))
print("Verts.: %i" % nv)
print("Edges.: %i" % ne)
print("Quads.: %i" % nq)
print("")

# try to create the output file
print("Creating '" + filename + "'...")
f = open(filename, "wt")
# write header
f.write("<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:2:2\">\n")
f.write("  <!-- Generated by the FEAT '" + myself + "' tool -->\n")
# write tool call as comment
f.write("  <!-- call: " + myself)
for i in range(1,8):
  f.write(" " + sys.argv[i])
f.write(" -->\n")
# write chart
f.write("  <Chart name=\"bnd\">\n")
f.write("    <Bezier dim=\"2\" size=\"5\" type=\"closed\">\n")
f.write("      <Points>\n")
f.write("        0 %g %g\n" % (x0, y0))
f.write("        0 %g %g\n" % (x1, y0))
f.write("        0 %g %g\n" % (x1, y1))
f.write("        0 %g %g\n" % (x0, y1))
f.write("        0 %g %g\n" % (x0, y0))
f.write("      </Points>\n")
f.write("      <Params>\n")
f.write("        0\n")
f.write("        1\n")
f.write("        2\n")
f.write("        3\n")
f.write("        4\n")
f.write("      </Params>\n")
f.write("    </Bezier>\n")
f.write("  </Chart>\n")
# write mesh
f.write("  <Mesh type=\"conformal:hypercube:2:2\" size=\"%i %i %i\">\n" % (nv, ne, nq))
f.write("    <Vertices>\n")
# write vertices in counter-clock-wise order
for i in range(0,n+1):
  s1 = 0.5*(1.0 + math.pow(s, i))
  s0 = 1.0 - s1
  f.write("      %g %g\n" % (x0+s0*(x1-x0), y0+s0*(y1-y0)))
  f.write("      %g %g\n" % (x0+s1*(x1-x0), y0+s0*(y1-y0)))
  f.write("      %g %g\n" % (x0+s1*(x1-x0), y0+s1*(y1-y0)))
  f.write("      %g %g\n" % (x0+s0*(x1-x0), y0+s1*(y1-y0)))
f.write("    </Vertices>\n")
# write boundary edges
f.write("    <Topology dim=\"1\">\n")
f.write("      0 1\n")
f.write("      1 2\n")
f.write("      2 3\n")
f.write("      3 0\n")
# write inner edges
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
# write onion layer quads
for i in range(0,n):
  f.write("      %i %i %i %i\n" % (4*i+0, 4*i+1, 4*i+4, 4*i+5))
  f.write("      %i %i %i %i\n" % (4*i+1, 4*i+2, 4*i+5, 4*i+6))
  f.write("      %i %i %i %i\n" % (4*i+2, 4*i+3, 4*i+6, 4*i+7))
  f.write("      %i %i %i %i\n" % (4*i+3, 4*i+0, 4*i+7, 4*i+4))
# write inner quad
f.write("      %i %i %i %i\n" % (4*n+0, 4*n+1, 4*n+3, 4*n+2))
f.write("    </Topology>\n")
f.write("  </Mesh>\n")
# write mesh parts
part_names = ["bnd:b", "bnd:r", "bnd:t", "bnd:l"]
for i in range(0,4):
  f.write("  <MeshPart name=\"" + part_names[i] + "\" parent=\"root\" chart=\"bnd\" topology=\"full\" size=\"2 1\">\n")
  f.write("    <Mapping dim=\"0\">\n")
  f.write("      %i\n" % (i))
  f.write("      %i\n" % ((i+1) % 4))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"1\">\n")
  f.write("      %i\n" % (i))
  f.write("    </Mapping>\n")
  f.write("    <Topology dim=\"1\">\n")
  f.write("      0 1\n")
  f.write("    </Topology>\n")
  f.write("    <Attribute dim=\"1\" name=\"param\">\n")
  f.write("      %i\n" % (i))
  f.write("      %i\n" % (i+1))
  f.write("    </Attribute>\n")
  f.write("  </MeshPart>\n");
f.write("</FeatMeshFile>\n")

print("Done!")