#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
################################################################################
# This tool generates a 3D cuboid mesh representing the domain
# [x0,x1]x[y0,y1]x[z0,z1] consisting of m-by-n-by-l congruent cuboids.
# The generated meshfile also contains six meshparts representing the
# six outer faces of the domain named:
# bnd:l     The left boundary face meshpart   (x=x0)
# bnd:r     The right boundary face meshpart  (x=x1)
# bnd:b     The bottom boundary face meshpart (y=y0)
# bnd:t     The top boundary face meshpart    (y=y1)
# bnd:f     The far boundary face meshpart    (z=z0)
# bnd:n     The near boundary face meshpart   (z=z1)
#
# USAGE: 3d_hexa_cube.py x0 x1 y0 y1 z0 z1 m n l filename
#
# Options:
# x0 x1     The X-range of the cuboid domain
# y0 y1     The Y-range of the cuboid domain
# z0 z1     The Z-range of the cuboid domain
# m         The number of hexahedra in X-direction
# n         The number of hexahedra in Y-direction
# l         The number of hexahedra in Z-direction
# filename  The name of the meshfile to be created
#
# Note: This script is compatible with both Python 2.x and Python 3.x
#
# \author Peter Zajac
#
import os.path
import sys

# generate YZ-plane meshpart (left, right)
def mp_yz(f, name, j, m, n, l):
  f.write("  <MeshPart name=\"" + name + "\" parent=\"root\" topology=\"none\" size=\"%i %i %i\">\n" % ((l+1)*(n+1), l*(n+1)+(l+1)*n, l*n))
  f.write("    <Mapping dim=\"0\">\n")
  for k in range(0, l+1):
    for i in range(0, n+1):
      f.write("      %i\n" % ((k*(n+1)+i)*(m+1)+j))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"1\">\n")
  for k in range(0, l+1):
    for i in range(0, n):
      f.write("      %i\n" % ((l+1)*(n+1)*m + (j*(l+1)+k)*n + i))
  for i in range(0, n+1):
    for k in range(0, l):
      f.write("      %i\n" % ((l+1)*((n+1)*m+(m+1)*n) + (i*(m+1)+j)*l + k))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"2\">\n")
  for k in range(0, l):
    for i in range(0, n):
      f.write("      %i\n" % ((l+1)*m*n + (j*l+k)*n + i))
  f.write("    </Mapping>\n")
  f.write("  </MeshPart>\n")

# generate XZ-plane meshpart (bottom, top)
def mp_xz(f, name, i, m, n, l):
  f.write("  <MeshPart name=\"" + name + "\" parent=\"root\" topology=\"none\" size=\"%i %i %i\">\n" % ((m+1)*(l+1), m*(l+1)+(m+1)*l, m*l))
  f.write("    <Mapping dim=\"0\">\n")
  for j in range(0, m+1):
    for k in range(0, l+1):
      f.write("      %i\n" % ((k*(n+1)+i)*(m+1)+j))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"1\">\n")
  for k in range(0, l+1):
    for j in range(0, m):
      f.write("      %i\n" % ((k*(n+1)+i)*m+j))
  for j in range(0, m+1):
    for k in range(0, l):
      f.write("      %i\n" % ((l+1)*((n+1)*m+(m+1)*n) + (i*(m+1)+j)*l + k))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"2\">\n")
  for j in range(0, m):
    for k in range(0, l):
      f.write("      %i\n" % (((l+1)*m+(m+1)*l)*n + (i*m+j)*l + k))
  f.write("    </Mapping>\n")
  f.write("  </MeshPart>\n")
  return

# generate XY-plane meshpart (near, far)
def mp_xy(f, name, k, m, n, l):
  f.write("  <MeshPart name=\"" + name + "\" parent=\"root\" topology=\"none\" size=\"%i %i %i\">\n" % ((m+1)*(n+1), m*(n+1)+(m+1)*n, m*n))
  f.write("    <Mapping dim=\"0\">\n")
  for i in range(0, n+1):
    for j in range(0, m+1):
      f.write("      %i\n" % ((k*(n+1)+i)*(m+1)+j))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"1\">\n")
  for i in range(0, n+1):
    for j in range(0, m):
      f.write("      %i\n" % ((k*(n+1)+i)*m+j))
  for j in range(0, m+1):
    for i in range(0, n):
      f.write("      %i\n" % ((l+1)*(n+1)*m + (j*(l+1)+k)*n + i))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"2\">\n")
  for i in range(0, n):
    for j in range(0, m):
      f.write("      %i\n" % ((k*n+i)*m+j))
  f.write("    </Mapping>\n")
  f.write("  </MeshPart>\n")


####################################################################################################
####################################################################################################
####################################################################################################

# get the script filename
myself = os.path.basename(sys.argv[0])

# do we have enough arguments?
if len(sys.argv) < 11:
  print("")
  print("USAGE: " + myself + " x0 x1 y0 y1 z0 z1 m n l filename")
  print("")
  print("Options:")
  print("x0 x1     The X-range of the cuboid domain")
  print("y0 y1     The Y-range of the cuboid domain")
  print("z0 z1     The Z-range of the cuboid domain")
  print("m         The number of hexahedra in X-direction")
  print("n         The number of hexahedra in Y-direction")
  print("l         The number of hexahedra in Z-direction")
  print("filename  The name of the meshfile to be created")
  print("")
  print("Information:")
  print("This tool generates a 3D cuboid mesh representing the domain")
  print("[x0,x1]x[y0,y1]x[z0,z1] consisting of m-by-n-by-l congruent hexahedra.")
  print("The generated meshfile also contains six meshparts representing the")
  print("six outer faces of the domain named:")
  print("bnd:l     The left boundary face meshpart   (x=x0)")
  print("bnd:r     The right boundary face meshpart  (x=x1)")
  print("bnd:b     The bottom boundary face meshpart (y=y0)")
  print("bnd:t     The top boundary face meshpart    (y=y1)")
  print("bnd:f     The far boundary face meshpart    (z=z0)")
  print("bnd:n     The near boundary face meshpart   (z=z1)")
  sys.exit(0)

# parse input arguments
x0 = float(sys.argv[1])
x1 = float(sys.argv[2])
y0 = float(sys.argv[3])
y1 = float(sys.argv[4])
z0 = float(sys.argv[5])
z1 = float(sys.argv[6])
m  = int(sys.argv[7])
n  = int(sys.argv[8])
l  = int(sys.argv[9])
filename = sys.argv[10]

# some basic sanity checks
if(x0 >= x1):
  print("ERROR: invalid X-range: x0 must be < x1")
  sys.exit(1)
if(y0 >= y1):
  print("ERROR: invalid Y-range: y0 must be < y1")
  sys.exit(1)
if(z0 >= z1):
  print("ERROR: invalid Z-range: z0 must be < z1")
  sys.exit(1)
if(m <= 0):
  print("ERROR: 'm' must be positive")
  sys.exit(1)
if(n <= 0):
  print("ERROR: 'n' must be positive")
  sys.exit(1)
if(l <= 0):
  print("ERROR: 'l' must be positive")
  sys.exit(1)

# compute vertex, edge, face and cell counts
nv = (m+1)*(n+1)*(l+1)
ne = m*(n+1)*(l+1) + (m+1)*n*(l+1) + (m+1)*(n+1)*l
nf = (m+1)*n*l + m*(n+1)*l + m*n*(l+1)
nc = m*n*l

# print some basic information
print("Domain: [%g , %g] x [%g , %g] x [%g , %g]" % (x0, x1, y0, y1, z0, z1))
print("Slices: %i x %i x %i" % (m, n, l))
print("Rects.: %g x %g x %g" % ((x1-x0)/float(m), (y1-y0)/float(n), (z1-z0)/float(l)))
print("Verts.: %i" % nv)
print("Edges.: %i" % ne)
print("Faces.: %i" % nf)
print("Hexas.: %i" % nc)
print("")

# try to create the output file
print("Creating '" + filename + "'...")
f = open(filename, "wt")
# write header
f.write("<FeatMeshFile version=\"1\" mesh=\"conformal:hypercube:3:3\">\n")
f.write("  <!-- Generated by the FEAT '" + myself + "' tool -->\n")
# write tool call as comment
f.write("  <!-- call: " + myself)
for i in range(1, 11):
  f.write(" " + sys.argv[i])
f.write(" -->\n")
# write mesh
f.write("  <Mesh type=\"conformal:hypercube:3:3\" size=\"%i %i %i %i\">\n" % (nv, ne, nf, nc))
f.write("    <Vertices>\n")
# write vertices in line-wise order
for k in range(0,l+1):
  z = z0 + (z1-z0)*(float(k)/float(l))
  for i in range(0,n+1):
    y = y0 + (y1-y0)*(float(i)/float(n))
    for j in range(0,m+1):
      x = x0 + (x1-x0)*(float(j)/float(m))
      f.write("      %g %g %g\n" % (x, y, z))
f.write("    </Vertices>\n")
f.write("    <Topology dim=\"1\">\n")
# write X-parallel edges
for k in range(0,l+1):
  for i in range(0,n+1):
    for j in range(0,m):
      f.write("      %i %i\n" % ((k*(n+1)+i)*(m+1)+j, (k*(n+1)+i)*(m+1)+j+1))
# write Y-parallel edges
for j in range(0,m+1):
  for k in range(0,l+1):
    for i in range(0,n):
      f.write("      %i %i\n" % ((k*(n+1)+i)*(m+1)+j, (k*(n+1)+i+1)*(m+1)+j))
# write Z-parallel edges
for i in range(0,n+1):
  for j in range(0,m+1):
    for k in range(0,l):
      f.write("      %i %i\n" % ((k*(n+1)+i)*(m+1)+j, ((k+1)*(n+1)+i)*(m+1)+j))
f.write("    </Topology>\n")
f.write("    <Topology dim=\"2\">\n")
# write XY-plane quads
for k in range(0,l+1):
  for i in range(0,n):
    for j in range(0,m):
      f.write("      %i %i %i %i\n" % ((k*(n+1)+i)*(m+1)+j, (k*(n+1)+i)*(m+1)+j+1, (k*(n+1)+i+1)*(m+1)+j, (k*(n+1)+i+1)*(m+1)+j+1))
# write YZ-plane quads
for j in range(0,m+1):
  for k in range(0,l):
    for i in range(0,n):
      f.write("      %i %i %i %i\n" % ((k*(n+1)+i)*(m+1)+j, (k*(n+1)+i+1)*(m+1)+j, ((k+1)*(n+1)+i)*(m+1)+j, ((k+1)*(n+1)+i+1)*(m+1)+j))
# write XZ-plane quads
for i in range(0,n+1):
  for j in range(0,m):
    for k in range(0,l):
      f.write("      %i %i %i %i\n" % ((k*(n+1)+i)*(m+1)+j, (k*(n+1)+i)*(m+1)+j+1, ((k+1)*(n+1)+i)*(m+1)+j, ((k+1)*(n+1)+i)*(m+1)+j+1))
f.write("    </Topology>\n")
# write hexahedra
f.write("    <Topology dim=\"3\">\n")
for k in range(0,l):
  for i in range(0,n):
    for j in range(0, m):
      f.write("      %i %i %i %i %i %i %i %i\n" % (
        (k*(n+1)+i)*(m+1)+j, (k*(n+1)+i)*(m+1)+j+1, (k*(n+1)+i+1)*(m+1)+j, (k*(n+1)+i+1)*(m+1)+j+1,
        ((k+1)*(n+1)+i)*(m+1)+j, ((k+1)*(n+1)+i)*(m+1)+j+1, ((k+1)*(n+1)+i+1)*(m+1)+j, ((k+1)*(n+1)+i+1)*(m+1)+j+1)
      )
f.write("    </Topology>\n")
f.write("  </Mesh>\n")
# write left YZ meshpart (x=x0)
mp_yz(f, "bnd:l", 0, m, n, l)
# write right YZ meshpart (x=x1)
mp_yz(f, "bnd:r", m, m, n, l)
# write bottom XZ meshpart (y=y0)
mp_xz(f, "bnd:b", 0, m, n, l)
# write top XZ meshpart (y=y1)
mp_xz(f, "bnd:t", n, m, n, l)
# write far XY-meshpart (z=z0)
mp_xy(f, "bnd:f", 0, m, n, l)
# write near XY-meshpart (z=z1)
mp_xy(f, "bnd:n", l, m, n, l)
f.write("</FeatMeshFile>\n")

print("Done!")
