#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
################################################################################
# This script generates a 2D triangular mesh representing the domain
# [x0,x1]x[y0,y1] consisting of m-by-n congruent or anisotropic triangle pairs.
# The generated meshfile also contains four meshparts representing the
# four outer edges of the domain named:
# bnd:l     The left boundary edge meshpart   (x=x0)
# bnd:r     The right boundary edge meshpart  (x=x1)
# bnd:b     The bottom boundary edge meshpart (y=y0)
# bnd:t     The top boundary edge meshpart    (y=y1)
#
# USAGE: 2d_tria_rect.py x0 x1 y0 y1 [m-list] [n-list] filename
#
# Options:
# x0 x1     The X-range of the rectangular domain
# y0 y1     The Y-range of the rectangular domain
# [m-list]  The rectangle list in X-direction, see below
# [n-list]  The rectangle list in Y-direction, see below
# filename  The name of the meshfile to be created
#
# Rectangle Lists:
# For both the X- and the Y-dimension, one may either specify the number of
# triangle pairs in that dimension, which will yield an equidistant discretisation,
# or alternatively a list of relative element sizes, which will yield an
# anisotropic discretisation.
#
# Simple Equidistant Example:
# The call
#      2d_tria_rect.py 0 2 0 1 4 2 mesh.xml
# will create a mesh discretising the domain [0,2]x[0,1] with 4x2 triangle pairs:
#
# +-----+-----+-----+-----+
# |  ./'|  ./'|  ./'|  ./'|
# |./'  |./'  |./'  |./'  |
# +-----+-----+-----+-----+
# |  ./'|  ./'|  ./'|  ./'|
# |./'  |./'  |./'  |./'  |
# +-----+-----+-----+-----+
#
#
# Mixed Equidistant/Anisotropic Example:
# The call
#      2d_tria_rect.py 0 5 0 2 [1 3 1] 2 mesh.xml
# will create a mesh discretising the domain [0,5]x[0,2] with 3x2 triangle pairs,
# where each element has the same Y-dimension and where the horizontally inner
# elements are 3 times as big in X-dimension as the horizontally outer elements:
#
# +------+------------------+------+
# |   ./'|         ...---'''|   ./'|
# |./'   |...---'''         |./'   |
# +------+------------------+------+
# |   ./'|         ...---'''|   ./'|
# |./'   |...---'''         |./'   |
# +------+------------------+------+
#
#
# Note: This script is compatible with both Python 2.x and Python 3.x
#
# \author Peter Zajac
#
import os.path
import sys

# generate horizontal meshpart (left, right)
def mp_horz(f, name, i, m, n):
  f.write("  <MeshPart name=\"" + name + "\" parent=\"root\" topology=\"none\" size=\"%i %i\">\n" % (m+1,m))
  f.write("    <Mapping dim=\"0\">\n")
  for j in range(0, m+1):
    f.write("      %i\n" % (i*(m+1) + j))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"1\">\n")
  for j in range(0, m):
    f.write("      %i\n" % (i*m + j))
  f.write("    </Mapping>\n")
  f.write("  </MeshPart>\n")

# generate vertical meshpart (bottom, top)
def mp_vert(f, name, j, m, n):
  f.write("  <MeshPart name=\"" + name + "\" parent=\"root\" topology=\"none\" size=\"%i %i\">\n" % (n+1,n))
  f.write("    <Mapping dim=\"0\">\n")
  for i in range(0, n+1):
    f.write("      %i\n" % (i*(m+1) + j))
  f.write("    </Mapping>\n")
  f.write("    <Mapping dim=\"1\">\n")
  for i in range(0, n):
    f.write("      %i\n" % ((n+1)*m + j*n + i))
  f.write("    </Mapping>\n")
  f.write("  </MeshPart>\n")

def pop_front(args):
  return (args[0], args[1:])

# parse anisotropic parameter list
def parse_aniso_list(args, mn):
  # empty list?
  if len(args) == 0:
    print("ERROR: " + mn + "-list is empty")
    sys.exit(1)

  # isotropic case ?
  a, args = (args[0], args[1:])
  if not a.startswith('['):
    return ([], int(a), args)

  # it's an anisotropy list
  lst = a
  while len(args) > 0:
    a, args = (args[0], args[1:])
    lst = lst + " " + a
    if a.endswith(']'):
      break
  if not lst.endswith(']'):
    print("ERROR: " + mn + "-list is incomplete")
    sys.exit(1)

  # remove braces and split by whitespaces
  lst = lst[1:-1].strip().split()

  # convert from string to float list
  l = []
  for t in lst:
    l = l + [float(t)]

  return (l, len(l), args)

# print aniso list
def print_aniso_list(il, lst, sum):
  s = ""
  for i in lst:
    s = s + " %f " % (il * i / sum)
  return "[" + s + "]"

####################################################################################################
####################################################################################################
####################################################################################################

# get the script filename
myself = os.path.basename(sys.argv[0])

# do we have enough arguments?
if len(sys.argv) < 8:
  print("")
  print("USAGE: " + myself + " x0 x1 y0 y1 [m-list] [n-list] filename")
  print("")
  print("Options:")
  print("x0 x1     The X-range of the rectangular domain")
  print("y0 y1     The Y-range of the rectangular domain")
  print("[m-list]  The rectangle list in X-direction")
  print("[n-list]  The rectangle list in Y-direction")
  print("filename  The name of the meshfile to be created")
  print("")
  print("Information:")
  print("This tool generates a 2D triangular mesh representing the domain")
  print("[x0,x1]x[y0,y1] consisting of m-by-n congruent or anisotropic triangle pairs.")
  print("The generated meshfile also contains four meshparts representing the")
  print("four outer edges of the domain named:")
  print("bnd:l     The left boundary edge meshpart   (x=x0)")
  print("bnd:r     The right boundary edge meshpart  (x=x1)")
  print("bnd:b     The bottom boundary edge meshpart (y=y0)")
  print("bnd:t     The top boundary edge meshpart    (y=y1)")
  sys.exit(0)

# parse input arguments
x0 = float(sys.argv[1])
x1 = float(sys.argv[2])
y0 = float(sys.argv[3])
y1 = float(sys.argv[4])
filename = sys.argv[-1]

# parse anisotropy lists
mn_args = sys.argv[5:-1]
m_list, m, mn_args = parse_aniso_list(mn_args, 'm')
n_list, n, mn_args = parse_aniso_list(mn_args, 'n')

# compute list offsets (anisotropic case only)
m_sum = 0.0
m_off = [0]
for i in m_list:
  m_sum = m_sum + i
  m_off = m_off + [m_sum]
n_sum = 0.0
n_off = [0]
for j in n_list:
  n_sum = n_sum + j
  n_off = n_off + [n_sum]

# some basic sanity checks
if(x0 >= x1):
  print("ERROR: invalid X-range: x0 must be < x1")
  sys.exit(1)
if(y0 >= y1):
  print("ERROR: invalid Y-range: y0 must be < y1")
  sys.exit(1)
if(m <= 0):
  print("ERROR: 'm' must be positive")
  sys.exit(1)
if(n <= 0):
  print("ERROR: 'n' must be positive")
  sys.exit(1)

# compute vertex, edge and quad counts
nv = (m+1)*(n+1)
ne = m*(n+1) + n*(m+1) + m*n # horizontal + vertical + diagonal
nq = 2*m*n

# print some basic information
print("Domain: [%g , %g] x [%g , %g]" % (x0, x1, y0, y1))
print("Slices: %i x %i" % (m, n))
if len(m_list) == 0:
  print("X-Size: isotropic: %g" % ((x1-x0)/float(m)))
else:
  print("X-Size: anisotropic: " + print_aniso_list(x1-x0, m_list, m_sum))
if len(n_list) == 0:
  print("Y-Size: isotropic: %g" % ((y1-y0)/float(n)))
else:
  print("Y-Size: anisotropic: " + print_aniso_list(y1-y0, n_list, n_sum))
print("Verts.: %i" % nv)
print("Edges.: %i" % ne)
print("Trias.: %i" % nq)
print("")

# try to create the output file
print("Creating '" + filename + "'...")
f = open(filename, "wt")
# write header
f.write("<FeatMeshFile version=\"1\" mesh=\"conformal:simplex:2:2\">\n")
f.write("  <!-- Generated by the FEAT '" + myself + "' tool -->\n")
# write tool call as comment
f.write("  <!-- call: " + myself)
for i in range(1,8):
  f.write(" " + sys.argv[i])
f.write(" -->\n")
# write mesh
f.write("  <Mesh type=\"conformal:simplex:2:2\" size=\"%i %i %i\">\n" % (nv, ne, nq))
f.write("    <Vertices>\n")
# write vertices in line-wise order
for i in range(0,n+1):
  # anisotropic in Y?
  if len(n_off) > 1:
    y = y0 + (y1-y0)*(n_off[i]/n_sum)
  else:
    y = y0 + (y1-y0)*(float(i)/float(n))
  for j in range(0,m+1):
    # anisotropic in X?
    if len(m_off) > 1:
      x = x0 + (x1-x0)*(m_off[j]/m_sum)
    else:
      x = x0 + (x1-x0)*(float(j)/float(m))
    f.write("      %g %g\n" % (x, y))
f.write("    </Vertices>\n")
f.write("    <Topology dim=\"1\">\n")
# write horizontal edges
for i in range(0,n+1):
  for j in range(0,m):
    f.write("      %i %i\n" % (i*(m+1)+j, i*(m+1)+j+1))
# write vertical edges
for j in range(0,m+1):
  for i in range(0,n):
    f.write("      %i %i\n" % (i*(m+1)+j, (i+1)*(m+1)+j))
# write diagonal edges
for i in range(0,n):
  for j in range(0,m):
    f.write("      %i %i\n" % (i*(m+1)+j, (i+1)*(m+1)+j+1))
f.write("    </Topology>\n")
f.write("    <Topology dim=\"2\">\n")
for i in range(0,n):
  for j in range(0,m):
    # lower left triangle
    f.write("      %i %i %i\n" % (i*(m+1)+j, i*(m+1)+j+1, (i+1)*(m+1)+j+1))
    # upper right triangle
    f.write("      %i %i %i\n" % ((i+1)*(m+1)+j+1, (i+1)*(m+1)+j, i*(m+1)+j))
f.write("    </Topology>\n")
f.write("  </Mesh>\n")
# write left meshpart
mp_vert(f, "bnd:l", 0, m, n)
# write right meshpart
mp_vert(f, "bnd:r", m, m, n)
# write bottom meshpart
mp_horz(f, "bnd:b", 0, m, n)
# write top meshpart
mp_horz(f, "bnd:t", n, m, n)
f.write("</FeatMeshFile>\n")

print("Done!")
