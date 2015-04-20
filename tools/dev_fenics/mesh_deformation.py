#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# vim: ft=python
import math
import dolfin
from dolfin import SubDomain, DOLFIN_EPS, RectangleMesh, BoxMesh

def gen_mesh(n_gridpoints,dim):
  if dim == 2:
    mesh = RectangleMesh(0,0,1,1,n_gridpoints,n_gridpoints,"left");
  else:
    mesh = BoxMesh(0,0,0,1,1,1,n_gridpoints,n_gridpoints,n_gridpoints);

  return mesh

class Top(SubDomain):
  def inside(self, x, on_boundary):
    tol = DOLFIN_EPS
    return on_boundary and  abs(x[1]-1.0) < tol

class NotTop(SubDomain):
  def inside(self, x, on_boundary):
    tol = DOLFIN_EPS
    return on_boundary and  abs(x[1]-1.0) > tol

class TestProblem:
  def __init__(self, bc0, bc1, dim=2):
    self.bc0 = bc0
    self.bc1 = bc1
    self.dim = dim

# This BC causes the nonlinear solver in RumpfSmoother to throw lots of NaNs.
# Just the presence in this file is enough, it does not have to be used.
# If an x[1] appears somewhere, everything is fine. Bug in FEniCS 1.4?
#bc0 = dolfin.Expression(("x[0]", "1 + t*sin(2*pi*x[0])"),t=0,pi=math.pi)
bc0 = dolfin.Expression(("x[0]", "x[1]"))
bc1 = dolfin.Expression(("x[0]", "x[1]"))

my_deform = TestProblem(bc0,bc1)
