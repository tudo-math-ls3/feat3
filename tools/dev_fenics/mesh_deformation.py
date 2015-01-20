#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# vim: ft=python
import math
import dolfin

class TestProblem:
	def __init__(self, bc0, bc1, bc2, dim=2):
		self.bc0 = bc0
		self.bc1 = bc1
		self.bc2 = bc2
		self.dim = dim

bc0 = dolfin.Expression(("x[0]", "x[1] + 0.5*t*x[1]*sin(2*pi*x[0])"),t=0,pi=math.pi)
bc1 = dolfin.Expression(("x[0]", "x[1] + 0.5*t*x[1]*sin(2*pi*x[0])"),t=0,pi=math.pi)
bc2 = dolfin.Expression(("x[0]", "x[1]"))

my_deform = TestProblem(bc0,bc1,bc2)
