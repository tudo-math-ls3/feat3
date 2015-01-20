#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# vim: ft=python
import math
import dolfin

# Class for analytic levelset function for the Rumpf smoother
class TestProblem:
	def __init__(self, displacement, scaling, dim=2):
		self.dim = dim
		self.value = value
		self.value.displacement = displacement
		self.value.scaling = scaling

# Circle with radius ... centered at (0.75, 0.5) for t = 0
value = dolfin.Expression("displacement + scaling*sqrt(  pow(x[0] - 0.25*(2. + cos(t)), 2) + pow( x[1] - 0.25*(2. + sin(3.*t)), 2))",displacement = 0., scaling = 1., t = 0.)
# Vertical line at displacement + 0.5
#value = dolfin.Expression("displacement + scaling*(x[0] - 0.5)",displacement = 0., scaling = 1., t = 0.)

my_lvlset_function = TestProblem(0.0, -1., 2)
