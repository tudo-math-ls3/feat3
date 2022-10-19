#!/usr/bin/env python
########################################################################################################################
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2022 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
########################################################################################################################
# Third-Party Library build script for Visual Studio 2022 (VC17)
# ----------------------------------------------------------------------------------------------------------------------
#
# \author Peter Zajac
########################################################################################################################
import os
import sys
import re
import shutil
import subprocess

sys.dont_write_bytecode = True

def modeline(bm, tag, mode):
  if bm.find(tag) >= 0:
    return '  <' + mode + '>true</' + mode + '>\n';
  else:
    return '  <' + mode + '>false</' + mode + '>\n';

# build tags
build_tags = [
  #'cuda',
  'omp',
  #'mkl',
  'mpi',
]

# enabled build modes
build_modes = ['dbg', 'opt', 'dbg-omp', 'opt-omp', 'dbg-mpi', 'opt-mpi']
platforms = ['x64']


# write build modes
fo = open("build-modes.xml", "wt")

print("Writing 'build-modes.xml'...")
fo.write('<?xml version="1.0" encoding="utf-8"?>\n')
fo.write('<Project DefaultTargets="Build" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">\n')
# write project configurations
fo.write('<!-- Project Configurations -->\n')
fo.write('<ItemGroup Label="ProjectConfigurations">\n')
for bm in build_modes:
  for pf in platforms:
    fo.write('<ProjectConfiguration Include="' + bm + '|' + pf + '">\n')
    fo.write('  <Configuration>' + bm + '</Configuration>\n')
    fo.write('  <Platform>' + pf + '</Platform>\n')
    fo.write('</ProjectConfiguration>\n')
fo.write('</ItemGroup>\n')
# write build mode mapping
fo.write('<!-- Configuration to Build-Mode mapping -->\n')
for bm in build_modes:
  fo.write('<PropertyGroup Label="BuildMode" Condition="\'$(Configuration)\'==\'' + bm + '\'">\n')
  fo.write(modeline(bm, 'dbg', 'DebugMode'))
  fo.write(modeline(bm, 'omp', 'EnableOMP'))
  fo.write(modeline(bm, 'mpi', 'EnableMPI'))
  fo.write(modeline(bm, 'mkl', 'EnableMKL'))
  fo.write(modeline(bm, 'cuda', 'EnableCUDA'))
  fo.write('</PropertyGroup>\n')
fo.write('</Project>\n')
fo.close()
print("done!")