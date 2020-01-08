# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import os
import subprocess
from distutils import spawn

def is_found(name):
  return spawn.find_executable(name) != None

def find_exe(name):
  return spawn.find_executable(name)

def remove_string(array, string):
  while array.count(string) > 0:
    array.remove(string)
  return array

def remove_substring(array, sub):
  for s in array:
    if sub in s:
      array.remove(s)
  return array

def get_output(command):
  pipe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
  pipe.stdin.close()
  output = pipe.stdout.read().decode().splitlines()
  return output

def get_output_utf8(command):
  pipe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
  pipe.stdin.close()
  output = pipe.stdout.read().splitlines()
  return output
