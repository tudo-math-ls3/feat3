import os
import subprocess


def find_exe(name):
  path = get_output("which " + name)
  if len(path) > 0:
    return path[0]
  else:
    return ""

def is_found(name):
  return find_exe(name) != ""

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
