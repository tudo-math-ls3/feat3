import os
import subprocess

def is_found(name):
  try:
    devnull = open(os.devnull)
    subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
  except OSError as e:
    if e.errno == os.errno.ENOENT:
      return False
  return True

def remove_string(array, string):
  while array.count(string) > 0:
    array.remove(string)
  return array
