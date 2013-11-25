import os
import subprocess

def is_found(name):
  try:
    with open(os.devnull, "w") as fp:
      subprocess.Popen([name], stdout=fp, stderr=subprocess.STDOUT).communicate()
  except OSError as e:
    if e.errno == os.errno.ENOENT:
      return False
  return True

def remove_string(array, string):
  while array.count(string) > 0:
    array.remove(string)
  return array

def get_output(command):
  #if "check_output" not in dir( subprocess ): #deactivated as its not available bevor python 2.7
  pipe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE)
  pipe.stdin.close()
  output = pipe.stdout.read().splitlines()
  #else:
  #  output = subprocess.check_output(command.split(), shell=True).splitlines()
  return output
