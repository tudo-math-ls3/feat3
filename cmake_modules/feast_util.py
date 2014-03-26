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
  pipe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
  pipe.stdin.close()
  output = pipe.stdout.read().decode().splitlines()
  return output

def get_output_utf8(command):
  pipe = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
  pipe.stdin.close()
  output = pipe.stdout.read().decode(encoding="UTF-8").splitlines()
  return output
