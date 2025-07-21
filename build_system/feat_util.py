# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import os
import sys
import subprocess
import json

#Stupid fix for the following problem:
#distutils will be removed in python 3.12, but shutil's which function is only
#implemented in the 3.2 Version, so for compatibility between python2 and 3
#we simply define which ourselves...
#adventage of which -> works on windows
def which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.

    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.
    """
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get("PATH", os.defpath)

    if not path:
        return None

    path = path.split(os.pathsep)
    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)
        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]
    seen = set()
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if _access_check(name, mode):
                    return name
    return None




def is_found(name):
  return which(name) != None

def find_exe(name):
  return which(name)

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

def preset_header():
  return {
    "version": 8,
    "cmakeMinimumRequired": {
      "major": 3,
      "minor": 28,
      "patch": 0,
    },
    "configurePresets": []
  }

def warning(*args):
  YELLOW = "\033[33;49;1m"
  RESET = "\033[0m"
  print(YELLOW + "Warning: " + " ".join(args) + RESET)

def error(*args):
  RED = "\033[31;49;1m"
  RESET = "\033[0m"
  print(RED + "Error: " + " ".join(args) + RESET)
  sys.exit(1)

def header(*args):
  BLUE = "\033[34;49m"
  RESET = "\033[0m"
  print(BLUE + " ".join(args) + RESET)

def append_preset_to_file(filename, configure_preset):
  presets = preset_header()

  if os.path.isfile(filename):
    with open(filename, "r") as in_file:
      content = in_file.read()
      if content:
        presets = json.loads(content)


  present = False
  for idx, preset in enumerate(presets["configurePresets"]):
    if preset["name"] == configure_preset["name"]:
      presets["configurePresets"][idx] = configure_preset
      present = True

  if not present:
    presets["configurePresets"].append(configure_preset)

  try:
    with open(filename, "w") as out:
      json.dump(presets, out, indent=2)
  except PermissionError:
    warning((f"Can not open file {filename} for writing.\n"
              "If you are trying to develop FEAT, move your repository to a location you have write permissions for.\n"
              "If you are trying to configure FEAT on a HPC-cluster and do not have write-access to your home directory,\n"
              "you can use the '--no-preset' flag to stop preset generation and supress this warning."))
