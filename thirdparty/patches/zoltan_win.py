import sys
import os
from common import patch_file

directory = sys.argv[1]
print("\nPatching Zoltan in " + directory + "...")

# patch "phg/phg_util.c"
patch_file(directory, os.path.join("src", "phg", "phg_util.c"), [
  [75, "{", "{/*"],
  [85, "}", "*/}"]
])
# patch "rcb/inertial2d.c"
patch_file(directory, os.path.join("src", "rcb", "inertial2d.c"), [
  [66, "#define max(a, b) ((a) < (b) ? (b) : (a))", "/*#define max(a, b) ((a) < (b) ? (b) : (a))*/"],
  [67, "#define min(a, b) ((a) > (b) ? (b) : (a))", "/*#define min(a, b) ((a) > (b) ? (b) : (a))*/"]
])
# patch "rcb/inertial3d.c"
patch_file(directory, os.path.join("src", "rcb", "inertial3d.c"), [
  [54, "#define max(a, b) ((a) < (b) ? (b) : (a))", "/*#define max(a, b) ((a) < (b) ? (b) : (a))*/"],
  [55, "#define min(a, b) ((a) > (b) ? (b) : (a))", "/*#define min(a, b) ((a) > (b) ? (b) : (a))*/"]
])
# patch "zz/zz_const.c"
patch_file(directory, os.path.join("src", "zz", "zz_const.h"), [
  [55, "#include <strings.h>", "/*#include <strings.h>*/"]
])
# patch "zz/zz_util.c"
patch_file(directory, os.path.join("src", "zz", "zz_util.c"), [
  [55, "#include <unistd.h>", "/*#include <unistd.h>*/"]
])
