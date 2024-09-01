import sys
import os
from common import patch_file

directory = sys.argv[1]
print("\nPatching SuperLU(dist) in " + directory + "...")

# don't use ParMETIS
patch_file(directory, os.path.join("SRC", "superlu_dist_config.h"), [
  [10, '#define HAVE_PARMETIS TRUE', '/* #undef HAVE_PARMETIS */'],
])
# don't include unistd.h header
patch_file(directory, os.path.join("SRC", "util.c"), [
  [25, '#include <unistd.h>', '/*#include <unistd.h>*/'],
])
