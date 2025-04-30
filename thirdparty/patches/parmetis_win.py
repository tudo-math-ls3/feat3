import sys
import os
from common import patch_file

directory = sys.argv[1]
print("\nPatching ParMETIS in " + directory + "...")

# remove some outdated MS-specific preprocessor mumbo-jumbo
patch_file(directory, os.path.join("GKlib", "gk_arch.h"), [
  [35, '#include "gk_ms_stdint.h"', '  #include <stdint.h>'],
  [36, '#include "gk_ms_inttypes.h"', '  #include <inttypes.h>'],
  [60, '', '#ifdef __MSC__\n#define __thread __declspec(thread)\n#endif'],
  [66, '#ifndef INFINITY', '#if 0']
])
# patch some broken pointer casts
patch_file(directory, os.path.join("GKlib", "gkregex.c"), [
  [5089, 'postorder (elem, mark_opt_subexp, (void *) (long) elem->token.opr.idx);',
     '    postorder (elem, mark_opt_subexp, (void *) (intptr_t) elem->token.opr.idx);'],
  [6301, 'int idx = (int) (long) extra;', '  int idx = (int) (intptr_t) extra;']
])
# patch POSIX file read
patch_file(directory, os.path.join("GKlib", "io.c"), [
  [18, '', '#include <io.h>'],
  [63, 'if ((rsize = read(fd, buf, tsize)) == -1)', 'if ((rsize = _read(fd, buf, tsize)) == -1)'],
  [84, 'if ((size = write(fd, buf, tsize)) == -1)', 'if ((size = _write(fd, buf, tsize)) == -1)']
])
# set index and real types
patch_file(directory, os.path.join("metis", "include", "metis.h"), [
  [33, '//#define IDXTYPEWIDTH 32', '#define IDXTYPEWIDTH 64'],
  [43, '//#define REALTYPEWIDTH 32', '#define REALTYPEWIDTH 64'],
  [75, '#define INT32_MIN    ((int32_t)_I32_MIN)', '//#define INT32_MIN    ((int32_t)_I32_MIN)'],
  [76, '#define INT32_MAX    _I32_MAX', '//#define INT32_MAX    _I32_MAX'],
  [77, '#define INT64_MIN    ((int64_t)_I64_MIN)', '//#define INT64_MIN    ((int64_t)_I64_MIN)'],
  [78, '#define INT64_MAX    _I64_MAX', '//#define INT64_MAX    _I64_MAX']
])
