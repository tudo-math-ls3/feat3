import sys
import os
from common import patch_file

directory = sys.argv[1]
print("\nPatching ParMETIS in " + directory + "...")

# remove some outdated MS-specific preprocessor mumbo-jumbo
patch_file(directory, os.path.join("GKlib", "gk_arch.h"), [
  #[35, '#include "ms_stdint.h"', '  #include <stdint.h>'],
  #[36, '#include "ms_inttypes.h"', '  #include <inttypes.h>'],
  #[48, '', '#if defined(_WIN32) && !defined(WIN32)\n  #define WIN32\n#endif'],
  #[61, '#ifdef __MSC__', '#if 0'],
  #[70, '', '#ifdef __MSC__\n#define __thread __declspec(thread)\n#endif']
])
# patch some broken pointer casts
patch_file(directory, os.path.join("GKlib", "gkregex.c"), [
  [5089, 'postorder (elem, mark_opt_subexp, (void *) (long) elem->token.opr.idx);',
     '    postorder (elem, mark_opt_subexp, (void *) (intptr_t) elem->token.opr.idx);'],
  [6301, 'int idx = (int) (long) extra;', '  int idx = (int) (intptr_t) extra;']
])
# patch more outdated MS-specific preprocessor mumbo-jumbo
patch_file(directory, os.path.join("METIS", "include", "metis.h"), [
  [66, '#ifdef COMPILER_MSC', '#if 0 /*def COMPILER_MSC*/']
])
