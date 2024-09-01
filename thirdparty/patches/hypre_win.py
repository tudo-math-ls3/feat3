import sys
import os
from common import patch_file

directory = sys.argv[1]
print("\nPatching Hypre in " + directory + "...")

# patch out non-existing header includes
patch_file(directory, os.path.join("src", "matrix_matrix", "HYPRE_ConvertParCSRMatrixToDistributedMatrix.c"), [
  [15, '#include <gmalloc.h>', '/*#include <gmalloc.h>*/']
])
patch_file(directory, os.path.join("src", "matrix_matrix", "HYPRE_ConvertPETScMatrixToDistributedMatrix.c"), [
  [15, '#include <gmalloc.h>', '/*#include <gmalloc.h>*/']
])

# patch 32 vs 64 bit int issues in utilities
patch_file(directory, os.path.join("src", "utilities", "_hypre_utilities.h"), [
  #[2171, 'hypre_ulongint h64 = HYPRE_XXH_PRIME64_5 + sizeof(input);', '  hypre_ulonglongint h64 = HYPRE_XXH_PRIME64_5 + sizeof(input);'],
  #[2173, 'hypre_ulongint k1 = input;', '  hypre_ulonglongint k1 = input;']
])
