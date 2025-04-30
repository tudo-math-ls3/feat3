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

# true and false are not ANSI-C keywords
patch_file(directory, os.path.join("src", "distributed_ls", "pilut", "parilut.c"), [
  [964, 'while (true) {', 'while (1) {'],
  [1048, 'while (true) {', 'while (1) {']
])
patch_file(directory, os.path.join("src", "distributed_ls", "pilut", "trifactor.c"), [
  [290, 'hypre_SetUpFactor( ddist, ldu, maxnz,   petotal, rind, imap, &maxsend,   true,', 'hypre_SetUpFactor( ddist, ldu, maxnz,   petotal, rind, imap, &maxsend, 1,'],
  [307, 'hypre_SetUpFactor( ddist, ldu, maxnz,   petotal, rind, imap, &maxsend,   false,', 'hypre_SetUpFactor( ddist, ldu, maxnz,   petotal, rind, imap, &maxsend, 0,']
])
