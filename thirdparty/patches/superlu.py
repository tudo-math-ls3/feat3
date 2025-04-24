import sys
import os
import shutil
from common import patch_file

directory = sys.argv[1]
print("\nPatching SuperLU(dist) in " + directory + "...")
# get the patches directory by searching for the script file
patches_directory = os.path.dirname(os.path.realpath(__file__))

# since we do extensive changes to the cmake lists, we simply copy
# a template into the directory
# first of all patch SRC/CMakeLists.txt
shutil.copyfile(os.path.join(patches_directory, "superlu_patches", "superlu_SRC_CMakeLists.in"),
                os.path.join(directory, "SRC", "CMakeLists.txt"))

# first of all patch CBLAS/CMakeLists.txt
shutil.copyfile(os.path.join(patches_directory, "superlu_patches", "superlu_CBLAS_CMakeLists.in"),
                os.path.join(directory, "CBLAS", "CMakeLists.txt"))