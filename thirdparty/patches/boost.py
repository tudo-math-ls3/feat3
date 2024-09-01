import sys
import os
from common import patch_file

directory = sys.argv[1]
print("\nPatching boost in " + directory + "...")

patch_file(directory, os.path.join("libs", "serialization", "CMake", "CMakeLists.txt"), [
  [3, 'cmake_minimum_required(VERSION 3.0)', 'cmake_minimum_required(VERSION 3.5...3.16)']
])
patch_file(directory, os.path.join("libs", "hana", "example", "cmake_integration", "CMakeLists.txt"), [
  [6, 'cmake_minimum_required(VERSION 3.0)', 'cmake_minimum_required(VERSION 3.5...3.16)']
])
patch_file(directory, os.path.join("libs", "hana", "test", "deploy", "CMakeLists.txt"), [
  [5, 'cmake_minimum_required(VERSION 3.0)', 'cmake_minimum_required(VERSION 3.5...3.16)']
])
patch_file(directory, os.path.join("libs", "predef", "test", "test_cmake", "CMakeLists.txt"), [
  [16, 'cmake_minimum_required( VERSION 3.0 )', 'cmake_minimum_required(VERSION 3.5...3.16)']
])
patch_file(directory, os.path.join("libs", "predef", "CMakeLists.txt"), [
  [23, 'cmake_minimum_required( VERSION 3.0 )', 'cmake_minimum_required(VERSION 3.5...3.16)']
])
patch_file(directory, os.path.join("libs", "callable_traits", "CMakeLists.txt"), [
  [26, 'cmake_minimum_required(VERSION 3.0)', 'cmake_minimum_required(VERSION 3.5...3.16)']
])
