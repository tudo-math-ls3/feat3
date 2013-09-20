# FEAST BUILD_ID mechanism:
#   big fat mapping of clang compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (clang)
#   CMAKE_CXX_COMPILER (clang)
#   FEAST_CXX_FLAGS_INTERNAL (string)
#
# Maintainer information:
# The clang compiler selection uses the -march flag to create architecture-aware
#   code. -march implies -mtune, and -mcpu is a deprecated synonym for -mtune.
#
# We must not use -march=native, because it prevents cross compilation. Instead,
#   compile a dummy program with "g++ -march=native --verbose" and take the options
#   it suggests as our settings.
#
# All settings below have been generated with clang 3.0
#
# IMPORTANT: If you want the compiler to be aware of specific cache sizes
#   and layout, you *must* use -march=native or pass things in manually,
#   because this build-system is too coarse-grained for this.
#
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# set compiler and compiler name
set (CMAKE_CXX_COMPILER "clang")
set (FEAST_COMPILER_NAME "clang")

# bail out for unsupported OSes
# see issue 00035
if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "SunOS"))
  message (STATUS "##############################################################")
  message (STATUS "clang compiler support is only available for Linux/SunOS        ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()

# verify that the given compiler is available
# cmake syntax lesson (for the find_program() command which is basically a wrapper
# around "which" on linux): quoting from cmake docs version 2.8.4
#   If the full path to a file is found the result is stored in the variable and
#   the search will not be repeated unless the variable is cleared. If nothing is
#   found, the result will be <VAR>-NOTFOUND, and the search will be attempted
#   again the next time find_file is invoked with the same variable.
# end quote
set (FINDPROG "FINDPROG-NOTFOUND")
find_program (FINDPROG NAMES ${CMAKE_CXX_COMPILER})
if (FINDPROG STREQUAL "FINDPROG-NOTFOUND")
   message (STATUS "##############################################################")
   message (STATUS "ERROR: Compiler \"${CMAKE_CXX_COMPILER}\" not found.          ")
   message (STATUS "       Please make sure that the compiler you selected is     ")
   message (STATUS "       available in your environment, for instance by loading ")
   message (STATUS "       the corresponding modules.                             ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()


# if compiler flags are not passed externally, determine our own
if (FEAST_CXX_FLAGS_INTERNAL STREQUAL "")

  # ensure that OS is 64-bit
  if ( (CMAKE_SYSTEM_NAME STREQUAL "Linux") )
    exec_program("uname -m" OUTPUT_VARIABLE UNAME_OUTPUT)
    if (NOT UNAME_OUTPUT STREQUAL "x86_64")
      message (STATUS "##############################################################")
      message (STATUS "Only 64-bit operating systems are supported by default.       ")
      message (STATUS "You might want to switch to manual mode.                       ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")
    endif ()
   else ()
     message (STATUS "TODO: Implement 64-bit OS detection for Non-Linux systems.")
  endif ()

  # generic settings independent of arch and optimisation level
  set (FEAST_CXX_FLAGS_INTERNAL "-pipe  -Wno-unused-parameter -m64")


  if (FEAST_DEBUG_MODE)
    # unoptimised settings for all archs
    # the following flag might be useful: -m128bit-long-double
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O0 -std=c++11 -Wall -Wextra -Wundef -ggdb -Wshorten-64-to-32 -Wconversion -Wstrict-aliasing=2 -Wunknown-pragmas -Wundef -Wno-unused-value")

  else ()
    # optimised settings for all currently supported archs
    # -ftree-vectorize is part of -O3 already
    # -ffast-math turns on additional non-IEEE-compliant optimisations
    # the other default flags are self-explanatory and are not included in -O3
    # the following flag might be useful: -m128bit-long-double
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O3 -ggdb -std=c++11")


  endif (FEAST_DEBUG_MODE)
endif (FEAST_CXX_FLAGS_INTERNAL STREQUAL "")

# check if compiler supports the given flags
set (COMPILER_SUPPORTS_FLAG ON)
include (CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("${FEAST_CXX_FLAGS_INTERNAL}" COMPILER_SUPPORTS_FLAG)
if (NOT COMPILER_SUPPORTS_FLAG)
  message (STATUS "##############################################################")
  message (STATUS "One of the flags ${FEAST_CXX_FLAGS_INTERNAL} is apparently    ")
  message (STATUS "unsupported by the clang Compiler.                            ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
