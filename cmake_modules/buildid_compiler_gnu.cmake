# FEAST BUILD_ID mechanism:
#   big fat mapping of GNU compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (gnu)
#   CMAKE_CXX_COMPILER (g++)
#   FEAST_CXX_FLAGS_INTERNAL (string)
#
# Maintainer information:
# The GNU compiler selection uses the -march flag to create architecture-aware
#   code. -march implies -mtune, and -mcpu is a deprecated synonym for -mtune.
#
# We must not use -march=native, because it prevents cross compilation. Instead,
#   compile a dummy program with "g++ -march=native --verbose" and take the options
#   it suggests as our settings.
#
# All settings below have been generated with gcc 4.4.3
#
# IMPORTANT: If you want the compiler to be aware of specific cache sizes
#   and layout, you *must* use -march=native or pass things in manually,
#   because this build-system is too coarse-grained for this.
#
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# set compiler and compiler name
set (CMAKE_CXX_COMPILER "g++")
set (FEAST_COMPILER_NAME "gnu")

# bail out for unsupported OSes
if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "SunOS"))
  message (STATUS "##############################################################")
  message (STATUS "GNU compiler support is only available for Linux/SunOS        ")
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

# ensure that GNU compiler version is at least 4.4
if ( (CMAKE_SYSTEM_NAME STREQUAL "Linux") OR
     (CMAKE_SYSTEM_NAME STREQUAL "SunOS") )
  set (GNU_VERSION_STRING "")
  set (GNU_VERSION_VALUE "")
  set (GNU_VERSION_OK ON)
  # ask the compiler about its version by executing it
  exec_program("${CMAKE_CXX_COMPILER} --version 2>&1 | head -n 1" OUTPUT_VARIABLE GNU_VERSION_STRING)
  # extract version from compiler output (because we know what the output looks like)
  string(REGEX REPLACE ".*([0-9]+[.]+[0-9]+[.]+[0-9]).*" "\\1" GNU_VERSION_VALUE "${GNU_VERSION_STRING}")
  # in case my regexp is not clever enough
  string(STRIP "${GNU_VERSION_VALUE}" GNU_VERSION_VALUE)
  message (STATUS "Found GNU compiler suite v.${GNU_VERSION_VALUE}")
  # and rely on string comparisons TODO convert to a number
  string(COMPARE GREATER "${GNU_VERSION_VALUE}" "4.4.0" GNU_VERSION_OK)
  if (NOT GNU_VERSION_OK)
    message (STATUS "##############################################################")
    message (STATUS "GNU Compiler Version ${GNU_VERSION_VALUE} is not supported,   ")
    message (STATUS "please upgrade to a version better than 4.4.                  ")
    message (STATUS "##############################################################")
    message (FATAL_ERROR "")
  endif()
else ()
  message (STATUS "TODO: Implement GNU version detection for Non-Linux systems.")
endif ()


# if compiler flags are not passed externally, determine our own
if (FEAST_CXX_FLAGS_INTERNAL STREQUAL "")

  # generic settings independent of arch and optimisation level
  set (FEAST_CXX_FLAGS_INTERNAL "-pipe")


  if (FEAST_DEBUG_MODE)
    # unoptimised settings for all archs
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O0 -std=c++0x -pedantic -Wall -Wextra -Wundef -g")

  else ()
    # optimised settings for all currently supported archs
    # -ftree-vectorize is part of -O3 already
    # -ffast-math turns on additional non-IEEE-compliant optimisations
    # -mfpmath=sse forces SSE math rather than FPU (i387) floating point math
    # the other default flags are self-explanatory and are not included in -O3
    set (FEAST_CXX_FLAGS_INTERNAL "-O3 -ffast-math -foptimize-register-move -fprefetch-loop-arrays -funroll-loops -mfpmath=sse")

    # please try to maintain the same order as in the buildid_arch module
    # Intel CPUs
    if (FEAST_CPU_TYPE STREQUAL "i486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=i486 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumpro")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentiumpro -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium2 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium3")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium3 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumm")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium-m -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4m")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium4m -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "coresolo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "coreduo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core2 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "penryn")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core2 -msse4.1 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "nehalem")
      string(COMPARE GREATER "${GNU_VERSION_VALUE}" "4.6.0" GNU_VERSION_OK)
      if (GNU_VERSION_OK)
        set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=corei7 -m64")
      else ()
        # note: older GCC versions do not support -match=corei7 , so emulate: core2+sse4.2=corei7
        set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core2 -msse4.2 -m64")
      endif ()
    elseif (FEAST_CPU_TYPE STREQUAL "westmere")
      string(COMPARE GREATER "${GNU_VERSION_VALUE}" "4.6.0" GNU_VERSION_OK)
      if (GNU_VERSION_OK)
        set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=corei7 -m64")
      else ()
        # note: older GCC versions do not support -match=corei7 , so emulate: core2+sse4.2=corei7
        set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core2 -msse4.2 -m64")
      endif ()
    elseif (FEAST_CPU_TYPE STREQUAL "itanium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=itanium")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium4m")
    elseif (FEAST_CPU_TYPE STREQUAL "nocona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=nocona -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=itanium2")

    # AMD CPUs
    elseif (FEAST_CPU_TYPE STREQUAL "amd486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "k5")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "k6")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=k6 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "athlonxp")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon-xp -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "opteron")
      # note that k8, opteron, athlon64 and athlon-fx are synonymous
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=k8 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon64")
      # note that k8, opteron, athlon64 and athlon-fx are synonymous
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=k8 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "opteronx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=k8-sse3 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "turionx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=k8-sse3 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "barcelona")
      # note that amdfam10 and barcelona are synonymous
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=barcelona -m64")

    # generic settings for all archs
    else ()
      message (STATUS "WARNING: Unsupported architecture/compiler combination found.")
      message (STATUS "Using generic optimisation flags. ")
    endif ()

  endif (FEAST_DEBUG_MODE)
endif (FEAST_CXX_FLAGS_INTERNAL STREQUAL "")

# check if compiler supports the given flags
set (COMPILER_SUPPORTS_FLAG ON)
include (CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("${FEAST_CXX_FLAGS_INTERNAL}" COMPILER_SUPPORTS_FLAG)
if (NOT COMPILER_SUPPORTS_FLAG)
  message (STATUS "##############################################################")
  message (STATUS "One of the flags ${FEAST_CXX_FLAGS_INTERNAL} is apparently             ")
  message (STATUS "unsupported by the GNU Compiler v.${GNU_VERSION_VALUE}.       ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
