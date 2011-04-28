# FEAST BUILD_ID mechanism:
#   big fat mapping of OPEN64 compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (open64)
#   CMAKE_CXX_COMPILER (openCC)
#   FEAST_CXX_FLAGS_INTERNAL (string)
#
# Maintainer information:
# The OPEN64 compiler selection uses the -march flag to create architecture-aware
#   code. -march implies -mtune, and -mcpu is a deprecated synonym for -mtune.
#   Generally, architecture-aware optimisations are not a strong point of OPEN64.
#
# We must not use -march=auto, because it prevents cross compilation. Instead,
#   compile a dummy program with "openCC -march=auto --verbose" and take the options
#   it suggests as our settings.
#
# All settings below have been generated with open64 4.2.4,
# see http://developer.amd.com/cpu/open64/onlinehelp/pages/x86_open64_help.htm
#
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# set compiler and compiler name
set (CMAKE_CXX_COMPILER "openCC")
set (FEAST_COMPILER_NAME "open64")

# bail out for unsupported OSes
# see issue 00035
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message (STATUS "##############################################################")
  message (STATUS "OPEN64 compiler support is only available for Linux           ")
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

# ensure that OPEN64 compiler version is at least 4.2.4
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set (OPEN64_VERSION_STRING "")
  set (OPEN64_VERSION_VALUE "")
  set (OPEN64_VERSION_OK ON)
  # ask the compiler about its version by executing it
  exec_program("${CMAKE_CXX_COMPILER} --version 2>&1 | head -n 1" OUTPUT_VARIABLE OPEN64_VERSION_STRING)
  # extract version from compiler output (because we know what the output looks like)
  # see isseu #34
  string(REGEX REPLACE ".*([0-9]+[.]+[0-9]+[.]+[0-9]).*" "\\1" OPEN64_VERSION_VALUE "${OPEN64_VERSION_STRING}")
  # in case my regexp is not clever enough
  string(STRIP "${OPEN64_VERSION_VALUE}" OPEN64_VERSION_VALUE)
  message (STATUS "Found OPEN64 compiler suite v.${OPEN64_VERSION_VALUE}")
  # and rely on string comparisons TODO convert to a number
  string(COMPARE GREATER "${OPEN64_VERSION_VALUE}" "4.2.3" OPEN64_VERSION_OK)
  if (NOT OPEN64_VERSION_OK)
    message (STATUS "##############################################################")
    message (STATUS "OPEN64 Compiler Version ${OPEN64_VERSION_VALUE} is not supported,")
    message (STATUS "please upgrade to a version better than 4.2.4.                ")
    message (STATUS "##############################################################")
    message (FATAL_ERROR "")
  endif()
else ()
  message (STATUS "TODO: Implement OPEN64 version detection for Non-Linux systems.")
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
  set (FEAST_CXX_FLAGS_INTERNAL "-pipe")


  if (FEAST_DEBUG_MODE)
    # unoptimised settings for all archs
    # Note: OPEN64 supports either 1998 ISO C++, or the same with GNU extensions (-std=gnu++98)
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O0 -std=gnu++98 -pedantic -Wall -Wextra -Wundef -g")

  else ()
    # optimised settings for all currently supported archs
    # -ftree-vectorize is ignored by OPEN64
    # -ffast-math turns on additional non-IEEE-compliant optimisations
    # -ffast-stdlib uses faster versions of some standard functions, when available
    # the other default flags are self-explanatory and are not included in -O3
    set (FEAST_CXX_FLAGS_INTERNAL "-O3 -ffast-math -ffast-stdlib -funroll-loops")

    # please try to maintain the same order as in the buildid_arch module
    # Intel CPUs
    if (FEAST_CPU_TYPE STREQUAL "i486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumpro")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium3")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumm")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4m")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=em64t")
    elseif (FEAST_CPU_TYPE STREQUAL "coresolo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core -msse3")
    elseif (FEAST_CPU_TYPE STREQUAL "coreduo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core -msse3")
    elseif (FEAST_CPU_TYPE STREQUAL "penryn")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core -msse4.1")
    elseif (FEAST_CPU_TYPE STREQUAL "nehalem")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core -msse4.2")
    elseif (FEAST_CPU_TYPE STREQUAL "westmere")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core -msse4.2")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=pentium4")
    elseif (FEAST_CPU_TYPE STREQUAL "nocona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=core -msse3")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")

    # AMD CPUs
    elseif (FEAST_CPU_TYPE STREQUAL "amd486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "k5")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "k6")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon")
    elseif (FEAST_CPU_TYPE STREQUAL "athlonxp")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon")
    elseif (FEAST_CPU_TYPE STREQUAL "opteron")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=opteron")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon64")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon64")
    elseif (FEAST_CPU_TYPE STREQUAL "opteronx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon64fx")
    elseif (FEAST_CPU_TYPE STREQUAL "turionx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=athlon64fx")
    elseif (FEAST_CPU_TYPE STREQUAL "barcelona")
      # note that amdfam10 and barcelona are synonymous
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -march=barcelona")

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
  message (STATUS "unsupported by the OPEN64 Compiler v.${OPEN64_VERSION_VALUE}. ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
