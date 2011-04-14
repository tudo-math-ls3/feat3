# FEAST BUILD_ID mechanism:
#   big fat mapping of INTEL compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (intel)
#   CMAKE_CXX_COMPILER (icpc)
#   FEAST_CXX_FLAGS (string)
#
#
# Maintainer information:
#
# We must not use -xHOST even through it generally gives the best code,
# because it prevents cross-compilation. Consequently, if -xHOST is considered
# absolutely necessary for one particular test, BUILD_ID=manual must be used.
#
# The Intel compiler suite follows the following conventions for creating
#   architecture-aware code, as outlined here:
#   http://software.intel.com/en-us/articles/performance-tools-for-software-developers-intel-compiler-options-for-sse-generation-and-processor-specific-optimizations/
#   or in the compiler man page.
#   Basically, there are two ways to set the target microarchitecture / processor category:
#     -mCODE generates code for a processor category identified by <CODE>,
#     -xCODE additionally generates optimisations specific to classes of Intel CPUs
#   <CODE> can be one of the following:
#      AVX (-x only)
#      SSE4.2 (-x only)
#      SSE4.1
#      SSE3
#      SSE2 (default value for -m)
#      IA32 (-m only)
# This consequently means that for any non-Intel processor, we must use
# the -m option. For all (not too ancient) Intel processors, we should use -x,
# and this is essentially everything that needs to be set for optimised,
# architecture-aware builds.
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# set compiler and compiler name
set (CMAKE_CXX_COMPILER "icpc")
set (FEAST_COMPILER_NAME "intel")


# bail out for unsupported OSes
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message (STATUS "##############################################################")
  message (STATUS "Intel compiler support is only available for Linux            ")
  message (FATAL_ERROR "##############################################################")
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
   message (FATAL_ERROR "##############################################################")
endif ()


# ensure that Intel compiler version is at least 11.1
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set (INTEL_VERSION_STRING "")
  set (INTEL_VERSION_VALUE "")
  set (INTEL_VERSION_OK ON)
  # ask the compiler about its version by executing it
  exec_program("${CMAKE_CXX_COMPILER} -V 2>&1 | head -n 1" OUTPUT_VARIABLE INTEL_VERSION_STRING)
  # extract version from compiler output (because we know what the output looks like)
  string(REGEX REPLACE ".*Version[ \t]+([0-9]+[.]+[0-9]).*" "\\1" INTEL_VERSION_VALUE "${INTEL_VERSION_STRING}")
  # in case my regexp is not clever enough
  string(STRIP "${INTEL_VERSION_VALUE}" INTEL_VERSION_VALUE)
  # and rely on string comparisons TODO convert to a number
  string(COMPARE GREATER "${INTEL_VERSION_VALUE}" "11.0" INTEL_VERSION_OK)
  message (STATUS "Found Intel compiler suite v.${INTEL_VERSION_VALUE}")
  if (NOT INTEL_VERSION_OK)
    message (STATUS "##############################################################")
    message (STATUS "Intel Compiler Version ${INTEL_VERSION_VALUE} is not supported,")
    message (STATUS "please upgrade to a version better than 11.1.                 ")
    message (FATAL_ERROR "##############################################################")
  endif()
else ()
  message (STATUS "TODO: Implement Intel version detection for Non-Linux systems.")
endif ()


# generic settings independent of arch and optimisation level
set (FEAST_CXX_FLAGS "")


if (FEAST_DEBUG_MODE)
  # unoptimised settings for all archs
  set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -O0 -std=c++0x -Wall -g")

else ()
  # optimised settings for all currently supported archs
  # please try to maintain the same order as in the buildid_arch module
  # Note: SSE2 is on by default, so we only have to specify what's better
  #       or worse
  set (FEAST_CXX_FLAGS "-O3")

  # Intel CPUs
  if (FEAST_CPU_TYPE STREQUAL "i486")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS}")
  elseif (FEAST_CPU_TYPE STREQUAL "pentium")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "pentiumpro")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "pentium2")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "pentium3")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "pentiumm")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE2")
  elseif (FEAST_CPU_TYPE STREQUAL "pentium4m")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE2")
  elseif (FEAST_CPU_TYPE STREQUAL "coresolo")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE3")
  elseif (FEAST_CPU_TYPE STREQUAL "coreduo")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE3")
  elseif (FEAST_CPU_TYPE STREQUAL "penryn")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE4.1")
  elseif (FEAST_CPU_TYPE STREQUAL "nehalem")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE4.2")
  elseif (FEAST_CPU_TYPE STREQUAL "westmere")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE4.2")
  elseif (FEAST_CPU_TYPE STREQUAL "itanium")
    # no setting necessary, the itanium version of the intel compiler
    # sets everything automatically
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS}")
  elseif (FEAST_CPU_TYPE STREQUAL "pentium4")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE2")
  elseif (FEAST_CPU_TYPE STREQUAL "nocona")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -xSSE3")
  elseif (FEAST_CPU_TYPE STREQUAL "itanium2")
    # no setting necessary, the itanium version of the intel compiler
    #sets everything automatically
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS}")

  # AMD CPUs
  elseif (FEAST_CPU_TYPE STREQUAL "amd486")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "k5")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "k6")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "athlon")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "athlonxp")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mIA32")
  elseif (FEAST_CPU_TYPE STREQUAL "opteron")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mSSE2")
  elseif (FEAST_CPU_TYPE STREQUAL "athlon64")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mSSE2")
  elseif (FEAST_CPU_TYPE STREQUAL "opteronx2")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mSSE3")
  elseif (FEAST_CPU_TYPE STREQUAL "turionx2")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mSSE3")
  elseif (FEAST_CPU_TYPE STREQUAL "barcelona")
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -mSSE4.1")

  # generic settings for all archs
  else ()
    message (STATUS "WARNING: Unsupported architecture/compiler combination found.")
    message (STATUS "Using generic optimisation flags. ")
  endif ()

  # check if compiler supports the given flags
  set (COMPILER_SUPPORTS_FLAG ON)
  include (CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("${FEAST_CXX_FLAGS}" COMPILER_SUPPORTS_FLAG)
  if (NOT COMPILER_SUPPORTS_FLAG)
    message (STATUS "##############################################################")
    message (STATUS "One of the flags ${FEAST_CXX_FLAGS} is apparently             ")
    message (STATUS "unsupported by the Intel Compiler v.${INTEL_VERSION_VALUE}.   ")
    message (FATAL_ERROR "##############################################################")
  endif ()

endif (FEAST_DEBUG_MODE)
