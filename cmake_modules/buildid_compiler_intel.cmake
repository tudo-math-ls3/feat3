# FEAST BUILD_ID mechanism:
#   big fat mapping of INTEL compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (intel)
#   CMAKE_CXX_COMPILER (icpc)
#   FEAST_CXX_FLAGS_INTERNAL (string)
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
# see issue 00035
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message (STATUS "##############################################################")
  message (STATUS "Intel compiler support is only available for Linux            ")
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


# ensure that Intel compiler version is at least 11.1
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set (INTEL_VERSION_STRING "")
  set (INTEL_VERSION_MAJOR 0)
  set (INTEL_VERSION_MINOR 0)
  # ask the compiler about its version by executing it
  exec_program("${CMAKE_CXX_COMPILER} -V 2>&1 | head -n 1" OUTPUT_VARIABLE INTEL_VERSION_STRING)
  # extract version from compiler output (because we know what the output looks like)
  string(REGEX REPLACE ".*Version[ \t]+([0-9]+)[.][0-9]+.*" "\\1" INTEL_VERSION_MAJOR "${INTEL_VERSION_STRING}")
  string(REGEX REPLACE ".*Version[ \t]+[0-9]+[.]([0-9]+).*" "\\1" INTEL_VERSION_MINOR "${INTEL_VERSION_STRING}")
  message (STATUS "Found Intel compiler suite v.${INTEL_VERSION_MAJOR}.${INTEL_VERSION_MINOR}")
  # and check for recent enough version
  if ("${INTEL_VERSION_MAJOR}.${INTEL_VERSION_MINOR}" VERSION_LESS "11.1")
    message (STATUS "##############################################################")
    message (STATUS "Intel Compiler Version ${INTEL_VERSION_MAJOR}.${INTEL_VERSION_MINOR} is not supported,")
    message (STATUS "please upgrade to a version better than 11.1.                 ")
    message (STATUS "##############################################################")
    message (FATAL_ERROR "")
  endif()
else ()
  message (STATUS "TODO: Implement Intel version detection for Non-Linux systems.")
endif ()

# if compiler flags are not passed externally, determine our own
if (FEAST_CXX_FLAGS_INTERNAL STREQUAL "")

  # ensure that OS is 64-bit
  if ( (CMAKE_SYSTEM_NAME STREQUAL "Linux") )
    exec_program("uname -m" OUTPUT_VARIABLE UNAME_OUTPUT)
    if (NOT UNAME_OUTPUT STREQUAL "x86_64")
      message (STATUS "##############################################################")
      message (STATUS "Only 64-bit operating systems are supported by default.       ")
      message (STATUS "You might want to switch to manual mode.                      ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")
    endif ()
   else ()
     message (STATUS "TODO: Implement 64-bit OS detection for Non-Linux systems.")
  endif ()

  # generic settings independent of arch and optimisation level
  set (FEAST_CXX_FLAGS_INTERNAL "-std=c++11")


  if (FEAST_DEBUG_MODE)
    # unoptimised settings for all archs
    # TODO: figure out how to set explicit 64-bit mode in this case
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O0 -Wall -g -Wp64 -mcmodel=large -Wshorten-64-to-32")
    # on Linux, use very very pedantic settings to detect (at run time)
    # bugs like out-of-bounds accesses. Warning: These settings result
    # in very slow execution, but they are occasionally worth it.
    # Since they already complain inside UMFPACK, leave them commented
    # out for the time being.
    # if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    #   set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -check all -debug -fp-stack-check -traceback -ftrapuv")
    # endif ()
    #
  else ()
    # optimised settings for all currently supported archs
    # please try to maintain the same order as in the buildid_arch module
    # Note: SSE2 is on by default, so we only have to specify what's better
    #       or worse
    # Note: 64 bit is enabled automatically by passing -xBLA
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O3 -mcmodel=large -g")

    # please try to maintain the same order as in the buildid_arch module
    # Intel CPUs
    if (FEAST_CPU_TYPE STREQUAL "i486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumpro")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium3")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumm")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse2")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4m")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse2")
    elseif (FEAST_CPU_TYPE STREQUAL "coresolo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse3")
    elseif (FEAST_CPU_TYPE STREQUAL "coreduo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse3")
    elseif (FEAST_CPU_TYPE STREQUAL "penryn")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse4.1")
    elseif (FEAST_CPU_TYPE STREQUAL "nehalem")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse4.2")
    elseif (FEAST_CPU_TYPE STREQUAL "westmere")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse4.2")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium")
      # no setting necessary, the itanium version of the intel compiler
      # sets everything automatically
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse2")
    elseif (FEAST_CPU_TYPE STREQUAL "nocona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xsse3")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium2")
      # no setting necessary, the itanium version of the intel compiler
      #sets everything automatically
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")

      # AMD CPUs
    elseif (FEAST_CPU_TYPE STREQUAL "amd486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "k5")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "k6")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "athlonxp")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -mia32")
    elseif (FEAST_CPU_TYPE STREQUAL "opteron")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -msse2")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon64")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -msse2")
    elseif (FEAST_CPU_TYPE STREQUAL "opteronx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -msse3")
    elseif (FEAST_CPU_TYPE STREQUAL "turionx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -msse3")
    elseif (FEAST_CPU_TYPE STREQUAL "barcelona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -msse4.1")

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
  message (STATUS "One of the flags ${FEAST_CXX_FLAGS_INTERNAL} is apparently    ")
  message (STATUS "unsupported by the Intel Compiler v.${INTEL_VERSION_MAJOR}.${INTEL_VERSION_MINOR}")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
