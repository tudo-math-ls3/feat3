# FEAST BUILD_ID mechanism:
#   big fat mapping of ORACLE compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (oracle)
#   CMAKE_CXX_COMPILER (sunCC)
#   FEAST_CXX_FLAGS_INTERNAL (string)
#
#
# Maintainer information:
#

#-dryrun
#-xtarget=native
# We must not use -xHOST even through it generally gives the best code,
# because it prevents cross-compilation. Consequently, if -xHOST is considered
# absolutely necessary for one particular test, BUILD_ID=manual must be used.
#
# The Oracle compiler suite follows the following conventions for creating
#   architecture-aware code, as outlined here:
#   http://software.intel.com/en-us/articles/performance-tools-for-software-developers-intel-compiler-options-for-sse-generation-and-processor-specific-optimizations/
#   or in the compiler man page.
#   Basically, there are two ways to set the target microarchitecture / processor category:
#     -mCODE generates code for a processor category identified by <CODE>,
#     -xCODE additionally generates optimisations specific to classes of Oracle CPUs
#   <CODE> can be one of the following:
#      AVX (-x only)
#      SSE4.2 (-x only)
#      SSE4.1
#      SSE3
#      SSE2 (default value for -m)
#      IA32 (-m only)
# This consequently means that for any non-Oracle processor, we must use
# the -m option. For all (not too ancient) Oracle processors, we should use -x,
# and this is essentially everything that needs to be set for optimised,
# architecture-aware builds.
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# set compiler and compiler name
set (CMAKE_CXX_COMPILER "sunCC")
set (FEAST_COMPILER_NAME "oracle")


# bail out for unsupported OSes
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message (STATUS "##############################################################")
  message (STATUS "Oracle compiler support is only available for Linux           ")
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


# ensure that Oracle compiler version is at least 12.1
message (STATUS "##############################################################")
message (STATUS "TODO: Implement Oracle version detection, because the output is")
message (STATUS "too weird for Dominik to understand (fundamentally different")
message (STATUS "versions between C, C++ and Fortran compilers in sunstudio 12.2!)")
message (STATUS "##############################################################")


# if compiler flags are not passed externally, determine our own
if (FEAST_CXX_FLAGS_INTERNAL STREQUAL "")

  # generic settings independent of arch and optimisation level
  set (FEAST_CXX_FLAGS_INTERNAL "")


  if (FEAST_DEBUG_MODE)
    # unoptimised settings for all archs
    # compat=5 is "according to the ANSI/ISO 1998 C++ standard as corrected in 2003"
    # which seems closest to chat we want
    # -verbose=template might be useful occasionally
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O0 -compat=5 -erroff=%none -g")

  else ()
    # optimised settings for all currently supported archs
    # please try to maintain the same order as in the buildid_arch module
    # Note: -fast expands pretty much everything we want, but unfortunately
    # also -xtarget=native. So we need to get rid of that again. Luckily,
    # this is easy: "You can override the values set by -fast by specifying
    # different values to the right of -fast on the command line."
    # -xlibmil inlines selected libmath functions (implied by -fast)
    # -xlibmopt uses optimised implemenation of libmath (implied by -fast)
    # -xprefetch=auto,explicit inserts prefetching instructions to improve overall instruction throughput
    set (FEAST_CXX_FLAGS_INTERNAL "-fast-xprefetch=auto,explicit")

    # please try to maintain the same order as in the buildid_arch module
    # Intel CPUs
    if (FEAST_CPU_TYPE STREQUAL "i486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=pentium -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumpro")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=pentium_pro -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium2")
      # strangely, pentium2 is not available
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium3")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=pentium3 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumm")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=pentium3 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4m")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=pentium4 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "coresolo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic64 -xarch=sse3 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "coreduo")
      # possible alternative: -xtarget=woodcrest?
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic64 -xchip=core2 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "penryn")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=penryn -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "nehalem")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=nehalem -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "westmere")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=nehalem -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium")
      message (STATUS "##############################################################")
      message (STATUS "Oracle compiler suite does not support itanium.               ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=pentium4")
    elseif (FEAST_CPU_TYPE STREQUAL "nocona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic64 -xchip=core2 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium2")
      message (STATUS "##############################################################")
      message (STATUS "Oracle compiler suite does not support itanium.               ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")

    # AMD CPUs
    elseif (FEAST_CPU_TYPE STREQUAL "amd486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "k5")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "k6")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "athlonxp")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "opteron")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=opteron -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon64")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=opteron -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "opteronx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=opteron -m64 -xarch=sse3")
    elseif (FEAST_CPU_TYPE STREQUAL "turionx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=opteron -m64 -xarch=sse3")
    elseif (FEAST_CPU_TYPE STREQUAL "barcelona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=barcelona -m64")

    # generic settings for all archs
    else ()
      message (STATUS "WARNING: Unsupported architecture/compiler combination found.")
      message (STATUS "Using generic optimisation flags. ")
      # filter out native part of -fast again
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -xtarget=generic64 -m64")

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
  message (STATUS "unsupported by the Oracle Compiler.                           ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
