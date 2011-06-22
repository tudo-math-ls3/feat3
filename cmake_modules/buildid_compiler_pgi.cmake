# FEAST BUILD_ID mechanism:
#   big fat mapping of PGI compiler to various modes and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (pgi)
#   CMAKE_CXX_COMPILER (pgCC)
#   FEAST_CXX_FLAGS_INTERNAL (string)
#
# Maintainer information:
# Run pgCC -help for lots of details
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# set compiler and compiler name
set (CMAKE_CXX_COMPILER "pgCC")
set (FEAST_COMPILER_NAME "pgi")

# bail out for unsupported OSes
# see issue 00035
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message (STATUS "##############################################################")
  message (STATUS "PGI compiler support is only available for Linux              ")
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

# ensure that PGI compiler version is at least
if ( CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set (PGI_VERSION_STRING "")
  set (PGI_VERSION_MAJOR "")
  set (PGI_VERSION_MINOR "")
  # ask the compiler about its version by executing it
  exec_program("${CMAKE_CXX_COMPILER} -V 2>&1 | head -n 2" OUTPUT_VARIABLE PGI_VERSION_STRING)
  # extract version from compiler output (because we know what the output looks like)
  # extract version from compiler output (because we know what the output looks like)
  string(REGEX REPLACE ".*Version[ \t]+([0-9]+)[.][0-9]+.*" "\\1" PGI_VERSION_MAJOR "${PGI_VERSION_STRING}")
  string(REGEX REPLACE ".*Version[ \t]+[0-9]+[.]([0-9]+).*" "\\1" PGI_VERSION_MINOR "${PGI_VERSION_STRING}")
  message (STATUS "Found Intel compiler suite v.${PGI_VERSION_MAJOR}.${PGI_VERSION_MINOR}")
  # and check for recent enough version
  if ("${PGI_VERSION_MAJOR}.${PGI_VERSION_MINOR}" VERSION_LESS "10.0")
    message (STATUS "##############################################################")
    message (STATUS "PGI Compiler Version ${PGI_VERSION_MAJOR}.${PGI_VERSION_MINOR} is not supported,")
    message (STATUS "please upgrade to a version better than 10.0.                 ")
    message (STATUS "##############################################################")
    message (FATAL_ERROR "")
  endif()
else ()
  message (STATUS "TODO: Implement PGI version detection for Non-Linux systems.")
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
  set (FEAST_CXX_FLAGS_INTERNAL "")


  if (FEAST_DEBUG_MODE)
    # unoptimised settings for all archs
    # -Mbounds: insert bound checking code
    # -traceback: run time crash tracebacks
    # -Minform=warn: -Wall
    # --remarks: let the compiler be a little bit verbose
    # -A: ANSI C++
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -O0 -g -Mbounds -traceback -Minform=warn --remarks -A")

  else ()
    # optimised settings for all currently supported archs
    set (FEAST_CXX_FLAGS_INTERNAL "-fast -O4")

    # please try to maintain the same order as in the buildid_arch module
    # Intel CPUs
    if (FEAST_CPU_TYPE STREQUAL "i486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=px -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=p5 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumpro")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=p6 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=p6 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium3")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=piii -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentiumm")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=p6 -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4m")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=piv -m32")
    elseif (FEAST_CPU_TYPE STREQUAL "coresolo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=p7 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "coreduo")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=core2 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "penryn")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=penryn -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "nehalem")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=nehalem -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "westmere")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=nehalem -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")
    elseif (FEAST_CPU_TYPE STREQUAL "pentium4")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=piv -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "nocona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=p7 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "itanium2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL}")

    # AMD CPUs
    elseif (FEAST_CPU_TYPE STREQUAL "amd486")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=px")
    elseif (FEAST_CPU_TYPE STREQUAL "k5")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=px")
    elseif (FEAST_CPU_TYPE STREQUAL "k6")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=px")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=athlon")
    elseif (FEAST_CPU_TYPE STREQUAL "athlonxp")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=athlonxp")
    elseif (FEAST_CPU_TYPE STREQUAL "opteron")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=k8 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "athlon64")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=k8 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "opteronx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=k8 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "turionx2")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=k8 -m64")
    elseif (FEAST_CPU_TYPE STREQUAL "barcelona")
      set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -tp=barcelona -m64")

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
  message (STATUS "unsupported by the PGI Compiler v.${PGI_VERSION_MAJOR}.${PGI_VERSION_MINOR}.")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
