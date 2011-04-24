# The cmake modules
#   buildid (this file)
#   buildid_help
#   buildid_mode
#   buildid_mpi
#   buildid_arch
#   buildid_compiler
#   buildid_compiler_gnu|intel|...
# implement FEAST's environment-aware build system. These modules are
# included appropriately by this file, and set the following cmake variables:
# see issue 00035
# TODO keep this list up to date
#   buildid (this file):
#     BUILD_ID (potentially corrected/reordered)
#   buildid_help:
#     nothing
#   buildid_mode:
#     FEAST_DEBUG_MODE (bool, cached)
#     FEAST_BUILD_MODE_NAME (string, debug|release|minsizerel|relwithdebinfo)
#     CMAKE_BUILD_TYPE (to cmake-equivalents of FEAST_BUILD_MODE_NAME)
#   buildid_mpi:
#     FEAST_MPI_ENV_NAME (string, serial|openmpi|mpich2)
#   buildid_arch:
#     FEAST_CPU_TYPE (string, lots of valid settings like nehalem)
#   buildid_compiler:
#     FEAST_COMPILER_NAME (string, gnu|intel|oracle|open64|pgi)
#     CMAKE_CXX_COMPILER (string)
#     FEAST_COMPILER_FLAGS (string, also passed automatically to add_definitions())
#
# Documentation on using the build system is available in the  "provide help"
# section below (actually implemented in the module buildid_help.cmake).
# This file implements the high-level logic to distinguish between the MANUAL,
# default and build-id-token based ways to configure FEAST.
#
# Some useful information for developers:
#   To get an idea of how things work, read the help file and this file, in particular
#      the logic it includes.
#   Adding a new MPI environment only requires changes to buildid_mpi.cmake,
#      and also to the code that overwrites the compiler with the MPI wrapper command
#      which we always use (in module buildid_compiler.cmake)
#   Adding a new compiler only requires changes to buildid_compiler.cmake and
#     creating a new buildid_compiler_foo.cmake file.
#   Adding a new target architecture requires changes in at least two files:
#     buildid_arch.cmake: proper detection code for the new arch (in the two macros
#                         provided there)
#     buildid_compiler_*.cmake: compiler flags for the new architecture (for as many
#                               compilers that have been tested)
#     The following procedure is highly recommended:
#       Compile for the predecessor architecture (e.g., for nehalem if the new
#       architecture is westmere), check if everything works, and then test additional
#       or replacement flags (e.g. for this example and the Intel compiler, switch from
#       SSE4.1 support to SSE4.2).
#
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock
#
# global TODO: cache everything that's necessary for proper ccmake support.



# initialise some variables

# don't cache BUILD_ID because otherwise, "cmake ." doesn't work and once
# a build-ID has been given, it must always be given, which is invonvenient.
set (BUILD_ID "" CACHE STRING "Use specific build id, or MANUAL, or empty for autodetection, or HELP to see valid settings")
set (DISPLAY_HELP_ONLY OFF)
set (FORCE_ERROR_MESSAGE OFF)


# display welcome message
message (STATUS "##############################################################")
message (STATUS "FEAST BUILD_ID mechanism")
message (STATUS "Help: \"cmake -D BUILD_ID=HELP /path/to/sources\"")


# provide help if requested
if ( (BUILD_ID STREQUAL "HELP") OR (BUILD_ID STREQUAL "help") )
  set (DISPLAY_HELP_ONLY ON)
  include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_help.cmake )
  message (STATUS "##############################################################")
  message (STATUS "Bailing out because \"help\" is not intended to actually      ")
  message (STATUS "configure anything.                                           ")
  message (STATUS "Ignore whatever cmake says below.                             ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")




# manual mode, i.e., skip build-ID mechanism completely
elseif ( (BUILD_ID STREQUAL "MANUAL") OR (BUILD_ID STREQUAL "manual") )
  message (STATUS "Manual cmake configure mode detected. Use at your own risk.")


# otherwise, either assume BUILD_ID really contains a build-ID,
# or set a default one in case BUILD_ID is not set;
# in any case, parse BUILD_ID
else ()

  # first of all, purge some settings, that cmake could grab from the shell environment
  set (CMAKE_CXX_COMPILER_ARG1 "")
  set (CMAKE_CXX_FLAGS "")
  set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

  # create a default build-ID
  if (BUILD_ID STREQUAL "")
    # tell the build system to guess the arch later on
    set (GUESS_ARCH ON)

    # and define the default ID for different operating systems
    # cmake-doc: the value of CMAKE_SYSTEM_NAME is the output of
    # uname -s on Unix-like systems (Linux, SunOS, Darwin, cygwin, ...)
    # http://cmake.org/cmake/help/cmake-2-8-docs.html#variable:CMAKE_SYSTEM_NAME
    if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
      # default build ID for Linux systems is
      #   - parallel build with OpenMPI
      #   - full compiler optimisations
      #   - GNU compiler suite
      #   - autodetection of underlying hardware
      # see issue 00035
      set (BUILD_ID opt-openmpi-gnu)

    elseif (CMAKE_SYSTEM_NAME STREQUAL "SunOS")
      # same default
      # see issue 00035
      set (BUILD_ID opt-openmpi-gnu)

    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      # TODO implement
      # see issue 00035
      message (STATUS "##############################################################")
      message (STATUS "ERROR: Default build-ID not implemented yet for Darwin.       ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")


    elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      # TODO implement
      message (STATUS "##############################################################")
      message (STATUS "ERROR: Default build-ID not implemented yet for Windows.      ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")


    else ()
      message (STATUS "##############################################################")
      message (STATUS "ERROR: Unsupported operating system found.                    ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")

    endif ()

  endif ()

  # parse "mode" token
  include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_mode.cmake )

  # parse "mpi" token
  include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_mpi.cmake )

  # parse arch token
  include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_arch.cmake )

  # parse compiler token
  include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler.cmake )

  # ok, all fine so far, so re-map BUILD_ID to our standard ordering of tokens
  # and print out a summary
  message (STATUS "##############################################################")
  message (STATUS "Supplied Build-ID: ${BUILD_ID}")
  set (BUILD_ID "${FEAST_BUILD_MODE_NAME}-${FEAST_MPI_ENV_NAME}-${FEAST_COMPILER_NAME}-${FEAST_CPU_TYPE}")
  message (STATUS "Debug mode         : ${FEAST_DEBUG_MODE}")
  if (NOT FEAST_MPI_ENV_NAME STREQUAL "serial")
    message (STATUS "MPI environment    : ${FEAST_MPI_ENV_NAME}")
  else ()
    message (STATUS "MPI environment    : not used, serial mode selected")
  endif ()
  message (STATUS "Architecture       : ${FEAST_CPU_TYPE}")
  message (STATUS "MPI wrapper command: ${CMAKE_CXX_COMPILER}")
  message (STATUS "Compiler           : ${FEAST_COMPILER_NAME}")
  message (STATUS "Flags              : ${FEAST_CXX_FLAGS}")
  message (STATUS "Using Build ID     : ${BUILD_ID}")
  message (STATUS "##############################################################")

endif ()
