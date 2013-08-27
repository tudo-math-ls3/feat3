# FEAST BUILD_ID mechanism:
#   big fat mapping of compiler to various modes, mpi envs and architectures
#
# This module sets the following variables:
#   FEAST_COMPILER_NAME (string, clang|gnu|intel|oracle|open64|pgi)
#   CMAKE_CXX_COMPILER (string, actual compiler to be used)
#   FEAST_CXX_FLAGS (list of compiler flags to be used, piped to cmake
#                    via add_definitions() later)
#
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# initialise local variables
set (FORCE_ERROR_MESSAGE OFF)


# determine and validate setting of "compiler" token,
# recall that this token is mandatory if tokens are used at all.
if (NOT DISPLAY_HELP_ONLY)

  # this is unfortunately needed to use command line values, because
  # in cmake, an unset variable is not strequal ""
  if (NOT FEAST_CXX_FLAGS STREQUAL "")
   set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS}")
  else ()
   set (FEAST_CXX_FLAGS_INTERNAL "")
  endif ()

  # map ambiguous compiler token to unique identifier,
  # this also checks if a supported compiler is requested,
  # and include the corresponding compiler file which in
  # turn contains the big fat mappings
  if ( (BUILD_ID MATCHES "^clang-.+|.+-clang-.+|.+-clang$") OR
       (BUILD_ID MATCHES "^llvm-.+|.+-llvm-.+|.+-llvm$") )
    include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler_clang.cmake )

  elseif ( (BUILD_ID MATCHES "^gnu-.+|.+-gnu-.+|.+-gnu$") OR
       (BUILD_ID MATCHES "^gcc-.+|.+-gcc-.+|.+-gcc$") )
    include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler_gnu.cmake )

  elseif ( (BUILD_ID MATCHES "^intel-.+|.+-intel-.+|.+-intel$") OR
       (BUILD_ID MATCHES "^icc-.+|.+-icc-.+|.+-icc$") )
    include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler_intel.cmake )

  elseif ( (BUILD_ID MATCHES "^sunstudio-.+|.+-sunstudio-.+|.+-sunstudio$") OR
           (BUILD_ID MATCHES "^oraclestudio-.+|.+-oraclestudio-.+|.+-oraclestudio$") OR
           (BUILD_ID MATCHES "^oracle-.+|.+-oracle-.+|.+-oracle$") )
    include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler_oracle.cmake )

  elseif (BUILD_ID MATCHES "^open64-.+|.+-open64-.+|.+-open64$")
    include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler_open64.cmake )

  elseif (BUILD_ID MATCHES "^pgi-.+|.+-pgi-.+|.+-pgi$")
    include ( ${FEAST_SOURCE_DIR}/cmake_modules/buildid_compiler_pgi.cmake )

  else ()
    # unsupported compiler found
    set (FORCE_ERROR_MESSAGE ON)
  endif ()


  # for parallel builds, use MPI wrapper commands instead of the compiler directly
  if (NOT FEAST_MPI_ENV_NAME STREQUAL "serial")
    # unfortunately, MPI wrapper commands are not standardised across implementations,
    # so try a few of the common ones and use the first one that matches.
    #
    # cmake syntax lesson (for the find_program() command which is basically a wrapper
    # around "which" on linux): quoting from cmake docs version 2.8.4
    #   If the full path to a file is found the result is stored in the variable and
    #   the search will not be repeated unless the variable is cleared. If nothing is
    #   found, the result will be <VAR>-NOTFOUND, and the search will be attempted
    #   again the next time find_file is invoked with the same variable.
    # end quote
    set (FINDPROG "FINDPROG-NOTFOUND")
    find_program (FINDPROG NAMES mpic++ mpiCC mpiCXX mpicxx)
    if (FINDPROG STREQUAL "FINDPROG-NOTFOUND")
      message (STATUS "##############################################################")
      message (STATUS "ERROR: MPI compiler wrapper not found.                        ")
      message (STATUS "       Please make sure that the compiler you selected is     ")
      message (STATUS "       available in your environment, for instance by loading ")
      message (STATUS "       the corresponding modules.                             ")
      message (STATUS "##############################################################")
      message (FATAL_ERROR "")
    else ()
      # strip path information
      set (PROGNAME)
      get_filename_component(PROGNAME ${FINDPROG} NAME)
      # and use the result as the compiler
      set (CMAKE_CXX_COMPILER ${PROGNAME})
    endif (FINDPROG STREQUAL "FINDPROG-NOTFOUND")

  endif (NOT FEAST_MPI_ENV_NAME STREQUAL "serial")


  # finally, pass all compiler flags to cmake
  if (FEAST_BUILD_MODE_NAME STREQUAL "debug")
    set (CMAKE_CXX_FLAGS_DEBUG "${FEAST_CXX_FLAGS_INTERNAL}")
  elseif (FEAST_BUILD_MODE_NAME STREQUAL "minsizerel")
    set (CMAKE_CXX_FLAGS_MINSIZEREL "${FEAST_CXX_FLAGS_INTERNAL}")
  elseif (FEAST_BUILD_MODE_NAME STREQUAL "opt")
    set (CMAKE_CXX_FLAGS_RELEASE "${FEAST_CXX_FLAGS_INTERNAL}")
  elseif (FEAST_BUILD_MODE_NAME STREQUAL "relwithdebinfo")
    # TODO This should really be done inside the compiler-specific module
    set (FEAST_CXX_FLAGS_INTERNAL "${FEAST_CXX_FLAGS_INTERNAL} -g")
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${FEAST_CXX_FLAGS_INTERNAL}")
  endif ()

  # and store them in our variable as well for proper screen output
  set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS_INTERNAL}")

endif (NOT DISPLAY_HELP_ONLY)


# in case something went wrong, or in case only a help message is requested,
# display some help. In order not to repeat ASCII art all over the place, some
# of the conditionals below are redundant.
if (FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "ERROR: Build-ID                                               ")
  message (STATUS "   \"${BUILD_ID}\"                                            ")
  message (STATUS "does not contain a supported compiler!                        ")
endif (FORCE_ERROR_MESSAGE)

if (DISPLAY_HELP_ONLY OR FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "Valid settings for token \"compiler\"                         ")
  message (STATUS "clang : LLVM/Clang c++ frontend                               ")
  message (STATUS "gnu   : GNU compiler suite                                    ")
  message (STATUS "intel : Intel compiler suite                                  ")
  message (STATUS "oracle: SunStudio / OracleStudio compiler suite               ")
  message (STATUS "open64: Open64 compiler suite                                 ")
  message (STATUS "pgi   : PGI compiler suite                                    ")
  message (STATUS "--------------------------------------------------------------")
  message (STATUS "Additional options:                                           ")
  message (STATUS "coverage: Enable code coverage                                ")
endif ()

if (FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "Something went wrong, please consult the output above.        ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")

endif ()
