# FEAST BUILD_ID mechanism:
#       detection of parallel/serial
#       detection of MPI environment to be used
#
# This module sets the following variables:
#   FEAST_MPI_ENV_NAME
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock


# initialise local variables
set (FORCE_ERROR_MESSAGE OFF)


# determine and validate setting of "mpi" token,
# recall that this token is mandatory if tokens are used at all.
# in case this big switch statement goes wrong, the variable FORCE_ERROR is set to "YES"
if (NOT DISPLAY_HELP_ONLY)

  # search for tokens and react accordingly
  # (for those not speaking reg-expr: the match means that a certain token
  #  appears at the beginning, in the middle or at the end of BUILD_ID)
  if ( (BUILD_ID MATCHES "^openmpi-.+|.+-openmpi-.+|.+-openmpi$") OR
       (BUILD_ID MATCHES "^ompi-.+|.+-ompi-.+|.+-ompi$") )
    set (FEAST_MPI_ENV_NAME "openmpi")

  elseif (BUILD_ID MATCHES "^mpich2-.+|.+-mpich2-.+|.+-mpich2$")
    set (FEAST_MPI_ENV_NAME "mpich2")

  elseif (BUILD_ID MATCHES "^serial-.+|.+-serial-.+|.+-serial$")
    message (FATAL_ERROR ERROR: "Build ID '${BUILD_ID}': Serial mode not supported yet!")
    set (FEAST_MPI_ENV_NAME "serial")

  else ()
    # an unsupported "mpi" token has been encountered, so bail out
    set (FORCE_ERROR_MESSAGE ON)
  endif ()

endif ()


# in case something went wrong, or in case only a help message is requested,
# display some help. In order not to repeat ASCII art all over the place, some
# of the conditionals below are redundant.
if (FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "ERROR: Build-ID                                               ")
  message (STATUS "   \"${BUILD_ID}\"                                            ")
  message (STATUS "does not contain a supported MPI environment!                 ")
endif ()

if (DISPLAY_HELP_ONLY OR FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "Valid settings for token \"mpi\"                              ")
  message (STATUS "openmpi: use Open-MPI                                         ")
  message (STATUS "mpich2 : use MPICH2                                           ")
  message (STATUS "serial : do not use MPI at all                                ")
endif ()

if (FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "Something went wrong, please consult the output above.        ")
  message (STATUS "##############################################################")
  message (FATAL_ERROR "")
endif ()
