# FEAST BUILD_ID mechanism: detection of additional backends
#
# This module sets the following variables:
#   FEAST_BACKEND_CUDA (bool, cached)
#
# \author Dirk Ribbrock


# initialise local variables
set (FORCE_ERROR_MESSAGE OFF)

# determine and validate setting of "mode" token,
# recall that this token is mandatory if tokens are used at all.
# in case this big switch statement goes wrong, the variable FORCE_ERROR is set to "YES"
if (NOT DISPLAY_HELP_ONLY)

  # search for tokens and react accordingly
  # (for those not speaking reg-expr: the match means that a certain token
  #  appears at the beginning, in the middle or at the end of BUILD_ID)
  if (BUILD_ID MATCHES "^cuda-.+|.+-cuda-.+|.+-cuda$")
    set (FEAST_BACKENDS_CUDA ON CACHE BOOL "" FORCE)

  endif ()

endif ()


# In case a help message is requested, display it.
if (DISPLAY_HELP_ONLY)
  message (STATUS "##############################################################")
  message (STATUS "Valid settings for token \"backend\"                             ")
  message (STATUS "cuda   : enable cuda support                     ")

endif ()
