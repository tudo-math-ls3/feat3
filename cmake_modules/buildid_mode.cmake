# FEAST BUILD_ID mechanism: detection of compile mode
#
# This module sets the following variables:
#   FEAST_DEBUG_MODE (bool, cached)
#   FEAST_BUILD_MODE_NAME (string, opt|debug)
#
# \author Dominik Goeddeke
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
  if ( (BUILD_ID MATCHES "^debug-.+|.+-debug-.+|.+-debug$") OR
       (BUILD_ID MATCHES "^noopt-.+|.+-noopt-.+|.+-noopt$") )
    set (FEAST_DEBUG_MODE ON CACHE BOOL "" FORCE)
    set (FEAST_BUILD_MODE_NAME "debug")

  elseif (BUILD_ID MATCHES "^opt-.+|.+-opt-.+|.+-opt$")
    set (FEAST_DEBUG_MODE OFF CACHE BOOL "" FORCE)
    set (FEAST_BUILD_MODE_NAME "opt")

  else ()
    # an unsupported "mode" token has been encountered, so bail out
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
  message (STATUS "does not contain a supported build mode!                      ")
endif (FORCE_ERROR_MESSAGE)

if ( (DISPLAY_HELP_ONLY) OR (FORCE_ERROR_MESSAGE) )
  message (STATUS "##############################################################")
  message (STATUS "Valid settings for token \"mode\"                             ")
  message (STATUS "opt        : enable full compiler optimisations               ")
  message (STATUS "debug|noopt: no compiler optimisations                        ")
endif ()

if (FORCE_ERROR_MESSAGE)
  message (STATUS "##############################################################")
  message (STATUS "Something went wrong, please consult the output above.        ")
  message (FATAL_ERROR "##############################################################")
endif ()
