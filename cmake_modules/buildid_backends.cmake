# FEAST BUILD_ID mechanism: detection of additional backends
#
# This module sets the following variables:
#   FEAST_BACKEND_CUDA (bool, cached)
#   FEAST_BACKEND_MKL  (bool, cached)
#   FEAST_GMP (bool, cached)
#   FEAST_BACKENDS_FLAGS (string, cached)
#
# \author Dirk Ribbrock


# initialise local variables
set (FORCE_ERROR_MESSAGE OFF)
set (FEAST_BACKENDS_FLAGS "none")

# determine and validate setting of "mode" token,
# recall that this token is mandatory if tokens are used at all.
# in case this big switch statement goes wrong, the variable FORCE_ERROR is set to "YES"
if (NOT DISPLAY_HELP_ONLY)

  # search for tokens and react accordingly
  # (for those not speaking reg-expr: the match means that a certain token
  #  appears at the beginning, in the middle or at the end of BUILD_ID)
  if (BUILD_ID MATCHES "^cuda-.+|.+-cuda-.+|.+-cuda$")
    set (FEAST_BACKENDS_CUDA ON CACHE BOOL "" FORCE)
    if (FEAST_BACKENDS_FLAGS MATCHES "none")
      set (FEAST_BACKENDS_FLAGS "cuda")
    else ()
      set (FEAST_BACKENDS_FLAGS "${FEAST_BACKENDS_FLAGS}-cuda")
    endif (FEAST_BACKENDS_FLAGS MATCHES "none")
  endif ()

  if (BUILD_ID MATCHES "^mkl-.+|.+-mkl-.+|.+-mkl$")
    set (FEAST_BACKENDS_MKL ON CACHE BOOL "" FORCE)
    if (FEAST_BACKENDS_FLAGS MATCHES "none")
      set (FEAST_BACKENDS_FLAGS "mkl")
    else ()
      set (FEAST_BACKENDS_FLAGS "${FEAST_BACKENDS_FLAGS}-mkl")
    endif (FEAST_BACKENDS_FLAGS MATCHES "none")
    # finally, pass all compiler flags to cmake
    if (FEAST_BUILD_MODE_NAME STREQUAL "debug")
      set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMKL_ILP64")
    elseif (FEAST_BUILD_MODE_NAME STREQUAL "opt")
      set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DMKL_ILP64")
    endif ()
    # and store them in our variable as well for proper screen output
    set (FEAST_CXX_FLAGS "${FEAST_CXX_FLAGS} -DMKL_ILP64")

  endif ()

  if (BUILD_ID MATCHES "^gmp-.+|.+-gmp-.+|.+-gmp$")
    set (FEAST_GMP ON CACHE BOOL "" FORCE)
    if (FEAST_BACKENDS_FLAGS MATCHES "none")
      set (FEAST_BACKENDS_FLAGS "gmp")
    else ()
      set (FEAST_BACKENDS_FLAGS "${FEAST_BACKENDS_FLAGS}-gmp")
    endif (FEAST_BACKENDS_FLAGS MATCHES "none")
  endif ()

endif ()


# In case a help message is requested, display it.
if (DISPLAY_HELP_ONLY)
  message (STATUS "##############################################################")
  message (STATUS "Valid settings for token \"backend\"                             ")
  message (STATUS "cuda   : enable cuda support                     ")
  message (STATUS "mkl    : enable mkl support                     ")
  message (STATUS "gmp    : enable gmp support                     ")

endif ()
