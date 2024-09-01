# - Find ALGLIB
# Find the native ALGLIB includes and library
#
#  Alglib_INCLUDES    - where to find Alglib includes
#  Alglib_LIBRARIES   - List of libraries when using Alglib.
#  Alglib_FOUND       - True if Alglib found.

if (Alglib_INCLUDES)
  # Already in cache, be silent
  set (Alglib_FIND_QUIETLY TRUE)
endif (Alglib_INCLUDES)

find_path (Alglib_INCLUDES
    alglibinternal.h
    alglibmisc.h
    ap.h
    dataanalysis.h
    diffequations.h
    fasttransforms.h
    integration.h
    interpolation.h
    linalg.h
    optimization.h
    solvers.h
    specialfunctions.h
    statistics.h
    stdafx.h
    PATHS
    /usr/include/alglib/
    /usr/include/libalglib/
    /usr/local/include/alglib3/
    /usr/local/include/libalglib/
    )

find_library (Alglib_LIBRARIES NAMES alglib alglib3)

# handle the QUIETLY and REQUIRED arguments and set ALGLIB_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Alglib REQUIRED_VARS Alglib_LIBRARIES Alglib_INCLUDES)

if(Alglib_FOUND)
  mark_as_advanced (Alglib_LIBRARIES Alglib_INCLUDES)
endif()

if(Alglib_FOUND AND NOT TARGET Alglib::Alglib)
  add_library(Alglib::Alglib UNKNOWN IMPORTED GLOBAL)
  set_target_properties(Alglib::Alglib PROPERTIES IMPORTED_LOCATION ${Alglib_LIBRARIES})
  target_include_directories(Alglib::Alglib INTERFACE ${Alglib_INCLUDES})
endif()
