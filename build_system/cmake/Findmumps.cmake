# Find-Module for FParser

include(FindPackageHandleStandardArgs)

find_library(mumps_LIBRARY NAMES dmumps)
find_path(mumps_INCLUDE_DIR NAMES dmumps_c.h dmumps_root.h dmumps_struc.h)

find_package_handle_standard_args(mumps REQUIRED_VARS mumps_LIBRARY mumps_INCLUDE_DIR)

if(mumps_FOUND)
  mark_as_advanced(mumps_LIBRARY)
  mark_as_advanced(mumps_INCLUDE_DIR)
endif()

if(mumps_FOUND AND NOT TARGET mumps::mumps)
  add_library(mumps::mumps UNKNOWN IMPORTED GLOBAL)
  set_target_properties(mumps::mumps PROPERTIES IMPORTED_LOCATION ${mumps_LIBRARY})
  target_include_directories(mumps::mumps INTERFACE ${mumps_INCLUDE_DIR})
endif()
