# Find-Module for Triangle

include(FindPackageHandleStandardArgs)

find_library(triangle_LIBRARY NAMES triangle)
find_path(triangle_INCLUDE_DIR NAMES triangle.h)

find_package_handle_standard_args(triangle REQUIRED_VARS triangle_LIBRARY triangle_INCLUDE_DIR)

if(triangle_FOUND)
  mark_as_advanced(triangle_LIBRARY)
  mark_as_advanced(triangle_INCLUDE_DIR)
endif()

if(triangle_FOUND AND NOT TARGET triangle::triangle)
  add_library(triangle::triangle UNKNOWN IMPORTED GLOBAL)
  set_target_properties(triangle::triangle PROPERTIES IMPORTED_LOCATION ${triangle_LIBRARY})
  target_include_directories(triangle::triangle INTERFACE ${triangle_INCLUDE_DIR})
endif()
