# Find-Module for zlib

include(FindPackageHandleStandardArgs)

find_library(zlib_LIBRARY NAMES z)
find_path(zlib_INCLUDE_DIR NAMES zlib.h)

find_package_handle_standard_args(zlib REQUIRED_VARS zlib_LIBRARY zlib_INCLUDE_DIR)

if(zlib_FOUND)
  mark_as_advanced(zlib_LIBRARY)
  mark_as_advanced(zlib_INCLUDE_DIR)
endif()

if(zlib_FOUND AND NOT TARGET zlib::zlib)
  add_library(zlib::zlib UNKNOWN IMPORTED GLOBAL)
  set_target_properties(zlib::zlib PROPERTIES IMPORTED_LOCATION ${zlib_LIBRARY})
  target_include_directories(zlib::zlib INTERFACE ${zlib_INCLUDE_DIR})
endif()
