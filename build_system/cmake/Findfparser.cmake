# Find-Module for FParser

include(FindPackageHandleStandardArgs)

find_library(fparser_LIBRARY NAMES fparser)
find_path(fparser_INCLUDE_DIR NAMES fparser.hh)

find_package_handle_standard_args(fparser REQUIRED_VARS fparser_LIBRARY fparser_INCLUDE_DIR)

if(fparser_FOUND)
  mark_as_advanced(fparser_LIBRARY)
  mark_as_advanced(fparser_INCLUDE_DIR)
endif()

if(fparser_FOUND AND NOT TARGET fparser::fparser)
  add_library(fparser::fparser UNKNOWN IMPORTED GLOBAL)
  set_target_properties(fparser::fparser PROPERTIES IMPORTED_LOCATION ${fparser_LIBRARY})
  target_include_directories(fparser::fparser INTERFACE ${fparser_INCLUDE_DIR})
endif()
