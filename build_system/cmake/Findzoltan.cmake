# Find-Module for zoltan

include(FindPackageHandleStandardArgs)

find_library(zoltan_LIBRARY NAMES zoltan)
find_path(zoltan_INCLUDE_DIR NAMES zoltan.h)

find_package_handle_standard_args(zoltan REQUIRED_VARS zoltan_LIBRARY zoltan_INCLUDE_DIR)

if(zoltan_FOUND)
  mark_as_advanced(zoltan_LIBRARY)
  mark_as_advanced(zoltan_INCLUDE_DIR)
endif()

if(zoltan_FOUND AND NOT TARGET Zoltan::Zoltan)
  add_library(Zoltan::Zoltan UNKNOWN IMPORTED GLOBAL)
  set_target_properties(Zoltan::Zoltan PROPERTIES IMPORTED_LOCATION ${zoltan_LIBRARY})
  target_include_directories(Zoltan::Zoltan INTERFACE ${zoltan_INCLUDE_DIR})
endif()
