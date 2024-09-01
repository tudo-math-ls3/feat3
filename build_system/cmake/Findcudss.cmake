# Find-Module for cudss

include(FindPackageHandleStandardArgs)

find_library(cudss_LIBRARY NAMES cudss)
find_path(cudss_INCLUDE_DIR NAMES cudss.h)

find_package_handle_standard_args(cudss REQUIRED_VARS cudss_LIBRARY cudss_INCLUDE_DIR)

if(cudss_FOUND)
  mark_as_advanced(cudss_LIBRARY)
  mark_as_advanced(cudss_INCLUDE_DIR)
endif()

if(cudss_FOUND AND NOT TARGET cuDSS::cuDSS)
add_library(cuDSS::cuDSS UNKNOWN IMPORTED GLOBAL)
set_target_properties(cuDSS::cuDSS PROPERTIES IMPORTED_LOCATION ${cudss_LIBRARY})
target_include_directories(cuDSS::cuDSS INTERFACE ${cudss_INCLUDE_DIR})
endif()
