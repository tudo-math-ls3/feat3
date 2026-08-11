# Find-Module for FParser

include(FindPackageHandleStandardArgs)

find_path(FloatX_INCLUDE_DIR NAMES floatx.hpp)

find_package_handle_standard_args(FloatX REQUIRED_VARS FloatX_INCLUDE_DIR)

if(FloatX_FOUND)
  mark_as_advanced(FloatX_INCLUDE_DIR)
endif()

if(FloatX_FOUND AND NOT TARGET FloatX::FloatX)
  add_library(FloatX::FloatX INTERFACE IMPORTED GLOBAL)
  target_include_directories(FloatX::FloatX INTERFACE ${FloatX_INCLUDE_DIR})
endif()
