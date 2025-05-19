# Find-Module for SuperLU

include(FindPackageHandleStandardArgs)

find_library(superlu_LIBRARY NAMES superlu_dist)
find_path(superlu_INCLUDE_DIR NAMES superlu_defs.h)

find_package_handle_standard_args(superlu REQUIRED_VARS superlu_LIBRARY superlu_INCLUDE_DIR)

if(superlu_FOUND)
  mark_as_advanced(superlu_LIBRARY)
  mark_as_advanced(superlu_INCLUDE_DIR)
endif()

if(superlu_FOUND AND NOT TARGET SuperLU::SuperLU)
  add_library(SuperLU::SuperLU UNKNOWN IMPORTED GLOBAL)
  set_target_properties(SuperLU::SuperLU PROPERTIES IMPORTED_LOCATION ${superlu_LIBRARY})
  target_include_directories(SuperLU::SuperLU INTERFACE ${superlu_INCLUDE_DIR})

  if(NOT TARGET OpenMP::OpenMP_C)
    find_package(OpenMP)
  endif()
  target_link_libraries(SuperLU::SuperLU INTERFACE OpenMP::OpenMP_C)

  if(NOT TARGET BLAS::BLAS)
    find_package(BLAS)
  endif()
  target_link_libraries(SuperLU::SuperLU INTERFACE BLAS::BLAS)
endif()
