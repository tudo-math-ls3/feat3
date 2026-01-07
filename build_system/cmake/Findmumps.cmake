# Find-Module for FParser

include(FindPackageHandleStandardArgs)

find_library(mumps_LIBRARY NAMES dmumps)
find_library(mumps_common_LIBRARY NAMES mumps_common)
find_library(pord_LIBRARY NAMES pord)
find_library(scalapack_LIBRARY NAMES scalapack)
find_path(mumps_INCLUDE_DIR NAMES dmumps_c.h dmumps_root.h dmumps_struc.h)

find_package_handle_standard_args(mumps REQUIRED_VARS mumps_LIBRARY mumps_common_LIBRARY pord_LIBRARY scalapack_LIBRARY mumps_INCLUDE_DIR)

if(mumps_FOUND)
  mark_as_advanced(mumps_LIBRARY)
  mark_as_advanced(mumps_INCLUDE_DIR)
endif()

if(mumps_FOUND AND NOT TARGET mumps::mumps)
  add_library(mumps::scalapack UNKNOWN IMPORTED GLOBAL)
  set_target_properties(mumps::scalapack PROPERTIES IMPORTED_LOCATION ${scalapack_LIBRARY})

  add_library(mumps::common UNKNOWN IMPORTED GLOBAL)
  set_target_properties(mumps::common PROPERTIES IMPORTED_LOCATION ${mumps_common_LIBRARY})

  add_library(mumps::pord UNKNOWN IMPORTED GLOBAL)
  set_target_properties(mumps::pord PROPERTIES IMPORTED_LOCATION ${pord_LIBRARY})

  add_library(mumps::mumps UNKNOWN IMPORTED GLOBAL)
  set_target_properties(mumps::mumps PROPERTIES IMPORTED_LOCATION ${mumps_LIBRARY})
  target_include_directories(mumps::mumps INTERFACE ${mumps_INCLUDE_DIR})

  target_link_libraries(mumps::mumps INTERFACE mumps::common mumps::scalapack gfortran)
  target_link_libraries(mumps::common INTERFACE mumps::pord)

  if(FEAT_HAVE_MKL)
    find_package(MKL CONFIG REQUIRED)
    target_link_libraries(mumps::mumps INTERFACE MKL::MKL)
  else()
    find_package(BLAS)
    target_link_libraries(mumps::mumps INTERFACE BLAS::BLAS)
  endif()

  if(FEAT_HAVE_MPI)
    # find_package(MPI COMPONENTS Fortran) depends on a working fortran compiler
    # We thus enable the language locally to avoid having to link against MPI manually
    enable_language(Fortran)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    target_link_libraries(mumps::mumps INTERFACE MPI::MPI_Fortran)
  endif()
endif()
