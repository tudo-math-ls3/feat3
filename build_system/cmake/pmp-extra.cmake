if(NOT TARGET pmp::pmp)
  FetchContent_GetProperties(pmp)

  file(GLOB pmp-sources ${pmp_SOURCE_DIR}/src/pmp/*.cpp ${pmp_SOURCE_DIR}/src/pmp/algorithms/*.cpp ${pmp_SOURCE_DIR}/src/pmp/io/*.cpp)
  file(GLOB pmp-headers ${pmp_SOURCE_DIR}/src/pmp/*.h ${pmp_SOURCE_DIR}/src/pmp/algorithms/*.h ${pmp_SOURCE_DIR}/src/pmp/io/*.h)

  add_library(feat-extern-pmp STATIC ${pmp-sources} ${pmp-headers})
  target_include_directories(feat-extern-pmp PUBLIC ${pmp_SOURCE_DIR}/src/)

  # NOTE(mmuegge): pmp includes the Eigen 3.4.0 headers, but does not link
  # against Eigen. This is intended.  pmp uses the type definitions of eigen to
  # support converting from Eigen vectors/matrices to pmp vectors/matrices, but
  # does not use the implementation of eigen in any way.
  target_include_directories(feat-extern-pmp PUBLIC ${pmp_SOURCE_DIR}/external/eigen-3.4.0)

  target_compile_features(feat-extern-pmp PRIVATE cxx_std_20)

  if(FEAT_HAVE_OMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(feat-extern-pmp PUBLIC OpenMP::OpenMP_CXX)
    endif()
  endif()

  add_library(pmp::pmp ALIAS feat-extern-pmp)
endif()
