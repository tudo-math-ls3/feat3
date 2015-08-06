cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

if (FEAST_HAVE_ALGLIB)
  file(GLOB alglib-list "${FEAST_SOURCE_DIR}/thirdparty/alglib/cpp/src/*.cpp")

  ADD_LIBRARY(thirdparty-alglib STATIC ${alglib-list})
  if (FEAST_COTIRE)
    cotire(thirdparty-alglib)
  endif (FEAST_COTIRE)
  SET_PROPERTY(TARGET thirdparty-alglib PROPERTY COMPILE_FLAGS " -w")
endif (FEAST_HAVE_ALGLIB)

if (FEAST_HAVE_FPARSER)
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/fparser/")
  file(GLOB fparser-list "${FEAST_SOURCE_DIR}/thirdparty/fparser/*.cc")

  # The 'fpoptimizer.cc' can not be compiled with clang++, so we have to disable
  # it by adding the corresponding pre-processor define
  if (${FEAST_COMPILER_ID} STREQUAL "clang")
    ADD_DEFINITIONS(-DFP_NO_SUPPORT_OPTIMIZER)
  endif (${FEAST_COMPILER_ID} STREQUAL "clang")

  ADD_LIBRARY(thirdparty-fparser STATIC ${fparser-list})
  SET_PROPERTY(TARGET thirdparty-fparser PROPERTY COMPILE_FLAGS " -w")
endif (FEAST_HAVE_FPARSER)

if (FEAST_HAVE_UMFPACK)
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/SuiteSparse_config/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/AMD/Include/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/UMFPACK/Include/")

  file(GLOB suite_sparse-list "${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/SuiteSparse_config/*.c")
  file(GLOB amd-list "${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse/AMD/Source/*.c")
  file(GLOB umfpack-list "${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse//UMFPACK/Source/*.c")

  ADD_DEFINITIONS(-DNCHOLMOD -DDLONG -DNBLAS -DNTIMER)
  ADD_LIBRARY(thirdparty-umfpack STATIC ${suite_sparse-list} ${amd-list} ${umfpack-list})
  if (FEAST_COTIRE)
    cotire(thirdparty-umfpack)
  endif (FEAST_COTIRE)
  SET_PROPERTY(TARGET thirdparty-umfpack PROPERTY COMPILE_FLAGS " -w")
endif (FEAST_HAVE_UMFPACK)
