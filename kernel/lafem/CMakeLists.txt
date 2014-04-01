cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

ADD_SUBDIRECTORY( arch )

# list of lafem source files
SET (kernel-lafem-list
  memory_pool.cpp
  )

if (FEAST_HAVE_UMFPACK)
  SET(kernel-lafem-list ${kernel-lafem-list} umfpack.cpp)
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/SuiteSparse_config/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/AMD/Include/")
  INCLUDE_DIRECTORIES("${FEAST_SOURCE_DIR}/thirdparty/UMFPACK/Include/")
endif (FEAST_HAVE_UMFPACK)

ADD_LIBRARY(kernel-lafem STATIC ${kernel-lafem-list})
TARGET_LINK_LIBRARIES(kernel-lafem kernel-lafem-arch kernel-util kernel-adjacency)

if (FEAST_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(kernel-lafem-cuda
    memory_pool.cu
    STATIC
  )
  CUDA_ADD_CUBLAS_TO_TARGET(kernel-lafem-cuda)
  TARGET_LINK_LIBRARIES(kernel-lafem kernel-lafem-cuda)
endif(FEAST_BACKENDS_CUDA)

# list of lafem tests
SET (test_list
  bicgstab-test
  block_to_scalar-test
  dense_matrix-test
  dense_vector-test
  dense_vector_blocked-test
  meta_filter-test
  meta_matrix-apply-test
  meta_mirror-test
  meta_vector-axpy-test
  meta_vector-comp_prod-test
  meta_vector-comp_invert-test
  meta_vector-dot-norm2-test
  meta_vector-scale-test
  pointstar_factory-test
  richardson-test
  sparse_matrix_conversion-test
  sparse_matrix_banded-test
  sparse_matrix_coo-test
  sparse_matrix_csr-test
  sparse_matrix_ell-test
  sparse_vector-test
  transposition-test
  umfpack-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feast test_system)

  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/lafem/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAST_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAST_VALGRIND)

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/lafem/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
  if (FEAST_VALGRIND)
    SET_PROPERTY(TEST ${test}_main PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_main PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAST_VALGRIND)

  ADD_TEST(${test}_cuda ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${CUDA_MEMCHECK_EXE} ${FEAST_BINARY_DIR}/kernel/lafem/${test} cuda cudadevicereset)
  SET_PROPERTY(TEST ${test}_cuda PROPERTY LABELS "cuda")
  if (FEAST_CUDAMEMCHECK AND FEAST_BACKENDS_CUDA)
    SET_PROPERTY(TEST ${test}_cuda PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAST_CUDAMEMCHECK AND FEAST_BACKENDS_CUDA)
ENDFOREACH(test)

# add all tests to lafem_tests
ADD_CUSTOM_TARGET(lafem_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests lafem_tests)
