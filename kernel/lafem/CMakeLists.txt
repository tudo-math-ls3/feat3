cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

ADD_SUBDIRECTORY( arch )

# list of lafem source files
# currently empty, but neccessary for sub library dependencies
# TODO remove dummy as soon as possible
SET (kernel-lafem-list
  empty.cpp
  )

ADD_LIBRARY(kernel-lafem STATIC ${kernel-lafem-list})
TARGET_LINK_LIBRARIES(kernel-lafem kernel-lafem-arch kernel-util kernel-adjacency)

if (FEAT_HAVE_CUDA)
  TARGET_LINK_LIBRARIES(kernel-lafem kernel-util-cuda)
endif(FEAT_HAVE_CUDA)

# list of lafem tests
SET (test_list
  bicgstab-test
  bicgstab-meta-test
  container_main_wrapper-test
  dense_matrix-test
  dense_vector-test
  dense_vector_blocked-test
  matrix_mult-test
  meta_filter-test
  meta_matrix-apply-test
  meta_matrix-io-test
  meta_mirror-test
  meta_to_scalar-test
  meta_vector-axpy-test
  meta_vector-comp_prod-test
  meta_vector-comp_invert-test
  meta_vector-dot-norm2-test
  meta_vector-io-test
  meta_vector-scale-test
  pointstar_factory-test
  richardson-test
  slip_filter-test
  sparse_matrix_conversion-test
  sparse_matrix_banded-test
  sparse_matrix_coo-test
  sparse_matrix_csr-test
  sparse_matrix_ell-test
  sparse_matrix_bcsr-test
  sparse_vector-test
  sparse_vector_blocked-test
  tuning-test
  unit_filter-test
  unit_filter_blocked-test
  vector_mirror-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)

  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/lafem/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/lafem/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_main PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_main PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  ADD_TEST(${test}_cuda ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${CUDA_MEMCHECK_EXE} ${FEAT_BINARY_DIR}/kernel/lafem/${test} cuda)
  SET_PROPERTY(TEST ${test}_cuda PROPERTY LABELS "cuda")
  if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    SET_PROPERTY(TEST ${test}_cuda PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    SET_PROPERTY(TEST ${test}_cuda PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
ENDFOREACH(test)

# add all tests to lafem_tests
ADD_CUSTOM_TARGET(lafem_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests lafem_tests)
