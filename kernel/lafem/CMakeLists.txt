# enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory( arch )

# list of lafem source files
set(kernel-lafem-list
  empty_lafem.cpp
  )

if (FEAT_EICKT)
  set(kernel-lafem-list
    ${kernel-lafem-list}
    dense_vector-eickt.cpp
    sparse_vector-eickt.cpp
    sparse_matrix_csr-eickt.cpp
    sparse_matrix_cscr-eickt.cpp
    )
endif (FEAT_EICKT)

add_library(kernel-lafem STATIC ${kernel-lafem-list})
if (FEAT_COTIRE)
  cotire(kernel-lafem)
endif (FEAT_COTIRE)

target_link_libraries(kernel-lafem PRIVATE kernel-lafem-arch kernel-util kernel-adjacency)

if (FEAT_HAVE_CUDA)
  target_link_libraries(kernel-lafem PRIVATE kernel-util-cuda)
endif (FEAT_HAVE_CUDA)

# list of lafem tests
set(test_list
  dense_matrix-test
  dense_vector-test
  dense_vector_blocked-test
  matrix_mirror_buffer-test
  matrix_mult-test
  mean_filter-test
  mean_filter_blocked-test
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
  slip_filter-test
  sparse_matrix_conversion-test
  sparse_matrix_banded-test
  sparse_matrix_csr-test
  sparse_matrix_cscr-test
  sparse_matrix_bcsr-test
  sparse_matrix_factory-test
  sparse_vector-test
  sparse_vector_blocked-test
  test_matrix_factory-test
  unit_filter-test
  unit_filter_blocked-test
  vector_mirror-test
)

# create all tests
foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  target_link_libraries(${test} feat test_system)

  add_test(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/lafem/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/lafem/${test} generic mkl)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

  if (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
    add_test(${test}_cuda_memcheck ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${CUDA_MEMCHECK_EXE} ${FEAT_BINARY_DIR}/kernel/lafem/${test} cuda)
    set_property(TEST ${test}_cuda_memcheck PROPERTY LABELS "cuda_memcheck")
    set_property(TEST ${test}_cuda_memcheck PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
    set_property(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
    set_property(TEST ${test}_cuda_memcheck PROPERTY FAIL_REGULAR_EXPRESSION "= Leaked")
  endif (FEAT_CUDAMEMCHECK AND FEAT_HAVE_CUDA)
endforeach (test)

# add all tests to lafem_tests
add_custom_target(lafem_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
add_dependencies(tests lafem_tests)
