cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

#list of hornet source files
ADD_LIBRARY(kernel-hornet STATIC memory_pool.cpp)

#list of hornet tests
SET (test_list dense_vector-test dense_matrix-test sparse_matrix_coo-test sparse_matrix_csr-test)

#create all tests
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
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_SOURCE_DIR}/test_system/run_mpi" "${FEAST_BINARY_DIR}/kernel/hornet/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "PASSED")

  ADD_TEST(${test}_cpu ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_SOURCE_DIR}/test_system/run_mpi" "${FEAST_BINARY_DIR}/kernel/hornet/${test}" cpu)
  SET_PROPERTY(TEST ${test}_cpu PROPERTY LABELS "cpu")
  SET_PROPERTY(TEST ${test}_cpu PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  SET_PROPERTY(TEST ${test}_cpu PROPERTY PASS_REGULAR_EXPRESSION "PASSED")
ENDFOREACH(test)

#add all tests to hornet_tests
ADD_CUSTOM_TARGET(hornet_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests hornet_tests)
