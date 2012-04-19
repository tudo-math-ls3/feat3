cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

ADD_LIBRARY(kernel-geometry STATIC)
TARGET_LINK_LIBRARIES (kernel-geometry kernel-util)

# list of util tests
SET (test_list shape-test standard_refinery-test-conf-quad standard_refinery-test-struct-quad standard_refinery-test-conf-hexa)

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
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_SOURCE_DIR}/test_system/run_mpi" "${FEAST_BINARY_DIR}/kernel/geometry/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "PASSED")
  SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
ENDFOREACH(test)

# add all tests to util_tests
ADD_CUSTOM_TARGET(geometry_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests geometry_tests)
