cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

# list of adjacency source files
ADD_LIBRARY(kernel-adjacency STATIC
  colouring.cpp
  cuthill_mckee.cpp
  graph.cpp
  permutation.cpp
)

# list of adjacency tests
SET (test_list
  colouring-test
  cuthill_mckee-test
  graph-test
  permutation-test
)

TARGET_LINK_LIBRARIES (kernel-adjacency kernel-util)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feast test_system)
  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${CURRENT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/adjacency/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "PASSED")
  SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
ENDFOREACH(test)

# add all tests to adjacency_tests
ADD_CUSTOM_TARGET(adjacency_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests adjacency_tests)
