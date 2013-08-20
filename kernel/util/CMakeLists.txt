cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

# list of util source files
ADD_LIBRARY(kernel-util STATIC
  colouring.cpp
  exception.cpp
  mesh_streamer.cpp
  param_section.cpp
)

# list of util tests
SET (test_list
  assertion-test
  binary_stream-test
  colouring-test
  cuthill_mckee-test
  graph-test
  inheritance-test
  linear_algebra-test
  math-test
  mesh_streamer-test
  meta_math-test
  param_section-test
  permutation-test
  string-test
  tiny_algebra-test
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
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/util/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "PASSED")
  SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
ENDFOREACH(test)

# add all tests to util_tests
ADD_CUSTOM_TARGET(util_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests util_tests)
