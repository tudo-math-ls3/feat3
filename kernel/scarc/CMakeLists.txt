cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})


ADD_LIBRARY(kernel-scarc STATIC scarc_error.cpp)
TARGET_LINK_LIBRARIES (kernel-scarc kernel-foundation kernel-util kernel-lafem)

#list of kernel tests
SET (test_list matrix-test vector-test solver_functor-test solver_pattern-test solver_expression-test)

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
    --test-command "${FEAST_BINARY_DIR}/kernel/scarc/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")

  ADD_TEST(${test}_cpu ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/scarc/${test}" cpu)
  SET_PROPERTY(TEST ${test}_cpu PROPERTY LABELS "cpu")
ENDFOREACH(test)

#add all tests to scarc_tests
ADD_CUSTOM_TARGET(scarc_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests scarc_tests)
