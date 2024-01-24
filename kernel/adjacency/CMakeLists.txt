# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# list of adjacency source files
ADD_LIBRARY(kernel-adjacency STATIC
  coloring.cpp
  cuthill_mckee.cpp
  graph.cpp
  permutation.cpp
)
if (FEAT_COTIRE)
  cotire(kernel-adjacency)
endif (FEAT_COTIRE)

# list of adjacency tests
SET (test_list
  coloring-test
  cuthill_mckee-test
  graph-test
  permutation-test
)

TARGET_LINK_LIBRARIES (kernel-adjacency kernel-util)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)

  ADD_TEST(${test}_generic ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/adjacency/${test} generic)
  SET_PROPERTY(TEST ${test}_generic PROPERTY LABELS "generic")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_generic PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_generic PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

ENDFOREACH(test)

# add all tests to adjacency_tests
ADD_CUSTOM_TARGET(adjacency_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests adjacency_tests)
