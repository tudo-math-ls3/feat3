# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#ADD_LIBRARY(kernel-assembly STATIC)
#TARGET_LINK_LIBRARIES (kernel-assembly kernel-util)

# list of assembly tests
SET (test_list
  bilinear_operator-test
  discrete_evaluator-test
  grid_transfer-test
  grid_transfer-mass-test
  grid_transfer-tip-test
  interpolator-test
  jump_stabil-test
  linear_functional-test
  mean_filter-test
  rew_projector-test
)

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
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/assembly/${test} generic)
  SET_PROPERTY(TEST ${test}_generic PROPERTY LABELS "generic")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_generic PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_generic PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

ENDFOREACH(test)

# add all tests to assembly_tests
ADD_CUSTOM_TARGET(assembly_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests assembly_tests)
