# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#ADD_LIBRARY(kernel-space STATIC)
#TARGET_LINK_LIBRARIES (kernel-space kernel-util)

# list of space tests
SET (test_list
  bogner_fox_schmit-test
  discontinuous-test
  element-regression-test
  lagrange1-test
  rannacher_turek-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system kernel-geometry-test-aux)

  ADD_TEST(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/kernel/space/${test})
  SET_PROPERTY(TEST ${test}_all PROPERTY LABELS "all")

  if (FEAT_VALGRIND)
    ADD_TEST(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/kernel/space/${test} generic)
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)

ENDFOREACH(test)

# add all tests to space_tests
ADD_CUSTOM_TARGET(space_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests space_tests)
