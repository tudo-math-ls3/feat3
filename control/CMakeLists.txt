# CmakeLists for directory "control"

#enable compiler output
set(CMAKE_VERBOSE_MAKEFILE ON)

#list of test_system tests
set(test_list checkpoint-test)

foreach (test ${test_list} )
  add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  target_link_libraries(${test} PRIVATE feat test_system)

  add_test(${test}_all ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${FEAT_CTEST_RUNNER} ${FEAT_BINARY_DIR}/control/${test})
  set_property(TEST ${test}_all PROPERTY LABELS "all")
  if (FEAT_VALGRIND)
    add_test(${test}_valgrind ${CMAKE_CTEST_COMMAND}
      --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
      --build-generator ${CMAKE_GENERATOR}
      --build-makeprogram ${CMAKE_MAKE_PROGRAM}
      --build-target ${test}
      --build-nocmake
      --build-noclean
      --test-command ${FEAT_CTEST_RUNNER} ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/control/${test} generic)
    set_property(TEST ${test}_valgrind PROPERTY LABELS "valgrind")
    set_property(TEST ${test}_valgrind PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    set_property(TEST ${test}_valgrind PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)
endforeach (test)

#add all tests to test_system_tests
add_custom_target(control_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
add_dependencies(tests control_tests)
