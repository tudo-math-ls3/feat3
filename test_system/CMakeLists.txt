# CmakeLists for directory "test_system"

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

ADD_LIBRARY(test_system STATIC test_system.cpp)

TARGET_LINK_LIBRARIES(test_system feat)

#list of test_system tests
SET ( test_list test_system-test)

#create all tests
FOREACH (test ${test_list})
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feat test_system)
  ADD_TEST(${test}_generic ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAT_SOURCE_DIR}" "${FEAT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --test-command ${VALGRIND_EXE} ${FEAT_BINARY_DIR}/test_system/${test} generic)
  SET_PROPERTY(TEST ${test}_generic PROPERTY LABELS "generic")
  if (FEAT_VALGRIND)
    SET_PROPERTY(TEST ${test}_generic PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_generic PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAT_VALGRIND)
ENDFOREACH(test)

#add all tests to test_system_tests
ADD_CUSTOM_TARGET(test_system_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests test_system_tests)
