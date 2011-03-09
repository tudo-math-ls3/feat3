cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

#our test_system basic library
ADD_LIBRARY( test_system SHARED test_system )

#list of test_system tests
SET ( test_list test_system-test )

#create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} test_system)
  ADD_TEST(${test}_nil ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_SOURCE_DIR}/test_system/run_mpi" "${FEAST_BINARY_DIR}/test_system/${test}" Nil)
  SET_PROPERTY(TEST ${test}_nil PROPERTY LABELS "Nil")
ENDFOREACH(test)

#add all tests to test_system_tests
ADD_CUSTOM_TARGET(test_system_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests test_system_tests)
