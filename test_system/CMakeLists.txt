cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

ADD_LIBRARY(test_system STATIC test_system.cpp)

if (FEAST_BACKENDS_CUDA)
  CUDA_ADD_LIBRARY(test_system-cuda STATIC
    cuda.cu
  )
  ADD_DEPENDENCIES(test_system test_system-cuda)
  TARGET_LINK_LIBRARIES(test_system test_system-cuda)
endif(FEAST_BACKENDS_CUDA)

#list of test_system tests
SET ( test_list test_system-test )

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
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/test_system/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")
  if (FEAST_VALGRIND)
    SET_PROPERTY(TEST ${test}_none PROPERTY PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors from")
    SET_PROPERTY(TEST ${test}_none PROPERTY FAIL_REGULAR_EXPRESSION "FAILED")
  endif (FEAST_VALGRIND)
ENDFOREACH(test)

#add all tests to test_system_tests
ADD_CUSTOM_TARGET(test_system_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests test_system_tests)
