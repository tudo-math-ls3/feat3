cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR})

ADD_SUBDIRECTORY( util )
ADD_SUBDIRECTORY( base_mesh )

#list of kernel tests
SET (test_list base_header-test logger-test work_group-test_1d work_group-test_2d_1 work_group-test_2d_2)

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
    --test-command "${FEAST_SOURCE_DIR}/test_system/run_mpi" "${FEAST_BINARY_DIR}/kernel/${test}" Nil)
  SET_PROPERTY(TEST ${test}_nil PROPERTY LABELS "Nil")
ENDFOREACH(test)

#add all tests to kernel_tests
ADD_CUSTOM_TARGET(kernel_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests kernel_tests)
