cmake_minimum_required (VERSION 2.8)

#enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

#top lvl include dir
INCLUDE_DIRECTORIES(${FEAST_SOURCE_DIR} ${FEAST_BINARY_DIR})

ADD_SUBDIRECTORY( util )
ADD_SUBDIRECTORY( lafem )
ADD_SUBDIRECTORY( geometry )
ADD_SUBDIRECTORY( foundation )

#list of kernel source files
ADD_LIBRARY(feast STATIC archs.cpp)
#list of libraries to include in libfeast
TARGET_LINK_LIBRARIES (feast kernel-util kernel-lafem kernel-geometry kernel-foundation)

#list of kernel tests
SET (test_list)

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
    --test-command "${FEAST_BINARY_DIR}/kernel/${test}" none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")

  ADD_TEST(${test}_cpu ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${FEAST_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command "${FEAST_BINARY_DIR}/kernel/${test}" cpu)
  SET_PROPERTY(TEST ${test}_cpu PROPERTY LABELS "cpu")
ENDFOREACH(test)

#add all tests to kernel_tests
ADD_CUSTOM_TARGET(kernel_tests DEPENDS ${test_list})
#build all tests through top lvl target tests
ADD_DEPENDENCIES(tests kernel_tests)
