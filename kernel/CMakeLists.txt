cmake_minimum_required (VERSION 2.8)

# enable compiler output
set (CMAKE_VERBOSE_MAKEFILE ON)

# top lvl include dir
INCLUDE_DIRECTORIES("${FEAST_BINARY_DIR}" "${FEAST_SOURCE_DIR}")

# include kernel subdirectories
ADD_SUBDIRECTORY( adjacency )
ADD_SUBDIRECTORY( assembly )
ADD_SUBDIRECTORY( cubature )
ADD_SUBDIRECTORY( foundation )
ADD_SUBDIRECTORY( scarc )
ADD_SUBDIRECTORY( geometry )
ADD_SUBDIRECTORY( lafem )
ADD_SUBDIRECTORY( space )
ADD_SUBDIRECTORY( trafo )
ADD_SUBDIRECTORY( util )

# list of kernel source files
ADD_LIBRARY(feast STATIC archs.cpp)

# list of libraries to include in feast lib
TARGET_LINK_LIBRARIES (feast
  kernel-adjacency
# kernel-assembly
# kernel-cubature
  kernel-foundation
  kernel-scarc
# kernel-geometry
  kernel-lafem
# kernel-space
# kernel-trafo
  kernel-util
)

if (FEAST_HAVE_ALGLIB)
  TARGET_LINK_LIBRARIES (feast thirdparty-alglib)
endif(FEAST_HAVE_ALGLIB)

if (FEAST_HAVE_UMFPACK)
  TARGET_LINK_LIBRARIES (feast thirdparty-umfpack)
endif(FEAST_HAVE_UMFPACK)

# list of kernel tests
SET (test_list
  shape-test
)

# create all tests
FOREACH (test ${test_list} )
  ADD_EXECUTABLE(${test} EXCLUDE_FROM_ALL ${test}.cpp)
  TARGET_LINK_LIBRARIES(${test} feast test_system)
  ADD_TEST(${test}_none ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${CURRENT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/${test} none)
  SET_PROPERTY(TEST ${test}_none PROPERTY LABELS "none")

  ADD_TEST(${test}_main ${CMAKE_CTEST_COMMAND}
    --build-and-test "${FEAST_SOURCE_DIR}" "${CURRENT_BINARY_DIR}"
    --build-generator ${CMAKE_GENERATOR}
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-target ${test}
    --build-nocmake
    --build-noclean
    --build-exe-dir "${CURRENT_BINARY_DIR}"
    --test-command ${VALGRIND_EXE} ${FEAST_BINARY_DIR}/kernel/${test} main)
  SET_PROPERTY(TEST ${test}_main PROPERTY LABELS "main")
ENDFOREACH(test)

# add all tests to kernel_tests
ADD_CUSTOM_TARGET(kernel_tests DEPENDS ${test_list})

# build all tests through top lvl target tests
ADD_DEPENDENCIES(tests kernel_tests)
